library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(prodlim)
library(pec)
library(lubridate)

# Create the dataframe with the provided data
data <- data.frame(
  Case = 1:32,
  Entry_Date = c("2/87", "2/87", "2/87", "2/87", "3/87", "3/87", "4/87", "4/87", "4/87", "5/87",
                 "6/87", "6/87", "6/87", "7/87", "8/87", "9/87", "9/87", "10/87", "12/87", "12/87",
                 "1/88", "2/88", "2/88", "3/88", "5/88", "5/88", "5/88", "6/88", "7/88", "7/88",
                 "7/88", "9/88"),
  Last_Observation_Date = c("2/89*", "2/89*", "12/88", "12/87", "6/87", "9/87", "2/89*", "5/87",
                            "2/89*", "2/89*", "9/87", "8/87", "8/89*", "10/87", "2/89*", "12/88*",
                            "11/88", "3/88", "8/88", "2/88", "10/88", "2/89*", "1/90*", "12/88", 
                            "5/89*", "4/89*", "11/89*", "12/88", "2/90*", "2/90*", "12/89*", "10/89*"),
  Treatment = c("Placebo", "AZT", "AZT", "AZT", "Placebo", "Placebo", "AZT", "Placebo", "Placebo",
                "AZT", "AZT", "Placebo", "AZT", "Placebo", "AZT", "Placebo", "Placebo", "AZT", 
                "Placebo", "Placebo", "Placebo", "AZT", "AZT", "Placebo", "Placebo", "AZT", 
                "Placebo", "Placebo", "AZT", "AZT", "AZT", "AZT")
)

# A. Identify censored observations
data$Status <- ifelse(grepl("\\*$", data$Last_Observation_Date), "Censored", "Complete")

# Convert dates to proper date format
data$Entry_Date <- as.Date(paste0("01/", data$Entry_Date), format = "%d/%m/%y")
data$Last_Observation_Date <- as.Date(paste0("01/", data$Last_Observation_Date), format = "%d/%m/%y")

# Calculate participation time in months for each patient
data$Participation_Time_Months <- round(as.numeric(difftime(data$Last_Observation_Date, data$Entry_Date, units = "days")) / 30)

# Display participation time and status (complete/censored) for each patient
participation_time <- data.frame(
  Case = data$Case,
  Participation_Time_Months = data$Participation_Time_Months, 
  Status = data$Status
)

participation_time

table(data$Status)

# Note: Cases 25, 27, and 32 were classified as censored even though they don't have an asterisk, 
# because their last observation date is after the end of the study in February 1989. 
# Since the event was not recorded before this date, they are considered censored.

# B. Life Table (by semesters)
# Calculate participation time in semesters
data$Participation_Time_Semesters <- ceiling(data$Participation_Time_Months / 6)

# Display participation time in semesters for each patient
participation_time_semesters <- data.frame(
  Case = data$Case,
  Participation_Time_Semesters = data$Participation_Time_Semesters,
  Status = data$Status
)

# Create the semester intervals
max_semesters <- max(data$Participation_Time_Semesters, na.rm = TRUE)
intervals <- 1:max_semesters

# Initialize the life table with N'_j
life_table <- data.frame(
  Interval = intervals,
  Nj = rep(NA, length(intervals)),   # Individuals at risk
  Dj = rep(NA, length(intervals)),   # Number of deaths
  Wj = rep(NA, length(intervals)),   # Number of censored cases
  N_prime_j = rep(NA, length(intervals)), # Adjusted risk set
  qj = rep(NA, length(intervals)),   # Probability of death in the interval
  pj = rep(NA, length(intervals)),   # Probability of survival in the interval
  Pj = rep(NA, length(intervals)),   # Cumulative survival probability
  Var_qj = rep(NA, length(intervals)), # Variance of qj
  Var_Pj = rep(NA, length(intervals))  # Variance of Pj
)

# Calculate values for the life table
for (j in intervals) {
  # Individuals at risk at the beginning of the interval
  life_table$Nj[j] <- sum(data$Participation_Time_Semesters >= j)
  
  # Number of deaths in the interval
  life_table$Dj[j] <- sum(data$Participation_Time_Semesters == j & data$Status == "Complete")
  
  # Number of censored cases in the interval
  life_table$Wj[j] <- sum(data$Participation_Time_Semesters == j & data$Status == "Censored")
  
  # Adjusted risk set N'_j
  life_table$N_prime_j[j] <- life_table$Nj[j] - life_table$Wj[j] / 2
  
  # Calculate the probability of death in the interval
  life_table$qj[j] <- ifelse(life_table$N_prime_j[j] > 0, life_table$Dj[j] / life_table$N_prime_j[j], 0)
  
  # Probability of surviving the interval
  life_table$pj[j] <- 1 - life_table$qj[j]
  
  # Cumulative survival probability
  life_table$Pj[j] <- ifelse(j == 1, life_table$pj[j], life_table$Pj[j-1] * life_table$pj[j])
  
  # Greenwood's variance for qj
  life_table$Var_qj[j] <- ifelse(life_table$N_prime_j[j] > 0, (life_table$qj[j] * life_table$pj[j]) / life_table$N_prime_j[j], NA)
  
  # Cumulative variance of Pj
  life_table$Var_Pj[j] <- ifelse(j == 1, life_table$Var_qj[j], life_table$Pj[j]^2 * sum(life_table$qj[1:j] / life_table$N_prime_j[1:j] / life_table$pj[1:j]))
}

# Display the life table 
life_table

# C. General Kaplan-Meier
# The Kaplan-Meier estimator is calculated only at times when events occur (failures)
# The censoring between failure times affects the respective jump at that time
# The graph and values remain constant between jumps

data$Status_num <- ifelse(data$Status == 'Complete', 1, 0)

# Fit the Kaplan-Meier model
km_fit <- survfit(Surv(Participation_Time_Months, Status_num) ~ 1, data = data)

# Summarize the Kaplan-Meier fit results
summary(km_fit)

plot(km_fit)

# Customized Kaplan-Meier plot
ggsurvplot(
  km_fit,
  data = data,
  conf.int = TRUE,           # Show confidence intervals
  palette = c("#E7B800"),    # Change the curve color
  xlab = "Time in months",   # X-axis label
  ylab = "Survival Probability",  # Y-axis label
  title = "Kaplan-Meier Survival Curve",  # Title
  legend = "none"            # Remove the legend
)

# D. Kaplan-Meier by Group
# Fit the Kaplan-Meier model by treatment group
km_fit2 <- survfit(Surv(Participation_Time_Months, Status_num) ~ Treatment, data = data)
km_fit2
summary(km_fit2)

# Plot survival curves by treatment group
ggsurvplot(
  km_fit2,
  data = data,
  pval = TRUE,
  conf.int = TRUE,
  palette = c("#E7B800", "#2E9FDF"),
  xlab = "Time (months)",
  ylab = "Survival Probability",
  title = "Survival Curves by Treatment"
)
