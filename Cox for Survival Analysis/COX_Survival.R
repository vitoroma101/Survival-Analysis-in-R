library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(prodlim)
library(pec)

# Create vectors for each column
unemployment_duration <- c(36, 25.32, 46.2, 10.32, 30.72, 9.6, 32.54, 60, 8.23, 60, 
                           2.47, 12, 2.59, 1.66, 1.82, 12.94, 4.75, 15.49, 10.08, 24,
                           30, 3.33, 54.21, 4.96, 11.36, 6, 5.46, 5.67, 19.03, 30, 3.7, 
                           35.5, 11.4, 24, 32.6, 60, 7.03, 3.28, 5.23, 11.84)
gender <- c("FEMALE", "FEMALE", "FEMALE", "FEMALE", "FEMALE", "FEMALE", "FEMALE", "FEMALE",
            "FEMALE", "FEMALE", "FEMALE", "FEMALE", "FEMALE", "FEMALE", "FEMALE", "FEMALE",
            "FEMALE", "FEMALE", "FEMALE", "FEMALE", "MALE", "MALE", "MALE", "MALE",
            "MALE", "MALE", "MALE", "MALE", "MALE", "MALE", "MALE", "MALE",
            "MALE", "MALE", "MALE", "MALE", "MALE", "MALE", "MALE", "MALE")
age <- c(40, 48, 50, 50, 39, 45, 45, 60, 60, 60, 18, 18, 20, 20, 30, 30, 30, 45, 45, 55,
         18, 20, 20, 25, 25, 25, 40, 40, 45, 50, 50, 50, 50, 60, 62, 62, 28, 30, 45, 45)
censorship <- c(1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1,
                0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0)
sector <- c("OTHER", "OTHER", "OTHER", "OTHER", "AGRICULTURE", "AGRICULTURE", "AGRICULTURE",
            "AGRICULTURE", "CONSTRUCTION", "CONSTRUCTION", "SERVICES", "SERVICES", "SERVICES",
            "SERVICES", "SERVICES", "SERVICES", "SERVICES", "SERVICES", "SERVICES", "SERVICES",
            "OTHER", "OTHER", "OTHER", "OTHER", "OTHER", "OTHER", "AGRICULTURE", "AGRICULTURE",
            "AGRICULTURE", "AGRICULTURE", "AGRICULTURE", "CONSTRUCTION", "CONSTRUCTION", "CONSTRUCTION",
            "CONSTRUCTION", "CONSTRUCTION", "SERVICES", "SERVICES", "SERVICES", "SERVICES")

# Create the dataframe with the created vectors
data <- data.frame(Unemployment_Duration = unemployment_duration,
                   Gender = gender,
                   Age = age,
                   Censorship = censorship,
                   Sector = sector)

# MULTIVARIATE COX REGRESSION

res.cox <- coxph(Surv(unemployment_duration, censorship) ~ gender + age + sector, data = data)
summary(res.cox)
