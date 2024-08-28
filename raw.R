# Load required libraries
library(tidyverse)
library(survival)
library(survminer)
library(dlookr)
library(gridExtra)

# Read in the data
data <- read.csv('data.csv', header = FALSE)
colnames(data) <- c('Age', 'Operation_year', 'Nb_pos_detected', 'Surv')

# Preview the data
head(data, 5)

# Profiling report
diagnose_numeric(data)
# Define time and event variables
T <- data$Age
E <- data$Surv

# Kaplan-Meier Fitter
km_fit <- survfit(Surv(T, E) ~ 1, data = data)

# Plot the Kaplan-Meier curve
ggsurvplot(km_fit, data = data, conf.int = FALSE, ggtheme = theme_minimal(), title = "Kaplan-Meier Estimate")

# Create cohorts based on the number of positive axillary nodes detected
groups <- data$Nb_pos_detected
i1 <- groups >= 1
i2 <- groups < 1

# Kaplan-Meier curve for cohort 1 (at least one positive axillary node detected)
km_fit1 <- survfit(Surv(T[i1], E[i1]) ~ 1, data = data[i1, ])
p1 <- ggsurvplot(km_fit1, data=data, conf.int = FALSE, title = "At least one positive axillary node detected")

# Kaplan-Meier curve for cohort 2 (no positive axillary nodes detected)
km_fit2 <- survfit(Surv(T[i2], E[i2]) ~ 1, data = data[i2, ])
p2 <- ggsurvplot(km_fit2, data=data, conf.int = FALSE, title = "No positive axillary nodes detected")

# Arrange both plots together
grid.arrange(p1$plot, p2$plot, ncol = 2)

# Cox Proportional Hazards Model
cph_fit <- coxph(Surv(Age, Surv) ~ Nb_pos_detected + Operation_year, data = data)

# Summary of the Cox model
summary(cph_fit)

# Plot the Cox model coefficients
ggforest(cph_fit, data = data)

# Predict survival curves for specific patients
patients <- c(8, 155, 288)
selected_patients <- data[patients, c('Nb_pos_detected', 'Operation_year')]

# Predict and plot survival curves for selected patients
surv_curves <- survfit(cph_fit, newdata = selected_patients)
ggsurvplot(surv_curves, data=data, ggtheme = theme_minimal(), title = "Survival Curves for Selected Patients")
