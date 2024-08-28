## Cox Model Survival Analysis

### About the Dataset
The dataset, `haberman.csv`, includes data from a study conducted between 1958 and 1970 at the University of Chicago's Billings Hospital. It focuses on the survival of breast cancer surgery patients.

- **Number of Instances:** 306
- **Number of Attributes:** 4
  - Age of patient at time of operation (numerical)
  - Patient's year of operation (year - 1900, numerical)
  - Number of positive axillary nodes detected (numerical)
  - Survival status (class attribute): 
    - `1` = Survived 5 years or longer
    - `2` = Died within 5 years
- **Missing Values:** None

### Files
- `haberman.csv`: The dataset used for survival analysis.
- `cox-survival.Rmd`: R Markdown file for performing survival analysis.

### Usability
- **License:** Unknown
- **Expected Update Frequency:** Not specified

### Analysis Overview
1. **Introduction:** Explains the dataset and the objectives of the survival analysis.
2. **Load Required Libraries:** Lists and loads the libraries needed for the analysis.
3. **Data Import and Preview:** Imports the dataset and provides an initial preview.
4. **Data Profiling:** Provides diagnostic statistics for numeric variables.
5. **Defining Time and Event Variables:** Sets up variables for survival analysis.
6. **Kaplan-Meier Estimator:** Fits and plots the Kaplan-Meier survival curve.
7. **Stratified Kaplan-Meier Curves:** Creates and plots survival curves for different patient cohorts based on positive axillary nodes detected.
8. **Cox Proportional Hazards Model:** Fits a Cox model to the data and summarizes the results.
9. **Visualizing Cox Model Coefficients:** Generates a forest plot of Cox model coefficients.
10. **Predicting Survival Curves for Specific Patients:** Predicts and plots survival curves for selected patients.

### Example Code
Below is an example of loading the dataset and performing a Kaplan-Meier survival analysis:

```r
library(tidyverse)
library(survival)
library(survminer)
library(dlookr)
library(gridExtra)

# Load the dataset
data <- read.csv('haberman.csv', header = FALSE)
colnames(data) <- c('Age', 'Operation_year', 'Nb_pos_detected', 'Surv')

# Kaplan-Meier Estimator
km_fit <- survfit(Surv(Age, Surv) ~ 1, data = data)
ggsurvplot(km_fit, data = data, conf.int = FALSE, ggtheme = theme_minimal(), title = "Kaplan-Meier Estimate")
```

### Conclusion
The analysis explores patient survival probabilities using Kaplan-Meier estimation and Cox Proportional Hazards modeling. It provides insights into how various factors influence survival outcomes for breast cancer patients.

---

Feel free to adjust any section or add additional details as needed!
