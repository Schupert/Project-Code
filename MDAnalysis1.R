### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-SimCSV")

### Set seed
set.seed(1)

### Import required libraries

### Import results
RE_Results_MD <- read.csv("RE_Results_MD.csv")

### Defining a change
Est_Change <- integer(length = RE_Results_MD$M_D)

for (i in 1:length(RE_Results_MD$M_D)) Est_Change[i] <- sum(RE_Results_MD$M_D[i] > 0)

Het_Change <- numeric(length = RE_Results_MD$Het_Change)

for (i in 1:length(RE_Results_MD$Het_Change)) Het_Change[i] <- sum(RE_Results_MD$M_D[i] > 0)

### Method using absolute difference in mean reuslts

Abs_Diff_Est <- RE_Results_MD$Init_Est - RE_Results_MD$Up_Est

### Method using relative difference in mean results

Abs_Diff_Est_byinitSE <- Abs_Diff_Est/RE_Results_MD$Init_Est_SE
Abs_Diff_Est_byupSE <- Abs_Diff_Est/RE_Results_MD$Up_Est_SE

### Method using change in variance

Abs_Diff_tau2 <- RE_Results_MD$Init_tau2 - RE_Results_MD$Up_tau2

### Method using relative change in variance

Abs_Diff_tau2_byinitSE <- Abs_Diff_tau2/RE_Results_MD$Init_tau2_SE
Abs_Diff_tau2_byupSE <- Abs_Diff_tau2/RE_Results_MD$Up_tau2_SE

### Method using changing width of CIs

### Difference in AIC, BIC