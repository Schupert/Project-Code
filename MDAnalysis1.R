### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-Results")

### Set seed
set.seed(1)

### Import required libraries
require(pROC)
require(data.table)

### Import results
RE_Results_MD_frame <- read.csv("RE_Results_MD.csv")
RE_Results_MD <- data.table(RE_Results_MD_frame)

### Defining a change
Est_Change <- integer(length = length(RE_Results_MD$M_D))

for (i in 1:length(RE_Results_MD$M_D)) Est_Change[i] <- sum(RE_Results_MD$M_D[i] > 0)

Het_Change_Lower <- numeric(length = length(RE_Results_MD$Het_New))

for (i in 1:length(RE_Results_MD$Het_New)) Het_Change_Lower[i] <- sum(RE_Results_MD$Het_New[i] > 0.3)

### Method using absolute difference in mean reuslts

Abs_Diff_Est <- RE_Results_MD$Init_Est - RE_Results_MD$Up_Est
ROC_Abs_Diff_Est <- roc(Est_Change ~ Abs_Diff_Est)
plot(ROC_Abs_Diff_Est)

### Method using relative difference in mean results

Abs_Diff_Est_byinitSE <- Abs_Diff_Est/RE_Results_MD$Init_Est_SE
ROC_Diff_Est_byinitSE <- roc(Est_Change ~ Abs_Diff_Est_byinitSE)
plot(ROC_Diff_Est_byinitSE)

Abs_Diff_Est_byupSE <- Abs_Diff_Est/RE_Results_MD$Up_Est_SE
ROC_Diff_Est_byupSE <- roc(Est_Change ~ Abs_Diff_Est_byupSE)
plot(ROC_Diff_Est_byupSE)

### Method using change in variance

Abs_Diff_tau2 <- RE_Results_MD$Init_tau2 - RE_Results_MD$Up_tau2
ROC_Diff_tau2 <- roc(Het_Change ~ Abs_Diff_tau2)
plot(ROC_Diff_tau2)

### Method using relative change in variance

Abs_Diff_tau2_byinitSE <- Abs_Diff_tau2/RE_Results_MD$Init_tau2_SE
Abs_Diff_tau2_byupSE <- Abs_Diff_tau2/RE_Results_MD$Up_tau2_SE

### Method using changing width of CIs

### Difference in AIC, BIC