### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-Results")

### Set seed
set.seed(1)

### Import required libraries
require(pROC)
require(data.table)
require(OptimalCutpoints)

### Import results
RE_Results_MD_frame <- read.csv("RE_Results_MD.csv")
RE_Results_MD <- data.table(RE_Results_MD_frame)

### Summary tables

MD_Means <- RE_Results_MD[, lapply(.SD, mean, na.rm=TRUE), by=list(M_D, Het_New, Num_Up)]
MD_SDs <- RE_Results_MD[, lapply(.SD, sd, na.rm=TRUE), by=list(M_D, Het_New, Num_Up)]
MD_SE <- 

### Write summary tables to csv

write.csv(MD_Means, file = "MD_Means.csv")
write.csv(MD_SDs, file = "MD_SDs.csv")


### Defining a change in Estimate and Heterogeneity (lower and higher)

Est_Change <- integer(length = length(RE_Results_MD$M_D))

for (i in 1:length(RE_Results_MD$M_D)) Est_Change[i] <- sum(RE_Results_MD$M_D[i] > 0)

Het_Change_Lower <- numeric(length = length(RE_Results_MD$Het_New))

for (i in 1:length(RE_Results_MD$Het_New)) Het_Change_Lower[i] <- sum(RE_Results_MD$Het_New[i] < 0.3)

Het_Change_Higher <- numeric(length = length(RE_Results_MD$Het_New))

for (i in 1:length(RE_Results_MD$Het_New)) Het_Change_Higher[i] <- sum(RE_Results_MD$Het_New[i] > 0.3)

### Method using absolute difference in mean reuslts

Abs_Diff_Est <- RE_Results_MD$Init_Est - RE_Results_MD$Up_Est
ROC_Abs_Diff_Est <- roc(Est_Change ~ Abs_Diff_Est)
plot(ROC_Abs_Diff_Est)
ROC_Abs_Diff_Est


### Playing with optimising sensitivity and specificity


### Graphics to adapt, requires a pre-calculation of mean and ci
y <- ggplot(subset(MD_Means, Het_New == 0.3), aes(x = M_D, y = Net_TE, colour = Num_Up))
y + geom_errorbar(aes(ymin=Net_TE-ci, ymax=Net_TE+ci), width=.1, position= "dodge") + geom_point()

product <- ROC_Abs_Diff_Est$sensitivities*ROC_Abs_Diff_Est$specificities
test_location <- match(max(product), product)
ROC_Abs_Diff_Est$sensitivities[test_location]
ROC_Abs_Diff_Est$specificities[test_location]
ROC_Abs_Diff_Est$thresholds[test_location]


## This can't handle the data size

data_test <- data.table(Abs_Diff_Est, Est_Change)
test <- optimal.cutpoints(X = Abs_Diff_Est~Est_Change, tag.healthy = 0, methods = "PROC01", data = data_test)

### Method using relative difference in mean results

Abs_Diff_Est_byinitSE <- Abs_Diff_Est/RE_Results_MD$Init_Est_SE
ROC_Diff_Est_byinitSE <- roc(Est_Change ~ Abs_Diff_Est_byinitSE)
plot(ROC_Diff_Est_byinitSE)

Abs_Diff_Est_byupSE <- Abs_Diff_Est/RE_Results_MD$Up_Est_SE
ROC_Diff_Est_byupSE <- roc(Est_Change ~ Abs_Diff_Est_byupSE)
plot(ROC_Diff_Est_byupSE)

### Method using change in variance

Abs_Diff_tau2 <- RE_Results_MD$Init_tau2 - RE_Results_MD$Up_tau2
ROC_Diff_tau2_Lower <- roc(Het_Change_Lower ~ Abs_Diff_tau2)
plot(ROC_Diff_tau2)

### Method using relative change in variance

Abs_Diff_tau2_byinitSE <- Abs_Diff_tau2/RE_Results_MD$Init_tau2_SE
Abs_Diff_tau2_byupSE <- Abs_Diff_tau2/RE_Results_MD$Up_tau2_SE

### Method using changing width of CIs

### Difference in AIC, BIC