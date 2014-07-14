### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-Results")

### Set seed
set.seed(1)

### Import required libraries
require(pROC)
require(data.table)
require(OptimalCutpoints)
require(ggplot2)

### Import results
RE_Results_OR_frame <- read.csv("RE_Results_OR.csv")
RE_Results_OR <- data.table(RE_Results_OR_frame)

### Introducing new variables

RE_Results_OR$I_U_Est_Diff <- RE_Results_OR$Up_Est - RE_Results_OR$Init_Est 

### Summary tables

OR_Means <- RE_Results_OR[, lapply(.SD, mean, na.rm=TRUE), by=list(OR_New, Het_New, Num_Up)]
OR_SDs <- RE_Results_OR[, lapply(.SD, sd, na.rm=TRUE), by=list(OR_New, Het_New, Num_Up)]
OR_Means_HetVary <- RE_Results_OR[, lapply(.SD, mean, na.rm=TRUE), by=list(OR_New, Num_Up)]


### Write summary tables to csv

write.csv(OR_Means, file = "OR_Means.csv")
write.csv(OR_SDs, file = "OR_SDs.csv")

### Graphical summary methods

## Updated meta-analysis
y <- ggplot(OR_Means, aes(x = OR_New, y = Up_Est, colour = as.factor(Num_Up)))
ci <- qnorm(0.975)* OR_SDs$Up_Est /sqrt(5000)
y + geom_errorbar(aes(ymin=Up_Est-ci, ymax=Up_Est+ci), width=.1) + geom_point() + scale_x_continuous(breaks= seq(1,2,0.2)) + xlab("Odds Ratio") + ylab("Updated Estimate")

y <- ggplot(OR_Means, aes(x = OR_New, y = Up_Est_p, colour = as.factor(Num_Up)))
ci <- qnorm(0.975)* OR_SDs$Up_Est_p /sqrt(5000)
y + geom_errorbar(aes(ymin=Up_Est_p-ci, ymax=Up_Est_p+ci), width=.1) + geom_point(aes(shape = factor(Het_New))) + scale_shape(solid = TRUE) + scale_x_continuous(breaks= seq(1,2,0.2)) + xlab("Odds Ratio") + ylab("Updated Estimate p")+ geom_abline(intercept = 0.05, slope = 0)

y <- ggplot(OR_Means, aes(x = OR_New, y = Up_tau2, colour = as.factor(Num_Up)))
ci <- qnorm(0.975)* OR_Means$Up_tau2_SE
y + geom_errorbar(aes(ymin=Up_tau2 - ci, ymax=Up_tau2 + ci), width=.1) + geom_point(aes(shape = factor(Het_New))) + scale_x_continuous(breaks= seq(1,2,0.2)) + xlab("Odds Ratio") + ylab("Updated Estimate")

## Simplified updated
y <- ggplot(OR_Means_HetVary, aes(x = OR_New, y = Up_Est, colour = as.factor(Num_Up)))
y + geom_errorbar(aes(ymin=Up_CI_lb, ymax=Up_CI_ub), width=.1) + geom_point() + scale_x_continuous(breaks= seq(1,2,0.2)) + xlab("Odds Ratio") + ylab("Updated Estimate") + geom_abline(intercept = 0, slope = 0)

## Meta-regression
y <- ggplot(OR_Means, aes(x = OR_New, y = MLM_Mod_Est, colour = as.factor(Num_Up)))
y + geom_errorbar(aes(ymin=MLM_Mod_CI_lb, ymax=MLM_Mod_CI_ub), width=.1) + geom_point() + scale_x_continuous(breaks= seq(1,2,0.2)) + xlab("Odds Ratio") + ylab("Meta-regression covairate")+ geom_abline(intercept = 0, slope = 0, alpha = 0.5)

y <- ggplot(OR_Means, aes(x = OR_New, y = MLM_tau2, colour = as.factor(Num_Up)))
ci <- qnorm(0.975)* OR_SDs$MLM_tau2 /sqrt(5000)
y + geom_errorbar(aes(ymin=MLM_tau2-ci, ymax=MLM_tau2+ci), width=.1) + geom_point(aes(shape = factor(Het_New))) + scale_shape(solid = TRUE) + scale_x_continuous(breaks= seq(1,2,0.2)) + xlab("Odds Ratio") + ylab("Meta-regression tau2")

y <- ggplot(OR_Means, aes(x = OR_New, y = MLM_Mod_p, colour = as.factor(Num_Up)))
ci <- qnorm(0.975)* OR_SDs$MLM_Mod_p /sqrt(5000)
y + geom_errorbar(aes(ymin=MLM_Mod_p-ci, ymax=MLM_Mod_p+ci), width=.1) + geom_point(aes(shape = factor(Het_New))) + scale_shape(solid = TRUE) + scale_x_continuous(breaks= seq(1,2,0.2)) + xlab("Odds Ratio") + ylab("Meta-regression covariate pval")+ geom_abline(intercept = 0.05, slope = 0)


## Network
y <- ggplot(OR_Means, aes(x = OR_New, y = Net_TE, colour = as.factor(Num_Up)))
y + geom_errorbar(aes(ymin=Net_TE_CI_lb, ymax=Net_TE_CI_ub), width=.1, alpha = 1/2) + geom_point() + scale_x_continuous(breaks= seq(1,2,0.2)) + xlab("Odds Ratio") + ylab("Network Estimate") + geom_abline(intercept = 0, slope = 0, alpha = 0.5)

y <- ggplot(OR_Means, aes(x = OR_New, y = Net_TE_tau2, colour = as.factor(Num_Up)))
ci <- qnorm(0.975)* OR_SDs$Net_TE_tau2 /sqrt(5000)
y + geom_errorbar(aes(ymin=Net_TE_tau2-ci, ymax=Net_TE_tau2+ci), width=.1) + geom_point(aes(shape = factor(Het_New))) + scale_shape(solid = TRUE) + scale_x_continuous(breaks= seq(1,2,0.2)) + xlab("Odds Ratio") + ylab("Network tau2")

y <- ggplot(OR_Means, aes(x = OR_New, y = Net_TE_pval, colour = as.factor(Num_Up)))
ci <- qnorm(0.975)* OR_SDs$Net_TE_pval /sqrt(5000)
y + geom_errorbar(aes(ymin=Net_TE_pval-ci, ymax=Net_TE_pval+ci), width=.1) + geom_point(aes(shape = factor(Het_New))) + scale_shape(solid = TRUE) + scale_x_continuous(breaks= seq(1,2,0.2)) + xlab("Odds Ratio") + ylab("Network pval")+ geom_abline(intercept = 0.05, slope = 0)


### Defining a change in Estimate and Heterogeneity (lower and higher)

Est_Change <- integer(length = length(RE_Results_OR$OR_New))

for (i in 1:length(RE_Results_OR$OR_New)) Est_Change[i] <- sum(RE_Results_OR$OR_New[i] > 0)

Het_Change_Lower <- numeric(length = length(RE_Results_OR$Het_New))

for (i in 1:length(RE_Results_OR$Het_New)) Het_Change_Lower[i] <- sum(RE_Results_OR$Het_New[i] < 0.3)

Het_Change_Higher <- numeric(length = length(RE_Results_OR$Het_New))

for (i in 1:length(RE_Results_OR$Het_New)) Het_Change_Higher[i] <- sum(RE_Results_OR$Het_New[i] > 0.3)

### Method using absolute difference in mean reuslts

Abs_Diff_Est <- RE_Results_OR$Init_Est - RE_Results_OR$Up_Est
ROC_Abs_Diff_Est <- roc(Est_Change ~ Abs_Diff_Est)
plot(ROC_Abs_Diff_Est)
ROC_Abs_Diff_Est




### Playing with optimising sensitivity and specificity

product <- ROC_Abs_Diff_Est$sensitivities*ROC_Abs_Diff_Est$specificities
test_location <- match(max(product), product)
ROC_Abs_Diff_Est$sensitivities[test_location]
ROC_Abs_Diff_Est$specificities[test_location]
ROC_Abs_Diff_Est$thresholds[test_location]

## This can't handle the data size
data_test <- data.table(Abs_Diff_Est, Est_Change)
test <- optimal.cutpoints(X = Abs_Diff_Est~Est_Change, tag.healthy = 0, methods = "PROC01", data = data_test)

### Method using relative difference in mean results

Abs_Diff_Est_byinitSE <- Abs_Diff_Est/RE_Results_OR$Init_Est_SE
ROC_Diff_Est_byinitSE <- roc(Est_Change ~ Abs_Diff_Est_byinitSE)
plot(ROC_Diff_Est_byinitSE)

Abs_Diff_Est_byupSE <- Abs_Diff_Est/RE_Results_OR$Up_Est_SE
ROC_Diff_Est_byupSE <- roc(Est_Change ~ Abs_Diff_Est_byupSE)
plot(ROC_Diff_Est_byupSE)

### Method using change in variance

Abs_Diff_tau2 <- RE_Results_OR$Init_tau2 - RE_Results_OR$Up_tau2
ROC_Diff_tau2_Lower <- roc(Het_Change_Lower ~ Abs_Diff_tau2)
plot(ROC_Diff_tau2)

### Method using relative change in variance

Abs_Diff_tau2_byinitSE <- Abs_Diff_tau2/RE_Results_OR$Init_tau2_SE
Abs_Diff_tau2_byupSE <- Abs_Diff_tau2/RE_Results_OR$Up_tau2_SE

### Method using changing width of CIs

### Difference in AIC, BIC