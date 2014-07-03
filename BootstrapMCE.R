### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-Results")

### Set seed
set.seed(1)

### Import required libraries
require(boot)
require(data.table)

### Import results file
MD_Results_frame <- read.csv("RE_Results_MD.csv")
MD_Results <- data.table(MD_Results_frame)
setkey(MD_Results, "M_D", "Het_New", "Num_Up", "Rep_Number")

### Select specific set of results
Working_copy1 <- MD_Results[J(0, 0.3, 20)]$Init_Est
Working_copy2 <- MD_Results[J(0, 0.3, 20)]$Up_Est
Working_copy3 <- MD_Results[J(0, 0.3, 20)]$Net_TE
Working_copy4 <- MD_Results[J(0, 0.3, 20)]$MLM_Est

Target_MCE <- 0.002
Reps <- 100
Mult <- 50

bootstrapMCE <- function(data, indices){
  d <- data[indices]
  return(mean(d))
}

### Function to bootstrap an estimate for MCE by finding multiple sds

estimate_MCE <- function(p, Boot_Repeats, Multiplier, data, indices){

  d <- data[indices]
  
  sd_values <- vector(length = p)
  predictor1 <- vector(length = p)

  for (j in 1:p){
    Rj <- j*Multiplier
    predictor1[j] <- 1/sqrt(Rj)
    X_Rj <- sample(d , size = Rj, replace = FALSE, prob = NULL)
    results <- boot(data = X_Rj, R = Boot_Repeats, statistic=bootstrapMCE)
    sd_values[j] <- sd(results$t)
  }

  model1 <- lm(sd_values~ 0 + predictor1)
  return(model1$coefficients[[1]])
}

### Actually calculates an estimate for MCE = 0.002
Total_Results1 <- boot(data = Working_copy1, R = Reps, statistic = estimate_MCE, p = 3, Boot_Repeats = 100, Multiplier = Mult)
R1 <- (mean(Total_Results1$t)/Target_MCE)^2
print(R1)

Total_Results2 <- boot(data = Working_copy2, R = Reps, statistic = estimate_MCE, p = 3, Boot_Repeats = 100, Multiplier = Mult)
R2 <- (mean(Total_Results2$t)/Target_MCE)^2
print(R2)

Total_Results3 <- boot(data = Working_copy3, R = Reps, statistic = estimate_MCE, p = 3, Boot_Repeats = 100, Multiplier = Mult)
R3 <- (mean(Total_Results3$t)/Target_MCE)^2
print(R3)

Total_Results4 <- boot(data = Working_copy3, R = Reps, statistic = estimate_MCE, p = 3, Boot_Repeats = 100, Multiplier = Mult)
R4 <- (mean(Total_Results4$t)/Target_MCE)^2
print(R4)

### For MCE = 0.001
Target_MCE <- 0.001

R1 <- (mean(Total_Results1$t)/Target_MCE)^2
print(R1)

R2 <- (mean(Total_Results2$t)/Target_MCE)^2
print(R2)

R3 <- (mean(Total_Results3$t)/Target_MCE)^2
print(R3)

R4 <- (mean(Total_Results4$t)/Target_MCE)^2
print(R4)