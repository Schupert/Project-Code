### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-SimCSV")

### Set seed
set.seed(1)

### Import required libraries
require(metafor)
require(data.table)
require(netmeta)
require(plyr)

### Import sample simulation, to determine length
MD_Results_frame <- read.csv("Total_Sim_MD.csv")
MD_Results <- data.table(MD_Results_frame)
setkey(MD_Results, "M_D", "Het_New", "Num_Up", "Rep_Number", "I_U")


MD_Analysis <- function(data_temp_all){
  data_temp_init <- data_temp_all[data_temp_all$I_U == "I",]
  NumUp <- mean(data_temp_all$Num_Up)
  MLM_var <- c(rep(0, 20), rep(1, NumUp))
  x <- rma(measure="MD", m1i = Group1Mean, sd1i = Group1SD, n1i = Group1Size, 
           m2i = Group2Mean, sd2i = Group2SD, n2i = Group2Size, 
           data=data_temp_init, method="SJ")
  y <- rma(measure="MD", m1i = Group1Mean, sd1i = Group1SD, n1i = Group1Size, 
           m2i = Group2Mean, sd2i = Group2SD, n2i = Group2Size, 
           data=data_temp_all, method="SJ")
  z <- rma(measure="MD", m1i = Group1Mean, sd1i = Group1SD, n1i = Group1Size, 
           m2i = Group2Mean, sd2i = Group2SD, n2i = Group2Size, 
           data=data_temp_all, method="SJ", mods = Init_Updated)
  Net_temp1 <- summary(escalc(measure="MD", m1i = Group1Mean, sd1i = Group1SD, n1i = Group1Size, 
                              m2i = Group2Mean, sd2i = Group2SD, n2i = Group2Size, 
                              data= data_temp_all, append = FALSE, replace = TRUE))
  
  Net_temp2 <- rbindlist(list(Net_temp1, Net_temp1[1:20,]))
  
  Net_treatlist1 <- c(rep("upd", 20 + k), rep("init", 20))
  
  Net_treatlist2 <- rep("plac", length(Net_treatlist1))
  
  zy <- netmeta(Net_temp2$yi, Net_temp2$sei, Net_treatlist2, Net_treatlist1, sm = "MD", warn = FALSE)
  Output <- c(
    Init_Est= x$b,
    Init_Est_SE= x$se,
    Init_Est_p= x$pval,
    Init_CI_lb = x$ci.lb,
    Init_CI_ub = x$ci.ub,
    Init_tau2 = x$tau2,
    Init_tau2_SE = x$se.tau2,
    Init_BIC = BIC.rma(x),
    Init_AIC = AIC.rma(x),
    Init_LogLik = logLik.rma(x),
    
    Up_Est = y$b,
    Up_Est_SE = y$se,
    Up_Est_p = y$pval,
    Up_CI_lb = y$ci.lb,
    Up_CI_ub = y$ci.ub,
    Up_tau2 = y$tau2,
    Up_tau2_SE = y$se.tau2,
    Up_BIC = BIC.rma(y),
    Up_AIC = AIC.rma(y),
    Up_LogLik = logLik.rma(y),
    
    MLM_Est = z$b[[1]],
    MLM_Est_SE = z$se[[1]], 
    MLM_Est_p = z$pval[[1]],
    MLM_CI_lb = z$ci.lb[[1]],
    MLM_CI_ub = z$ci.ub[[1]],
    MLM_Mod_Est = z$b[[2]],
    MLM_Mod_Est_SE = z$se[[2]],
    MLM_Mod_p = z$pval[[2]],
    MLM_Mod_CI_lb = z$ci.lb[[2]],
    MLM_Mod_CI_ub = z$ci.ub[[2]],
    MLM_tau2 = z$tau2,
    MLM_tau2_SE = z$se.tau2,
    MLM_BIC = BIC.rma(z),
    MLM_AIC = AIC.rma(z),
    MLM_LogLik = logLik.rma(z),
    
    Net_TE = zy$TE.random[[3]],
    Net_TE_SE = zy$seTE.random[[3]],
    Net_TE_CI_lb = zy$lower.random[[3]],
    Net_TE_CI_ub = zy$upper.random[[3]],
    Net_TE_pval = zy$pval.random[[3]],
    Net_TE_tau2 = zy$tau^2,
    Net_TE_I2 = zy$I2
    )
  return(Output)
}

system.time(
test1 <- ddply(MD_Results[Rep_Number == 1], .(M_D, Het_New, Num_Up, Rep_Number), MD_Analysis, 
                 .progress = "text")
)

a <- MD_Results[Rep_Number == 1]
names(a)

### Set working directory to results
setwd("D:\\Biostats and Epidemiology\\Project\\Project-Results")

### Write whole simulation to .csv
write.csv(RE_Results_MD, file = "RE_Results_MD.csv")

### Checking

sum(is.na(RE_Results_MD)==TRUE)
sum(RE_Results_MD$Up_tau2_se == 0)

which(is.na(RE_Results_MD) == TRUE)