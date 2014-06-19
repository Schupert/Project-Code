### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-SimCSV")

### Set seed
set.seed(1)

### Import required libraries
require(metafor)

### Import sample simulation results, to determine length
OR_Results <- read.csv("OR_10_1_2.csv")

### Set variables
O_R_New <- c(1, 1.2, 1.4, 1.6, 1.8, 2)
Hetero_New <- c(0.1, 0.3, 0.5)
Updated_Studies <- c(2, 5, 10, 15, 20)
Repeats <- max(OR_Results$Rep_Number)
Database_Length <- length(O_R_New)*length(Hetero_New)*length(Updated_Studies)*Repeats

### Working with meta-analysis

RE_Results_OR <- data.frame(Rep_Number = integer(length = Database_Length),
                            OR_New = numeric(length = Database_Length), 
                            Het_New = numeric(length = Database_Length), 
                            Num_Up = integer(length = Database_Length),
                            Init_Est = numeric(length = Database_Length),
                            Init_Est_SE = numeric(length = Database_Length),
                            Init_Est_p = numeric(length = Database_Length),
                            Init_CI_lb = numeric(length = Database_Length),
                            Init_CI_ub = numeric(length = Database_Length),
                            Init_tau2 = numeric(length = Database_Length),
                            Init_tau2_SE = numeric(length = Database_Length),
                            Init_BIC = numeric(length = Database_Length),
                            Init_AIC = numeric(length = Database_Length),
                            Init_LogLik = numeric(length = Database_Length),
                            Up_Est = numeric(length = Database_Length),
                            Up_Est_SE = numeric(length = Database_Length),
                            Up_Est_p = numeric(length = Database_Length),
                            Up_CI_lb = numeric(length = Database_Length),
                            Up_CI_ub = numeric(length = Database_Length),
                            Up_tau2 = numeric(length = Database_Length),
                            Up_tau2_SE = numeric(length = Database_Length),
                            Up_BIC = numeric(length = Database_Length),
                            Up_AIC = numeric(length = Database_Length),
                            Up_LogLik = numeric(length = Database_Length),
                            MLM_Est = numeric(length = Database_Length),
                            MLM_Est_SE = numeric(length = Database_Length),
                            MLM_Est_p = numeric(length = Database_Length),
                            MLM_CI_lb = numeric(length = Database_Length),
                            MLM_CI_ub = numeric(length = Database_Length),
                            MLM_Mod_Est = numeric(length = Database_Length),
                            MLM_Mod_Est_SE = numeric(length = Database_Length),
                            MLM_Mod_p = numeric(length = Database_Length),
                            MLM_Mod_CI_lb = numeric(length = Database_Length),
                            MLM_Mod_CI_ub = numeric(length = Database_Length),
                            MLM_tau2 = numeric(length = Database_Length),
                            MLM_tau2_SE = numeric(length = Database_Length),
                            MLM_BIC = numeric(length = Database_Length),
                            MLM_AIC = numeric(length = Database_Length),
                            MLM_LogLik = numeric(length = Database_Length),
                            stringsAsFactors=FALSE)


## Set loops
counter <- 1

for (m in O_R_New){
  
  for (l in Hetero_New){
    
    for (k in Updated_Studies){
      
      nam <- paste("OR", as.integer(m*10), as.integer(l*10), k, sep = "_")
      imported_data <- read.csv(file = paste(nam, ".csv", sep = ""))
      
      ## Dummy variable for MLM
      
      Init_Updated <- c(rep(0, 20), rep(1, k))
      
      for(n in 1:Repeats){
        
        ## Set temporary dataset
        data_temp <- imported_data[imported_data$Rep_Number == n,]
        
        ## Initial meta-analysis (random-effects)
        x <- rma(measure="OR", ai = Group1Outcome1, bi = Group1Outcome2, n1i = Group1Size, 
                 ci = Group2Outcome1, di = Group2Outcome2, n2i = Group2Size, 
                 data=data_temp[data_temp$I_U == "I",], method="SJ")
        
        ## Input values
        
        RE_Results_OR$Init_Est[counter] <- x$b
        RE_Results_OR$Init_Est_SE[counter] <- x$se
        RE_Results_OR$Init_Est_p[counter] <- x$pval
        RE_Results_OR$Init_CI_lb[counter] <- x$ci.lb
        RE_Results_OR$Init_CI_ub[counter] <- x$ci.ub
        RE_Results_OR$Init_tau2[counter] <- x$tau2
        RE_Results_OR$Init_tau2_SE[counter] <- x$se.tau2
        RE_Results_OR$Init_BIC[counter] <- BIC.rma(x)
        RE_Results_OR$Init_AIC[counter] <- AIC.rma(x)
        RE_Results_OR$Init_LogLik[counter] <- logLik.rma(x)
        
        ## Updated meta-analysis
        y <- rma(measure="OR", ai = Group1Outcome1, bi = Group1Outcome2, n1i = Group1Size, 
                 ci = Group2Outcome1, di = Group2Outcome2, n2i = Group2Size, 
                 data=data_temp, method="SJ")
        
        ## Input values
        
        RE_Results_OR$Up_Est[counter] <- y$b
        RE_Results_OR$Up_Est_SE[counter] <- y$se
        RE_Results_OR$Up_Est_p[counter] <- y$pval
        RE_Results_OR$Up_CI_lb[counter] <- y$ci.lb
        RE_Results_OR$Up_CI_ub[counter] <- y$ci.ub
        RE_Results_OR$Up_tau2[counter] <- y$tau2
        RE_Results_OR$Up_tau2_SE[counter] <- y$se.tau2
        RE_Results_OR$Up_BIC[counter] <- BIC.rma(y)
        RE_Results_OR$Up_AIC[counter] <- AIC.rma(y)
        RE_Results_OR$Up_LogLik[counter] <- logLik.rma(y)
        
        ## Multi-level meta-analysis
        
        z <- rma(measure="OR", ai = Group1Outcome1, bi = Group1Outcome2, n1i = Group1Size, 
                 ci = Group2Outcome1, di = Group2Outcome2, n2i = Group2Size, 
                 data=data_temp, method="SJ", mods = Init_Updated)
        
        ## Input values
        
        RE_Results_OR$MLM_Est[counter] <- z$b[[1]]
        RE_Results_OR$MLM_Est_SE[counter] <- z$se[[1]]
        RE_Results_OR$MLM_Est_p[counter] <- z$pval[[1]]
        RE_Results_OR$MLM_CI_lb[counter] <- z$ci.lb[[1]]
        RE_Results_OR$MLM_CI_ub[counter] <- z$ci.ub[[1]]
        RE_Results_OR$MLM_Mod_Est[counter] <- z$b[[2]]
        RE_Results_OR$MLM_Mod_Est_SE[counter] <-z$se[[2]]
        RE_Results_OR$MLM_Mod_p[counter] <- z$pval[[2]]
        RE_Results_OR$MLM_Mod_CI_lb[counter] <- z$ci.lb[[2]]
        RE_Results_OR$MLM_Mod_CI_ub[counter] <- z$ci.ub[[2]]
        RE_Results_OR$MLM_tau2[counter] <- z$tau2
        RE_Results_OR$MLM_tau2_SE[counter] <- z$se.tau2
        RE_Results_OR$MLM_BIC[counter] <- BIC.rma(z)
        RE_Results_OR$MLM_AIC[counter] <- AIC.rma(z)
        RE_Results_OR$MLM_LogLik[counter] <- logLik.rma(z)
        
        ## Input static values
        RE_Results_OR$Rep_Number[counter] <- n
        RE_Results_OR$OR_New[counter] <- m
        RE_Results_OR$Het_New[counter] <- l
        RE_Results_OR$Num_Up[counter] <- k
        
        
        counter <- counter + 1
        
      }
    }
  }
}


### Set working directory to results
setwd("D:\\Biostats and Epidemiology\\Project\\Project-Results")

### Write whole simulation to .csv
write.csv(RE_Results_OR, file = "RE_Results_OR.csv")

sum(is.na(RE_Results_OR)==TRUE)
sum(RE_Results_OR$Up_tau2_se == 0)
head(RE_Results_OR)