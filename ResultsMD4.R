### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-SimCSV")

### Set seed
set.seed(1)

### Import required libraries
require(metafor)
require(data.table)
require(netmeta)

### Import sample simulation, to determine length
MD_Results_frame <- read.csv("Total_Sim_MD.csv")
MD_Results <- data.table(MD_Results_frame)
setkey(MD_Results, "M_D", "Het_New", "Num_Up", "Rep_Number", "I_U")

### Set variables
Mean_Difference_Change <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
Hetero_New <- c(0.1, 0.3, 0.5)
Updated_Studies <- c(2, 5, 10, 15, 20)
Repeats <- max(MD_Results$Rep_Number)
Database_Length <- length(Mean_Difference_Change)*length(Hetero_New)*length(Updated_Studies)*Repeats

### Working with meta-analysis

RE_Results_MD <- data.table(Rep_Number = integer(length = Database_Length),
                            M_D = numeric(length = Database_Length), 
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
                            Net_TE = numeric(length = Database_Length),
                            Net_TE_SE = numeric(length = Database_Length),
                            Net_TE_CI_lb = numeric(length = Database_Length),
                            Net_TE_CI_ub = numeric(length = Database_Length),
                            Net_TE_pval = numeric(length = Database_Length),
                            Net_TE_tau2 = numeric(length = Database_Length),
                            Net_TE_I2 = numeric(length = Database_Length)
                            )


## Set loops

counter <- 1

for (m in Mean_Difference_Change){
  
  for (l in Hetero_New){
    
    for (k in Updated_Studies){
      
      ## Dummy variable for MLM
      
      Init_Updated <- c(rep(0, 20), rep(1, k))
      
      for(n in 1:Repeats){
        
        ## Set temporary dataset
        data_temp_init <- MD_Results[J(m, l, k, n, "I")]
        data_temp_all <- MD_Results[J(m, l, k, n)]
        
        ## Initial meta-analysis (random-effects)
        x <- rma(measure="MD", m1i = Group1Mean, sd1i = Group1SD, n1i = Group1Size, 
                 m2i = Group2Mean, sd2i = Group2SD, n2i = Group2Size, 
                 data=data_temp_init, method="SJ")
        
        ## Input values
        
        RE_Results_MD[counter, Init_Est:= x$b]
        RE_Results_MD[counter, Init_Est_SE:= x$se]
        RE_Results_MD[counter, Init_Est_p:= x$pval] 
        RE_Results_MD[counter, Init_CI_lb := x$ci.lb]
        RE_Results_MD[counter, Init_CI_ub := x$ci.ub] 
        RE_Results_MD[counter, Init_tau2 := x$tau2] 
        RE_Results_MD[counter, Init_tau2_SE := x$se.tau2]
        RE_Results_MD[counter, Init_BIC := BIC.rma(x)] 
        RE_Results_MD[counter, Init_AIC := AIC.rma(x)]
        RE_Results_MD[counter, Init_LogLik := logLik.rma(x)]
        
        ## Updated meta-analysis
        y <- rma(measure="MD", m1i = Group1Mean, sd1i = Group1SD, n1i = Group1Size, 
                 m2i = Group2Mean, sd2i = Group2SD, n2i = Group2Size, 
                 data=data_temp_all, method="SJ")
        
        ## Input values
        
        RE_Results_MD[counter, Up_Est := y$b] 
        RE_Results_MD[counter, Up_Est_SE := y$se] 
        RE_Results_MD[counter, Up_Est_p := y$pval]
        RE_Results_MD[counter, Up_CI_lb := y$ci.lb]
        RE_Results_MD[counter, Up_CI_ub := y$ci.ub] 
        RE_Results_MD[counter, Up_tau2 := y$tau2]
        RE_Results_MD[counter, Up_tau2_SE := y$se.tau2] 
        RE_Results_MD[counter, Up_BIC := BIC.rma(y)] 
        RE_Results_MD[counter, Up_AIC := AIC.rma(y)] 
        RE_Results_MD[counter, Up_LogLik := logLik.rma(y)]
        
        ## Multi-level meta-analysis
        
        z <- rma(measure="MD", m1i = Group1Mean, sd1i = Group1SD, n1i = Group1Size, 
                 m2i = Group2Mean, sd2i = Group2SD, n2i = Group2Size, 
                 data=data_temp_all, method="SJ", mods = Init_Updated)
        
        ## Input values
        
        RE_Results_MD[counter, MLM_Est := z$b[[1]] ] 
        RE_Results_MD[counter, MLM_Est_SE := z$se[[1]] ]
        RE_Results_MD[counter, MLM_Est_p := z$pval[[1]]]  
        RE_Results_MD[counter, MLM_CI_lb := z$ci.lb[[1]]]
        RE_Results_MD[counter, MLM_CI_ub := z$ci.ub[[1]]] 
        RE_Results_MD[counter, MLM_Mod_Est := z$b[[2]]] 
        RE_Results_MD[counter, MLM_Mod_Est_SE := z$se[[2]]]
        RE_Results_MD[counter, MLM_Mod_p := z$pval[[2]]]  
        RE_Results_MD[counter, MLM_Mod_CI_lb := z$ci.lb[[2]]] 
        RE_Results_MD[counter, MLM_Mod_CI_ub := z$ci.ub[[2]]] 
        RE_Results_MD[counter, MLM_tau2 := z$tau2] 
        RE_Results_MD[counter, MLM_tau2_SE := z$se.tau2] 
        RE_Results_MD[counter, MLM_BIC := BIC.rma(z)] 
        RE_Results_MD[counter, MLM_AIC := AIC.rma(z)]
        RE_Results_MD[counter, MLM_LogLik := logLik.rma(z)] 
        
        ## Network Meta Analysis
        
        Net_temp1 <- summary(escalc(measure="MD", m1i = Group1Mean, sd1i = Group1SD, n1i = Group1Size, 
                            m2i = Group2Mean, sd2i = Group2SD, n2i = Group2Size, 
                            data= data_temp_all, append = FALSE, replace = TRUE))
        
        Net_temp2 <- rbindlist(list(Net_temp1, Net_temp1[1:20,]))
        
        Net_treatlist1 <- c(rep("upd", 20 + k), rep("init", 20))
        
        Net_treatlist2 <- rep("plac", length(Net_treatlist1))
        
        zy <- netmeta(Net_temp2$yi, Net_temp2$sei, Net_treatlist2, Net_treatlist1, sm = "MD", warn = FALSE)
        
        ## Input values
        
        RE_Results_MD[counter, Net_TE := zy$TE.random[[3]] ]
        RE_Results_MD[counter, Net_TE_SE := zy$seTE.random[[3]] ]
        RE_Results_MD[counter, Net_TE_CI_lb := zy$lower.random[[3]] ]
        RE_Results_MD[counter, Net_TE_CI_ub := zy$upper.random[[3]] ]
        RE_Results_MD[counter, Net_TE_pval := zy$pval.random[[3]] ]
        RE_Results_MD[counter, Net_TE_tau2 := zy$tau^2 ]
        RE_Results_MD[counter, Net_TE_I2 := zy$I2 ]
        
        ## Input static values
        RE_Results_MD[counter, Rep_Number := n] 
        RE_Results_MD[counter, M_D := m]
        RE_Results_MD[counter, Het_New := l]
        RE_Results_MD[counter, Num_Up := k]
        
        counter <- counter + 1
        
      }
    }
  }
}


### Set working directory to results
setwd("D:\\Biostats and Epidemiology\\Project\\Project-Results")

### Write whole simulation to .csv
write.csv(RE_Results_MD, file = "RE_Results_MD.csv")

### Checking

sum(is.na(RE_Results_MD)==TRUE)
sum(RE_Results_MD$Up_tau2_se == 0)

which(is.na(RE_Results_MD) == TRUE)