### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-SimCSV")

### Set seed
set.seed(1)

### Import required libraries
require(metafor)
require(data.table)

### Import sample simulation, to determine length
OR_Results_frame <- read.csv("Total_Sim_OR.csv")
OR_Results <- data.table(OR_Results_frame)
setkey(OR_Results, "Rep_O_R", "Het_New", "Num_Up", "Rep_Number", "I_U")

### Set variables
O_R_New <- c(1, 1.2, 1.4, 1.6, 1.8, 2)
Hetero_New <- c(0.1, 0.3, 0.5)
Updated_Studies <- c(2, 5, 10, 15, 20)
Repeats <- max(OR_Results$Rep_Number)
Database_Length <- length(O_R_New)*length(Hetero_New)*length(Updated_Studies)*Repeats

### Working with meta-analysis

RE_Results_OR <- data.table(Rep_Number = integer(length = Database_Length),
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
      
      ## Dummy variable for MLM
      
      Init_Updated <- c(rep(0, 20), rep(1, k))
      
      for(n in 1:Repeats){
        
        ## Set temporary dataset
        data_temp_init <- OR_Results[J(m, l, k, n, "I")]
        data_temp_all <- OR_Results[J(m, l, k, n)]
        
        ## Initial meta-analysis (random-effects)
        x <- rma(measure="OR", ai = Group1Outcome1, bi = Group1Outcome2, n1i = Group1Size, 
                 ci = Group2Outcome1, di = Group2Outcome2, n2i = Group2Size, 
                 data=data_temp_init, method="SJ")
        
        ## Input values
        
        RE_Results_OR[counter, Init_Est:= x$b]
        RE_Results_OR[counter, Init_Est_SE:= x$se]
        RE_Results_OR[counter, Init_Est_p:= x$pval] 
        RE_Results_OR[counter, Init_CI_lb := x$ci.lb]
        RE_Results_OR[counter, Init_CI_ub := x$ci.ub] 
        RE_Results_OR[counter, Init_tau2 := x$tau2] 
        RE_Results_OR[counter, Init_tau2_SE := x$se.tau2]
        RE_Results_OR[counter, Init_BIC := BIC.rma(x)] 
        RE_Results_OR[counter, Init_AIC := AIC.rma(x)]
        RE_Results_OR[counter, Init_LogLik := logLik.rma(x)]
        
        ## Updated meta-analysis
        y <- rma(measure="OR", ai = Group1Outcome1, bi = Group1Outcome2, n1i = Group1Size, 
                 ci = Group2Outcome1, di = Group2Outcome2, n2i = Group2Size, 
                 data=data_temp_all, method="SJ")
        
        ## Input values
        
        RE_Results_OR[counter, Up_Est := y$b] 
        RE_Results_OR[counter, Up_Est_SE := y$se] 
        RE_Results_OR[counter, Up_Est_p := y$pval]
        RE_Results_OR[counter, Up_CI_lb := y$ci.lb]
        RE_Results_OR[counter, Up_CI_ub := y$ci.ub] 
        RE_Results_OR[counter, Up_tau2 := y$tau2]
        RE_Results_OR[counter, Up_tau2_SE := y$se.tau2] 
        RE_Results_OR[counter, Up_BIC := BIC.rma(y)] 
        RE_Results_OR[counter, Up_AIC := AIC.rma(y)] 
        RE_Results_OR[counter, Up_LogLik := logLik.rma(y)]
        
        ## Multi-level meta-analysis
        
        z <- rma(measure="OR", ai = Group1Outcome1, bi = Group1Outcome2, n1i = Group1Size, 
                 ci = Group2Outcome1, di = Group2Outcome2, n2i = Group2Size, 
                 data=data_temp_all, method="SJ", mods = Init_Updated)
        
        ## Input values
        
        RE_Results_OR[counter, MLM_Est := z$b[[1]] ] 
        RE_Results_OR[counter, MLM_Est_SE := z$se[[1]] ]
        RE_Results_OR[counter, MLM_Est_p := z$pval[[1]]]  
        RE_Results_OR[counter, MLM_CI_lb := z$ci.lb[[1]]]
        RE_Results_OR[counter, MLM_CI_ub := z$ci.ub[[1]]] 
        RE_Results_OR[counter, MLM_Mod_Est := z$b[[2]]] 
        RE_Results_OR[counter, MLM_Mod_Est_SE := z$se[[2]]]
        RE_Results_OR[counter, MLM_Mod_p := z$pval[[2]]]  
        RE_Results_OR[counter, MLM_Mod_CI_lb := z$ci.lb[[2]]] 
        RE_Results_OR[counter, MLM_Mod_CI_ub := z$ci.ub[[2]]] 
        RE_Results_OR[counter, MLM_tau2 := z$tau2] 
        RE_Results_OR[counter, MLM_tau2_SE := z$se.tau2] 
        RE_Results_OR[counter, MLM_BIC := BIC.rma(z)] 
        RE_Results_OR[counter, MLM_AIC := AIC.rma(z)]
        RE_Results_OR[counter, MLM_LogLik := logLik.rma(z)] 
        
        
        ## Input static values
        RE_Results_OR[counter, Rep_Number := n] 
        RE_Results_OR[counter, OR_New := m]
        RE_Results_OR[counter, Het_New := l]
        RE_Results_OR[counter, Num_Up := k]
        
        counter <- counter + 1
        
      }
    }
  }
}


### Set working directory to results
setwd("D:\\Biostats and Epidemiology\\Project\\Project-Results")

### Write whole simulation to .csv
write.csv(RE_Results_OR, file = "RE_Results_OR.csv")

###

sum(is.na(RE_Results_OR)==TRUE)
sum(RE_Results_OR$Up_tau2_se == 0)