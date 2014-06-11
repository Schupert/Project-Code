### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-SimCSV")

### Set seed
set.seed(1)

### Import required libraries
require(metafor)

### Import sample simulation, to determine length
MD_Results <- read.csv("MD_10_1_2.csv")

### Set variables
Mean_Difference_Change <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
Hetero_New <- c(0.1, 0.3, 0.5)
Updated_Studies <- c(2, 5, 10, 15, 20)
Repeats <- max(MD_Results$Rep_Number)
Database_Length <- length(Mean_Difference_Change)*length(Hetero_New)*length(Updated_Studies)*Repeats

### Working with meta-analysis

RE_Results_MD <- data.frame(Rep_Number = integer(length = Database_Length),
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
                            stringsAsFactors=FALSE)


## Set loops

counter <- 1

for (m in Mean_Difference_Change){
  
  for (l in Hetero_New){
    
    for (k in Updated_Studies){
      
      nam <- paste("MD", as.integer(m*10), as.integer(l*10), k, sep = "_")
      imported_data <- read.csv(file = paste(nam, ".csv", sep = ""))
      
      for(n in 1:Repeats){
        
        ## Set temporary dataset
        data_temp <- imported_data[imported_data$Rep_Number == n,]
        
        ## Initial meta-analysis (random-effects)
        x <- rma(measure="MD", m1i = Group1Mean, sd1i = Group1SD, n1i = Group1Size, 
                 m2i = Group2Mean, sd2i = Group2SD, n2i = Group2Size, 
                 data=data_temp[data_temp$I_U == "I",], method="SJ")
        
        ## Input values
        
        RE_Results_MD$Init_Est[counter] <- x$b
        RE_Results_MD$Init_Est_SE[counter] <- x$se
        RE_Results_MD$Init_Est_p[counter] <- x$pval
        RE_Results_MD$Init_CI_lb[counter] <- x$ci.lb
        RE_Results_MD$Init_CI_ub[counter] <- x$ci.ub
        RE_Results_MD$Init_tau2[counter] <- x$tau2
        RE_Results_MD$Init_tau2_SE[counter] <- x$se.tau2
        RE_Results_MD$Init_BIC[counter] <- BIC.rma(x)
        RE_Results_MD$Init_AIC[counter] <- AIC.rma(x)
        RE_Results_MD$Init_LogLik[counter] <- logLik.rma(x)
        
        ## Updated meta-analysis
        y <- rma(measure="MD", m1i = Group1Mean, sd1i = Group1SD, n1i = Group1Size, 
                 m2i = Group2Mean, sd2i = Group2SD, n2i = Group2Size, 
                 data=data_temp, method="SJ")
        
        ## Input values
        
        RE_Results_MD$Up_Est[counter] <- y$b
        RE_Results_MD$Up_Est_SE[counter] <- y$se
        RE_Results_MD$Up_Est_p[counter] <- y$pval
        RE_Results_MD$Up_CI_lb[counter] <- y$ci.lb
        RE_Results_MD$Up_CI_ub[counter] <- y$ci.ub
        RE_Results_MD$Up_tau2[counter] <- y$tau2
        RE_Results_MD$Up_tau2_SE[counter] <- y$se.tau2
        RE_Results_MD$Up_BIC[counter] <- BIC.rma(y)
        RE_Results_MD$Up_AIC[counter] <- AIC.rma(y)
        RE_Results_MD$Up_LogLik[counter] <- logLik.rma(y)
        
        ## Input static values
        RE_Results_MD$Rep_Number[counter] <- n
        RE_Results_MD$M_D[counter] <- m
        RE_Results_MD$Het_New[counter] <- l
        RE_Results_MD$Num_Up[counter] <- k
        
        counter <- counter + 1
        
      }
    }
  }
}


### Set working directory to results
setwd("D:\\Biostats and Epidemiology\\Project\\Project-Results")

### Write whole simulation to .csv
write.csv(RE_Results_MD, file = "RE_Results_MD.csv")
