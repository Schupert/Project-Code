### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-SimCSV")

### Set seed
set.seed(1)

### Import required libraries

### Determine Study Sizes

Small_Study_Sizes <- as.integer(runif(40, 20, 50))
Medium_Study_Sizes <- as.integer(runif(40, 51, 200))
Large_Study_Sizes <- as.integer(runif(40, 201, 1000))
Study_Size_Allocation_Dummy <- runif(40, 0, 1)
Sample_Sizes_All <- vector()

for (i in 1:40){
  if (Study_Size_Allocation_Dummy[i] <0.5){Sample_Sizes_All[i] <- Small_Study_Sizes[i]}
  else if (Study_Size_Allocation_Dummy[i] >0.9){Sample_Sizes_All[i] <- Large_Study_Sizes[i]}
  else {Sample_Sizes_All[i] <- Medium_Study_Sizes[i]}
}

### Function to simulate a single observational study. Control_prop represents the proportion of 
### outcome 1 (deaths) that are controls. Pic is the proportion of controls that are expected to die

Log_Odds_Ratio <- function(StudySize, O_R, Heterogeneity, Control_Prop){
  StudyLogOR <- rnorm(1, log(O_R), Heterogeneity)
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- as.integer(StudySize - Group1Size)
  Pic <- runif(1, 0.05, 0.65)
  Group1Out1 <- as.integer(rbinom(1, Group1Size, Pic))
  Group1Out2 <- as.integer(Group1Size - Group1Out1)
  Pit <- (Pic * exp(StudyLogOR))/(1 - Pic + (Pic * exp(StudyLogOR)))
  Group2Out1 <- as.integer(rbinom(1, Group2Size, Pit))
  Group2Out2 <- as.integer(Group2Size - Group2Out1)
  if (Group1Out1 == 0){Group1Out1 <- 0.5}
  if (Group2Out1 == 0){Group2Out1 <- 0.5}
  if (Group1Out2 == 0){Group1Out2 <- 0.5}
  if (Group2Out2 == 0){Group2Out2 <- 0.5}
  return(c(Group1Out1, Group1Out2, Group2Out1, Group2Out2, Group1Size, Group2Size))
}

### Set simulation values

O_R_New <- c(1, 1.2, 1.4, 1.6, 1.8, 2)
Hetero_New <- c(0.1, 0.3, 0.5)
Updated_Studies <- c(2, 5, 10, 15, 20)
Initial_Studies <- 20
Control_Prop <- 0.5
Repeats <- 10

### Run Simulation

Total_Simulation_OR <- data.frame(Rep_Number = integer(),
                                  Study_ID = integer(), 
                                  I_U = factor(), 
                                  Group1Out1 = integer(), 
                                  Group1Out2 = integer(), 
                                  Group2Out1 = integer(), 
                                  Group2Out2 = integer(), 
                                  Group1Size = integer(), 
                                  Group2Size = integer(), 
                                  Study_O_R = numeric(), 
                                  Het_New = numeric(), 
                                  Num_Up = integer(), 
                                  stringsAsFactors=FALSE)


for (m in O_R_New){
  
  for (l in Hetero_New){
    
    for (k in Updated_Studies){
      
      ### Set vectors for correct length, taking account of changin numbers of studies
      
      Rep_Number <- vector(length = (Initial_Studies + k) * Repeats)
      Study_ID <- vector(length = (Initial_Studies + k) * Repeats)
      I_U <- vector(length = (Initial_Studies + k) * Repeats)
      Group1Outcome1 <- vector(length = (Initial_Studies + k) * Repeats)
      Group1Outcome2 <- vector(length = (Initial_Studies + k) * Repeats)
      Group2Outcome1 <- vector(length = (Initial_Studies + k) * Repeats)
      Group2Outcome2 <- vector(length = (Initial_Studies + k) * Repeats)
      Group1Size <- vector(length = (Initial_Studies + k) * Repeats)
      Group2Size <- vector(length = (Initial_Studies + k) * Repeats)
      
      
      ## Set counters
      Study_ID_Counter <- 1
      Repeat_Counter <- 1
      Total_Counter <- 1
      
      
      for (j in 1:Repeats){
        
        Study_ID_Counter <- 1
        
        ## Run initial studies
        for (i in 1:(Initial_Studies)){ 
          x <- Log_Odds_Ratio(Sample_Sizes_All[Study_ID_Counter], 1, 0.3, Control_Prop)
          Rep_Number[Total_Counter] <- Repeat_Counter
          Study_ID[Total_Counter] <- Study_ID_Counter
          I_U[Total_Counter] <- "I"
          Group1Outcome1[Total_Counter] <- x[1]
          Group1Outcome2[Total_Counter] <- x[2]
          Group2Outcome1[Total_Counter] <- x[3]
          Group2Outcome2[Total_Counter] <- x[4]
          Group1Size[Total_Counter] <- x[5]
          Group2Size[Total_Counter] <- x[6]
          
          Total_Counter <- Total_Counter + 1
          Study_ID_Counter <- Study_ID_Counter + 1
        } 
        
        ## Run updated studies
        for (i in (Initial_Studies + 1):(Initial_Studies + k)){ 
          x <- Log_Odds_Ratio(Sample_Sizes_All[Study_ID_Counter], m, l, Control_Prop)
          Rep_Number[Total_Counter] <- Repeat_Counter
          Study_ID[Total_Counter] <- Study_ID_Counter
          I_U[Total_Counter] <- "U"
          Group1Outcome1[Total_Counter] <- x[1]
          Group1Outcome2[Total_Counter] <- x[2]
          Group2Outcome1[Total_Counter] <- x[3]
          Group2Outcome2[Total_Counter] <- x[4]
          Group1Size[Total_Counter] <- x[5]
          Group2Size[Total_Counter] <- x[6]
          
          Total_Counter <- Total_Counter + 1
          Study_ID_Counter <- Study_ID_Counter + 1
        } 
        Repeat_Counter <- Repeat_Counter + 1
      }
      
      Rep_O_R <- rep(m, (Initial_Studies + k) * Repeats)
      Het_New <- rep(l, (Initial_Studies + k) * Repeats)
      Num_Up <- rep(k, (Initial_Studies + k) * Repeats) 
      
      ### Assign set of repeats a name, then write to csv
      
      nam <- paste("OR", as.integer(m*10), as.integer(l*10), k, sep = "_")
      assign(nam, data.frame(Rep_Number, Study_ID, I_U, Group1Outcome1, Group1Outcome2, Group2Outcome1, Group2Outcome2, Group1Size, Group2Size))
      write.csv(get(nam), file = paste(nam, ".csv", sep = ""))
      Temp_DF <- data.frame(Rep_Number, Study_ID, I_U, Group1Outcome1, Group1Outcome2, Group2Outcome1, Group2Outcome2, Group1Size, Group2Size, Rep_O_R, Het_New, Num_Up)
      Total_Simulation_OR <- rbind (Temp_DF, Total_Simulation_OR)
      
    }
  }
}

### Write whole simulation to .csv
write.csv(Total_Simulation_OR, file = "Total_Sim_OR.csv")

### Check no negative values

for (i in Total_Simulation_OR){ print(sum( i< 0)) }