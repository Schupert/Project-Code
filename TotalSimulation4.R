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

### Simulating mean difference function

Mean_Difference <- function(StudySize, MeanDifference, Heterogeneity, BioVariation){
  StudyMeanDifference <- rnorm(1, MeanDifference, sqrt(Heterogeneity))
  GroupSize <- rbinom(1, StudySize, 0.5)
  if (GroupSize <= 1){GroupSize <- 2}
  if ((StudySize - GroupSize) <= 1){GroupSize <- GroupSize - 2}
  Group1 <- rnorm(GroupSize, 10, BioVariation)
  Group2 <- rnorm((StudySize - GroupSize), (10 + StudyMeanDifference), BioVariation)
  Group1Mean <- mean(Group1)
  Group2Mean <- mean(Group2)
  Group1SD <- sd(Group1)
  Group2SD <- sd(Group2)
  return(c(Group1Mean, Group2Mean, Group1SD, Group2SD, GroupSize, (StudySize-GroupSize)))
}

### Set up Variables which change across simulations

Mean_Difference_Change <- seq(0, 1, 0.2)
Hetero_New <- c(0.1, 0.3, 0.5)
Hetero_Init <- 0.3
Updated_Studies <- c(2, 5, 10, 15, 20)
Initial_Studies <- 20
Variance <- 0.5
Repeats <- 5000

### Run Simulation

Total_Simulation_MD <- data.frame(Rep_Number = integer(),
                                  Study_ID = integer(), 
                                  I_U = factor(), 
                                  Group1Mean = numeric(), 
                                  Group2Mean = numeric(), 
                                  Group1SD = numeric(), 
                                  Group2SD = numeric(), 
                                  Group1Size = integer(), 
                                  Group2Size = integer(), 
                                  M_D = numeric(), 
                                  Het_New = numeric(), 
                                  Num_Up = integer(), 
                                  stringsAsFactors=FALSE)


for (m in Mean_Difference_Change){
  
  for (l in Hetero_New){
    
    for (k in Updated_Studies){
      
      ### Set vectors for correct length, taking account of changin numbers of studies
      
      Rep_Number <- vector(length = (Initial_Studies + k) * Repeats)
      Study_ID <- vector(length = (Initial_Studies + k) * Repeats)
      I_U <- vector(length = (Initial_Studies + k) * Repeats)
      Group1Mean <- vector(length = (Initial_Studies + k) * Repeats)
      Group2Mean <- vector(length = (Initial_Studies + k) * Repeats)
      Group1SD <- vector(length = (Initial_Studies + k) * Repeats)
      Group2SD <- vector(length = (Initial_Studies + k) * Repeats)
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
          x <- Mean_Difference(Sample_Sizes_All[Study_ID_Counter], 0, Hetero_Init, Variance)
          Rep_Number[Total_Counter] <- Repeat_Counter
          Study_ID[Total_Counter] <- Study_ID_Counter
          I_U[Total_Counter] <- "I"
          Group1Mean[Total_Counter] <- x[1]
          Group2Mean[Total_Counter] <- x[2]
          Group1SD[Total_Counter] <- x[3]
          Group2SD[Total_Counter] <- x[4]
          Group1Size[Total_Counter] <- x[5]
          Group2Size[Total_Counter] <- x[6]
          
          Total_Counter <- Total_Counter + 1
          Study_ID_Counter <- Study_ID_Counter + 1
        } 
        
        ## Run updated studies
        for (i in (Initial_Studies + 1):(Initial_Studies + k)){ 
          x <- Mean_Difference(Sample_Sizes_All[Study_ID_Counter], m, l, Variance)
          Rep_Number[Total_Counter] <- Repeat_Counter
          Study_ID[Total_Counter] <- Study_ID_Counter
          I_U[Total_Counter] <- "U"
          Group1Mean[Total_Counter] <- x[1]
          Group2Mean[Total_Counter] <- x[2]
          Group1SD[Total_Counter] <- x[3]
          Group2SD[Total_Counter] <- x[4]
          Group1Size[Total_Counter] <- x[5]
          Group2Size[Total_Counter] <- x[6]
          
          Total_Counter <- Total_Counter + 1
          Study_ID_Counter <- Study_ID_Counter + 1
        } 
        Repeat_Counter <- Repeat_Counter + 1
      }
      
      M_D <- rep(m, (Initial_Studies + k) * Repeats)
      Het_New <- rep(l, (Initial_Studies + k) * Repeats)
      Num_Up <- rep(k, (Initial_Studies + k) * Repeats) 
      
      ### Assign set of repeats a name, then write to csv
      
      nam <- paste("MD", as.integer(m*10), as.integer(l*10), k, sep = "_")
      assign(nam, data.frame(Rep_Number, Study_ID, I_U, Group1Mean, Group2Mean, Group1SD, Group2SD, Group1Size, Group2Size))
      write.csv(get(nam), file = paste(nam, ".csv", sep = ""))
      Temp_DF <- data.frame(Rep_Number, Study_ID, I_U, Group1Mean, Group2Mean, Group1SD, Group2SD, Group1Size, Group2Size, M_D, Het_New, Num_Up)
      Total_Simulation_MD <- rbind (Temp_DF, Total_Simulation_MD)
      
    }
  }
}


### Write whole simulation to .csv
write.csv(Total_Simulation_MD, file = "Total_Sim_MD.csv")

### Check no negative values

for (i in Total_Simulation_MD){ print(sum( i< 0)) }
sum(is.na(Total_Simulation_MD)==TRUE)

### Function to simulate a single observational study. Control_prop represents the proportion of 
### outcome 1 (deaths) that are controls. Pic is the proportion of controls that are expected to die

Log_Odds_Ratio <- function(StudySize, O_R, Heterogeneity, Control_Prop){
  StudyLogOR <- rnorm(1, log(O_R), sqrt(Heterogeneity))
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
Control_Prop <- 0.5

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
          x <- Log_Odds_Ratio(Sample_Sizes_All[Study_ID_Counter], 1, Hetero_Init, Control_Prop)
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

### Function to simulate a single observational study. Control_prop represents the proportion of deaths 
### (outcome 1) in group 1 (controls). Proportions of the study in each group are also set.
### Need rewrite to make slope heterogeneous, not ORs on the slope

Log_Odds_Ratio_Four <- function(StudySize, log_OR_Slope, Heterogeneity, 
                                Control_Prop, Group1_Prop, Group2_Prop, Group3_Prop){
  
  Study_logOR_Slope <- rnorm(1, log_OR_Slope, sqrt(Heterogeneity))
  if (Study_logOR_Slope < -0.32) Study_logOR_Slope <- runif(1, -0.32, -0.25)
  StudyOR1 <- exp(Study_logOR_Slope)
  StudyOR2 <- exp(Study_logOR_Slope*2)
  StudyOR3 <- exp(Study_logOR_Slope*3)
  
  Group1Size <- as.integer(StudySize*Group1_Prop)
  Group2Size <- as.integer(StudySize*Group2_Prop)
  Group3Size <- as.integer(StudySize*Group3_Prop)
  Group4Size <- StudySize - Group1Size - Group2Size - Group3Size
  
  Pic <- runif(1, 0.05, 0.65)
  Group1Out1 <- as.integer(rbinom(1, Group1Size, Pic))
  Group1Out2 <- as.integer(Group1Size - Group1Out1)
  
  Pit1 <- (Pic * StudyOR1)/(1 - Pic + (Pic * StudyOR1))
  Group2Out1 <- as.integer(rbinom(1, Group2Size, Pit1))
  Group2Out2 <- as.integer(Group2Size - Group2Out1)
  
  Pit2 <- (Pic * StudyOR2)/(1 - Pic + (Pic * StudyOR2))
  Group3Out1 <- as.integer(rbinom(1, Group3Size, Pit2))
  Group3Out2 <- as.integer(Group3Size - Group3Out1)
  
  Pit3 <- (Pic * StudyOR3)/(1 - Pic + (Pic * StudyOR3))
  Group4Out1 <- as.integer(rbinom(1, Group4Size, Pit3))
  Group4Out2 <- as.integer(Group4Size - Group4Out1)
  
  if (Group1Out1 == 0){Group1Out1 <- 0.5}
  if (Group2Out1 == 0){Group2Out1 <- 0.5}
  if (Group1Out2 == 0){Group1Out2 <- 0.5}
  if (Group2Out2 == 0){Group2Out2 <- 0.5}
  if (Group3Out1 == 0){Group1Out1 <- 0.5}
  if (Group4Out1 == 0){Group2Out1 <- 0.5}
  if (Group3Out2 == 0){Group1Out2 <- 0.5}
  if (Group4Out2 == 0){Group2Out2 <- 0.5}
  
  return(c(Group1Out1, Group1Out2, Group2Out1, Group2Out2, Group3Out1, Group3Out2, Group4Out1, Group4Out2, 
           Group1Size, Group2Size, Group3Size, Group4Size))
}


### Set simulation values

O_R_New <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

### Run Simulation

Total_Simulation_OR_Four <- data.frame(Rep_Number = integer(),
                                       Study_ID = integer(), 
                                       I_U = factor(), 
                                       Group1Out1 = integer(), 
                                       Group1Out2 = integer(), 
                                       Group2Out1 = integer(), 
                                       Group2Out2 = integer(), 
                                       Group3Out1 = integer(), 
                                       Group3Out2 = integer(),
                                       Group4Out1 = integer(), 
                                       Group4Out2 = integer(),
                                       Group1Size = integer(), 
                                       Group2Size = integer(), 
                                       Group3Size = integer(), 
                                       Group4Size = integer(), 
                                       Study_O_R = list(), 
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
      Group3Outcome1 <- vector(length = (Initial_Studies + k) * Repeats)
      Group3Outcome2 <- vector(length = (Initial_Studies + k) * Repeats)
      Group4Outcome1 <- vector(length = (Initial_Studies + k) * Repeats)
      Group4Outcome2 <- vector(length = (Initial_Studies + k) * Repeats)
      Group1Size <- vector(length = (Initial_Studies + k) * Repeats)
      Group2Size <- vector(length = (Initial_Studies + k) * Repeats)
      Group3Size <- vector(length = (Initial_Studies + k) * Repeats)
      Group4Size <- vector(length = (Initial_Studies + k) * Repeats)
      
      
      ## Set counters
      Study_ID_Counter <- 1
      Repeat_Counter <- 1
      Total_Counter <- 1
      
      
      for (j in 1:Repeats){
        
        Study_ID_Counter <- 1
        
        ## Run initial studies
        for (i in 1:(Initial_Studies)){ 
          x <- Log_Odds_Ratio_Four(Sample_Sizes_All[Study_ID_Counter], 0, Hetero_Init, Control_Prop, 0.4, 0.2, 0.2)
          Rep_Number[Total_Counter] <- Repeat_Counter
          Study_ID[Total_Counter] <- Study_ID_Counter
          I_U[Total_Counter] <- "I"
          Group1Outcome1[Total_Counter] <- x[1]
          Group1Outcome2[Total_Counter] <- x[2]
          Group2Outcome1[Total_Counter] <- x[3]
          Group2Outcome2[Total_Counter] <- x[4]
          Group3Outcome1[Total_Counter] <- x[5]
          Group3Outcome2[Total_Counter] <- x[6]
          Group4Outcome1[Total_Counter] <- x[7]
          Group4Outcome2[Total_Counter] <- x[8]
          Group1Size[Total_Counter] <- x[9]
          Group2Size[Total_Counter] <- x[10]
          Group3Size[Total_Counter] <- x[11]
          Group4Size[Total_Counter] <- x[12]
          
          Total_Counter <- Total_Counter + 1
          Study_ID_Counter <- Study_ID_Counter + 1
        } 
        
        ## Run updated studies
        for (i in (Initial_Studies + 1):(Initial_Studies + k)){ 
          x <- Log_Odds_Ratio_Four(Sample_Sizes_All[Study_ID_Counter], m, l, Control_Prop, 0.4, 0.2, 0.2)
          Rep_Number[Total_Counter] <- Repeat_Counter
          Study_ID[Total_Counter] <- Study_ID_Counter
          I_U[Total_Counter] <- "U"
          Group1Outcome1[Total_Counter] <- x[1]
          Group1Outcome2[Total_Counter] <- x[2]
          Group2Outcome1[Total_Counter] <- x[3]
          Group2Outcome2[Total_Counter] <- x[4]
          Group3Outcome1[Total_Counter] <- x[5]
          Group3Outcome2[Total_Counter] <- x[6]
          Group4Outcome1[Total_Counter] <- x[7]
          Group4Outcome2[Total_Counter] <- x[8]
          Group1Size[Total_Counter] <- x[9]
          Group2Size[Total_Counter] <- x[10]
          Group3Size[Total_Counter] <- x[11]
          Group4Size[Total_Counter] <- x[12]
          
          Total_Counter <- Total_Counter + 1
          Study_ID_Counter <- Study_ID_Counter + 1
        } 
        Repeat_Counter <- Repeat_Counter + 1
      }
      
      Rep_O_R <- rep(m, (Initial_Studies + k) * Repeats)
      Het_New <- rep(l, (Initial_Studies + k) * Repeats)
      Num_Up <- rep(k, (Initial_Studies + k) * Repeats) 
      
      ### Assign set of repeats a name, then write to csv
      
      nam <- paste("ORFour", as.integer(m*10), as.integer(l*10), k, sep = "_")
      assign(nam, data.frame(Rep_Number, Study_ID, I_U, Group1Outcome1, Group1Outcome2, Group2Outcome1, 
                             Group2Outcome2, Group3Outcome1, Group3Outcome2, Group4Outcome1, Group4Outcome2,
                             Group1Size, Group2Size, Group3Size, Group4Size))
      write.csv(get(nam), file = paste(nam, ".csv", sep = ""))
      Temp_DF <- data.frame(Rep_Number, Study_ID, I_U, Group1Outcome1, Group1Outcome2, Group2Outcome1, 
                            Group2Outcome2, Group3Outcome1, Group3Outcome2, Group4Outcome1, Group4Outcome2,
                            Group1Size, Group2Size, Group3Size, Group4Size)
      Total_Simulation_OR_Four <- rbind (Temp_DF, Total_Simulation_OR_Four)
      
    }
  }
}

### Write whole simulation to .csv
write.csv(Total_Simulation_OR_Four, file = "Total_Sim_OR_Four.csv")

### Check no negative values

for (i in Total_Simulation_OR_Four){ print(sum( i< 0)) }