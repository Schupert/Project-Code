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
  StudyMeanDifference <- rnorm(1, MeanDifference, Heterogeneity)
  GroupSize <- rbinom(1, StudySize, 0.5)
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
Updated_Studies <- c(2, 5, 10, 15, 20)
Initial_Studies <- 20
Variance <- 0.5
Repeats <- 1000

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
  x <- Mean_Difference(Sample_Sizes_All[Study_ID_Counter], 0, 0.3, Variance)
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
which(is.na(Total_Simulation_MD)==TRUE)