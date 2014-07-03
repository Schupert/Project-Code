### Set working directory
setwd("D:\\Biostats and Epidemiology\\Project\\Project-SimCSV")

### Set seed
set.seed(1)

### Import required libraries
require(metafor)
require(data.table)
require(netmeta)

###### Import sample simulation, to determine length
MD_Results_frame <- read.csv("MD_10_3_5.csv")
MD_Results <- data.table(MD_Results_frame)



x <- summary(escalc(measure="MD", m1i = Group1Mean, sd1i = Group1SD, n1i = Group1Size, 
                    m2i = Group2Mean, sd2i = Group2SD, n2i = Group2Size, 
                    data=MD_Results[Rep_Number == 1], append = FALSE, replace = TRUE))
x2 <- rbind(x,x[1:20,])

a <- c(rep("upd", 20 + 5), rep("init", 20))
b <- rep("plac", length(a))

y1 <- netmeta(x2$yi, x2$sei, a, b, sm = "MD")



OR_Results_Frame <- read.csv("OR_10_3_5.csv")
OR_Results <- data.table(OR_Results_Frame)

x <- summary(escalc(measure="OR", ai = Group1Outcome1, bi = Group1Outcome2, n1i = Group1Size, 
         ci = Group2Outcome1, di = Group2Outcome2, n2i = Group2Size, 
         data=MD_Results[Rep_Number == 1], append = FALSE, replace = TRUE))

x2 <- rbind(x,x[1:20,])

y2 <- netmeta(x2$yi, x2$sei, a, b, sm = "OR")
