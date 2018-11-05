

df<-read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/BLUPsDay6_thinned_Oct19_geneticdistance.csv")
distance  <- as.matrix(dist(scale(df)))
write.csv(distance,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/BigDataBLUPs_Genetic_Distance_Day6.csv", row.names = T)

df<-read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/BLUPsDay9_thinned_Oct19_geneticdistance.csv")
distance  <- as.matrix(dist(scale(df)))
write.csv(distance,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/BigDataBLUPs_Genetic_Distance.csv", row.names = T)

df<-read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/BLUPsDay12_thinned_Oct19_geneticdistance.csv")
distance  <- as.matrix(dist(scale(df)))
write.csv(distance,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/BigDataBLUPs_Genetic_Distance_Day12.csv", row.names = T)
