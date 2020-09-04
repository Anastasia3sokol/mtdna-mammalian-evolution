rm(list=ls(all=TRUE))

Grantham = read.table("/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Grantham_1.csv", sep='\t', header = TRUE) 

    string.counter<-function(strings, pattern){  
  counts<-NULL
  for(i in 1:length(strings)){
    counts[i]<-length(attr(gregexpr(pattern,strings[i])[[1]], "match.length")[attr(gregexpr(pattern,strings[i])[[1]], "match.length")>0])
  }
  return(counts)
}
################################### 
  for(i in 1:nrow(Grantham)){
  vec = unlist(strsplit(Grantham$NewDistance[i],';'))
  N <- c()
  Output <- c()
  
  if (length(vec)==0){N=N+1;next
  }else{
      for (j in 1:length(vec)){
     if (string.counter(vec[j], pattern=">")) {Output = c(Output, vec[j])
    }   }
  Grantham[i, "only_distances"] <- paste(Output, collapse=";")
  } }
#################################### 

    for(i in 1:nrow(Grantham)){
      vec = unlist(strsplit(Grantham$only_distances[i],';'))
        if (length(vec)==0){N=N+1;next
      }else{
          Distances <- c()
      for (j in 1:length(vec)){
    AA_with_distance = vec[j]
    Distance = gsub("(.*)\\:",'',AA_with_distance)# оставляет все после знака :
    Distances = c(Distances, Distance)
 }
    Grantham[i, "distances"] <- paste(Distances, collapse=";") 
    }    }
################################################
    for(i in 1:nrow(Grantham)){
      
      vec = unlist(strsplit(Grantham$NewDistance[i],';'))
      Number <- c()
      if (length(vec)==0){N=N+1;next
      }else{
             for (j in 1:length(vec)){
            if (string.counter(vec[j], pattern="Non_Synon:")) {
         Number = gsub("(.*)\\:",'',vec[j])
          }
        Grantham[i, "All_Non_Synon"] <- paste(Number, collapse=";")
      }
     }
    }
    for(i in 1:nrow(Grantham)){
      Grantham$All_Synon[i] <-as.numeric(Grantham$TotalDiv[i]) - as.numeric(Grantham$All_Non_Synon[i])
    }
    for(i in 1:nrow(Grantham)){
           if (string.counter(Grantham$All_Synon[i], pattern=0)) {next} else
      {Grantham$KnKs[i] = as.numeric(Grantham$All_Non_Synon[i]) / as.numeric(Grantham$All_Synon[i])
      } 
    }
   
 for (i in 1:nrow(Grantham)){
   Grantham$KnKs[i][is.infinite(Grantham$KnKs[i])] = NA
      }
###################################################
        for(i in 1:nrow(Grantham)){
      vec = unlist(strsplit(Grantham$distances[i],';'))
      if (length(vec)==0){N=N+1;next
      }else{
        
      Distances_mean <- 0
      for (j in 1:length(vec)){
        num = as.numeric(vec[j])
        Distances_mean = (Distances_mean + num)
        
      }
      Grantham[i, "distances_mean"] <- paste((as.numeric(Distances_mean/length(vec))), collapse=";") 
    }
    }

#Grantham_1= Grantham[, -8:-9]

KnKs <-aggregate(as.numeric(Grantham$KnKs),list(Grantham$Species),mean,na.rm=TRUE)
colnames(KnKs) <- c("Species","KnKs")

Distance_for_each_species <- aggregate(as.numeric(Grantham$distances_mean),list(Grantham$Species),mean,na.rm=TRUE)
colnames(Distance_for_each_species) <- c("Species","Distance")
#Distance_for_each_species$Species = sub("_", " ", data$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

df <- merge(Distance_for_each_species,
            KnKs,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)


S <- (aggregate((as.numeric(Grantham$All_Synon)),list(Grantham$Species),mean,na.rm = TRUE))
colnames(S) <- c("Species","Synon")
N <- (aggregate((as.numeric(Grantham$All_Non_Synon)),list(Grantham$Species),mean,na.rm = TRUE))
colnames(N) <- c("Species","Non_Synon")

df1 <- merge(S,N,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)
data <-merge(df,df1,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

for(i in 1:nrow(data)){
 data$proportion_of_S[i] = (data$Synon[i])/(as.numeric(data$Synon[i]) + as.numeric(data$Non_Synon[i]))
 data$proportion_of_NS[i] = (data$Non_Synon[i])/(as.numeric(data$Synon[i]) + as.numeric(data$Non_Synon[i]))
}
write.table(data,file = "/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Distances&KnKs.csv",quote = F, row.names = FALSE,sep = '\t')


library(ape)
library(gdata)
library(ggplot2)

data = read.table("/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Distances&KnKs.csv", sep='\t', header = TRUE) 
GenLength <- read.xls("/home/anastasia/mtdna-mammalian-evolution/Body/1Raw/GenerationLengthForMammals.xlsx")# табличка с продолжительностью жизни от Алины

data$Species = sub("_", " ", data$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

NewData = merge(data,GenLength, by.x = "Species", by.y = "Scientific_name",all = FALSE,no.dups = TRUE,)
Data <- NewData[,-13:-19]
Data <- Data[,-14]
Data <- Data[,-8]

#GenerationLength_Distances

cor.test(x = Data$GenerationLength_d, y = Data$Distance, method = "spearm")
GenerationLength_Distances_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$Distance)

plot(x = Data$GenerationLength_d, y = Data$Distance)

ggplot(Data, aes(x = log(GenerationLength_d), y = Distance,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(Data, aes(x = log(GenerationLength_d), y = Distance, fill = Order)) + 
  geom_boxplot()

#GenerationLength_KnKs

cor.test(x = Data$GenerationLength_d, y = Data$KnKs, method = "spearm")
GenerationLength_KnKs_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$KnKs)

plot(x = Data$GenerationLength_d, y = Data$KnKs)

ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs)) + 
  geom_point() +
  geom_smooth()

ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs, col = factor(Order))) + 
  geom_point()

ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs, fill = Order)) + 
  geom_boxplot()

#GenerationLength_proportion_of_S

cor.test(x = Data$GenerationLength_d, y = Data$proportion_of_S, method = "spearm")
GenerationLength_proportion_of_S_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$proportion_of_S)

plot(x = Data$GenerationLength_d, y = Data$proportion_of_S)

ggplot(Data, aes(x = log(GenerationLength_d), y = proportion_of_S)) + 
  geom_point() +
  geom_smooth()

ggplot(Data, aes(x = log(GenerationLength_d), y = proportion_of_S, col = factor(Order))) + 
  geom_point()

ggplot(Data, aes(x = log(GenerationLength_d), y = proportion_of_S, fill = Order)) + 
  geom_boxplot()

#GenerationLength_proportion_of_NS

cor.test(x = Data$GenerationLength_d, y = Data$proportion_of_NS, method = "spearm")
GenerationLength_proportion_of_NS_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$proportion_of_NS)

plot(x = Data$GenerationLength_d, y = Data$proportion_of_NS)

ggplot(Data, aes(x = log(GenerationLength_d), y = proportion_of_NS)) + 
  geom_point() +
  geom_smooth()

ggplot(Data, aes(x = log(GenerationLength_d), y = proportion_of_NS, col = factor(Order))) + 
  geom_point()

ggplot(Data, aes(x = log(GenerationLength_d), y = proportion_of_NS, fill = Order)) + 
  geom_boxplot()

