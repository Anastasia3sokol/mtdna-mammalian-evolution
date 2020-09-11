rm(list=ls(all=TRUE))

Grantham = read.table("/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Grantham.csv", sep='\t', header = TRUE) 
#Grantham = Grantham1[1:10,]

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
  Grantham[i, "only_distances_temp"] <- paste(Output, collapse=";")
  } }
#################################### 

    for(i in 1:nrow(Grantham)){
      vec = unlist(strsplit(Grantham$only_distances_temp[i],';'))
        if (length(vec)==0){N=N+1;next
      }else{
          Distances <- c()
      for (j in 1:length(vec)){
    AA_with_distance = vec[j]
    Distance = gsub("(.*)\\:",'',AA_with_distance)# оставляет все после знака :
    Distances = c(Distances, Distance)
 }
    Grantham[i, "distances_temp"] <- paste(Distances, collapse=";") 
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
           #if (string.counter(Grantham$All_Synon[i], pattern=0)) {next} else
     # {
        Grantham$KnKs[i] = as.numeric(Grantham$All_Non_Synon[i]) / as.numeric(Grantham$All_Synon[i])
      #} 
    }
   
 for (i in 1:nrow(Grantham)){
   Grantham$KnKs[i][is.infinite(Grantham$KnKs[i])] = NA
 }
   
for(i in 1:nrow(Grantham)){
  Grantham$FractionOfSyn[i] =  as.numeric(Grantham$All_Synon[i])/ as.numeric(Grantham$TotalDiv)
  Grantham$FractionOfNonsyn[i] =  as.numeric(Grantham$All_Non_Synon[i])/ as.numeric(Grantham$TotalDiv)
    }
###################################################
          for(i in 1:nrow(Grantham)){
      vec = unlist(strsplit(Grantham$distances_temp[i],';'))
      if (length(vec)==0){N=N+1;next
      }else{
        
        SummOfAllGrantham <- 0
      for (j in 1:length(vec)){
        num = as.numeric(vec[j])
        NS = Grantham$All_Non_Synon[j]
        SummOfAllGrantham = (SummOfAllGrantham + num)
        
      }
      Grantham[i, "SummOfAllGrantham"] <- paste((as.numeric(SummOfAllGrantham)), collapse=";") 
      Grantham[i,"AverageGrantham"]<- paste(as.numeric(SummOfAllGrantham)/as.numeric(Grantham$All_Non_Synon[i]))
    }
    }

KnKs <-aggregate(as.numeric(Grantham$KnKs),list(Grantham$Species),mean,na.rm=TRUE)
colnames(KnKs) <- c("Species","KnKs")

SummOfAllGrantham <- aggregate(as.numeric(Grantham$SummOfAllGrantham),list(Grantham$Species),mean,na.rm=TRUE)
colnames(SummOfAllGrantham) <- c("Species","SummOfAllGrantham")

AverageGrantham <- aggregate(as.numeric(Grantham$AverageGrantham),list(Grantham$Species),mean,na.rm=TRUE)
colnames(AverageGrantham) <- c("Species","AverageGrantham")

#Distance_for_each_species$Species = sub("_", " ", data$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

df <- merge(AverageGrantham,
            SummOfAllGrantham,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

df1 <- merge(df,
            KnKs,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)





MeanOfAllSyn <- (aggregate((as.numeric(Grantham$All_Synon)),list(Grantham$Species),mean,na.rm = TRUE))
colnames(MeanOfAllSyn) <- c("Species","MeanOfAllSyn")
MeanOfAllNonsyn <- (aggregate((as.numeric(Grantham$All_Non_Synon)),list(Grantham$Species),mean,na.rm = TRUE))
colnames(MeanOfAllNonsyn) <- c("Species","MeanOfAllNonsyn")

df2 <- merge(MeanOfAllSyn,MeanOfAllNonsyn,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)
Data1 <-merge(df1,df2,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

FractionOfSyn <-aggregate(as.numeric(Grantham$FractionOfSyn),list(Grantham$Species),mean,na.rm=TRUE)
colnames(FractionOfSyn) <- c("Species","FractionOfSyn")

FractionOfNonsyn <-aggregate(as.numeric(Grantham$FractionOfNonsyn),list(Grantham$Species),mean,na.rm=TRUE)
colnames(FractionOfNonsyn) <- c("Species","FractionOfNonsyn")

Data2 <-merge(FractionOfSyn,FractionOfNonsyn,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

Data <-merge(Data1,Data2,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

write.table(Data,file = "/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Distances_KnKs.csv",quote = F, row.names = FALSE,sep = '\t')



library(ape)
library(gdata)
library(ggplot2)

data = read.table("/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Distances_KnKs.csv", sep='\t', header = TRUE) 
GenLength <- read.xls("/home/anastasia/mtdna-mammalian-evolution/Body/1Raw/GenerationLengthForMammals.xlsx")# табличка с продолжительностью жизни от Алины

data$Species = sub("_", " ", data$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

Data = merge(data,GenLength, by.x = "Species", by.y = "Scientific_name",all = FALSE,no.dups = TRUE,)
#Data <- NewData[,-13:-19]
#Data <- Data[,-14]
#Data <- Data[,-8]

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

#GenerationLength and AverageGrantham

cor.test(x = Data$GenerationLength_d, y = Data$AverageGrantham, method = "spearm", na.action = "na.exclude")
GenerationLengthAverageGrantham_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$AverageGrantham)

plot(x = Data$GenerationLength_d, y = Data$AverageGrantham)

ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGrantham)) + 
  geom_point() +
  geom_smooth()

ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGrantham, col = factor(Order))) + 
  geom_point()

ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGrantham, fill = Order)) + 
  geom_boxplot()



