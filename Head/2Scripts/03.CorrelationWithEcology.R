rm(list=ls(all=TRUE))

Grantham = read.table("/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Grantham.csv", sep='\t', header = TRUE) 
#Grantham = Grantham1[1:5000,]

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
  Grantham$FractionOfSyn[i] =  as.numeric(Grantham$All_Synon[i])/ as.numeric(Grantham$TotalDiv[i])
  Grantham$FractionOfNonsyn[i] =  as.numeric(Grantham$All_Non_Synon[i])/ as.numeric(Grantham$TotalDiv[i])
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
      Grantham[i,"SummOfAllGrantham"] <- paste((as.numeric(SummOfAllGrantham)), collapse=";") 
      Grantham[i,"AverageGrantham"]<- paste(as.numeric(SummOfAllGrantham)/as.numeric(Grantham$All_Non_Synon[i]))
    }
    }

KnKs <-aggregate(as.numeric(Grantham$KnKs),list(Grantham$Species),median,na.rm=TRUE)
colnames(KnKs) <- c("Species","KnKs")

SummOfAllGrantham <- aggregate(as.numeric(Grantham$SummOfAllGrantham),list(Grantham$Species),median,na.rm=TRUE)
colnames(SummOfAllGrantham) <- c("Species","SummOfAllGrantham")

AverageGrantham <- aggregate(as.numeric(Grantham$AverageGrantham),list(Grantham$Species),median,na.rm=TRUE)
colnames(AverageGrantham) <- c("Species","AverageGrantham")

MedianOfAllSyn <- (aggregate((as.numeric(Grantham$All_Synon)),list(Grantham$Species),median,na.rm = TRUE))
colnames(MedianOfAllSyn) <- c("Species","MedianOfAllSyn")
MedianOfAllNonsyn <- (aggregate((as.numeric(Grantham$All_Non_Synon)),list(Grantham$Species),median,na.rm = TRUE))
colnames(MedianOfAllNonsyn) <- c("Species","MedianOfAllNonsyn")

FractionOfSyn <-aggregate(as.numeric(Grantham$FractionOfSyn),list(Grantham$Species),median,na.rm=TRUE)
colnames(FractionOfSyn) <- c("Species","FractionOfSyn")
FractionOfNonsyn <-aggregate(as.numeric(Grantham$FractionOfNonsyn),list(Grantham$Species),median,na.rm=TRUE)
colnames(FractionOfNonsyn) <- c("Species","FractionOfNonsyn")


df <- merge(AverageGrantham, SummOfAllGrantham,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

df1 <- merge(df,KnKs,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

df2 <- merge(MedianOfAllSyn,MedianOfAllNonsyn,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

df3 <- merge(FractionOfSyn,FractionOfNonsyn,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

df4 <- merge(df1,df2,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

df5 <- merge(df4,df3,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

write.table(df5,file = "/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Distances_KnKs_with_W.csv",quote = F, row.names = FALSE,sep = '\t')


library(ape)
library(gdata)
library(ggplot2)

Data = read.table("/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Distances_KnKs_with_W.csv", sep='\t', header = TRUE) 
GenLength<- read.xls("/home/anastasia/mtdna-mammalian-evolution/Body/1Raw/GenerationLengthForMammals.xlsx")# табличка с продолжительностью жизни от Алины
#GenLength <- GenLength1[1:1,]
Data$Species = sub("_", " ", Data$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

Data = merge(Data,GenLength, by.x = "Species", by.y = "Scientific_name",all = FALSE,no.dups = TRUE,)
#Data <- NewData[,-13:-19]
#Data <- Data[,-14]
#Data <- Data[,-8]

#GenerationLength_Distances

cor.test(x = Data$GenerationLength_d, y = Data$SummOfAllGrantham , method = "spearm")
GenerationLength_SummOfAllGrantham_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$SummOfAllGrantham )

ggplot(Data, aes(x = log(GenerationLength_d), y = SummOfAllGrantham ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(Data, aes(x = log(GenerationLength_d), y = SummOfAllGrantham , fill = Order)) + 
  geom_boxplot()

#GenerationLength_KnKs

cor.test(x = Data$GenerationLength_d, y = Data$KnKs, method = "spearm")
GenerationLength_KnKs_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$KnKs)

ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs)) + 
  geom_point() +
  geom_smooth()

ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs, col = factor(Order))) + 
  geom_point()

ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs, fill = Order)) + 
  geom_boxplot()

#GenerationLength_proportion_of_S

cor.test(x = Data$GenerationLength_d, y = Data$FractionOfSyn, method = "spearm")
GenerationLength_FractionOfSyn_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$FractionOfSyn)

ggplot(Data, aes(x = log(GenerationLength_d), y = FractionOfSyn)) + 
  geom_point() +
  geom_smooth()

ggplot(Data, aes(x = log(GenerationLength_d), y = FractionOfSyn, col = factor(Order))) + 
  geom_point()

ggplot(Data, aes(x = log(GenerationLength_d), y = FractionOfSyn, fill = Order)) + 
  geom_boxplot()

#GenerationLength_proportion_of_NS

cor.test(x = Data$GenerationLength_d, y = Data$FractionOfNonsyn, method = "spearm")
GenerationLength_FractionOfNonsyn_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$FractionOfNonsyn)

ggplot(Data, aes(x = log(GenerationLength_d), y = FractionOfNonsyn)) + 
  geom_point() +
  geom_smooth()

ggplot(Data, aes(x = log(GenerationLength_d), y = FractionOfNonsyn, col = factor(Order))) + 
  geom_point()

ggplot(Data, aes(x = log(GenerationLength_d), y = FractionOfNonsyn, fill = Order)) + 
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

#SummOfAllNonsyn ~ +GenerLength 

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












#(i) FractionOfSyn ~ -GenerLength; 

cor.test(x = Data$GenerationLength_d, y = Data$FractionOfSyn , method = "spearm")
GenerationLength_FractionOfSyn_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$FractionOfSyn)

ggplot(Data, aes(x = log(GenerationLength_d), y = FractionOfSyn ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(Data, aes(x = log(GenerationLength_d), y = FractionOfSyn , fill = Order)) + 
  geom_boxplot()


#(ii) FractionOfNonsyn ~ +GenerLength; 

cor.test(x = Data$GenerationLength_d, y = Data$FractionOfNonsyn , method = "spearm")
GenerationLength_FractionOfNonsyn_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$FractionOfNonsyn)

ggplot(Data, aes(x = log(GenerationLength_d), y = FractionOfNonsyn ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(Data, aes(x = log(GenerationLength_d), y = FractionOfNonsyn , fill = Order)) + 
  geom_boxplot()

#(iii) SummOfAllGrantham ~ +GenerLength; 
cor.test(x = Data$GenerationLength_d, y = Data$SummOfAllGrantham , method = "spearm")
GenerationLength_SummOfAllGrantham_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$SummOfAllGrantham)

ggplot(Data, aes(x = log(GenerationLength_d), y = SummOfAllGrantham ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(Data, aes(x = log(GenerationLength_d), y = SummOfAllGrantham , fill = Order)) + 
  geom_boxplot()

#(iv) SummOfAllNonsyn ~ +GenerLength 
cor.test(x = Data$GenerationLength_d, y = Data$MeanOfAllNonsyn , method = "spearm")
GenerationLength_MeanOfAllNonsyn_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$MeanOfAllNonsyn)

ggplot(Data, aes(x = log(GenerationLength_d), y = MeanOfAllNonsyn ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(Data, aes(x = log(GenerationLength_d), y =MeanOfAllNonsyn , fill = Order)) + 
  geom_boxplot()

#(v) SummOfAllGrantham ~+SummOfAllNonsyn (проверка на вшивость)
cor.test(x = Data$SummOfAllGrantham, y = Data$MeanOfAllNonsyn, method = "spearm")
GenerationLength_MeanOfAllNonsyn_fit  <- cor.test(x = Data$SummOfAllGrantham, y = Data$MeanOfAllNonsyn)

ggplot(Data, aes(x = log(SummOfAllGrantham), y = MeanOfAllNonsyn ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(Data, aes(x = log(SummOfAllGrantham), y = MeanOfAllNonsyn , fill = Order)) + 
  geom_boxplot()

#(vi) AverageGranthamN ~ +GenerLength; 
cor.test(x = Data$GenerationLength_d, y = Data$AverageGranthamN , method = "spearm")
GenerationLength_AverageGranthamN_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$AverageGranthamN)

ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGranthamN ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGranthamN, fill = Order)) + 
  geom_boxplot()

#(vii) AverageGranthamNS ~ +GenerLength; 
cor.test(x = Data$GenerationLength_d, y = Data$AverageGranthamNS , method = "spearm")
GenerationLength_AverageGranthamNS_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$AverageGranthamNS)

ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGranthamNS ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGranthamNS, fill = Order)) + 
  geom_boxplot()

#(vii) KnKs ~ +GenerLength; 
cor.test(x = Data$GenerationLength_d, y = Data$KnKs , method = "spearm")
GenerationLength_AverageGranthamNS_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$KnKs)

ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs, fill = Order)) + 
  geom_boxplot()
