  rm(list=ls(all=TRUE))
  
  #wd = getwd()
  #wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
  #setwd(wd)
  
  #CytB 42975
  #ATP6 8415
  #ATP8 1890
  #ND1 5355
  #ND2 13770
  #ND3 2745
  #ND4 4500
  #ND4L 1800
  #ND5 4095
  #ND6 225
  #COX1 3780
  #COX2 3420
  #COX3 2880  
  
  Grantham = read.csv("../../Body/2Derived/COX3.csv", sep='\t', header = T) 
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
  
  write.table(df5,file = "../../Body/2Derived/COX3_Distances.csv",quote = F, row.names = FALSE,sep = '\t')
  
  
  
  library(ape)
  library(gdata)
  library(ggplot2)
  
  Data = read.table("../../Body/2Derived/Distances_KnKs_RG.csv", sep='\t', header = TRUE) 
  GenLength<- read.xls("../../Body/1Raw/GenerationLengthForMammals.xlsx")# табличка с продолжительностью жизни от Алины
  #GenLength <- GenLength1[1:1,]
  Data$Species = sub("_", " ", Data$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
  
  Data = merge(Data,GenLength, by.x = "Species", by.y = "Scientific_name",all = FALSE,no.dups = TRUE,)
  
  Data <- Data[,-11:-20]
  Data <- Data[,-12]
  Data <- Data[,-9]
  write.table(Data,file = "../../Body/2Derived/CytB_Distances.csv",quote = F, row.names = FALSE,sep = '\t')
  
  #############################Табличка с грантхэмом, кнкс и экологией (в экологии будут пропуски) для Александра Купцова
  
  Data = read.table("../../Body/2Derived/Distances_KnKs_RG.csv", sep='\t', header = TRUE) 
  GenLength <- read.xls("../../Body/1Raw/GenerationLengthForMammals.xlsx")# табличка с продолжительностью жизни от Алины
  Taxa <- read.csv("../../Body/1Raw/TaxaWithClasses.txt", sep = " ")# таксономия от Али
  
  setdiff(Data$Species,Taxa$Species)# В таксономии Али не хватает 33 вида
  
  Data = merge(Data,Taxa, by.x = "Species", by.y = "Species",all.x =TRUE,no.dups = TRUE,)
  
  #	Canis_lupus_familiaris - домашняя собака
  # Ctenomys_sp - какой то грызун
  # Cynopterus_sp - летучая мышь
  # Eospalax_baileyi - крот
  # Eospalax_cansus - еще крот
  # Gorilla_gorilla_gorilla 
  # Mammuthus_sp - маммонт
  # Mastomys_sp - грызун
  # Muntiacus_sp - 
  # Mus_musculus_castaneus
  # Mus_musculus_domesticus
  # Mus_musculus_molossinus
  # Mus_musculus_musculus
  # Myotis_sp - еще летучая мышь
  # Niviventer_sp
  # Ochotona_sp
  # Peromyscus_sp
  # Pteronotus_sp - еще летучая мышь
  # Sus_scrofa_domesticus - домашняя свинья, надо удалять(
  
  #Data[Data$Species == "Canis_lupus_familiaris",]$Class = "Mammalia"
  Data[Data$Species == "Ctenomys_sp",]$Class = "Mammalia"
  Data[Data$Species == "Cynopterus_sp",]$Class = "Mammalia"
  Data[Data$Species == "Eospalax_baileyi",]$Class = "Mammalia"
  Data[Data$Species == "Eospalax_cansus",]$Class = "Mammalia"
  Data[Data$Species == "Gorilla_gorilla_gorilla",]$Class = "Mammalia"
  Data[Data$Species == "Mammuthus_sp",]$Class = "Mammalia"
  Data[Data$Species == "Mastomys_sp",]$Class = "Mammalia"
  Data[Data$Species == "Muntiacus_sp",]$Class = "Mammalia"
  Data[Data$Species == "Mus_musculus_castaneus",]$Class = "Mammalia"
  #Data[Data$Species == "Mus_musculus_domesticus",]$Class = "Mammalia"
  Data[Data$Species == "Mus_musculus_molossinus",]$Class = "Mammalia"
  Data[Data$Species == "Mus_musculus_musculus",]$Class = "Mammalia"
  Data[Data$Species == "Myotis_sp",]$Class = "Mammalia"
  Data[Data$Species == "Niviventer_sp",]$Class = "Mammalia"
  Data[Data$Species == "Ochotona_sp",]$Class = "Mammalia"
  Data[Data$Species == "Peromyscus_sp",]$Class = "Mammalia"
  Data[Data$Species == "Pteronotus_sp",]$Class = "Mammalia"
  #Data[Data$Species == "Sus_scrofa_domesticus",]$Class = "Mammalia"
  
  Data = Data[Data$Species == "Canis_lupus_familiaris",] # попытка удалить домашних
  
  Data = Data[Data$Class == "Mammalia",]
  Data1 = Data[is.na(Data$AverageGrantham)==FALSE,]
  Data1$Species = sub("_", " ", Data1$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
  
  NewData = merge(Data1,GenLength, by.x = "Species", by.y = "Scientific_name",all.x  = TRUE,no.dups = TRUE,)
  
  NewData <- NewData[,-12:-21]
  NewData <- NewData[,-13]
  NewData <- NewData[,-10]
  NewData <- NewData[,-9]

  IUCN = read.csv("../../Body/1Raw/Red_book/IUCN.csv", sep=';', header = TRUE) #табличка по красной книге от Алины https://github.com/mitoclub/red-book/blob/master/Body/1Raw/IUCN.csv
  names(IUCN)[3] <- "Species"
  
  Dist_with_IUCN <- merge(NewData, IUCN,by.x = "Species", by.y = "Species",all.x  = T,no.dups = TRUE,)
  
  Dist_with_IUCN = Dist_with_IUCN [-11:-21]
  Dist_with_IUCN = Dist_with_IUCN [,-11]
  Dist_with_IUCN = Dist_with_IUCN [-12:-28]# тут обрезается куча данных, наверное, ненужных. В переменной category информация о том, в каком состоянии вид
  
  
  write.table(Dist_with_IUCN,file = "../../Body/2Derived/Grantham_RedBook.csv",quote = F, row.names = FALSE,sep = '\t')

  ##################### Корреляции
  
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
  cor.test(x = Data$GenerationLength_d, y = Data$MedianOfAllNonsyn , method = "spearm")
  GenerationLength_MedianOfAllNonsyn_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$MedianOfAllNonsyn)
  
  ggplot(Data, aes(x = log(GenerationLength_d), y = MedianOfAllNonsyn ,col = factor(Order) ))+
    geom_point(size = 2)
  
  ggplot(Data, aes(x = log(GenerationLength_d), y =MedianOfAllNonsyn , fill = Order)) + 
    geom_boxplot()
  
  #(v) SummOfAllGrantham ~+SummOfAllNonsyn (проверка на вшивость)
  cor.test(x = Data$SummOfAllGrantham, y = Data$MedianOfAllNonsyn, method = "spearm")
  GenerationLength_MedianOfAllNonsyn_fit  <- cor.test(x = Data$SummOfAllGrantham, y = Data$MedianOfAllNonsyn)
  
  ggplot(Data, aes(x = log(SummOfAllGrantham), y = MedianOfAllNonsyn ,col = factor(Order) ))+
    geom_point(size = 2)
  
  ggplot(Data, aes(x = log(SummOfAllGrantham), y = MedianOfAllNonsyn , fill = Order)) + 
    geom_boxplot()
  
  #(vi) AverageGranthamN ~ +GenerLength; 
  cor.test(x = Data$GenerationLength_d, y = Data$AverageGrantham , method = "spearm")
  GenerationLength_AverageGrantham_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$AverageGrantham)
  
  ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGrantham ,col = factor(Order) ))+
    geom_point(size = 2)
  
  ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGrantham, fill = Order)) + 
    geom_boxplot()
  
  #(vii) AverageGrantham ~ +GenerLength; 
  cor.test(x = Data$GenerationLength_d, y = Data$AverageGrantham , method = "spearm")
  GenerationLength_AverageGrantham_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$AverageGrantham)
  
  ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGrantham ,col = factor(Order) ))+
    geom_point(size = 2)
  
  ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGrantham, fill = Order)) + 
    geom_boxplot()
  
  #(viii) KnKs ~ +GenerLength; 
  cor.test(x = Data$GenerationLength_d, y = Data$KnKs , method = "spearm")
  GenerationLength_KnKs_fit  <- cor.test(x = Data$GenerationLength_d, y = Data$KnKs)
  
  ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs ,col = factor(Order) ))+
    geom_point(size = 2)
  
  ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs, fill = Order)) + 
    geom_boxplot()
####################
  
quantile(Data$GenerationLength_d, seq(from=0, to=1, by=0.25),na.rm = TRUE)  

boxplot(Data[Data$GenerationLength_d < 637.3329,]$AverageGrantham, +
          Data[Data$GenerationLength_d < 1825.0000 & Data$GenerationLength_d > 637.3329,]$AverageGrantham, +
          Data[Data$GenerationLength_d < 2869.6933 & Data$GenerationLength_d > 1825.0000,]$AverageGrantham, +
          Data[Data$GenerationLength_d > 2869.6933,]$AverageGrantham, names = c('Q1','Q2', 'Q3','Q4'), outline = FALSE, xlab = "GenLenght", ylab = "AverageGrantham")

boxplot(Data[Data$GenerationLength_d < 637.3329,]$KnKs, +
          Data[Data$GenerationLength_d < 1825.0000 & Data$GenerationLength_d > 637.3329,]$KnKs, +
          Data[Data$GenerationLength_d < 2869.6933 & Data$GenerationLength_d > 1825.0000,]$KnKs, +
          Data[Data$GenerationLength_d > 2869.6933,]$KnKs, names = c('Q1','Q2', 'Q3','Q4'), outline = FALSE, xlab = "GenLenght", ylab = "KnKs")

boxplot(Data[Data$Order == "Eulipotyphla",]$KnKs,+
          Data[Data$Order == "Rodentia",]$KnKs, +
          Data[Data$Order == "Chiroptera",]$KnKs, +
          Data[Data$Order == "Carnivora",]$KnKs, +
          Data[Data$Order == "Primates",]$KnKs, +
          Data[Data$Order =="Cetartiodactyla",]$KnKs, names = c('Eulipotyphla',"Rodentia",'Chiroptera','Carnivora',"Primates",'Cetartiodactyla'), outline = FALSE, xlab = "GenLenght", ylab = "KnKs")

TempData = subset(Data, Data$Order == 'Eulipotyphla'|Data$Order =="Rodentia"|Data$Order =='Chiroptera'|Data$Order =='Carnivora'|Data$Order =="Primates"|Data$Order =='Cetartiodactyla')

ggplot(TempData, aes(x=Order, y=KnKs,fill=Order)) + 
     geom_violin(trim=FALSE)+ 
  stat_summary(fun=median, geom="crossbar", size=1, color="red")+
  labs(title="KnKs~GenLen",x="Order", y = "KnKs")+
  theme_classic()

boxplot(Data[Data$Order == "Eulipotyphla",]$AverageGrantham,+
          Data[Data$Order == "Rodentia",]$AverageGrantham, +
          Data[Data$Order == "Chiroptera",]$AverageGrantham, +
          Data[Data$Order == "Carnivora",]$AverageGrantham, +
          Data[Data$Order == "Primates",]$AverageGrantham, +
          Data[Data$Order =="Cetartiodactyla",]$AverageGrantham, names = c('Eulipotyphla',"Rodentia",'Chiroptera','Carnivora',"Primates",'Cetartiodactyla'), outline = FALSE, xlab = "GenLenght", ylab = "AverageGrantham")

ggplot(TempData, aes(x=Order, y=AverageGrantham,fill=Order)) + 
  geom_violin(trim=FALSE)+ 
  stat_summary(fun=median, geom="crossbar", size=1, color="red")+
  labs(title="AverageGrantham~GenLen",x="Order", y = "AverageGrantham")+
  theme_classic()

boxplot(Data[Data$Order == "Carnivora",]$GenerationLength_d,+
          Data[Data$Order == "Cetartiodactyla",]$GenerationLength_d, +
          Data[Data$Order == "Chiroptera",]$GenerationLength_d, +
          Data[Data$Order == "Eulipotyphla",]$GenerationLength_d, +
          Data[Data$Order == "Primates",]$GenerationLength_d, +
          Data[Data$Order =="Rodentia",]$GenerationLength_d, names = c('Eulipotyphla',"Rodentia",'Chiroptera','Carnivora',"Primates",'Cetartiodactyla'), outline = FALSE, xlab = "GenLenght", ylab = "G")

median(Data[Data$Order == "Primates",]$GenerationLength_d)
median(Data[Data$Order =="Cetartiodactyla",]$GenerationLength_d)# продолжительность жизни китопарнокопытных чуть больше, чем у приматов
# Carnivora Cetartiodactyla Chiroptera Eulipotyphla Primates Rodentia 


