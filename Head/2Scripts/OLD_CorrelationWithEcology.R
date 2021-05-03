# Скрипт не нужен?

  rm(list=ls(all=TRUE))
  
  wd = getwd()
  wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
  setwd(wd)
  
  Grantham = read.csv("../../Body/2Derived/CytB.csv", sep='\t', header = T) 
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
  
  df5 = df5[is.na(df5$AverageGrantham)==F,]
  
  write.table(df5,file = "../../Body/2Derived/CytB_Distances.csv",quote = F, row.names = FALSE,sep = '\t')

  
  library(ape)
  library(gdata)
  library(gtools)
  
  Distances = read.table("../../Body/2Derived/CytB_Distances.csv", sep='\t', header = TRUE) 
  GenLength<- read.xls("../../Body/1Raw/GenerationLengthForMammals.xlsx")# табличка с продолжительностью жизни от Алины
  
  Distances$Species = sub("_", " ", Distances$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
  Data = merge(Distances,GenLength, by.x = "Species", by.y = "Scientific_name",all = FALSE,no.dups = TRUE,)
  
  RB = merge(Distances,GenLength, by.x = "Species", by.y = "Scientific_name",all.x = T,no.dups = TRUE,)
  Data <- Data[,-12:-20]
  Data <- Data[,-13]
  Data <- Data[,-9]
  Data = Data[is.na(Data$AverageGrantham)==F,]
 
  Data = Data[-220,]#дублирующийся Neophocaena phocaenoides
  
  RB <- RB[,-12:-20]
  RB <- RB[,-13]
  RB <- RB[,-9]
  RB = RB[is.na(Data$AverageGrantham)==F,]
  
  dnds = read.csv("../../Body/2Derived/DnDs.txt", sep = " ")
  dnds$Species = sub("_", " ", dnds$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
  
  Data1 = merge(Data,dnds, by.x = "Species", by.y = "Species",all.x = T,no.dups = TRUE,)
  names(Data1) = c("1.Species","2.AverageGrantham","SummOfAllGrantham","3.KnKs","MedianOfAllSyn","MedianOfAllNonsyn","FractionOfSyn","FractionOfNonsyn","Order",'family',"5.GenerationLength_d","Gene","4.DnDs" )
  Data = Data1[c(1, mixedorder(names(Data1)[-1]) + 1)];names(Data) = c("Species","AverageGrantham","KnKs","DnDs","GenerationLength_d",'family',"FractionOfNonsyn","FractionOfSyn","Gene","MedianOfAllNonsyn","MedianOfAllSyn","Order","SummOfAllGrantham" )
  
  Data2 = merge(RB,dnds, by.x = "Species", by.y = "Species",all=F,no.dups = TRUE,)
  names(Data2) = c("1.Species","2.AverageGrantham","SummOfAllGrantham","3.KnKs","MedianOfAllSyn","MedianOfAllNonsyn","FractionOfSyn","FractionOfNonsyn","Order",'family',"5.GenerationLength_d","Gene","4.DnDs" )
  Data2 = Data2[c(1, mixedorder(names(Data2)[-1]) + 1)];names(Data2) = c("Species","AverageGrantham","KnKs","DnDs","GenerationLength_d",'family',"FractionOfNonsyn","FractionOfSyn","Gene","MedianOfAllNonsyn","MedianOfAllSyn","Order","SummOfAllGrantham" )
  Data = Data[is.na(Data$DnDs)==F,]
  Data2 = Data2[is.na(Data2$DnDs)==F,]
  write.table(Data,file = "../../Body/2Derived/CytB_with_ecology.csv",quote = F, row.names = FALSE,sep = '\t')
  write.table(Data2,file = "../../Body/2Derived/CytB_with_ecologyForRedBook.csv",quote = F, row.names = FALSE,sep = '\t')
 