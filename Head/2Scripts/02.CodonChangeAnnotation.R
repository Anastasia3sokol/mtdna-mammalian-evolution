rm(list=ls(all=TRUE))

wd = getwd()
wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
setwd(wd)


#### libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("Biostrings")

library(seqinr)
library(Biostrings)
library(dplyr)
library(RGenetics)
SGC1 <- getGeneticCode("SGC1")  # Vertebrate Mitochondrial code for translate function

Codons2 = read.table("../../Body/2Derived/PolymorphicPairwiseCodons.txt") #датасет, созданный в скрипте 01,содержит вектор замещений SubstVec
Codons = Codons2[Codons2$Gene == 'CytB',]#42975
Codons = Codons[1:1000,]
#Функция, которая считает общее количество замещений
TotalDiv <- function(x) {Div = length(unlist(strsplit(x,';')))}
Codons$TotalDiv <- apply(as.matrix(Codons$SubstVec),1,FUN = TotalDiv) #посчитанные замещения прикрепляются в колонку TotalDiv в датасет Codons

Gr = read.table("../../Body/1Raw/Grantham - Sheet1.csv", sep=',', header = F)# таблица с дистанциями из статьи

Synon <- 0
Non_Synon <- 0
Composition <- ''
Polarity<- ""
Volume <- ""
Grantham <- ""

#x = "3:TGC>AGC;7:ACA>ACC;8:CAT>CAC;9:CCC>CCT;14:GCT>GCG;17:GCG>ACG;26:AAC>AGC;369:GCA>ACC;374:AAC>AAT;377:TTA>ATA;380:GCC>GCT;"

Translation <- function(x) 
{
  
  VecOfCodons = unlist(strsplit(x,';'));
  
  AminoSubs=c();
  Output <- c()
  Grantham = c()
  VecOfDistances = c()
  Indel <- 0
  N <- 0
  Check <- c("AAA","AAT","AAC","AAG","ATA","ATT","ATC","ATG","ACA","ACT","ACC","ACG","AGA","AGT","AGC","AGG","TAA","TAT","TAC","TAG","TTA","TTT","TTC","TTG","TCA","TCT","TCC","TCG","TGA","TGT","TGC","TGG","CAA","CAT","CAC","CAG","CTA","CTT","CTC","CTG","CCA","CCT","CCC","CCG","CGA","CGT","CGC","CGG","GAA","GAT","GAC","GAG","GTA","GTT","GTC","GTG","GCA","GCT","GCC","GCG","GGA","GGT","GGC","GGG")
  for (j in (1:length(VecOfCodons))) 
  {
    ##j =1
    CodonSubst = VecOfCodons[j]
    CodonSubst = gsub("(.*)\\:",'',CodonSubst)
    CodonSubst1 = gsub(">(.*)",'',CodonSubst)
    CodonSubst2 = gsub("(.*)>",'',CodonSubst)
    
    if (CodonSubst1 %in% Check & CodonSubst2 %in% Check ) {
      
      Codon1 <- DNAString(CodonSubst1)
      Codon2  <- DNAString(CodonSubst2)
      
      Codon1.Character = as.character(Codon1)
      Codon2.Character = as.character(Codon2)
      
      
      A1 = as.character(Biostrings::translate(Codon1, genetic.code=SGC1))
      A2 = as.character(Biostrings::translate(Codon2, genetic.code=SGC1))
      
      Subs = paste(A1,'>',A2,sep = '')
      
      AminoSubs = c(AminoSubs, Subs)  
      
      if (A1 == A2) {Synon <- Synon + 1} else {Non_Synon <- Non_Synon +1}
     
      GranthamNew <- Gr[which(Gr[,1]== A1), which(Gr[1,] == A2)]

      Grantham <<-c(Grantham,paste(GranthamNew))
      
      AminoGrantham = paste(Subs, ':', GranthamNew,';','Synon:', Synon,';','Non_Synon:', Non_Synon, sep='')
      
      VecOfDistances = c(VecOfDistances, AminoGrantham)
      
    }else{N=N+1;break}
  }
  return(VecOfDistances)
}

Codons$Distances = lapply(as.character(Codons$SubstVec), Translation)

for(i in 1:nrow(Codons)){
  Codons[i, 'NewDistance'] = paste(unlist(Codons$Distance[i]), collapse = ',')
}

#newCodons = Codons[, -7]
#write.table(newCodons,file = "../../Body/2Derived/Grantham.csv",quote = F, row.names = FALSE,sep = '\t')

###########################################################################################################################################################


#wd = getwd()
#wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
#setwd(wd)

#Grantham = read.csv("../../Body/2Derived/CytB.csv", sep='\t', header = T) 
Grantham = Codons

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

Distances = df5[is.na(df5$AverageGrantham)==F,]

write.table(Data,file = "../../Body/2Derived/CytB_Distances.csv",quote = F, row.names = FALSE,sep = '\t')


library(ape)
library(gdata)
library(gtools)

#Distances = read.table("../../Body/2Derived/CytB_Distances.csv", sep='\t', header = TRUE) 
GenLength<- read.xls("../../Body/1Raw/GenerationLengthForMammals.xlsx")# табличка с продолжительностью жизни от Алины

Distances$Species = sub("_", " ", Distances$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
DataGL = merge(Distances,GenLength, by.x = "Species", by.y = "Scientific_name",all = FALSE,no.dups = TRUE,)
DataGL = DataGL[-220,]#дублирующийся Neophocaena phocaenoides

DataRB = merge(Distances,GenLength, by.x = "Species", by.y = "Scientific_name",all.x = T,no.dups = TRUE,)#Здесь надо проследить чтобы виды не обрезались
DataGL <- DataGL[,-11:-20]
DataGL <- DataGL[,-12]
DataGL <- DataGL[,-9]

DataRB <- DataRB[,-11:-20]
DataRB <- DataRB[,-12]
DataRB <- DataRB[,-9]

dnds = read.csv("../../Body/2Derived/DnDs.txt", sep = " ")
dnds$Species = sub("_", " ", dnds$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

DataGL = merge(DataGL,dnds, by.x = "Species", by.y = "Species",all.x = T,no.dups = TRUE,)
names(DataGL) = c("1.Species","2.AverageGrantham","SummOfAllGrantham","3.KnKs","MedianOfAllSyn","MedianOfAllNonsyn","FractionOfSyn","FractionOfNonsyn","Order","5.GenerationLength_d","Gene","4.DnDs" )
Data = DataGL[c(1, mixedorder(names(DataGL)[-1]) + 1)];names(Data) = c("Species","AverageGrantham","KnKs","DnDs","GenerationLength_d","FractionOfNonsyn","FractionOfSyn","Gene","MedianOfAllNonsyn","MedianOfAllSyn","Order","SummOfAllGrantham" )

DataRB = merge(DataRB,dnds, by.x = "Species", by.y = "Species",all=F,no.dups = TRUE,)
names(DataRB) = c("1.Species","2.AverageGrantham","SummOfAllGrantham","3.KnKs","MedianOfAllSyn","MedianOfAllNonsyn","FractionOfSyn","FractionOfNonsyn","Order","5.GenerationLength_d","Gene","4.DnDs" )
DataRB = DataRB[c(1, mixedorder(names(DataRB)[-1]) + 1)];names(DataRB) = c("Species","AverageGrantham","KnKs","DnDs","GenerationLength_d","FractionOfNonsyn","FractionOfSyn","Gene","MedianOfAllNonsyn","MedianOfAllSyn","Order","SummOfAllGrantham" )
Data = Data[is.na(Data$DnDs)==F,]
DataRB = DataRB[is.na(DataRB$DnDs)==F,]
write.table(Data,file = "../../Body/2Derived/CytB_with_ecology.csv",quote = F, row.names = FALSE,sep = '\t')
write.table(DataRB,file = "../../Body/2Derived/CytB_with_ecologyForRedBook.csv",quote = F, row.names = FALSE,sep = '\t')

