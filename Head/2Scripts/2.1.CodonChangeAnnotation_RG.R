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

Codons = read.table("../../Body/2Derived/PolymorphicPairwiseCodons.txt") #датасет, созданный в скрипте 01,содержит вектор замещений SubstVec
#Codons <- Codons2[1:10,]
Codons_CytB = Codons[Codons$Gene == 'CytB',]#42975
Codons_ATP6 = Codons[Codons$Gene == 'ATP6',]#8415
Codons_ATP8 = Codons[Codons$Gene == 'ATP8',]#1890
Codons_ND1 = Codons[Codons$Gene == 'ND1',]#5355
Codons_ND2 = Codons[Codons$Gene == 'ND2',]#13770
Codons_ND3 = Codons[Codons$Gene == 'ND3',]#2745
Codons_ND4 = Codons[Codons$Gene == 'ND4',]#4500
Codons_ND4L = Codons[Codons$Gene == 'ND4L',]#1800
Codons_ND5 = Codons[Codons$Gene == 'ND5',]#4095
Codons_ND6 = Codons[Codons$Gene == 'ND6',]#225
Codons_COX1 = Codons[Codons$Gene == 'COX1',]#3780
Codons_COX2 = Codons[Codons$Gene == 'COX2',]#3420
Codons_COX3 = Codons[Codons$Gene == 'COX3',]#2880

####################
#Codons = Codons[Codons$Gene == 'COX3',]
####################

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

newCodons = Codons[, -7]
write.table(newCodons,file = "../../Body/2Derived/Grantham.csv",quote = F, row.names = FALSE,sep = '\t')
