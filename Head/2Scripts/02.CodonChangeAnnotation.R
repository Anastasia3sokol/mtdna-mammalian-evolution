  rm(list=ls(all=TRUE))
  
  wd = getwd()
  wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
  setwd(wd)
  
  #### libraries
  # if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
  # BiocManager::install("Biostrings")
  
  library(seqinr)
  library(Biostrings)
  library(dplyr)
  library(RGenetics)
  SGC1 <- getGeneticCode("SGC1")  # Vertebrate Mitochondrial code for translate function
  
  Codons = read.table("../../Body/2Derived/PolymorphicPairwiseCodons.txt") #датасет, созданный в скрипте 01,содержит вектор замещений SubstVec
  #Codons <- Codons2[1:100,]
  
  ### functions which take example line as input and return different metrics: 
  # number of all codon substitutions,
  # number of of syn substitutions,
  # number of nonsyn substitutions, 
  # within nonsyn substitutions - average grantham distance between the first and the second AA
  
  #Функция, которая считает общее количество замещений
  TotalDiv <- function(x) {Div = length(unlist(strsplit(x,';')))}
  Codons$TotalDiv <- apply(as.matrix(Codons$SubstVec),1,FUN = TotalDiv) #посчитанные замещения прикрепляются в колонку TotalDiv в датасет Codons
  
  # А дальше считаем дистанции
  data(aaindex)
  names(aaindex) # 544 codes related to different amino-acid properties
  aaindex$GRAR740101 # composition from Grantham 
  aaindex$GRAR740102 # polarity from Grantham 
  aaindex$GRAR740103 # volume from Grantham 
  
  GranthamDistance = data.frame(cbind(aaindex$GRAR740101$I,aaindex$GRAR740102$I,aaindex$GRAR740103$I)) # создается таблица с данными по структуре, полярности и обьему каждой аминокислоты
  names(GranthamDistance)=c("Composition","Polarity","Volume")
  GranthamDistance$Aaa = row.names(GranthamDistance) # создается колонка с трехбуквенными названиями аминокислот
  GranthamDistance$A = a(GranthamDistance$Aaa) # создается колонка с однобуквенными названиями аминокислот
  
  Synon <- 0
  Non_Synon <- 0
  Composition <- ''
  Polarity<- ""
  Volume <- ""
  Grantham <- ""
  
  ###########333
  
  one_line = c()
  
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
      ##j =11
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
        
        SquareOfDiffInComposition =  (GranthamDistance[GranthamDistance$A == A1,]$Composition - GranthamDistance[GranthamDistance$A == A2,]$Composition)^2
        # Composition <<- c(Composition,paste(SquareOfDiffInComposition)) 
        
        SquareOfDiffInPolarity =    (GranthamDistance[GranthamDistance$A == A1,]$Polarity - GranthamDistance[GranthamDistance$A == A2,]$Polarity)^2
        # Polarity <<- c(Polarity,paste(SquareOfDiffInPolarity))
        
        SquareOfDiffInVolume   =    (GranthamDistance[GranthamDistance$A == A1,]$Volume  - GranthamDistance[GranthamDistance$A == A2,]$Volume)^2
        # Volume <<- c(Volume,paste(SquareOfDiffInVolume))
        
        GranthamNew <- (1.833*SquareOfDiffInComposition + 0.1018*SquareOfDiffInPolarity + 0.000399*SquareOfDiffInVolume)^0.5
        Grantham <<-c(Grantham,paste(GranthamNew))
        
        # Output=paste(Synon,Non_Synon,Grantham,sep=';')
        # one_line = rbind(one_line, c(AminoSubs, Grantham))
        
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

write.table(newCodons,file = "../../Body/2Derived/Grantham.csv",quote = F, row.names = FALSE,
              sep = '\t')
  
