rm(list=ls(all=TRUE))
  
  # if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
  # BiocManager::install("Biostrings")
  
  library(seqinr)

  library(Biostrings)
  SGC1 <- getGeneticCode("SGC1")  # Vertebrate Mitochondrial code for translate function
  
  #### aaindex from Seqinr package #  List of 544 physicochemical and biological properties for the 20 amino-acids
  data(aaindex)
  names(aaindex) # 544 codes related to different amino-acid properties!!! A lot to play!!!
  aaindex$GRAR740101 # composition from Grantham  = compared with numbers in Science = the same!!
  aaindex$GRAR740102 # polarity from Grantham 
  aaindex$GRAR740103 # volume from Grantham 
  
  GranthamDistance = data.frame(cbind(aaindex$GRAR740101$I,aaindex$GRAR740102$I,aaindex$GRAR740103$I))
  names(GranthamDistance)=c("Composition","Polarity","Volume")
  GranthamDistance$Aaa = row.names(GranthamDistance)
  GranthamDistance$A = a(GranthamDistance$Aaa)
  str(GranthamDistance)
  
  Trios = read.table("../../Body/2Derived/DifferenceBetweenSpeciesFamilies1CytB.txt", header = TRUE)
  
  AllUniqueSubstTogether = c(as.character(Trios$AllCodonSubstBetweenSp1AndSp2),as.character(Trios$AllCodonSubstBetweenSp1AndSp3),as.character(Trios$AllCodonSubstBetweenSp2AndSp3))
  length(AllUniqueSubstTogether) # 131481 / 3 = 43827
  AllUniqueSubstTogether = unique(AllUniqueSubstTogether)
  length(AllUniqueSubstTogether) # 3587!!! che - che
  
  # for (step in 1:length(AllUniqueSubstTogether))
  #{ # step = 1
  
  #temp = as.character(AllUniqueSubstTogether[step])
  
  ComparisonOfTwoSpecies <- function(x) {
  #x = temp
  Substitutions = unlist(strsplit(x,';')) # 430:GCC>GCT;
  NumbOfSubst = length(Substitutions)
  Nons1=0
  Nons2=0
  Nons3=0
  NonsAll = 0
  Indel = 0
  Synon = 0 
  Grantham1 = 0
  Grantham2 = 0
  Grantham3 = 0
  
  for (i in 1:length(Substitutions))
    { # i = 1
    TwoCodons = Substitutions[i] # TwoCodons = '430:GCC>GCT'
    TwoCodons = gsub(".*:",'',TwoCodons) # GCC>GCT
    TwoCodons = unlist(strsplit(TwoCodons,'>'))
    Codon1 <- DNAString(TwoCodons[1])
    Codon2  <- DNAString(TwoCodons[2])
    Codon1.Character = as.character(Codon1)
    Codon2.Character = as.character(Codon2)
    if  (grepl('-',Codon1.Character) | grepl('-',Codon2.Character)) {Indel=Indel+1; break}
    if  (!grepl('-',Codon1.Character) & !grepl('-',Codon2.Character)) 
    {
      A1 = as.character(Biostrings::translate(Codon1, genetic.code=SGC1))
      A2 = as.character(Biostrings::translate(Codon2, genetic.code=SGC1))
      if (A1 == A2) {Synon=Synon+1}
    }
    if  (A1 != A2) 
    { # estimate Grantham distance 
      # from https://en.wikipedia.org/wiki/Amino_acid_replacement:
      # Grantham's distance depends on 3 properties: composition, polarity and molecular volume.[4]
      # Distance difference D for each pair of amino acid i and j is calculated as: {\displaystyle D_{ij}=[\alpha (c_{i}-c_{j})^{2}+\beta (p_{i}-p_{j})^{2}+\gamma (v_{i}-v_{j})^{2}]}{\displaystyle D_{ij}=[\alpha (c_{i}-c_{j})^{2}+\beta (p_{i}-p_{j})^{2}+\gamma (v_{i}-v_{j})^{2}]}
      # where c = composition, p = polarity, and v = molecular volume; and are constants of squares of the inverses of the mean distance for each property, respectively equal to 1.833, 0.1018, 0.000399. According to Grantham's distance, most similar amino acids are leucine and isoleucine and the most distant are cysteine and tryptophan.
      
      diff = 0
      if (unlist(strsplit(Codon1.Character,""))[1] != unlist(strsplit(Codon2.Character, ""))[1]) {diff = diff+1}
      if (unlist(strsplit(Codon1.Character,""))[2] != unlist(strsplit(Codon2.Character, ""))[2]) {diff = diff+1}
      if (unlist(strsplit(Codon1.Character,""))[3] != unlist(strsplit(Codon2.Character, ""))[3]) {diff = diff+1}
      if (diff == 1) {Nons1 = Nons1+1}
      if (diff == 2) {Nons2 = Nons2+1}
      if (diff == 3) {Nons3 = Nons3+1}
      NonsAll = NonsAll + 1
      
  #    SquareOfDiffInComposition =  1.833*(GranthamDistance[GranthamDistance$A == A1,]$Composition - GranthamDistance[GranthamDistance$A == A2,]$Composition)^2
  #    SquareOfDiffInPolarity =    0.1018*(GranthamDistance[GranthamDistance$A == A1,]$Polarity - GranthamDistance[GranthamDistance$A == A2,]$Polarity)^2
  #    SquareOfDiffInVolume   =    0.000399*(GranthamDistance[GranthamDistance$A == A1,]$Volume  - GranthamDistance[GranthamDistance$A == A2,]$Volume)^2
       SquareOfDiffInComposition =  (GranthamDistance[GranthamDistance$A == A1,]$Composition - GranthamDistance[GranthamDistance$A == A2,]$Composition)^2
       SquareOfDiffInPolarity =    (GranthamDistance[GranthamDistance$A == A1,]$Polarity - GranthamDistance[GranthamDistance$A == A2,]$Polarity)^2
       SquareOfDiffInVolume   =    (GranthamDistance[GranthamDistance$A == A1,]$Volume  - GranthamDistance[GranthamDistance$A == A2,]$Volume)^2
      
          GranthamNew = (SquareOfDiffInComposition + SquareOfDiffInPolarity + SquareOfDiffInVolume)^0.5
      if (diff == 1) {Grantham1 = Grantham1+GranthamNew}
      if (diff == 2) {Grantham2 = Grantham2+GranthamNew}
      if (diff == 3) {Grantham3 = Grantham3+GranthamNew}
      
     }
   }  
  Output=paste(NumbOfSubst,Synon,NonsAll,Nons1,Nons2,Nons3,Indel,Grantham1,Grantham2,Grantham3,sep=',')
  # Trios$DistanceSp1AndSp2[step] = Output
  }
  
DistancesForAllUniqueSubstTogether = apply(as.matrix(AllUniqueSubstTogether),1,ComparisonOfTwoSpecies)
ResDataFrame = data.frame(AllUniqueSubstTogether,DistancesForAllUniqueSubstTogether)

# 1 
names(ResDataFrame) = c('AllUniqueSubstTogether','DistancesAllCodonSubstBetweenSp1AndSp2')
Trios = merge(Trios,ResDataFrame, by.x = 'AllCodonSubstBetweenSp1AndSp2', by.y = 'AllUniqueSubstTogether')

names(ResDataFrame) = c('AllUniqueSubstTogether','DistancesAllCodonSubstBetweenSp1AndSp3')
Trios = merge(Trios,ResDataFrame, by.x = 'AllCodonSubstBetweenSp1AndSp3', by.y = 'AllUniqueSubstTogether')

names(ResDataFrame) = c('AllUniqueSubstTogether','DistancesAllCodonSubstBetweenSp2AndSp3')
Trios = merge(Trios,ResDataFrame, by.x = 'AllCodonSubstBetweenSp2AndSp3', by.y = 'AllUniqueSubstTogether')

#Trios$DistanceSp1AndSp2 = apply(as.matrix(Trios$AllCodonSubstBetweenSp1AndSp2),1,ComparisonOfTwoSpecies)
#Trios$DistanceSp1AndSp3 = apply(as.matrix(Trios$AllCodonSubstBetweenSp1AndSp3),1,ComparisonOfTwoSpecies)
#Trios$DistanceSp2AndSp3 = apply(as.matrix(Trios$AllCodonSubstBetweenSp2AndSp3),1,ComparisonOfTwoSpecies)

Trios$AllCodonSubstBetweenSp1AndSp2<-NULL
Trios$AllCodonSubstBetweenSp1AndSp3<-NULL
Trios$AllCodonSubstBetweenSp2AndSp3<-NULL

write.table(Trios,"../../Body/2Derived/TrioTracerA.PairWiseAnnotation.txt")
