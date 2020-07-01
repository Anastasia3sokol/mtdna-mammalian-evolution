rm(list=ls(all=TRUE))

#### libraries
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("Biostrings")

library(seqinr)
library(Biostrings)
SGC1 <- getGeneticCode("SGC1")  # Vertebrate Mitochondrial code for translate function


#### read data:
Codons = read.table("../../Body/2Derived/PolymorphicPairwiseCodons.txt")
dim(Codons)

#  "3:TGC>AGC;7:ACA>ACC;8:CAT>CAC;9:CCC>CCT;14:GCT>GCG;17:GCG>ACG;26:AAC>AGC;369:GCA>ACC;374:AAC>AAT;377:TTA>ATA;380:GCC>GCT;"

### functions which take example line as input and return different metrics: 
# number of all codon substitutions,
# number of of syn substitutions,
# number of nonsyn substitutions, 
# within nonsyn substitutions - average grantham distance between the first and the second AA

#### 1 - count number of all codon changes
# x = Codons$SubstVec[1] 
# "3:TGC>AGC;7:ACA>ACC;8:CAT>CAC;9:CCC>CCT;14:GCT>GCG;17:GCG>ACG;26:AAC>AGC;369:GCA>ACC;374:AAC>AAT;377:TTA>ATA;380:GCC>GCT;"

TotalDiv <- function(x) {Div = length(unlist(strsplit(x,';')))}
Codons$TotalDiv = apply(as.matrix(Codons$SubstVec),1,FUN = TotalDiv)
summary(Codons$TotalDiv)

#### 2 - count number of syn and nons changes
# x = Codons$SubstVec[1] 
# "3:TGC>AGC;7:ACA>ACC;8:CAT>CAC;9:CCC>CCT;14:GCT>GCG;17:GCG>ACG;26:AAC>AGC;369:GCA>ACC;374:AAC>AAT;377:TTA>ATA;380:GCC>GCT;"

NumberSynNons <- function(x) 
{
  VecOfCodons = unlist(strsplit(x,';'));
  for (i in 1:length(VecOfCodons))
  { # i = 1
    CodonSubst = VecOfCodons[i]
    CodonSubst = gsub("(.*)\\:",'',CodonSubst)
    CodonSubst1 = gsub(">(.*)",'',CodonSubst)
    CodonSubst2 = gsub("(.*)>",'',CodonSubst)
    
    Codon1 <- DNAString(CodonSubst1)
    Codon2  <- DNAString(CodonSubst2)
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
      
    
return  paste(number of Syn subst, number of Nons subst, list of all AA subst: C>S,A>V...)
    
  "5;3;C>S,A>V,A>V"
  
  apply  
    
    
    
        
    
  }
  
Div = length(unlist(strsplit(x,';')))
  
}


#### from DIMA: TracerA

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
  # Grantham's distance depends on 3 p


