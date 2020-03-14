rm(list=ls(all=TRUE))
  
Trios = read.table("../../Body/2Derived/TrioTracerA.PairWiseAnnotation.txt")  

# NumbOfSubst,Synon,NonsAll,Nons1,Nons2,Nons3,Indel,Grantham1,Grantham2,Grantham3

#Trios$KnKsSp1Sp2 = unlist(strsplit(as.character(Trios$DistancesAllCodonSubstBetweenSp1AndSp2),','))[2]

Divergence <- function(x) {
  # x = as.character(Trios$DistancesAllCodonSubstBetweenSp1AndSp2)
  Div = as.numeric(unlist(strsplit(x,','))[1])}

KnKs <- function(x) {
  # x = as.character(Trios$DistancesAllCodonSubstBetweenSp1AndSp2)
  Ks = as.numeric(unlist(strsplit(x,','))[2])
  Kn = as.numeric(unlist(strsplit(x,','))[3])
  KnKs = Kn/Ks }

Ks <- function(x) {
  # x = as.character(Trios$DistancesAllCodonSubstBetweenSp1AndSp2)
  Ks = as.numeric(unlist(strsplit(x,','))[2])   }

AverageGrantham <- function(x) {
  # x = as.character(Trios$DistancesAllCodonSubstBetweenSp1AndSp2)
  Nons = as.numeric(unlist(strsplit(x,','))[4])
  Grantham = as.numeric(unlist(strsplit(x,','))[8])
  AverageGrantham = Grantham/Nons }

Trios$DivergenceSp1Sp2 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp2),1,Divergence)
Trios$KnKsSp1Sp2 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp2),1,KnKs)
Trios$KsSp1Sp2 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp2),1,Ks)
Trios$AverageGranthamSp1Sp2 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp2),1,AverageGrantham)

Trios$DivergenceSp1Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp3),1,Divergence)
Trios$KnKsSp1Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp3),1,KnKs)
Trios$KsSp1Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp3),1,Ks)
Trios$AverageGranthamSp1Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp3),1,AverageGrantham)

Trios$DivergenceSp2Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp2AndSp3),1,Divergence)
Trios$KnKsSp2Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp2AndSp3),1,KnKs)
Trios$KsSp2Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp2AndSp3),1,Ks)
Trios$AverageGranthamSp2Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp2AndSp3),1,AverageGrantham)

### GENETIC EXPECTATIONS IN SIMPLE PAIRWISE ANALYSES:
### high Kn/Ks should be associated with high average Grantham
cor.test(Trios$KnKsSp1Sp2,Trios$AverageGranthamSp1Sp2, method = 'spearman') # - again negative!
cor.test(Trios$KnKsSp1Sp3,Trios$AverageGranthamSp1Sp3, method = 'spearman') # - again negative!
cor.test(Trios$KnKsSp2Sp3,Trios$AverageGranthamSp2Sp3, method = 'spearman') # - again negative!

### outgroup is 1: 2 and 3 are ingroup species. 
# If 2 is more long-lived than 3 => Trios$GenLenSp2/Trios$GenLenSp3 > 1 &  Trios$KnKsSp1Sp2/Trios$KnKsSp1Sp3 > 1
# if 2 is more short lived than 3 => Trios$GenLenSp2/Trios$GenLenSp3 < 1 &  Trios$KnKsSp1Sp2/Trios$KnKsSp1Sp3 < 1
# FUCK, it is opposite

# elephant has higher Kn/Ks: opposite, DONT UNDERSTAND
cor.test(Trios$GenLenSp2/Trios$GenLenSp3,Trios$KnKsSp1Sp2/Trios$KnKsSp1Sp3, method = 'spearman') # out 1
cor.test(Trios$GenLenSp1/Trios$GenLenSp3,Trios$KnKsSp1Sp2/Trios$KnKsSp2Sp3, method = 'spearman') # out 2
cor.test(Trios$GenLenSp1/Trios$GenLenSp2,Trios$KnKsSp1Sp3/Trios$KnKsSp2Sp3, method = 'spearman') # out 3

# elephant has shorter Ks branch: nothing
cor.test(Trios$GenLenSp2/Trios$GenLenSp3,Trios$KsSp1Sp2/Trios$KsSp1Sp3, method = 'spearman') # out 1
cor.test(Trios$GenLenSp1/Trios$GenLenSp3,Trios$KsSp1Sp2/Trios$KsSp2Sp3, method = 'spearman') # out 2
cor.test(Trios$GenLenSp1/Trios$GenLenSp2,Trios$KsSp1Sp3/Trios$KsSp2Sp3, method = 'spearman') # out 3

# elephant has shorter branch: something
cor.test(Trios$GenLenSp2/Trios$GenLenSp3,Trios$DivergenceSp1Sp2/Trios$DivergenceSp1Sp3, method = 'spearman') # out 1
cor.test(Trios$GenLenSp1/Trios$GenLenSp3,Trios$DivergenceSp1Sp2/Trios$DivergenceSp2Sp3, method = 'spearman') # out 2
cor.test(Trios$GenLenSp1/Trios$GenLenSp2,Trios$DivergenceSp1Sp3/Trios$DivergenceSp2Sp3, method = 'spearman') # out 3

# elephant has higher Grhantham: something
cor.test(Trios$GenLenSp2/Trios$GenLenSp3,Trios$AverageGranthamSp1Sp2/Trios$AverageGranthamSp1Sp3, method = 'spearman') # out 1
cor.test(Trios$GenLenSp1/Trios$GenLenSp3,Trios$AverageGranthamSp1Sp2/Trios$AverageGranthamSp2Sp3, method = 'spearman') # out 2
cor.test(Trios$GenLenSp1/Trios$GenLenSp2,Trios$AverageGranthamSp1Sp3/Trios$AverageGranthamSp2Sp3, method = 'spearman') # out 3




