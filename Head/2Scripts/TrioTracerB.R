rm(list=ls(all=TRUE))
  
Trios = read.table("../../Body/2Derived/TrioTracerA.PairWiseAnnotation.txt")  

# NumbOfSubst,Synon,NonsAll,Nons1,Nons2,Nons3,Indel,Grantham1,Grantham2,Grantham3

Trios$KnKsSp1Sp2 = unlist(strsplit(as.character(Trios$DistancesAllCodonSubstBetweenSp1AndSp2),','))[2]

Divergence <- function(x) {
  # x = as.character(Trios$DistancesAllCodonSubstBetweenSp1AndSp2)
  Div = as.numeric(unlist(strsplit(x,','))[1])}


KnKs <- function(x) {
  # x = as.character(Trios$DistancesAllCodonSubstBetweenSp1AndSp2)
  Ks = as.numeric(unlist(strsplit(x,','))[2])
  Kn = as.numeric(unlist(strsplit(x,','))[3])
  KnKs = Kn/Ks }

AverageGrantham <- function(x) {
  # x = as.character(Trios$DistancesAllCodonSubstBetweenSp1AndSp2)
  Nons = as.numeric(unlist(strsplit(x,','))[4])
  Grantham = as.numeric(unlist(strsplit(x,','))[8])
  AverageGrantham = Grantham/Nons }

Trios$DivergenceSp1Sp2 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp2),1,Divergence)
Trios$KnKsSp1Sp2 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp2),1,KnKs)
Trios$AverageGranthamSp1Sp2 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp2),1,AverageGrantham)

Trios$DivergenceSp1Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp3),1,Divergence)
Trios$KnKsSp1Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp3),1,KnKs)
Trios$AverageGranthamSp1Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp3),1,AverageGrantham)

Trios$DivergenceSp2Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp2AndSp3),1,Divergence)
Trios$KnKsSp2Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp2AndSp3),1,KnKs)
Trios$AverageGranthamSp2Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp2AndSp3),1,AverageGrantham)

# control tests
DIvergenceVec = c(Trios$DivergenceSp1Sp2,Trios$DivergenceSp1Sp3,Trios$DivergenceSp2Sp3)
KnKsVec = c(Trios$KnKsSp1Sp2,Trios$KnKsSp1Sp3,Trios$KnKsSp2Sp3)
AverageGranthamVec = c(Trios$AverageGranthamSp1Sp2,Trios$AverageGranthamSp1Sp3,Trios$AverageGranthamSp2Sp3)

cor.test(KnKsVec,DIvergenceVec, method = 'spearman') # negative
cor.test(KnKsVec,AverageGranthamVec, method = 'spearman') # negative

cor.test(Trios$KnKsSp1Sp2,Trios$DivergenceSp1Sp2, method = 'spearman') # negative
cor.test(Trios$KnKsSp1Sp3,Trios$AverageGranthamSp1Sp3, method = 'spearman') # negative
cor.test(Trios$KnKsSp2Sp3,Trios$AverageGranthamSp2Sp3, method = 'spearman') # negative

cor.test(Trios$KnKsSp1Sp2,Trios$AverageGranthamSp1Sp2, method = 'spearman') # negative

plot(Trios$KnKsSp1Sp2,Trios$AverageGranthamSp1Sp2)
plot(Trios$KnKsSp1Sp3,Trios$AverageGranthamSp1Sp3)
plot(Trios$KnKsSp2Sp3,Trios$AverageGranthamSp2Sp3)


# find outgroup as the most distant species in the trio: