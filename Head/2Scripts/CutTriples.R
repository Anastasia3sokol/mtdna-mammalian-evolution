rm(list=ls(all=TRUE))

library(ggplot2)


Trios = read.table('../../Body/2Derived/TrioTracerA.PairWiseAnnotation.txt') #### read or triples

TriosTree = read.table('../../Body/1Raw/RRT/cytb.threesomes.neighbours4.RRT.txt', header = FALSE, skip = 3) ### Using tree triples
names(TriosTree) = c('Family','Ingroup0','Ingroup1','Outgroup','Gene','X1','X2','X3','X4','X5','X6','X7','X8',
                 'X9','X10','X11','X12','X13','X14','X15','X16','X17','X18','X19','X20')


GenerationL = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t') ### read Ecology for tree
GenerationL = GenerationL[ ,c(2,11)] # Take information about name and GenLength
GenerationL$Scientific_name = gsub(' ','_',GenerationL$Scientific_name) 

TriosTree1 = merge(TriosTree, GenerationL, by.x = 'Outgroup', by.y = 'Scientific_name')
TriosTree1 = merge(TriosTree1, GenerationL, by.x = 'Ingroup0', by.y = 'Scientific_name' ) # in1
TriosTree1 = merge(TriosTree1, GenerationL, by.x = 'Ingroup1', by.y = 'Scientific_name' ) # in2

names(TriosTree1)
names(TriosTree1)[c(26,27,28)] =c('GenLenOut','GenLenIn2','GenLenIn1') 

TriosTree1$GenLenDif = TriosTree1$GenLenIn1-TriosTree1$GenLenIn2


Kn <- function(x) {
  # x = as.character(Trios$DistancesAllCodonSubstBetweenSp1AndSp2)
  Kn = as.numeric(unlist(strsplit(x,','))[3])
}

Ks <- function(x) {
  # x = as.character(Trios$DistancesAllCodonSubstBetweenSp1AndSp2)
  Ks = as.numeric(unlist(strsplit(x,','))[2])   }


Trios$KsSp1Sp2 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp2),1,Ks)
Trios$KnSp1Sp2 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp2),1,Kn)


Trios$KsSp1Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp3),1,Ks)
Trios$KnSp1Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp1AndSp3),1,Kn)


Trios$KsSp2Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp2AndSp3),1,Ks)
Trios$KnSp2Sp3 = apply(as.matrix(Trios$DistancesAllCodonSubstBetweenSp2AndSp3),1,Kn)

Trios1 = data.frame()
for (i in 1:nrow(Trios)){
  # i = 1
  if (Trios$KnSp1Sp2[i] <= Trios$KnSp1Sp3[i] & Trios$KnSp1Sp2[i] <= Trios$KnSp2Sp3[i])
    {
    Trios1 = rbind(Trios[i,],Trios1) 
    } 
}

Trios1$DifGenLenSp1Sp2 = Trios1$GenLenSp1 - Trios1$GenLenSp2

cor.test(Trios1$DifGenLenSp1Sp2,Trios1$KsSp1Sp2, method = 'spearman')
cor.test(Trios1$DifGenLenSp1Sp2,Trios1$Contrast, method = 'spearman')

cor.test(TriosTree1$GenLenDif,TriosTree1$X4, method = 'spearman')
cor.test(TriosTree1$GenLenDif,TriosTree1$X19, method = 'spearman')


nrow(Trios1[Trios1$Species1 == 'Acomys_cahirinus',]) ### 37 times 1 animal

nrow(TriosTree1[TriosTree1$Ingroup1 == 'Acomys_cahirinus',]) ### 5 times 1 the same animal

Trios1$ContrastInKs = Trios1$KsSp1Sp3 - Trios1$KsSp2Sp3 
Trios1$ContrastInKn = Trios1$KnSp1Sp3 - Trios1$KnSp2Sp3



ggplot(data = Trios1, aes(x = DifGenLenSp1Sp2, y = ContrastInKn))+
  geom_bin2d()+
  scale_fill_continuous(type = "viridis")+
  stat_smooth(method="lm", se=FALSE)+
  theme_bw()


ggplot(data = Trios2, aes(x = DifGenLenSp1Sp2, y = ContrastInKs))+
  geom_bin2d()+
  scale_fill_continuous(type = "viridis")+
  theme_bw()
  
