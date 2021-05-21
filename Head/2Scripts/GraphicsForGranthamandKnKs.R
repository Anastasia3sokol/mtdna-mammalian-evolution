rm(list=ls(all=TRUE))

library(ggplot2)
library(nortest)

gr = read.table('/home/diliushchenko/mtdna-mammalian-evolution/Body/2Derived/TrioTracerA.PairWiseAnnotation.txt')

countavggrall = function(x)
{
  NonsAll = as.numeric(unlist(strsplit(x,','))[3])
  Grantham1 = as.numeric(unlist(strsplit(x,','))[8])
  Grantham2 = as.numeric(unlist(strsplit(x,','))[9])
  Grantham3 = as.numeric(unlist(strsplit(x,','))[10])
  AverageGrantham = (Grantham1+Grantham2+Grantham3)/NonsAll
  return(AverageGrantham)
}

countavggr1 = function(x)
{
  Nons1 = as.numeric(unlist(strsplit(x,','))[4])
  Grantham1 = as.numeric(unlist(strsplit(x,','))[8])
  AverageGrantham = Grantham1/Nons1
  return(AverageGrantham)
}

gr$avggrallIn0In1 =  apply(as.matrix(gr$DistanceIn0In1),1,countavggrall)
gr$avggrallIn0Out =  apply(as.matrix(gr$DistanceIn0Out),1,countavggrall)
gr$avggrallIn1Out =  apply(as.matrix(gr$DistanceIn1Out),1,countavggrall)


gr$avggr1In0In1 =  apply(as.matrix(gr$DistanceIn0In1),1,countavggr1)
gr$avggr1In0Out =  apply(as.matrix(gr$DistanceIn0Out),1,countavggr1)
gr$avggr1In1Out =  apply(as.matrix(gr$DistanceIn1Out),1,countavggr1)


gr = gr[,-c(5,6,7)]


genlen = read.table('/home/diliushchenko/mtdna-mammalian-evolution/Body/1Raw/GenerationLenghtforMammals.xlsx.txt',header = TRUE, sep = '\t')
genlen = genlen[ ,c(2,11)] # Take information about name and GenLength
genlen$Scientific_name = gsub(' ','_',genlen$Scientific_name) 


grf = merge(gr, genlen, by.x = 'Outgroup', by.y = 'Scientific_name')   # GenLenOut
colnames(grf)[colnames(grf) == 'GenerationLength_d'] <- 'GenLenOut'

grf = merge(grf, genlen, by.x = 'Ingroup0', by.y = 'Scientific_name' ) # GenLenIn0
colnames(grf)[colnames(grf) == 'GenerationLength_d'] <- 'GenLenIn0'

grf = merge(grf, genlen, by.x = 'Ingroup1', by.y = 'Scientific_name' ) # GenLenIn1
colnames(grf)[colnames(grf) == 'GenerationLength_d'] <- 'GenLenIn1'
names(grf)


#### calculate different values

grf$GenLenDif = grf$GenLenIn0/grf$GenLenIn1
grf$GenLenMinus = grf$GenLenIn0 - grf$GenLenIn1
grf$GenLenPlusOut = (grf$GenLenIn0 + grf$GenLenOut) / (grf$GenLenIn1 + grf$GenLenOut)

grf$distgallIn0In1Out = grf$avggrallIn0Out/grf$avggrallIn1Out
grf$distg1In0In1Out = grf$avggr1In0Out/grf$avggr1In1Out


grf$minusdistgallIn0In1Out = grf$avggrallIn0Out-grf$avggrallIn1Out
grf$minusdistg1In0In1Out = grf$avggr1In0Out-grf$avggr1In1Out

grf = grf[!is.na(grf$distg1In0In1Out),]


## CONTROL CHECK:
lillie.test(grf$distgallIn0In1Out)
lillie.test(grf$GenLenMinus) ## both have p-value < 0.05


cor.test(grf$GenLenMinus, grf$distgallIn0In1Out, method = 'spearman') # positive and good
cor.test(grf$GenLenMinus, grf$distg1In0In1Out, method = 'spearman') # positive and good

cor.test(grf$GenLenMinus, grf$minusdistgallIn0In1Out, method = 'spearman') # positive and good
cor.test(grf$GenLenMinus, grf$minusdistg1In0In1Out, method = 'spearman')# positive and good

cor.test(grf$GenLenDif, grf$distgallIn0In1Out, method = 'spearman')# positive and good
cor.test(grf$GenLenDif, grf$distg1In0In1Out, method = 'spearman')# positive and good

cor.test(grf$GenLenDif, grf$minusdistgallIn0In1Out, method = 'spearman')# positive and good
cor.test(grf$GenLenDif, grf$minusdistg1In0In1Out, method = 'spearman')# positive and good

cor.test(grf$GenLenPlusOut, grf$distg1In0In1Out, method = 'spearman') # the same corr, genlen for out not important
summary(lm(formula = distgallIn0In1Out ~  GenLenDif + GenLenOut, data = grf)) # from this result GenLenOut not important for Grandham

cor.test(grf$GenLenDif, grf$avggrallIn0In1, method = 'spearman')# rho 0.76
cor.test(grf$GenLenDif, grf$avggr1In0In1, method = 'spearman')# rho 0.0323 , cor just 1 nucleotide Grantham not enough 

##### Draw plots

ggplot(data = grf, aes(x = GenLenMinus, y = distgallIn0In1Out))+
  geom_bin2d()+
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95, col='red')

ggplot(data = grf, aes(x = GenLenMinus, y = distg1In0In1Out))+
  geom_bin2d()+
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95, col='red')



#### make an with median and quant
mm = median(grf$GenLenMinus)

mmm = function(x)
{
  if (x < mm){return(0)}
  else{return(1)}
}

grf$median = apply(as.matrix(grf$GenLenMinus),1,mmm)

ggplot(data= grf, aes(x = GenLenMinus, y = distg1In0In1Out, group = median))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = quantile(grf$distg1In0In1Out, c(0.1, 0.9)))

d = quantile(grf$GenLenMinus, names = FALSE)
quant = function(x)
{
  a = 0
  if (x < d[2]) {a = 1}
  if (x >= d[2] & x < d[3]) {a = 2}
  if (x >= d[3] & x < d[4]){a = 3}
  if (x >= d[4]){a = 4}
  return(a)
}

grf$quant = apply(as.matrix(grf$GenLenMinus),1,quant)

ggplot(data= grf, aes(x = GenLenMinus, y = distg1In0In1Out, group = quant))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = quantile(grf$distg1In0In1Out, c(0.1, 0.9)))


wilcox.test(grf$distg1In0In1Out~grf$median, paired = FALSE) ### p-value < 0.05
wilcox.test(grf$distgallIn0In1Out~grf$median, paired = FALSE)

summary(aov(distgallIn0In1Out ~ quant, data = grf))
summary(aov(distg1In0In1Out ~ quant, data = grf))


aggregate(grf$distg1In0In1Out, by = list(grf$quant), FUN = mean)



#####Start read KnKs

KnKs = read.table('/home/diliushchenko/mtdna-mammalian-evolution/Body/1Raw/RRT/cytb.threesomes.neighbours5.RRT.txt', header = FALSE, skip = 3)
names(KnKs) = c('Family','Ingroup1','Ingroup2','Outgroup','Gene','X1','X2','X3','X4','X5','X6','X7','X8',
                'X9','X10','X11','X12','X13','X14','X15','X16','X17','X18','X19','X20')

colnames(KnKs)[colnames(KnKs) == 'X4'] <- 'NonsNondegIn1MinusIn2'
colnames(KnKs)[colnames(KnKs) == 'X9'] <- 'Nons2FDegenIn1MinusIn2'
summary(KnKs$NonsNondegIn1MinusIn2)
summary(KnKs$Nons2FDegenIn1MinusIn2)

colnames(KnKs)[colnames(KnKs) == 'X14'] <-'Syn2FDegenIn1MinusIn2'
colnames(KnKs)[colnames(KnKs) == 'X19'] <- 'Syn4FDegenIn1MinusIn2'
summary(KnKs$Syn2FDegenIn1MinusIn2)
summary(KnKs$Syn4FDegenIn1MinusIn2)

FinalKnKs = merge(KnKs, genlen, by.x = 'Outgroup', by.y = 'Scientific_name')   # GenLenOut
colnames(FinalKnKs)[colnames(FinalKnKs) == 'GenerationLength_d'] <- 'GenLenOut'

FinalKnKs = merge(FinalKnKs, genlen, by.x = 'Ingroup1', by.y = 'Scientific_name' ) # GenLenIn1
colnames(FinalKnKs)[colnames(FinalKnKs) == 'GenerationLength_d'] <- 'GenLenIn0'

FinalKnKs = merge(FinalKnKs, genlen, by.x = 'Ingroup2', by.y = 'Scientific_name' ) # GenLenIn2
colnames(FinalKnKs)[colnames(FinalKnKs) == 'GenerationLength_d'] <- 'GenLenIn1'

FinalKnKs$GenLenIn0MinusIn1 = FinalKnKs$GenLenIn0-FinalKnKs$GenLenIn1
summary(FinalKnKs$GenLenIn0MinusIn1)

#### CONTROL CHECK

lillie.test(FinalKnKs$NonsNondegIn1MinusIn2)
lillie.test(FinalKnKs$GenLenIn0MinusIn1)

## the longer the GenLength the higher the Nons => elephants vs mice have more Kn (NonsNondegIn1MinusIn2)
cor.test(FinalKnKs$GenLenIn0MinusIn1,FinalKnKs$NonsNondegIn1MinusIn2, method = 'spearman')
cor.test(FinalKnKs$GenLenIn0MinusIn1,FinalKnKs$Nons2FDegenIn1MinusIn2, method = 'spearman')
plot(FinalKnKs$GenLenIn0MinusIn1,FinalKnKs$NonsNondegIn1MinusIn2)

## the longer the GenLength the lower the Syn => elephants vs mice have less Ks (Syn4FDegenIn1MinusIn2)
cor.test(FinalKnKs$GenLenIn0MinusIn1,FinalKnKs$Syn4FDegenIn1MinusIn2, method = 'spearman')
cor.test(FinalKnKs$GenLenIn0MinusIn1,FinalKnKs$Syn2FDegenIn1MinusIn2, method = 'spearman')
plot(FinalKnKs$GenLenIn0MinusIn1,FinalKnKs$Syn4FDegenIn1MinusIn2)



ggplot(data = FinalKnKs, aes(x = GenLenIn0MinusIn1, y = NonsNondegIn1MinusIn2))+
  geom_bin2d()+
  geom_smooth(method = lm)

ggplot(data = FinalKnKs, aes(x = GenLenIn0MinusIn1, y = Syn4FDegenIn1MinusIn2))+
  geom_bin2d()+
  geom_smooth(method = lm)




##### TO DO ###
# 1. Ask Konstantin to calculate diff in KnKs Between In0-Out and In1-Out, in older table we saw just difference between Ingroups
