rm(list=ls(all=TRUE))

##### READ AND MERGE GENERATION LENGTH AND TRIO DATA
GenerationL = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GenerationL = GenerationL[ ,c(2,11)] # Take information about name and GenLength
GenerationL$Scientific_name = gsub(' ','_',GenerationL$Scientific_name) 

Trios = read.table('../../Body/1Raw/RRT/cytb.threesomes.neighbours4.RRT.txt', header = FALSE, skip = 3) # neighbours2/3/4
names(Trios) = c('Family','Ingroup1','Ingroup2','Outgroup','Gene','X1','X2','X3','X4','X5','X6','X7','X8',
                 'X9','X10','X11','X12','X13','X14','X15','X16','X17','X18','X19','X20')
# Column names from table with 16 columns (but in Kostya's table there are 20 columns!&): 
# Proc.Nati.Acad.Sci.USA Vol.82, pp.1741-1745, March1985, Table 1
colnames(Trios)[colnames(Trios) == 'X4'] <- 'NonsNondegIn1MinusIn2'
colnames(Trios)[colnames(Trios) == 'X9'] <- 'Nons2FDegenIn1MinusIn2'
summary(Trios$NonsNondegIn1MinusIn2)
summary(Trios$Nons2FDegenIn1MinusIn2)
colnames(Trios)[colnames(Trios) == 'X14'] <-'Syn2FDegenIn1MinusIn2'
colnames(Trios)[colnames(Trios) == 'X19'] <- 'Syn4FDegenIn1MinusIn2'
summary(Trios$Syn2FDegenIn1MinusIn2)
summary(Trios$Syn4FDegenIn1MinusIn2)

Trios1 = merge(Trios, GenerationL, by.x = 'Outgroup', by.y = 'Scientific_name')   # GenLenOut
colnames(Trios1)[colnames(Trios1) == 'GenerationLength_d'] <- 'GenLenOut'

Trios1 = merge(Trios1, GenerationL, by.x = 'Ingroup1', by.y = 'Scientific_name' ) # GenLenIn1
colnames(Trios1)[colnames(Trios1) == 'GenerationLength_d'] <- 'GenLenIn1'

Trios1 = merge(Trios1, GenerationL, by.x = 'Ingroup2', by.y = 'Scientific_name' ) # GenLenIn2
colnames(Trios1)[colnames(Trios1) == 'GenerationLength_d'] <- 'GenLenIn2'
names(Trios1)

##### RUN RELATIVE RATIO TEST ANALYSES

Trios1$GenLenIn1MinusIn2 = Trios1$GenLenIn1-Trios1$GenLenIn2
summary(Trios1$GenLenIn1MinusIn2)


## CONTROL CHECK:
cor.test(Trios1$NonsNondegIn1MinusIn2,Trios1$Nons2FDegenIn1MinusIn2, method = 'spearman') # positive and good
cor.test(Trios1$Syn2FDegenIn1MinusIn2,Trios1$Syn4FDegenIn1MinusIn2, method = 'spearman') # positive and good
cor.test(Trios1$NonsNondegIn1MinusIn2,Trios1$Syn4FDegenIn1MinusIn2, method = 'spearman') # positive but not so strong

Trios1$NonsToSyn = Trios1$NonsNondegIn1MinusIn2 / Trios1$Syn4FDegenIn1MinusIn2
summary(Trios1$NonsToSyn); nrow(Trios1)
Trios1 = Trios1[abs(Trios1$NonsToSyn) < Inf,]; nrow(Trios1)
cor.test(Trios1$NonsToSyn,Trios1$Syn4FDegenIn1MinusIn2, method = 'spearman') # nothing
cor.test(Trios1$NonsToSyn,Trios1$NonsNondegIn1MinusIn2, method = 'spearman') # a bit negative

## the longer the GenLength the higher the Nons => elephants vs mice have more Kn (NonsNondegIn1MinusIn2)
cor.test(Trios1$GenLenIn1MinusIn2,Trios1$NonsNondegIn1MinusIn2, method = 'spearman')
cor.test(Trios1$GenLenIn1MinusIn2,Trios1$Nons2FDegenIn1MinusIn2, method = 'spearman')
plot(Trios1$GenLenIn1MinusIn2,Trios1$NonsNondegIn1MinusIn2)

## the longer the GenLength the lower the Syn => elephants vs mice have less Ks (Syn4FDegenIn1MinusIn2)
cor.test(Trios1$GenLenIn1MinusIn2,Trios1$Syn4FDegenIn1MinusIn2, method = 'spearman')
cor.test(Trios1$GenLenIn1MinusIn2,Trios1$Syn2FDegenIn1MinusIn2, method = 'spearman')
plot(Trios1$GenLenIn1MinusIn2,Trios1$Syn4FDegenIn1MinusIn2)

## 
cor.test(Trios1$GenLenIn1MinusIn2,Trios1$NonsToSyn, method = 'spearman')
plot(Trios1$GenLenIn1MinusIn2,Trios1$Syn4FDegenIn1MinusIn2)


##### by families:
table(Trios1$Family)

cor.test(Trios1[Trios1$Family == 'Cercopithecidae',]$GenLenDif,Trios1[Trios1$Family == 'Cercopithecidae',]$X4, method = 'spearman')
cor.test(Trios1[Trios1$Family == 'Cercopithecidae',]$GenLenDif,Trios1[Trios1$Family == 'Cercopithecidae',]$X9, method = 'spearman')
plot(Trios1[Trios1$Family == 'Cercopithecidae',]$GenLenDif,Trios1[Trios1$Family == 'Cercopithecidae',]$X4)

cor.test(Trios1[Trios1$Family == 'Cercopithecidae',]$GenLenDif,Trios1[Trios1$Family == 'Cercopithecidae',]$X19, method = 'spearman')
cor.test(Trios1[Trios1$Family == 'Cercopithecidae',]$GenLenDif,Trios1[Trios1$Family == 'Cercopithecidae',]$X14, method = 'spearman')
plot(Trios1$GenLenDif,Trios1$X19)


