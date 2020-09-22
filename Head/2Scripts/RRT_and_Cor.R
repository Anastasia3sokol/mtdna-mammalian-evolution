rm(list=ls(all=TRUE))

#wd = getwd()
#wd = paste(wd, '/mtdna-mammalian-evolution/Body/1Raw',sep='')
#setwd(wd)

GenerationL = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
# GenerationL = read.table("GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GenerationL = GenerationL[ ,c(2,11)] # Take information about name and GenLength
GenerationL$Scientific_name = gsub(' ','_',GenerationL$Scientific_name) 

Trios = read.table('../../Body/1Raw/RRT/cytb.threesomes.neighbours4.RRT.txt', header = FALSE, skip = 3) # neighbours2/3/4
# Trios = read.table('RRT/cytb.threesomes.neighbours2.RRT.txt', header = FALSE, skip = 3)
names(Trios) = c('Family','Ingroup0','Ingroup1','Outgroup','Gene','X1','X2','X3','X4','X5','X6','X7','X8',
                 'X9','X10','X11','X12','X13','X14','X15','X16','X17','X18','X19','X20')

Trios1 = merge(Trios, GenerationL, by.x = 'Outgroup', by.y = 'Scientific_name')
Trios1 = merge(Trios1, GenerationL, by.x = 'Ingroup0', by.y = 'Scientific_name' ) # in1
Trios1 = merge(Trios1, GenerationL, by.x = 'Ingroup1', by.y = 'Scientific_name' ) # in2

names(Trios1)
##### names(Trios1)[c(26,27,28)] =c('GenLenOut','GenLenIn1','GenLenIn2') #### HERE IS ERROR, NOT NORMAL INDEX FOR GENLEN 
names(Trios1)[c(26,27,28)] =c('GenLenOut','GenLenIn2','GenLenIn1')
#wd = gsub('1Raw','2Derived',wd)
#setwd(wd)
#write.table(Trios1, 'threesomes_nb2_RRT_wGL.txt')

Trios1$GenLenDif = Trios1$GenLenIn1-Trios1$GenLenIn2

## elephants have more Kn

cor.test(Trios1$GenLenDif,Trios1$X4, method = 'spearman')
cor.test(Trios1$GenLenDif,Trios1$X9, method = 'spearman')
plot(Trios1$GenLenDif,Trios1$X4)

## elephants have less Ks

cor.test(Trios1$GenLenDif,Trios1$X19, method = 'spearman')
cor.test(Trios1$GenLenDif,Trios1$X14, method = 'spearman')
plot(Trios1$GenLenDif,Trios1$X19)


##### by families:
table(Trios1$Family)

cor.test(Trios1[Trios1$Family == 'Cercopithecidae',]$GenLenDif,Trios1[Trios1$Family == 'Cercopithecidae',]$X4, method = 'spearman')
cor.test(Trios1[Trios1$Family == 'Cercopithecidae',]$GenLenDif,Trios1[Trios1$Family == 'Cercopithecidae',]$X9, method = 'spearman')
plot(Trios1[Trios1$Family == 'Cercopithecidae',]$GenLenDif,Trios1[Trios1$Family == 'Cercopithecidae',]$X4)

cor.test(Trios1[Trios1$Family == 'Cercopithecidae',]$GenLenDif,Trios1[Trios1$Family == 'Cercopithecidae',]$X19, method = 'spearman')
cor.test(Trios1[Trios1$Family == 'Cercopithecidae',]$GenLenDif,Trios1[Trios1$Family == 'Cercopithecidae',]$X14, method = 'spearman')
plot(Trios1$GenLenDif,Trios1$X19)


