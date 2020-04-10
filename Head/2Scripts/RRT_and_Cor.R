rm(list=ls(all=TRUE))

wd = getwd()
wd = paste(wd, '/mtdna-mammalian-evolution/Body/1Raw',sep='')
setwd(wd)

GenerationL = read.table("GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GenerationL = GenerationL[ ,c(2,11)] # Take information about name and GenLength
GenerationL$Scientific_name = gsub(' ','_',GenerationL$Scientific_name) 

Trios = read.table('RRT/cytb.threesomes.neighbours2.RRT.txt', header = FALSE, skip = 3)
names(Trios) = c('Family','Ingroup0','Ingroup1','Outgroup','Gene','Kn[0]','Kn[0]','Kn[0]','Kn[0]','Kn[0]','Kn[2]','Kn[2]','Kn[2]',
                 'Kn[2]','Kn[2]','Ks[2]','Ks[2]','Ks[2]','Ks[2]','Ks[2]','Ks[4]','Ks[4]','Ks[4]','Ks[4]','Ks[4]')

Trios1 = merge(Trios1, GenerationL, by.x = 'Ingroup0', by.y = 'Scientific_name' )
Trios1 = merge(Trios1, GenerationL, by.x = 'Ingroup1', by.y = 'Scientific_name' )
Trios1 = merge(Trios, GenerationL, by.x = 'Outgroup', by.y = 'Scientific_name' )
names(Trios1)[c(26,27,28)] =c('GenLenIn0','GenLenIn1','GenLenOut') 

wd = gsub('1Raw','2Derived',wd)
setwd(wd)
write.table(Trios1, 'threesomes_nb2_RRT_wGL.txt')

GenLenDif = Trios1$GenLenIn0 - Trios1$GenLenIn1

cor.test(GenLenDif,Trios1[,9], method = 'spearman')
plot(GenLenDif,Trios1[,9], main="GenLen~K13-K23",
     xlab="GenerationLength", ylab="AminoContrast", pch=19)


cor.test(GenLenDif,Trios1[,24], method = 'spearman')
plot(GenLenDif,Trios1[,24], main="GenLen~K13-K23",
     xlab="GenerationLength", ylab="KsContrast", pch=19)


cor.test(Trios1[,9],Trios1[,24], method = 'spearman')
plot(Trios1[,24],Trios1[,24], main="AC~Ks",
     xlab="AminoContrast", ylab="KsContrast", pch=19)

