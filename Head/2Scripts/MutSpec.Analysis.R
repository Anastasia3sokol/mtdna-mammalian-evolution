rm(list=ls(all=TRUE))

library(seqinr)
library(dplyr)

mut = read.table('../../Body/3Results/MutSpec.Final.WithGenLen.txt', header = TRUE)
mut$TrioId = paste(mut$Ingroup1,mut$Ingroup0,mut$Outgroup,sep='.')
nrow(mut) # 37108/9277 = 4 => each trio has exactly 4 mutations!!!!!!!!!!!!????????????????
length(unique(mut$TrioId)) # 9277
freq=as.data.frame(table(mut$TrioId))
names(freq)=c('TrioId','NumberOfMut')
summary(freq$NumberOfMut)
temp = mut[mut$TrioId == freq$TrioId[1],] # and all these mutations are on codons 6, 23, 8 and 10.   
