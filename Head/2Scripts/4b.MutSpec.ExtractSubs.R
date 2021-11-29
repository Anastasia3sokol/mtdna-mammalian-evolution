rm(list=ls(all=TRUE))

mut = read.csv('../../Body/2Derived/MutSpec.DerivedSubs.csv')
mut$MutSpec = as.character(mut$MutSpec)

new_mut = data.frame()


### extract information
for ( i in 1:nrow(mut))
{#i = 1
  subs = unlist(strsplit(mut$MutSpec[i],';'))
  for (c in 2:length(subs))
  { #c = 2
    sub = unlist(strsplit(susb[c],'\\.'))
    df = data.frame(mut$Outgroup[i],mut$Ingroup0[i],mut$Ingroup1[i],sub[1],sub[2],sub[3])
    new_mut = rbind(new_mut, df)
  }
}

names(new_mut) = c('Outgroup','Ingroup0','Ingroup1','CodonPos','Substitution','From>To')

write.table(new_mut,file = '../../Body/2Derived/MutSpec.ExtractedSubs.txt')

genlen= read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
genlen = genlen[ ,c(2,11)] # Take information about name and GenLength
genlen$Scientific_name = gsub(' ','_',genlen$Scientific_name) 


finaldf = merge(new_mut, genlen, by.x = 'Ingroup0', by.y = 'Scientific_name' ) # GenLenIn0
colnames(finaldf)[colnames(finaldf) == 'GenerationLength_d'] <- 'GenLenIn0'

finaldf = merge(finaldf, genlen, by.x = 'Ingroup1', by.y = 'Scientific_name' ) # GenLenIn1
colnames(finaldf)[colnames(finaldf) == 'GenerationLength_d'] <- 'GenLenIn1'
names(finaldf)

write.table(finaldf, file = '../../Body/3Result/MutSpec.Final.WithGenLen.txt')



