rm(list=ls(all=TRUE))

mut = read.table('../../Body/2Derived/MutSpec.DerivedSubs.txt')
mut$MutSpec = as.character(mut$MutSpec)

new_mut = data.frame()


### extract information
for ( i in 1:nrow(mut))
{#i = 1
  subs = unlist(strsplit(mut$MutSpec[i],';'))
  for (c in 2:length(lenn))
  { #c = 2
    sub = unlist(strsplit(lenn[c],'\\.'))
    df = data.frame(mut$Outgroup[i],mut$Ingroup0[i],mut$Ingroup1[i],sub[1],sub[2],sub[3])
    new_mut = rbind(new_mut, df)
  }
}

names(new_mut) = c('Outgroup','Ingroup0','Ingroup1','CodonPos','Substitution','From>To')

write.table(new_mut,file = '../../Body/2Derived/MutSpec.ExtractedSubs.txt')
