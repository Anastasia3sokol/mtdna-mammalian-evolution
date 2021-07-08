rm(list=ls(all=TRUE))

library(seqinr)
library(dplyr)

mut = read.csv('4Fold_MutSpec_Trio.csv')

mut = mut[,c(2,3,4,5,1,6)]

names(mut) = c('Ingroup0','Ingroup1','Outgroup','Family','Substitution','Rel_Of_Subs')

In1Out = subset(mut, Rel_Of_Subs =='subsIn1Out')
In1Out = select(In1Out, -Ingroup0)
names(In1Out) = c('Species1','Species2','Family','Substitution','Rel_Of_Subs')

In0Out = subset(mut, Rel_Of_Subs =='subsIn0Out')
In0Out = select(In0Out, -Ingroup1)
names(In0Out) = c('Species1','Species2','Family','Substitution','Rel_Of_Subs')

In0In1 = subset(mut, Rel_Of_Subs =='subsIn0In1')
In0In1 = select(In0In1, -Outgroup)
names(In0In1) = c('Species1','Species2','Family','Substitution','Rel_Of_Subs')

final_mut = rbind(In0In1,In0Out,In1Out)

write.table(final_mut,'Mut_Spec_4Fold_Without_Normalization.txt')
