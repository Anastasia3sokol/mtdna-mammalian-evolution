#write.table(Final,"../../Body/2Derived/DnDs255.txt")
rm(list=ls(all=TRUE))

library(ape)
library(seqinr)

List = list.files("../../Body/1Raw/MutSpecTerminalsNucFa/", pattern=".*\\.terminals.nuc.fa")# функция, которая выбирает из папки файлы, имеющие в названии ".terminals.nuc.fa"

path="../../Body/1Raw/MutSpecTerminalsNucFa/"

Final=data.frame('Gene','Species','DnDs'); names(Final)=c('Gene','Species','DnDs')
Final=Final[-1,]

for (file in List) # length(List)
{# file = 'Anguilla_rostrata.ATP6.terminals.nuc.fa'
  Sp <- read.alignment(paste(path,file,sep=''), format = "fasta")
  Species=unlist(strsplit(file,'\\.'))[1]
  Gene=unlist(strsplit(file,'\\.'))[2]
  
  name_of_d = Sp$nam[duplicated(Sp$seq)]
  
  temp = unique.matrix(as.matrix(Sp))

  #Sp$seq = unique(Sp$seq)
  res = as.data.frame(as.matrix(dnds(as.DNAbin(temp),code = 2 ,codonstart = 1, quiet = F)))
  
  a = as.data.frame(res[[1]])
  DnDs = median(a[-1,])
  
  if (Gene == '') {Gene=unlist(strsplit(file,'\\.'))[3]}
  
  SpSeq=as.data.frame(unlist(Sp$seq))# создается таблица из последовательностей
  SpName=as.data.frame(Sp$nam) ; names(SpName) = "Name"# создается таблица из организмов 
  SpName = as.data.frame(SpName[SpName$Name!=name_of_d,]) 
  #res = as.data.frame(res[colnames(res)!=name_of_d,]) 
  
  sp=as.data.frame(cbind(res[,1])); names(sp)=c('DnDs')# таблицы с именами и последовательностями соединяются
  #sp$Name=gsub('>','',sp$Name); sp=sp[sp$Name != 'OUTGRP',]#убираются значки >
  
  OneLine=data.frame(Gene,Species,DnDs)
  Final=rbind(Final,OneLine)
}          

length(unique(Final$Species[Final$Gene == 'CytB'])) # 

Final1 = Final[Final$Gene == 'CytB',]
for (i in 1:nrow(Final1)){
  Final1$DnDs[i][is.infinite(Final1$DnDs[i])] = NA
}
Final1 = Final1[is.na(Final$DnDs)==F,]

write.table(Final,"../../Body/2Derived/DnDs.txt")

