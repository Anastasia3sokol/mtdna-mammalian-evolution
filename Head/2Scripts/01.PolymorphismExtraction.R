rm(list=ls(all=TRUE))

wd = getwd()
wd = paste(wd, '/mtdna-mammalian-evolution/Body/1Raw',sep='')
setwd(wd)

# List = list.files("../../Body/1Raw/PolymorphismsFromMutSpec/CYTB terminals/", pattern=".*\\.terminals.nuc.fa")
List = list.files("../../Body/1Raw/MutSpecTerminalsNucFa/", pattern=".*\\.terminals.nuc.fa")# функция, которая выбирает из папки файлы, имеющие в названии ".terminals.nuc.fa"

path="../../Body/1Raw/MutSpecTerminalsNucFa/"

Final=data.frame('Gene','Species','FirstName','SecondName','SubstVec'); names(Final)=c('Gene','Species','FirstName','SecondName','SubstVec')
Final=Final[-1,]

for (file in List) # length(List)
{# file = "Abbottina_obtusirostris.CytB.terminals.nuc.fa"
  Sp <- read.table(paste(path,file,sep=''), header = FALSE)
  Species=unlist(strsplit(file,'\\.'))[1]
  Gene=unlist(strsplit(file,'\\.'))[2]
  
  # if Species == 'XXXXX_sp.' => in this case Gene = ''; for example in Mastomys_sp..COX1.terminals.nuc.fa
  # to solve it we can do the next:
  if (Gene == '') {Gene=unlist(strsplit(file,'\\.'))[3]}

  Odd = seq(2,nrow(Sp),2)#обращается к четным строчкам файла (они содержат последовательности)
  NonOdd = seq(1,nrow(Sp),2)#обращается к нечетным строчкам файла (они содержат имена последовательностей - организмы)
  SpSeq=data.frame(Sp[Odd,])# создается таблица из последовательностей
  SpName=data.frame(Sp[NonOdd,])# создается таблица из организмов
  Sp=cbind(SpName,SpSeq); names(Sp)=c('Name','Seq')# таблицы с именами и последовательностями соединяются
  Sp$Name=gsub('>','',Sp$Name); Sp=Sp[Sp$Name != 'OUTGRP',]#убираются значки >
  
  if (nrow(Sp)>=10) 
  {
    Sp=Sp[sample(seq(1,nrow(Sp),1),10),]
    {
    for (i in 1:(nrow(Sp)-1))
      { # i = 1
      FirstName = Sp$Name[i]  
      FirstSeq = as.character(Sp$Seq[i])  
      FirstSeq = unlist(strsplit(FirstSeq,''))
      
      for (j in (i+1):(nrow(Sp)))
       { # j = 2
          SecondName = Sp$Name[j]  
          SecondSeq = as.character(Sp$Seq[j])  
          SecondSeq = unlist(strsplit(SecondSeq,''))
          
          SubstVec=''
          
          for (codon in 1:(length(FirstSeq)/3)) # should be divided by 3 without the rest
          { # codon=1
            FirstCodon = paste(FirstSeq[codon*3-2],FirstSeq[codon*3-1],FirstSeq[codon*3],sep='')
            SecondCodon = paste(SecondSeq[codon*3-2],SecondSeq[codon*3-1],SecondSeq[codon*3],sep='')
            if (FirstCodon != SecondCodon) {SubstVec = paste(SubstVec,codon,':',FirstCodon,'>',SecondCodon,';',sep='')}
          }
          OneLine=data.frame(Gene,Species,FirstName,SecondName,SubstVec)
          Final=rbind(Final,OneLine)
        }
      }
    }
  }
}          
 
length(unique(Final$Species)) # 
         
wd = gsub('/1Raw', '/2Derived', wd)
setwd(wd)
#write.table(DifferenceBetweenSpecies, file = "DifferenceBetweenSpeciesFamilies.txt")

write.table(Final,file = "PolymorphicPairwiseCodons.txt")
#write.table(Final,"../../Body/2Derived/PolymorphicPairwiseCodons.txt")
