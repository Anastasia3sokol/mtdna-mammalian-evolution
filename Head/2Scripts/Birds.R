rm(list=ls(all=TRUE))

wd = getwd()
wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
setwd(wd)

Birds = read.table("../../Body/1Raw/Birds/ALL_PHENOTYPES.txt", sep=',', header = TRUE) # Группы птиц от Вали
Data = read.table("../../Body/2Derived/Distances_KnKs_RG.csv", sep='\t', header = TRUE) # всякие дистанции по всем группам 
Data2 = read.table("../../Body/2Derived/Grantham_KnKs_RedBook.csv", sep='\t', header = TRUE)
Dist_with_Birds <- merge(Data, Birds ,by.x = "Species", by.y = "Tetrastes_sewerzowi",all = FALSE, no.dups = TRUE)# Всего лишь 42 вида птиц Вали совпадают с нашими
a = setdiff(Data$Species,Birds$Tetrastes_sewerzowi)

Taxons = read.table("../../Body/1Raw/TaxaWithClasses.txt", sep=' ', header = TRUE)# таксономия от Али


Data_Taxa <- merge(Data, Taxons ,by.x = "Species", by.y = "Species",all = FALSE, no.dups = TRUE)
My_birds = Data_Taxa[which(Data_Taxa$Class == "Aves"),]# Всего у нас 189 птиц


a = setdiff(My_birds$Species,Dist_with_Birds$Species)
#Data$Species = sub("_", " ", Data$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
