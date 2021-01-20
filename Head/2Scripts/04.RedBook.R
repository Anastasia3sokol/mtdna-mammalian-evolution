rm(list=ls(all=TRUE))

wd = getwd()
wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
setwd(wd)

IUCN = read.csv("../../Body/1Raw/Red_book/IUCN.csv", sep=';', header = TRUE) #табличка по красной книге от Алины https://github.com/mitoclub/red-book/blob/master/Body/1Raw/IUCN.csv
Data = read.table("../../Body/2Derived/Distances_KnKs_Ecology_RG.csv", sep='\t', header = TRUE)# табличка из скрипта 03 с дистанциями, KnKS и экологией

names(IUCN)[3] <- "Species"

Dist_with_IUCN <- merge(Data, IUCN,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

Dist_with_IUCN = Dist_with_IUCN [-10:-21]
Dist_with_IUCN = Dist_with_IUCN [-11:-27]# тут обрезается куча данных, наверное, ненужных. В переменной category информация о том, в каком состоянии вид

#LC - Least Concern,таксоны, вызывающие наименьшие опасения
#VU - Vulnerable, Уязвимые таксоны
#EN - Endangered, Вымирающие таксоны 
#NT - Near Threatened, Таксоны, близкие к уязвимому положению
#DD - Data Deficient, Таксоны, для оценки угрозы которым недостаточно данных
#CR - Critically Endangered, Таксоны, находящиеся на грани полного исчезновения

ggplot(Dist_with_IUCN, aes(x = log(GenerationLength_d), y = KnKs, col = factor(category))) + 
  geom_point()

write.table(Dist_with_IUCN,file = "../../Body/2Derived/Dist_with_IUCN.csv",quote = F, row.names = FALSE,sep = '\t')
