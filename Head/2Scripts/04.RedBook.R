rm(list=ls(all=TRUE))

wd = getwd()
wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
setwd(wd)

library(ggplot2)

IUCN = read.csv("../../Body/1Raw/Red_book/IUCN.csv", sep=';', header = TRUE) #табличка по красной книге от Алины https://github.com/mitoclub/red-book/blob/master/Body/1Raw/IUCN.csv
Data = read.table("../../Body/2Derived/Distances_KnKs_Ecology_RG.csv", sep='\t', header = TRUE)# табличка из скрипта 03 с дистанциями, KnKS и экологией

names(IUCN)[3] <- "Species"

Dist_with_IUCN <- merge(Data, IUCN,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

Dist_with_IUCN = Dist_with_IUCN [-11:-21]
Dist_with_IUCN = Dist_with_IUCN [,-11]
Dist_with_IUCN = Dist_with_IUCN [-12:-28]# тут обрезается куча данных, наверное, ненужных. В переменной category информация о том, в каком состоянии вид

#CR - Critically Endangered, Таксоны, находящиеся на грани полного исчезновения (1)
#EN - Endangered, Вымирающие таксоны (2)
#VU - Vulnerable, Уязвимые таксоны (3)
#NT - Near Threatened, Таксоны, близкие к уязвимому положению (4)
#LC - Least Concern,таксоны, вызывающие наименьшие опасения (5)
#DD - Data Deficient, Таксоны, для оценки угрозы которым недостаточно данных - выкинуть

ggplot(Dist_with_IUCN, aes(x = log(GenerationLength_d), y = KnKs, col = factor(category))) + 
  geom_point()

ggplot(Dist_with_IUCN, aes(x = log(GenerationLength_d), y = AverageGrantham, col = factor(category))) + 
  geom_point()

write.table(Dist_with_IUCN,file = "../../Body/2Derived/Dist_with_IUCN.csv",quote = F, row.names = FALSE,sep = '\t')

Data = subset(Dist_with_IUCN, Dist_with_IUCN[,11] != "DD") # обрезаем Data Deficient тк ничего про них не знаем, их всего лишь 10
Data$category = sub("NT","LC",Data$category) # категория NT присоединина к категории LC

Data = subset(Dist_with_IUCN, Dist_with_IUCN[,11] != "DD") # обрезаем Data Deficient тк ничего про них не знаем, их всего лишь 10

Data$CategoryPoisson = Data$category
Data$CategoryPoisson = sub("LC",0,Data$CategoryPoisson)
Data$CategoryPoisson = sub("NT",1,Data$CategoryPoisson)
Data$CategoryPoisson = sub("VU",2,Data$CategoryPoisson)
Data$CategoryPoisson = sub("EN",3,Data$CategoryPoisson)
Data$CategoryPoisson = sub("CR",4,Data$CategoryPoisson)
Data$CategoryPoisson = as.numeric(Data$CategoryPoisson) # from character to numeric
summary(Data$CategoryPoisson)
table(Data$CategoryPoisson)

Data$CategoryBinomial = as.numeric(Data$CategoryPoisson)
Data$CategoryBinomial = as.numeric(sub("1",0,Data$CategoryBinomial))
Data$CategoryBinomial = as.numeric(sub("2|3|4",1,Data$CategoryBinomial))
table(Data$CategoryBinomial)

#model1 = glm(KnKs ~ log2(GenerationLength_d) + category, family = 'poisson', data = Data); 
#summary(model1)

model1KP.KnKs_Poisson = glm(CategoryPoisson ~ scale(log2(GenerationLength_d)) + scale(KnKs), family = 'poisson', data = Data); 
summary(model1KP.KnKs_Poisson)

model1KPKnKs_Binomial = glm(CategoryBinomial ~ scale(log2(GenerationLength_d)) + scale(KnKs), family = 'binomial', data = Data); 
summary(model1KPKnKs_Binomial)

ggplot(Data,aes(GenerationLength_d, KnKs))+
  geom_point(size = 1)+
  geom_smooth(method = "lm")+
  facet_grid(.~category)

#model3 = glm(AverageGrantham ~ log2(GenerationLength_d) + category, family = 'poisson', data = Data); # model1 = glm(Data$KnKs ~ log(Data$GenerationLength_d) + IucnRanks, family = 'poisson', data = Data); 
#summary(model3)

model3KP.AverageGrantham_Poisson = glm(CategoryPoisson ~ scale(log2(GenerationLength_d)) + scale(AverageGrantham), family = 'poisson', data = Data); # model1 = glm(Data$KnKs ~ log(Data$GenerationLength_d) + IucnRanks, family = 'poisson', data = Data); 
summary(model3KP.AverageGrantham_Poisson)

model3KP.AverageGrantham_Binomial = glm(CategoryBinomial ~ scale(log2(GenerationLength_d)) + scale(AverageGrantham), family = 'binomial', data = Data); # model1 = glm(Data$KnKs ~ log(Data$GenerationLength_d) + IucnRanks, family = 'poisson', data = Data); 
summary(model3KP.AverageGrantham_Binomial)

ggplot(Data,aes(GenerationLength_d, AverageGrantham))
  geom_point(size = 2)+
  geom_smooth(method = "lm")+
  facet_grid(.~category)


#########
# чем отличается log и log2?

Data2 = Data

Data2$category = sub("2",'1',Data2$category)
Data2$category = sub("3",'1',Data2$category)
Data2$category = sub("4",'1',Data2$category)
Data2$category = sub("5",'0',Data2$category) #Least Concern

model2 = glm(KnKs ~ log2(GenerationLength_d) + category, family = 'binomial', data = Data2); 
summary(model2)

ggplot(Data2,aes(GenerationLength_d, KnKs))+
  geom_point(size = 2)+
  geom_smooth(method = "lm")+
  facet_grid(.~category)

model4 = glm(AverageGrantham ~ log2(GenerationLength_d) + category, family = 'poisson', data = Data2); # model1 = glm(Data$KnKs ~ log(Data$GenerationLength_d) + IucnRanks, family = 'poisson', data = Data); 
summary(model4)

ggplot(Data2,aes(GenerationLength_d, AverageGrantham))+
  geom_point(size = 2)+
  geom_smooth(method = "lm")+
  facet_grid(.~category)

### PICs 

library(ape)
library(geiger)
library(caper)

cor.test(pic(Data$AverageGrantham, Data$), pic(data_tree$GenerationLength_d, tree_w), method = 'spearman')
# rho = 0.07723414, p-value = 0.04939


#Primates - приматы
#Rodentia - грызуны
#Carnivora - хищные
#Afrosoricida - Афросорици́ды
#Eulipotyphla - насекомоядные
#Chiroptera - рукокрылые
#Cetartiodactyla - китопарнокопытные
#Proboscidea - хоботные
#Didelphimorphia - опоссумы
#Lagomorpha - зайцеобразные
#Monotremata - однопроходные 
#Scandentia - тупаи