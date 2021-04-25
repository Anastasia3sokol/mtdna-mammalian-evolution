rm(list=ls(all=TRUE))

library(shades)

wd = getwd()
wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
setwd(wd)

library(ggplot2)

IUCN = read.csv("../../Body/1Raw/Red_book/IUCN.csv", sep=';', header = TRUE) #табличка по красной книге от Алины https://github.com/mitoclub/red-book/blob/master/Body/1Raw/IUCN.csv
Data = read.table("../../Body/2Derived/CytB_with_ecologyForRedBook.csv", sep='\t', header = TRUE)# табличка из скрипта 03 с дистанциями, KnKS и экологией

names(IUCN)[3] <- "Species"

Dist_with_IUCN <- merge(Data, IUCN,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

Dist_with_IUCN = Dist_with_IUCN [-26:-42]
Dist_with_IUCN = Dist_with_IUCN [-20:-24]
Dist_with_IUCN = Dist_with_IUCN [-13:-17]# тут обрезается куча данных, наверное, ненужных. В переменной category информация о том, в каком состоянии вид
Dist_with_IUCN = Dist_with_IUCN [-8];Dist_with_IUCN = Dist_with_IUCN [-10]

#CR - Critically Endangered, Таксоны, находящиеся на грани полного исчезновения (1)
#EN - Endangered, Вымирающие таксоны (2)
#VU - Vulnerable, Уязвимые таксоны (3)
#NT - Near Threatened, Таксоны, близкие к уязвимому положению (4)
#LC - Least Concern,таксоны, вызывающие наименьшие опасения (5)
#DD - Data Deficient, Таксоны, для оценки угрозы которым недостаточно данных - выкинуть


Data = subset(Dist_with_IUCN, Dist_with_IUCN$category != "DD") # обрезаем Data Deficient тк ничего про них не знаем, их всего лишь 10
Data$category = sub("NT","LC",Data$category) # категория NT присоединина к категории LC

ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs, col = factor(category),hsv(0.5, 0.7, 0.9, maxColorValue=255, alpha=0)))  +
  
  geom_point() 



ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGrantham, col = factor(category))) + 
  geom_point()

write.table(Data,file = "../../Body/2Derived/CytB_RedBook.csv",quote = F, row.names = FALSE,sep = '\t')

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

write.table(Data,file = "../../Body/2Derived/CytB_RedBook.csv",quote = F, row.names = FALSE,sep = '\t')

model1KP.KnKs_Poisson = glm(CategoryPoisson ~ scale(log2(GenerationLength_d)) + scale(KnKs), family = 'poisson', data = Data); 
summary(model1KP.KnKs_Poisson)

model1KPKnKs_Binomial = glm(CategoryBinomial ~ scale(log2(GenerationLength_d)) + scale(KnKs), family = 'binomial', data = Data); 
summary(model1KPKnKs_Binomial)

ggplot(Data,aes(GenerationLength_d, KnKs), color = RGB)+
  geom_point(size = 1)+
  geom_smooth(method = "lm")+
  facet_grid(.~category)

model3KP.AverageGrantham_Poisson = glm(CategoryPoisson ~ scale(log2(GenerationLength_d)) + scale(AverageGrantham), family = 'poisson', data = Data); # model1 = glm(Data$KnKs ~ log(Data$GenerationLength_d) + IucnRanks, family = 'poisson', data = Data); 
summary(model3KP.AverageGrantham_Poisson)

model3KP.AverageGrantham_Binomial = glm(CategoryBinomial ~ scale(log2(GenerationLength_d)) + scale(AverageGrantham), family = 'binomial', data = Data); # model1 = glm(Data$KnKs ~ log(Data$GenerationLength_d) + IucnRanks, family = 'poisson', data = Data); 
summary(model3KP.AverageGrantham_Binomial)

ggplot(Data,aes(GenerationLength_d, AverageGrantham))+
  geom_point(size = 2)+
  geom_smooth(method = "lm")+
  facet_grid(.~category)


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