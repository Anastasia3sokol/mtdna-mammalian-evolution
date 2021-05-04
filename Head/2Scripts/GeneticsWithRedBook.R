rm(list=ls(all=TRUE))

library(shades)

#wd = getwd()
#wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
#setwd(wd)

library(ggplot2)

IUCN = read.csv("../../Body/1Raw/Red_book/IUCN.csv", sep=';', header = TRUE) #табличка по красной книге от Алины https://github.com/mitoclub/red-book/blob/master/Body/1Raw/IUCN.csv
Data = read.table("../../Body/2Derived/CytB_with_ecologyForRedBook.csv", sep='\t', header = TRUE)
names(IUCN)[3] <- "Species"

Dist_with_IUCN <- merge(Data, IUCN,by.x = "Species", by.y = "Species",all = FALSE,no.dups = TRUE,)

Dist_with_IUCN = Dist_with_IUCN[Dist_with_IUCN$class == 'MAMMALIA' ,]

Dist_with_IUCN = Dist_with_IUCN [-27:-43]
Dist_with_IUCN = Dist_with_IUCN [-21:-25]
Dist_with_IUCN = Dist_with_IUCN [-14:-18]# тут обрезается куча данных, наверное, ненужных. В переменной category информация о том, в каком состоянии вид

#CR - Critically Endangered, Таксоны, находящиеся на грани полного исчезновения (1)
#EN - Endangered, Вымирающие таксоны (2)
#VU - Vulnerable, Уязвимые таксоны (3)
#NT - Near Threatened, Таксоны, близкие к уязвимому положению (4)
#LC - Least Concern,таксоны, вызывающие наименьшие опасения (5)
#DD - Data Deficient, Таксоны, для оценки угрозы которым недостаточно данных - выкинуть
#EX - вымершие

Data = subset(Dist_with_IUCN, Dist_with_IUCN$category != "DD") # обрезаем Data Deficient тк ничего про них не знаем, их всего лишь 10
Data$category = sub("NT","LC",Data$category) # категория NT присоединина к категории LC
Data$category = sub("LR/cd","LC",Data$category) # категория LR/cd присоединина к категории LC

ggplot(Data, aes(x = log(GenerationLength_d), y = KnKs, col = factor(category),hsv(0.5, 0.7, 0.9, maxColorValue=255, alpha=0)))  +
  geom_point() 

ggplot(Data, aes(x = log(GenerationLength_d), y = AverageGrantham, col = factor(category))) + 
  geom_point()

Data$CategoryPoisson = Data$category
Data$CategoryPoisson = sub("LC",0,Data$CategoryPoisson)
Data$CategoryPoisson = sub("VU",1,Data$CategoryPoisson)
Data$CategoryPoisson = sub("EN",2,Data$CategoryPoisson)
Data$CategoryPoisson = sub("CR",3,Data$CategoryPoisson)
Data$CategoryPoisson = sub("EX",4,Data$CategoryPoisson)
Data$CategoryPoisson = as.numeric(Data$CategoryPoisson) # from character to numeric
summary(Data$CategoryPoisson)
table(Data$CategoryPoisson)

Data$CategoryBinomial = as.numeric(Data$CategoryPoisson)
Data$CategoryBinomial = as.numeric(sub("1",0,Data$CategoryBinomial))
Data$CategoryBinomial = as.numeric(sub("1|2|3|4",1,Data$CategoryBinomial))
table(Data$CategoryBinomial)

write.table(Data,file = "../../Body/2Derived/CytB_RedBook.csv",quote = F, row.names = FALSE,sep = '\t')

######## КОРРЕЛЯЦИИ

IUCN_AverageGrantham_Binomial = cor.test(Data$AverageGrantham,Data$CategoryBinomial, method = 'spearman')
IUCN_AverageGrantham_Poisson = cor.test(Data$AverageGrantham,Data$CategoryPoisson, method = 'spearman')
IUCN_DnDs_Binomial = cor.test(Data$DnDs,Data$CategoryBinomial, method = 'spearman')
IUCN_DnDs_Poisson = cor.test(Data$DnDs,Data$CategoryPoisson, method = 'spearman')
IUCN_KnKs_Binomial = cor.test(Data$KnKs,Data$CategoryBinomial, method = 'spearman')
IUCN_KnKs_Poisson = cor.test(Data$KnKs,Data$CategoryPoisson, method = 'spearman')

estimate = data.frame(IUCN_AverageGrantham_Binomial$estimate,IUCN_AverageGrantham_Poisson$estimate,IUCN_DnDs_Binomial$estimate,IUCN_DnDs_Poisson$estimate,
                   IUCN_KnKs_Binomial$estimate,IUCN_KnKs_Poisson$estimate);names(estimate) = c('GranthamBin','GranthamPoi','DnDsBin','DnDsPoi','KnKsBin','KnKsPoi')
pvalue = data.frame(IUCN_AverageGrantham_Binomial$p.value,IUCN_AverageGrantham_Poisson$p.value,
             IUCN_DnDs_Binomial$p.value,IUCN_DnDs_Poisson$p.value,IUCN_KnKs_Binomial$p.value,
             IUCN_KnKs_Poisson$p.value);names(pvalue) = c('GranthamBin','GranthamPoi','DnDsBin','DnDsPoi','KnKsBin','KnKsPoi')

fit = rbind (estimate,pvalue); rownames(fit)= c('estimate', 'p-value')

########### GLM

Binomial = glm(CategoryBinomial ~ scale(AverageGrantham) + scale(KnKs)+scale(DnDs), family = 'binomial', data = Data); 
summary(Binomial)

Poisson = glm(CategoryPoisson ~ scale(AverageGrantham) + scale(KnKs) + scale(DnDs), family = 'poisson', data = Data); 
summary(Poisson)


########## ПИКИ

library(ape)
library(geiger)
library(caper)
library(stringr)
library(phytools)
library(picante)

Data = Data[Data$Species!='Neophocaena phocaenoides',] 

tree = read.tree("../../Body/1Raw/mammalia_species.nwk") 
Data$Species = sub(" ", "_", Data$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)


row.names(Data) = Data$Species

nrow(Data) # 347
tree_w = treedata(tree, Data[, c('Species', 'AverageGrantham', 'KnKs','DnDs' ,"CategoryPoisson",'CategoryBinomial')], 
                  sort=T, warnings=T)$phy
data_tree = as.data.frame(treedata(tree_w, Data[, c('Species', 'AverageGrantham', 'KnKs','DnDs',"CategoryPoisson",'CategoryBinomial')], 
                                   sort=T, warnings=T)$data)

setdiff(Data$Species, data_tree$Species)

data_tree$AverageGrantham = as.numeric(as.character(data_tree$AverageGrantham))
data_tree$KnKs = as.numeric(as.character(data_tree$KnKs))
data_tree$CategoryPoisson = as.numeric(as.character(data_tree$CategoryPoisson))
data_tree$CategoryBinomial = as.numeric(as.character(data_tree$CategoryBinomial))

IUCN_AverageGrantham_Binomial = cor.test(pic(data_tree$AverageGrantham, tree_w), pic((data_tree$CategoryBinomial), tree_w), method = 'spearman')
IUCN_AverageGrantham_Poisson = cor.test(pic(data_tree$AverageGrantham, tree_w), pic((data_tree$CategoryPoisson), tree_w), method = 'spearman')
IUCN_DnDs_Binomial = cor.test(pic(data_tree$DnDs, tree_w), pic((data_tree$CategoryBinomial), tree_w), method = 'spearman')
IUCN_DnDs_Poisson = cor.test(pic(data_tree$DnDs, tree_w), pic((data_tree$CategoryPoisson), tree_w), method = 'spearman')
IUCN_KnKs_Binomial = cor.test(pic(data_tree$KnKs, tree_w), pic((data_tree$CategoryBinomial), tree_w), method = 'spearman')
IUCN_KnKs_Poisson = cor.test(pic(data_tree$KnKs, tree_w), pic((data_tree$CategoryPoisson), tree_w), method = 'spearman')

estimate = data.frame(IUCN_AverageGrantham_Binomial$estimate,IUCN_AverageGrantham_Poisson$estimate,IUCN_DnDs_Binomial$estimate,IUCN_DnDs_Poisson$estimate,
                      IUCN_KnKs_Binomial$estimate,IUCN_KnKs_Poisson$estimate);names(estimate) = c('GranthamBin','GranthamPoi','DnDsBin','DnDsPoi','KnKsBin','KnKsPoi')
pvalue = data.frame(IUCN_AverageGrantham_Binomial$p.value,IUCN_AverageGrantham_Poisson$p.value,
                    IUCN_DnDs_Binomial$p.value,IUCN_DnDs_Poisson$p.value,IUCN_KnKs_Binomial$p.value,
                    IUCN_KnKs_Poisson$p.value);names(pvalue) = c('GranthamBin','GranthamPoi','DnDsBin','DnDsPoi','KnKsBin','KnKsPoi')

fit_PICs = rbind (estimate,pvalue); rownames(fit_PICs)= c('estimate', 'p-value')

############# МАНН-УИТНИ

wilcox.test(Data$AverageGrantham~Data$CategoryBinomial)

wilcox.test(Data$KnKs~Data$CategoryBinomial)

wilcox.test(Data$DnDs~Data$CategoryBinomial)

############ ЛЯМБДА


CategoryBinomial <- data_tree$CategoryBinomial
names(CategoryBinomial) <- rownames(data_tree)

CategoryPoisson <- data_tree$CategoryPoisson
names(CategoryPoisson) <- rownames(data_tree)

Grantham <- data_tree$AverageGrantham
names(Grantham) <- rownames(data_tree)

DnDs <- as.matrix(as.numeric(data_tree$DnDs))
rownames(DnDs) <- rownames(data_tree)

phylosig(tree_w, CategoryBinomial, method = "lambda", test = TRUE) 
phylosig(tree_w, CategoryPoisson, method = "lambda", test = TRUE) 
phylosig(tree_w, Grantham, method = "lambda", test = TRUE)
phylosig(tree_w, DnDs, method = "lambda", test = TRUE)

############## PGLS

MutComp = comparative.data(tree_w, Data, Species, vcv=TRUE)

summary(pgls(scale(AverageGrantham) ~ scale(CategoryBinomial), MutComp, lambda="ML"))
summary(pgls(scale(AverageGrantham) ~ scale(CategoryPoisson), MutComp, lambda="ML"))

summary(pgls(scale(DnDs) ~ scale(CategoryBinomial), MutComp, lambda="ML"))
summary(pgls(scale(DnDs) ~ scale(CategoryPoisson), MutComp, lambda="ML"))

summary(pgls(scale(KnKs) ~ scale(CategoryBinomial), MutComp, lambda="ML"))
summary(pgls(scale(KnKs) ~ scale(CategoryPoisson), MutComp, lambda="ML"))

############ ДЕЛЕНИЕ ПО СЕМЕЙСТВАМ
unique(Data$specieOrder) # [1] "RODENTIA","CARNIVORA","PRIMATES","EULIPOTYPHLA","CHIROPTERA","CETARTIODACTYLA","CINGULATA","PROBOSCIDEA","DIDELPHIMORPHIA"
#"LAGOMORPHA","MONOTREMATA","DASYUROMORPHIA","SCANDENTIA"

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

length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'RODENTIA',]$Species)) # 3
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'CARNIVORA',]$Species)) # 2
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'PRIMATES',]$Species)) # 15
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'EULIPOTYPHLA',]$Species)) # 0
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'CHIROPTERA',]$Species)) # 0
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'CETARTIODACTYLA',]$Species)) # 0
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'CINGULATA',]$Species)) # 0
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'PROBOSCIDEA',]$Species)) # 1
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'DIDELPHIMORPHIA',]$Species)) # 0
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'LAGOMORPHA',]$Species)) # 0
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'MONOTREMATA',]$Species)) # 0
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'DASYUROMORPHIA',]$Species)) # 1
length(unique(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'SCANDENTIA',]$Species)) # 0


boxplot(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'RODENTIA',]$AverageGrantham,+
          Data[Data$CategoryBinomial == 0 & Data$specieOrder == 'RODENTIA',]$AverageGrantham,+
          Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'CARNIVORA',]$AverageGrantham,+
          Data[Data$CategoryBinomial == 0 & Data$specieOrder == 'CARNIVORA',]$AverageGrantham,+
          Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'PRIMATES',]$AverageGrantham,+
          Data[Data$CategoryBinomial == 0 & Data$specieOrder == 'PRIMATES',]$AverageGrantham, names = c("Threatened Rodentia",'LC Rodentia',"Threatened Carnivora",'LC Carnivora',"Threatened Primates",'LC Primates'), outline = FALSE, xlab = "GenLenght", ylab = "AverageGrantham",notch = TRUE)

boxplot(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'RODENTIA',]$DnDs,+
          Data[Data$CategoryBinomial == 0 & Data$specieOrder == 'RODENTIA',]$DnDs,+
          Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'CARNIVORA',]$DnDs,+
          Data[Data$CategoryBinomial == 0 & Data$specieOrder == 'CARNIVORA',]$DnDs,+
          Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'PRIMATES',]$DnDs,+
          Data[Data$CategoryBinomial == 0 & Data$specieOrder == 'PRIMATES',]$DnDs, names = c("Threatened Rodentia",'LC Rodentia',"Threatened Carnivora",'LC Carnivora',"Threatened Primates",'LC Primates'), outline = FALSE, xlab = "GenLenght", ylab = "DnDs",notch = TRUE)

boxplot(Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'RODENTIA',]$KnKs,+
          Data[Data$CategoryBinomial == 0 & Data$specieOrder == 'RODENTIA',]$KnKs,+
          Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'CARNIVORA',]$KnKs,+
          Data[Data$CategoryBinomial == 0 & Data$specieOrder == 'CARNIVORA',]$KnKs,+
          Data[Data$CategoryBinomial == 1 & Data$specieOrder == 'PRIMATES',]$KnKs,+
          Data[Data$CategoryBinomial == 0 & Data$specieOrder == 'PRIMATES',]$KnKs, names = c("Threatened Rodentia",'LC Rodentia',"Threatened Carnivora",'LC Carnivora',"Threatened Primates",'LC Primates'), outline = FALSE, xlab = "GenLenght", ylab = "KnKs",notch = TRUE)


########################################### Red book ~ Longivity

Data = Data[is.na(Data$GenerationLength_d)==F,]

MutComp = comparative.data(tree_w, Data, Species, vcv=TRUE)

model3KP.AverageGrantham_Poisson = glm(CategoryPoisson ~ scale(log2(GenerationLength_d)) + scale(AverageGrantham), family = 'poisson', data = Data); # model1 = glm(Data$KnKs ~ log(Data$GenerationLength_d) + IucnRanks, family = 'poisson', data = Data); 
summary(model3KP.AverageGrantham_Poisson)

model3KP.AverageGrantham_Binomial = glm(CategoryBinomial ~ scale(log2(GenerationLength_d)) + scale(AverageGrantham), family = 'binomial', data = Data); # model1 = glm(Data$KnKs ~ log(Data$GenerationLength_d) + IucnRanks, family = 'poisson', data = Data); 
summary(model3KP.AverageGrantham_Binomial)


summary(pgls(CategoryPoisson ~ log2(GenerationLength_d)+scale(AverageGrantham), MutComp, lambda="ML"))
summary(pgls(CategoryBinomial ~ 0 + log2(GenerationLength_d)+scale(AverageGrantham), MutComp, lambda="ML"))

summary(pgls(scale(AverageGrantham) ~ log2(GenerationLength_d)+(CategoryBinomial), MutComp, lambda="ML"))
summary(pgls(scale(AverageGrantham) ~ log2(GenerationLength_d)+(CategoryPoisson), MutComp, lambda="ML"))

summary(pgls(scale(AverageGrantham) ~ log2(GenerationLength_d)*(CategoryBinomial), MutComp, lambda="ML"))
summary(pgls(scale(AverageGrantham) ~ log2(GenerationLength_d)*(CategoryPoisson), MutComp, lambda="ML"))

