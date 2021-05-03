rm(list=ls(all=TRUE)) 

library(ape)
library(geiger)
library(caper)
library(stringr)
library(phytools)
library(picante)

wd = getwd()
wd = paste(wd, '/mtdna-mammalian-evolution/Body/1Raw',sep='')
setwd(wd)

tree = read.tree("../../Body/1Raw/mammalia_species.nwk") 
data_redbook = read.csv ("../../Body/2Derived/CytB_RedBook.csv", sep = "\t")# нет пропусков в генетике и КК, но есть в систематике и продолжительности жизни. Нужно смеорджить с еще систематикой!

data = data_redbook[is.na(data_redbook$GenerationLength_d)==F,]
data_redbook = data_redbook[,-5]

wilcox.test(data$AverageGrantham~data$CategoryBinomial)
#p-value = 0.002681
wilcox.test(data$KnKs~data$CategoryBinomial)
#p-value = 0.0001416
wilcox.test(data$DnDs~data$CategoryBinomial)
#p-value = 9.494e-06

data$Species = sub(" ", "_", data$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)


#data = data[data$Species!='Neophocaena_phocaenoides',] 277 и 278 строки
data = data[-277,]


row.names(data) = data$Species

nrow(data) # 344
tree_w = treedata(tree, data[, c('Species', 'AverageGrantham', 'KnKs','DnDs' ,'GenerationLength_d',"CategoryPoisson",'CategoryBinomial')], 
                  sort=T, warnings=T)$phy
data_tree = as.data.frame(treedata(tree_w, data[, c('Species', 'AverageGrantham', 'KnKs','DnDs' ,'GenerationLength_d',"CategoryPoisson",'CategoryBinomial')], 
                                   sort=T, warnings=T)$data)

setdiff(data$Species, data_tree$Species)

data_tree$AverageGrantham = as.numeric(as.character(data_tree$AverageGrantham))
data_tree$GenerationLength_d = as.numeric(as.character(data_tree$GenerationLength_d))
data_tree$KnKs = as.numeric(as.character(data_tree$KnKs))

#GenLen & Grantham
cor.test(pic(data_tree$AverageGrantham, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w), method = 'spearman')
#rho = 0.07973827; p-value = 0.1459; 
cor.test(data_tree$AverageGrantham,data_tree$GenerationLength_d, method = 'spearman')
#rho = 0.3490021; p-value = 4.975e-11;

#GenLen & KnKs
cor.test(pic(data_tree$KnKs, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w), method = 'spearman')
#rho = 0.06296045; p-value = 0.2512;
cor.test(data_tree$KnKs,data_tree$GenerationLength_d, method = 'spearman')
#rho = 0.3860619; p-value = 2.378e-13;

#GenLen & DnDs
cor.test(pic(data_tree$DnDs, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w), method = 'spearman')
#rho = 0.1271365; p-value = 0.02011;
cor.test(data_tree$KnKs,data_tree$GenerationLength_d, method = 'spearman')
#rho = 0.3860619; p-value = 2.378e-13;

#RedBook & Grantham
MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)

summary(pgls(scale(AverageGrantham) ~ scale(GenerationLength_d), MutComp, lambda="ML"))
summary(pgls(scale(KnKs) ~ log2(GenerationLength_d), MutComp, lambda="ML"))
summary(pgls(scale(DnDs) ~ log2(GenerationLength_d), MutComp, lambda="ML"))

summary(pgls(CategoryPoisson ~ log2(GenerationLength_d)+scale(AverageGrantham), MutComp, lambda="ML"))
summary(pgls(CategoryBinomial ~ 0 + log2(GenerationLength_d)+scale(AverageGrantham), MutComp, lambda="ML"))

summary(pgls(scale(AverageGrantham) ~ log2(GenerationLength_d)+(CategoryBinomial), MutComp, lambda="ML"))
summary(pgls(scale(AverageGrantham) ~ log2(GenerationLength_d)+(CategoryPoisson), MutComp, lambda="ML"))

summary(pgls(scale(AverageGrantham) ~ log2(GenerationLength_d)*(CategoryBinomial), MutComp, lambda="ML"))
summary(pgls(scale(AverageGrantham) ~ log2(GenerationLength_d)*(CategoryPoisson), MutComp, lambda="ML"))

#tree0 <- rescale(tree_w, model = "lambda", 1)


GenLen <- data_tree$GenerationLength_d
names(GenLen) <- rownames(data_tree)

Grantham <- data_tree$AverageGrantham
names(Grantham) <- rownames(data_tree)

DnDs <- as.matrix(as.numeric(data_tree$DnDs))
rownames(DnDs) <- rownames(data_tree)

phylosig(tree_w, GenLen, method = "lambda", test = TRUE) 
phylosig(tree_w, Grantham, method = "lambda", test = TRUE)
phylosig(tree_w, DnDs, method = "lambda", test = TRUE)


hist(data$GenerationLength_d,breaks=10)

##### boxplots by quartiles
boxplot(data[data$GenerationLength_d<=quantile(data$GenerationLength_d,0.25),]$AverageGrantham,
        data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.25) & data$GenerationLength_d<=quantile(data$GenerationLength_d,0.5),]$AverageGrantham,
        data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.5) & data$GenerationLength_d<=quantile(data$GenerationLength_d,0.75),]$AverageGrantham,
        data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.75),]$AverageGrantham,
        names=c('1','2','3','4'), outline = FALSE, notch = TRUE)
boxplot(data[data$GenerationLength_d<=quantile(data$GenerationLength_d,0.25),]$DnDs,
        data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.25) & data$GenerationLength_d<=quantile(data$GenerationLength_d,0.5),]$DnDs,
        data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.5) & data$GenerationLength_d<=quantile(data$GenerationLength_d,0.75),]$DnDs,
        data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.75),]$DnDs,
        names=c('1','2','3','4'), outline = FALSE, notch = TRUE)

wilcox.test(
  data[data$GenerationLength_d<=quantile(data$GenerationLength_d,0.25),]$AverageGrantham,
  data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.25) & data$GenerationLength_d<=quantile(data$GenerationLength_d,0.5),]$AverageGrantham)
wilcox.test(
  data[data$GenerationLength_d<=quantile(data$GenerationLength_d,0.25),]$AverageGrantham,
  data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.5) & data$GenerationLength_d<=quantile(data$GenerationLength_d,0.75),]$AverageGrantham)
wilcox.test(
  data[data$GenerationLength_d<=quantile(data$GenerationLength_d,0.25),]$AverageGrantham,
  data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.75),]$AverageGrantham)
wilcox.test(
  data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.25) & data$GenerationLength_d<=quantile(data$GenerationLength_d,0.5),]$AverageGrantham,
  data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.5) & data$GenerationLength_d<=quantile(data$GenerationLength_d,0.75),]$AverageGrantham)
wilcox.test(
  data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.25) & data$GenerationLength_d<=quantile(data$GenerationLength_d,0.5),]$AverageGrantham,
  data[data$GenerationLength_d>quantile(data$GenerationLength_d,0.75),]$AverageGrantham)


OrderFreq = data.frame(table(data_redbook$specieOrder));
FrequentOrders = OrderFreq[OrderFreq$Freq >= 3,]$Var1; length(FrequentOrders) # 3!!!
OrderFreq = OrderFreq[OrderFreq$Var1 %in% FrequentOrders,]
names(OrderFreq)=c('Order','NumberOfSpecies')

Mammalia = data[data$specieOrder %in% FrequentOrders,]
agg = aggregate(list(Mammalia$AverageGrantham,Mammalia$GenerationLength_d,Mammalia$DnDs,Mammalia$KnKs), by = list(Mammalia$specieOrder), FUN = median)
names(agg) = c('Order','AverageGrantham','GenerationLength_d','DnDs','KnKs')

test1 = cor.test(agg$AverageGrantham,agg$GenerationLength_d,method = 'spearman') ### 0.008503, Rho = 0.8434347 !!!#не работает с пропусками в продолжительности жизни
plot(agg$GenerationLength_d,agg$AverageGrantham)
test2 = cor.test(agg$DnDs,agg$GenerationLength_d,method = 'spearman') ### 0.01071, Rho = 0.8571429 !!!
plot(agg$DnDs,agg$AverageGrantham)
test3 = cor.test(agg$KnKs,agg$GenerationLength_d,method = 'spearman') ### 0.004563, Rho = 0.9047619 !!!
plot(agg$KnKs,agg$AverageGrantham)

OrderGrantham_GenLen = c(test1$p.value,test1$estimate)
OrderDnDs_GenLen = c(test2$p.value,test2$estimate)
OrderKnKs_GenLen = c(test3$p.value,test3$estimate)

OrderFinal = data.frame(OrderGrantham_GenLen,OrderDnDs_GenLen,OrderKnKs_GenLen); rownames(OrderFinal) = c('p-value','rho')

FamFreq = data.frame(table(data_redbook$family));
FrequentFamilies = FamFreq[FamFreq$Freq >= 3,]$Var1; length(FrequentFamilies) # 3!!!
FamFreq = FamFreq[FamFreq$Var1 %in% FrequentFamilies,]
names(FamFreq)=c('Family','NumberOfSpecies')

# Family
Mammalia1 = data[data$family %in% FrequentFamilies,]
agg1 = aggregate(list(Mammalia1$AverageGrantham,Mammalia1$GenerationLength_d,Mammalia1$DnDs,Mammalia1$KnKs), by = list(Mammalia1$family), FUN = median)
names(agg1) = c('Family','AverageGrantham','GenerationLength_d','DnDs','KnKs')

test4 = cor.test(agg1$AverageGrantham,agg1$GenerationLength_d,method = 'spearman') ### 0.00281, Rho = 0.5108606 !!!
plot(agg1$GenerationLength_d,agg1$AverageGrantham)
test5 = cor.test(agg1$DnDs,agg1$GenerationLength_d,method = 'spearman') ### 0.003045, Rho = 0.5128299 !!!
plot(agg1$DnDs,agg1$AverageGrantham)
test6 = cor.test(agg1$KnKs,agg1$GenerationLength_d,method = 'spearman') ### 0.01945, Rho = 0.4109991  !!!
plot(agg1$KnKs,agg1$AverageGrantham)

FamilyGrantham_GenLen = c(test4$p.value,test4$estimate)
FamilyDnDs_GenLen = c(test5$p.value,test5$estimate)
FamilyKnKs_GenLen = c(test6$p.value,test6$estimate)

FamilyFinal = data.frame(FamilyGrantham_GenLen,FamilyDnDs_GenLen,FamilyKnKs_GenLen); rownames(FamilyFinal) = c('p-value','rho')

agg = merge(agg,OrderFreq)
agg = agg[order(agg$GenerationLength_d),]##"Eulipotyphla","Rodentia","Didelphimorphia","Lagomorpha","Chiroptera","Carnivora","Cetartiodactyla","Primates"
agg$AverageGrantham = round(agg$AverageGrantham,2)
agg$GenerationLength_d = round(agg$GenerationLength_d,0)

boxplot(Mammalia[Mammalia$specieOrder == 'EULIPOTYPHLA',]$AverageGrantham,
        Mammalia[Mammalia$specieOrder == 'RODENTIA',]$AverageGrantham,
        Mammalia[Mammalia$specieOrder == 'DIDELPHIMORPHIA',]$AverageGrantham,
        Mammalia[Mammalia$specieOrder == 'LAGOMORPHA',]$AverageGrantham,
        Mammalia[Mammalia$specieOrder == 'CHIROPTERA',]$AverageGrantham,
        Mammalia[Mammalia$specieOrder == 'CARNIVORA',]$AverageGrantham,
        Mammalia[Mammalia$specieOrder == 'CETARTIODACTYLA',]$AverageGrantham,
        Mammalia[Mammalia$specieOrder == 'PRIMATES',]$AverageGrantham,
        names=c("Eulipotyphla","Rodentia","Didelphimorphia","Lagomorpha","Chiroptera","Carnivora","Cetartiodactyla","Primates"), outline = FALSE, notch = TRUE,  xlab = "GenLenght", ylab = "AverageGrantham")

boxplot(Mammalia[Mammalia$specieOrder == 'EULIPOTYPHLA',]$DnDs,
        Mammalia[Mammalia$specieOrder == 'RODENTIA',]$DnDs,
        Mammalia[Mammalia$specieOrder == 'DIDELPHIMORPHIA',]$DnDs,
        Mammalia[Mammalia$specieOrder == 'LAGOMORPHA',]$DnDs,
        Mammalia[Mammalia$specieOrder == 'CHIROPTERA',]$DnDs,
        Mammalia[Mammalia$specieOrder == 'CARNIVORA',]$DnDs,
        Mammalia[Mammalia$specieOrder == 'CETARTIODACTYLA',]$DnDs,
        Mammalia[Mammalia$specieOrder == 'PRIMATES',]$DnDs,
        names=c("Eulipotyphla","Rodentia","Didelphimorphia","Lagomorpha","Chiroptera","Carnivora","Cetartiodactyla","Primates"), outline = FALSE, notch = TRUE, xlab = "GenLenght", ylab = "DnDs")

boxplot(Mammalia[Mammalia$specieOrder == 'EULIPOTYPHLA',]$KnKs,  
        Mammalia[Mammalia$specieOrder == 'RODENTIA',]$KnKs,
        Mammalia[Mammalia$specieOrder == 'DIDELPHIMORPHIA',]$KnKs,
        Mammalia[Mammalia$specieOrder == 'LAGOMORPHA',]$KnKs,
        Mammalia[Mammalia$specieOrder == 'CHIROPTERA',]$KnKs,
        Mammalia[Mammalia$specieOrder == 'CARNIVORA',]$KnKs,
        Mammalia[Mammalia$specieOrder == 'CETARTIODACTYLA',]$KnKs,
        Mammalia[Mammalia$specieOrder == 'PRIMATES',]$KnKs,
        names=c("Eulipotyphla","Rodentia","Didelphimorphia","Lagomorpha","Chiroptera","Carnivora","Cetartiodactyla","Primates"), outline = FALSE, notch = TRUE, xlab = "GenLenght", ylab = "KnKs")

Rodenta = Mammalia[Mammalia$specieOrder == 'Rodentia',]# надо проверить на домашних


# Eulipotyphla 34 и 1
# Rodentia 98 LC и 3 
# Chiroptera 62 LC и 2 
# Carnivora 30 LC и 2 
# Cetartiodactyla 19 и 13 
# Primates  18 и 28 

###### red book попытка коррелировать красную книгу по таксонам

MammaliaRB = data_redbook[data_redbook$specieOrder %in% FrequentOrders,]
aggRB = aggregate(list(MammaliaRB$AverageGrantham,MammaliaRB$DnDs,MammaliaRB$KnKs), by = list(MammaliaRB$specieOrder), FUN = median)
names(aggRB) = c('Order','AverageGrantham','DnDs','KnKs')

NuOfSpecies0 = c(1,2,3,4,5,6,7,8,9,10)
NuOfSpecies1 = c(1:10)
NuOfSpecies = as.data.frame(NuOfSpecies0,NuOfSpecies1)

#for (i in 1:length(unique(MammaliaRB$specieOrder))){
  #i = 4
#NuOfSpecies$NuOfSpecies0[i] = as.numeric(nrow(Mammalia[MammaliaRB$specieOrder == MammaliaRB$specieOrder[i]& MammaliaRB$CategoryBinomial == 0,]))
#NuOfSpecies$NuOfSpecies1[i] = as.numeric(nrow(Mammalia[MammaliaRB$specieOrder == MammaliaRB$specieOrder[i]& MammaliaRB$CategoryBinomial == 1,]))
#}

nrow(Mammalia[MammaliaRB$CategoryBinomial == 1,])#67 краснокнижных из 334

RODENTIA0 = nrow(Mammalia[MammaliaRB$specieOrder == 'RODENTIA'& MammaliaRB$CategoryBinomial == 0,])# 99
RODENTIA1 = nrow(Mammalia[MammaliaRB$specieOrder == 'RODENTIA'& MammaliaRB$CategoryBinomial == 1,])# 4
DIDELPHIMORPHIA0 = nrow(Mammalia[MammaliaRB$specieOrder == 'DIDELPHIMORPHIA'& MammaliaRB$CategoryBinomial == 0,])# 9
DIDELPHIMORPHIA1 = nrow(Mammalia[MammaliaRB$specieOrder == 'DIDELPHIMORPHIA'& MammaliaRB$CategoryBinomial == 1,])# 0
LAGOMORPHA0 = nrow(Mammalia[MammaliaRB$specieOrder == 'LAGOMORPHA'& MammaliaRB$CategoryBinomial == 0,])# 15
LAGOMORPHA1 = nrow(Mammalia[MammaliaRB$specieOrder == 'LAGOMORPHA'& MammaliaRB$CategoryBinomial == 1,])# 0
CHIROPTERA0 = nrow(Mammalia[MammaliaRB$specieOrder == 'CHIROPTERA'& MammaliaRB$CategoryBinomial ==0,])# 62
CHIROPTERA1 = nrow(Mammalia[MammaliaRB$specieOrder == 'CHIROPTERA'& MammaliaRB$CategoryBinomial == 1,])# 2
CARNIVORA0 = nrow(Mammalia[MammaliaRB$specieOrder == 'CARNIVORA'& MammaliaRB$CategoryBinomial == 0,])# 30
CARNIVORA1 = nrow(Mammalia[MammaliaRB$specieOrder == 'CARNIVORA'& MammaliaRB$CategoryBinomial == 1,])# 9
CETARTIODACTYLA0 = nrow(Mammalia[MammaliaRB$specieOrder == 'CETARTIODACTYLA'& MammaliaRB$CategoryBinomial == 0,])# 19
CETARTIODACTYLA1 = nrow(Mammalia[MammaliaRB$specieOrder == 'CETARTIODACTYLA'& MammaliaRB$CategoryBinomial == 1,])# 13
PRIMATES0 = nrow(Mammalia[MammaliaRB$specieOrder == 'PRIMATES'& MammaliaRB$CategoryBinomial == 0,])# 18
PRIMATES1 = nrow(Mammalia[MammaliaRB$specieOrder == 'PRIMATES'& MammaliaRB$CategoryBinomial == 1,])# 28

# Оставляю только RODENTIA,CHIROPTERA,CARNIVORA,CETARTIODACTYLA,PRIMATES тк у других нет видов под угрозой

AverageGrantham_Poisson_OrderRodentia = glm(CategoryPoisson ~ scale(AverageGrantham) + scale(DnDs), family = 'poisson', data = data_redbook[data_redbook$Order == 'Rodentia'& data_redbook$CategoryBinomial == 0,]); 
summary(AverageGrantham_Poisson_OrderRodentia)                                    


AverageGrantham_Poisson_OrderRodentia = glm(CategoryPoisson ~ scale(AverageGrantham) + scale(DnDs) + scale(KnKs), family = 'poisson', data = Mammalia[Mammalia$Order == 'Rodentia',]); 
summary(AverageGrantham_Poisson_OrderRodentia)                                    
#(Intercept)             -2.9199     0.4439  -6.578 4.76e-11 ***
#scale(AverageGrantham)   0.3904     0.2796   1.396  0.16261    
#scale(DnDs)             -0.2204     0.4895  -0.450  0.65257    
#scale(KnKs)              0.4172     0.1607   2.595  0.00945 **  

AverageGrantham_Poisson_OrderEulipotyphla = glm(CategoryPoisson ~ scale(AverageGrantham)+ scale(DnDs) + scale(KnKs), family = 'poisson', data = Mammalia[Mammalia$Order == 'Eulipotyphla',]); 
summary(AverageGrantham_Poisson_OrderEulipotyphla)
#(Intercept)                       -4.378      1.778  -2.463   0.0138 *
#scale(log2(GenerationLength_d))    1.211      1.119   1.082   0.2791  
#scale(AverageGrantham)            -0.489      1.442  -0.339   0.7346 


AverageGrantham_Poisson_OrderChiroptera = glm(CategoryPoisson ~ scale(AverageGrantham)+ scale(DnDs) + scale(KnKs), family = 'poisson', data = Mammalia[Mammalia$Order == 'Chiroptera',]); 
summary(AverageGrantham_Poisson_OrderChiroptera)
#(Intercept)            -4.13573    1.10429  -3.745  0.00018 ***
#scale(AverageGrantham)  0.27841    0.68179   0.408  0.68302    
#scale(DnDs)             0.77312    0.34619   2.233  0.02554 *  
#scale(KnKs)            -0.04689    0.80572  -0.058  0.95360 

AverageGrantham_Poisson_OrderCarnivora = glm(CategoryPoisson ~ scale(AverageGrantham)+ scale(DnDs) + scale(KnKs), family = 'poisson', data = Mammalia[Mammalia$Order == 'Carnivora',]); 
summary(AverageGrantham_Poisson_OrderCarnivora)
#(Intercept)             -1.2306     0.3227  -3.814 0.000137 ***
#scale(AverageGrantham)   0.0408     0.3486   0.117 0.906842    
#scale(DnDs)              0.2777     0.1958   1.418 0.156061    
#scale(KnKs)              0.3413     0.3748   0.911 0.362512 

AverageGrantham_Poisson_OrderCetartiodactyla = glm(CategoryPoisson ~ scale(AverageGrantham)+ scale(DnDs) + scale(KnKs), family = 'poisson', data = Mammalia[Mammalia$Order == 'Cetartiodactyla',]); 
summary(AverageGrantham_Poisson_OrderCetartiodactyla)
#(Intercept)            -0.820115   0.285253  -2.875  0.00404 **
#scale(AverageGrantham)  0.009696   0.290208   0.033  0.97335   
#scale(DnDs)             0.267461   0.207082   1.292  0.19651   
#scale(KnKs)             0.083856   0.411350   0.204  0.83847

AverageGrantham_Poisson_OrderPrimates = glm(CategoryPoisson ~ scale(AverageGrantham)+ scale(DnDs) + scale(KnKs), family = 'poisson', data = Mammalia[Mammalia$Order == 'Primates',]); 
summary(AverageGrantham_Poisson_OrderPrimates)
#(Intercept)             -0.1045     0.1699  -0.615    0.538  
#scale(AverageGrantham)   0.1521     0.1287   1.182    0.237  
#scale(DnDs)             -0.6157     0.2422  -2.542    0.011 *
#scale(KnKs)              0.1968     0.1677   1.174    0.240



AverageGrantham_Binomial_OrderRodentia = glm(CategoryBinomial ~ scale(AverageGrantham)+ scale(DnDs) + scale(KnKs), family = 'binomial', data = Mammalia[Mammalia$Order == 'Rodentia',]); 
summary(AverageGrantham_Binomial_OrderRodentia)                                    
#(Intercept)             -3.9008     0.7488  -5.210 1.89e-07 ***
#scale(AverageGrantham)   0.5796     0.4095   1.415   0.1570    
#scale(DnDs)             -0.2758     0.7802  -0.354   0.7237    
#scale(KnKs)              0.5596     0.3068   1.824   0.0682  

AverageGrantham_Binomial_OrderEulipotyphla = glm(CategoryBinomial ~ scale(AverageGrantham)+ scale(DnDs) + scale(KnKs), family = 'binomial', data = Mammalia[Mammalia$Order == 'Eulipotyphla',]); 
summary(AverageGrantham_Binomial_OrderEulipotyphla)
#(Intercept)            -2.557e+01  3.651e+04  -0.001    0.999
#scale(AverageGrantham)  3.057e-15  4.007e+04   0.000    1.000
#scale(DnDs)            -1.057e-14  3.829e+04   0.000    1.000
#scale(KnKs)            -4.455e-16  3.901e+04   0.000    1.000


AverageGrantham_Binomial_OrderChiroptera = glm(CategoryBinomial ~ scale(AverageGrantham)+ scale(DnDs) + scale(KnKs), family = 'binomial', data = Mammalia[Mammalia$Order == 'Chiroptera',]); 
summary(AverageGrantham_Binomial_OrderChiroptera)
#(Intercept)            -2.657e+01  4.638e+04  -0.001        1
#scale(AverageGrantham) -1.501e-14  4.704e+04   0.000        1
#scale(DnDs)            -8.452e-15  4.686e+04   0.000        1
#scale(KnKs)            -7.554e-15  4.727e+04   0.000        1

AverageGrantham_Binomial_OrderCarnivora = glm(CategoryBinomial ~ scale(AverageGrantham)+ scale(DnDs) + scale(KnKs), family = 'binomial', data = Mammalia[Mammalia$Order == 'Carnivora',]); 
summary(AverageGrantham_Binomial_OrderCarnivora)
#(Intercept)             -2.9063     0.8048  -3.611 0.000305 ***
#scale(AverageGrantham)   0.5294     0.8417   0.629 0.529369    
#scale(DnDs)              0.1738     0.5327   0.326 0.744213    
#scale(KnKs)              0.4604     0.8986   0.512 0.608388 

AverageGrantham_Binomial_OrderCetartiodactyla = glm(CategoryBinomial ~ scale(AverageGrantham)+ scale(DnDs) + scale(KnKs), family ='binomial', data = Mammalia[Mammalia$Order == 'Cetartiodactyla',]); 
summary(AverageGrantham_Binomial_OrderCetartiodactyla)
#(Intercept)            -3.52688    1.25132  -2.819  0.00482 **
#scale(AverageGrantham) -0.65952    1.27739  -0.516  0.60564   
#scale(DnDs)            -0.02166    0.94350  -0.023  0.98168   
#scale(KnKs)             0.73670    1.48506   0.496  0.61984

AverageGrantham_Binomial_OrderPrimates = glm(CategoryBinomial ~ scale(AverageGrantham)+ scale(DnDs) + scale(KnKs), family = 'binomial', data = Mammalia[Mammalia$Order == 'Primates',]); 
summary(AverageGrantham_Binomial_OrderPrimates)
#(Intercept)             -0.9540     0.3972  -2.402   0.0163 *
#scale(AverageGrantham)   0.4140     0.3731   1.110   0.2671  
#scale(DnDs)             -1.3843     0.6614  -2.093   0.0364 *
#scale(KnKs)              0.7904     0.4160   1.900   0.0574 




model1KP.KnKs_Binomial = glm(CategoryBinomial ~ scale(log2(GenerationLength_d)) + scale(KnKs) , family = 'binomial', data = Mammalia); 
summary(model1KP.KnKs_Binomial)
#(Intercept)                      -3.4742     0.3868  -8.981  < 2e-16 ***
#scale(log2(GenerationLength_d))   1.1842     0.3121   3.794 0.000148 ***
#scale(KnKs)                       0.4957     0.2062   2.404 0.016222 

model1KP.DnDs_Poisson = glm(CategoryPoisson ~ scale(log2(GenerationLength_d)) + scale(DnDs), family = 'poisson', data = Mammalia); 
summary(model1KP.DnDs_Poisson)
#(Intercept)                     -1.89801    0.16210 -11.709  < 2e-16 ***
#scale(log2(GenerationLength_d))  0.96497    0.12850   7.510 5.93e-14 ***
#scale(DnDs)                      0.09587    0.08278   1.158    0.247

model1KP.DnDs_Binomial = glm(CategoryBinomial ~ scale(log2(GenerationLength_d)) + scale(DnDs), family = 'binomial', data = Mammalia); 
summary(model1KP.DnDs_Binomial)

model3KP.AverageGrantham_Poisson = glm(CategoryPoisson ~ scale(log2(GenerationLength_d)) + scale(AverageGrantham), family = 'poisson', data = Mammalia); # model1 = glm(Data$KnKs ~ log(Data$GenerationLength_d) + IucnRanks, family = 'poisson', data = Data); 
summary(model3KP.AverageGrantham_Poisson)
#(Intercept)                     -3.40938    0.37750  -9.031  < 2e-16 ***
#scale(log2(GenerationLength_d))  1.36380    0.31965   4.267 1.99e-05 ***
#scale(DnDs)                      0.02696    0.19494   0.138     0.89 

model3KP.AverageGrantham_Binomial = glm(CategoryBinomial ~ scale(log2(GenerationLength_d)) + scale(AverageGrantham), family = 'binomial', data = Mammalia); # model1 = glm(Data$KnKs ~ log(Data$GenerationLength_d) + IucnRanks, family = 'poisson', data = Data); 
summary(model3KP.AverageGrantham_Binomial)


### попытка смэтчить не вручную
data("lalonde")
m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde)

m.out1
summary(m.out1)
####


data = data[is.na(data$CategoryBinomial)==F,]
data_m = matchit(CategoryBinomial ~ AverageGrantham + DnDs, data = data)
summary(data_m)

Mammalia = Mammalia[is.na(Mammalia$CategoryBinomial)==F,]
Mammalia_m = matchit(CategoryBinomial ~ AverageGrantham + DnDs + family, data = Mammalia)
summary(Mammalia_m)

library(boxplotdbl) # install.packages('boxplotdbl')
X = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$AverageGrantham)
Y = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$GenerationLength_d)
par(mar = c(4, 4, 4, 4))
boxplotdou(Y,X,ylim = c(0,100), xlim = c(0,10000), name.on.axis = FALSE, cex = 1, pch = 0, cex.lab = 1, cex.axis = 1, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

X = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$DnDs)
Y = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$GenerationLength_d)
par(mar = c(4, 4, 4, 4))
boxplotdou(Y,X, xlim = c(0,8500), name.on.axis = FALSE, cex = 1, pch = 0, cex.lab = 1, cex.axis = 1, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

# we can see that the strong increasing effect is starting from the animals with > 1000 days of generation length
dev.off()
