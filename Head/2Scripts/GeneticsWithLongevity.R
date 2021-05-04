#4 Генетика ~ Продолжительность жизни

rm(list=ls(all=TRUE))

#wd = getwd()
#wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
#setwd(wd)

library(ggplot2)
data = read.table("../../Body/2Derived/CytB_with_ecology.csv", sep='\t', header = TRUE) 

###### КОРЕЛЯЦИИ

#(i) FractionOfSyn ~ -GenerLength; 

cor.test(x = data$GenerationLength_d, y = data$FractionOfSyn , method = "spearm")
GenerationLength_FractionOfSyn_fit  <- cor.test(x = data$GenerationLength_d, y = data$FractionOfSyn, method = "spearm")

ggplot(data, aes(x = log(GenerationLength_d), y = FractionOfSyn ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(data, aes(x = log(GenerationLength_d), y = FractionOfSyn , fill = Order)) + 
  geom_boxplot()

#(ii) FractionOfNonsyn ~ +GenerLength; 

cor.test(x = data$GenerationLength_d, y = data$FractionOfNonsyn , method = "spearm")
GenerationLength_FractionOfNonsyn_fit  <- cor.test(x = data$GenerationLength_d, y = data$FractionOfNonsyn, method = "spearm")

ggplot(data, aes(x = log(GenerationLength_d), y = FractionOfNonsyn ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(data, aes(x = log(GenerationLength_d), y = FractionOfNonsyn , fill = Order)) + 
  geom_boxplot()

#(iii) SummOfAllGrantham ~ +GenerLength; 
cor.test(x = data$GenerationLength_d, y = data$SummOfAllGrantham , method = "spearm")
GenerationLength_SummOfAllGrantham_fit  <- cor.test(x = data$GenerationLength_d, y = data$SummOfAllGrantham, method = "spearm")

ggplot(data, aes(x = log(GenerationLength_d), y = SummOfAllGrantham ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(data, aes(x = log(GenerationLength_d), y = SummOfAllGrantham , fill = Order)) + 
  geom_boxplot()

#(iv) MedianOfAllNonsyn ~ +GenerLength 
cor.test(x = data$GenerationLength_d, y = data$MedianOfAllNonsyn , method = "spearm")
GenerationLength_MedianOfAllNonsyn_fit  <- cor.test(x = data$GenerationLength_d, y = data$MedianOfAllNonsyn, method = "spearm")

ggplot(data, aes(x = log(GenerationLength_d), y = MedianOfAllNonsyn ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(data, aes(x = log(GenerationLength_d), y =MedianOfAllNonsyn , fill = Order)) + 
  geom_boxplot()

#(v) SummOfAllGrantham ~+SummOfAllNonsyn (проверка на вшивость)
cor.test(x = data$SummOfAllGrantham, y = data$MedianOfAllNonsyn, method = "spearm")
GenerationLength_MedianOfAllNonsyn_Grantham_fit  <- cor.test(x = data$SummOfAllGrantham, y = data$MedianOfAllNonsyn, method = "spearm")

ggplot(data, aes(x = log(SummOfAllGrantham), y = MedianOfAllNonsyn ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(data, aes(x = log(SummOfAllGrantham), y = MedianOfAllNonsyn , fill = Order)) + 
  geom_boxplot()

#(vii) AverageGrantham ~ +GenerLength; 
cor.test(x = data$GenerationLength_d, y = data$AverageGrantham , method = "spearm")
GenerationLength_AverageGrantham_fit  <- cor.test(x = data$GenerationLength_d, y = data$AverageGrantham, method = "spearm")

ggplot(data, aes(x = log(GenerationLength_d), y = AverageGrantham ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(data, aes(x = log(GenerationLength_d), y = AverageGrantham, fill = Order)) + 
  geom_boxplot()

#(viii) KnKs ~ +GenerLength; 
cor.test(x = data$GenerationLength_d, y = data$KnKs , method = "spearm")
GenerationLength_KnKs_fit  <- cor.test(x = data$GenerationLength_d, y = data$KnKs, method = "spearm")

ggplot(data, aes(x = log(GenerationLength_d), y = KnKs ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(data, aes(x = log(GenerationLength_d), y = KnKs, fill = Order)) + 
  geom_boxplot()

#(ix) DnDs ~ +GenerLength;
cor.test(x = data$GenerationLength_d, y = data$DnDs , method = "spearm")
GenerationLength_DnDs_fit  <- cor.test(x = data$GenerationLength_d, y = data$DnDs, method = "spearm")

ggplot(data, aes(x = log(GenerationLength_d), y = DnDs ,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(data, aes(x = log(GenerationLength_d), y = DnDs, fill = Order)) + 
  geom_boxplot()

estimate = data.frame(GenerationLength_FractionOfSyn_fit$estimate,GenerationLength_FractionOfNonsyn_fit$estimate,GenerationLength_SummOfAllGrantham_fit$estimate,
                      GenerationLength_MedianOfAllNonsyn_fit$estimate,GenerationLength_MedianOfAllNonsyn_Grantham_fit$estimate,GenerationLength_AverageGrantham_fit$estimate,
                      GenerationLength_KnKs_fit$estimate,GenerationLength_DnDs_fit$estimate);names(estimate) = c('FractionOfSyn','FractionOfNonsyni','SummOfAllGrantham','MedianOfAllNonsyn','MedianOfAllNonsyn_Grantham','AverageGrantham','KnKs','DnDs')
pvalue = data.frame(GenerationLength_FractionOfSyn_fit$p.value,GenerationLength_FractionOfNonsyn_fit$p.value,GenerationLength_SummOfAllGrantham_fit$p.value,
                    GenerationLength_MedianOfAllNonsyn_fit$p.value,GenerationLength_MedianOfAllNonsyn_fit$p.value,GenerationLength_AverageGrantham_fit$p.value,
                    GenerationLength_KnKs_fit$p.value,GenerationLength_DnDs_fit$p.value);names(pvalue) = c('FractionOfSyn','FractionOfNonsyni','SummOfAllGrantham','MedianOfAllNonsyn','MedianOfAllNonsyn_Grantham','AverageGrantham','KnKs','DnDs')
  
fit = rbind (estimate,pvalue); rownames(fit)= c('estimate', 'p-value') ########## Табличка с результатами корреляций

##### boxplots by quartiles

Data = data
quantile(Data$GenerationLength_d, seq(from=0, to=1, by=0.25),na.rm = TRUE)  

boxplot(Data[Data$GenerationLength_d < 637.3329,]$AverageGrantham, +
          Data[Data$GenerationLength_d < 1825.0000 & Data$GenerationLength_d > 637.3329,]$AverageGrantham, +
          Data[Data$GenerationLength_d < 2869.6933 & Data$GenerationLength_d > 1825.0000,]$AverageGrantham, +
          Data[Data$GenerationLength_d > 2869.6933,]$AverageGrantham, names = c('Q1','Q2', 'Q3','Q4'), outline = FALSE, xlab = "GenLenght", ylab = "AverageGrantham",notch = TRUE)

boxplot(Data[Data$GenerationLength_d < 637.3329,]$KnKs, +
          Data[Data$GenerationLength_d < 1825.0000 & Data$GenerationLength_d > 637.3329,]$KnKs, +
          Data[Data$GenerationLength_d < 2869.6933 & Data$GenerationLength_d > 1825.0000,]$KnKs, +
          Data[Data$GenerationLength_d > 2869.6933,]$KnKs, names = c('Q1','Q2', 'Q3','Q4'), outline = FALSE, xlab = "GenLenght", ylab = "KnKs",notch = TRUE)


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

##################### ПИКИ

library(ape)
library(geiger)
library(caper)
library(stringr)
library(phytools)
library(picante)

tree = read.tree("../../Body/1Raw/mammalia_species.nwk") 
data$Species = sub(" ", "_", data$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

row.names(data) = data$Species

nrow(data) # 344
tree_w = treedata(tree, data[, c('Species', 'AverageGrantham', 'KnKs','DnDs' ,'GenerationLength_d')], 
                  sort=T, warnings=T)$phy
data_tree = as.data.frame(treedata(tree_w, data[, c('Species', 'AverageGrantham', 'KnKs','DnDs' ,'GenerationLength_d')], 
                                   sort=T, warnings=T)$data)

setdiff(data$Species, data_tree$Species)# 4 вида отвалились при присоединении дерева

data_tree$AverageGrantham = as.numeric(as.character(data_tree$AverageGrantham))
data_tree$GenerationLength_d = as.numeric(as.character(data_tree$GenerationLength_d))
data_tree$KnKs = as.numeric(as.character(data_tree$KnKs))

#GenLen & Grantham
GranthamPIC = cor.test(pic(data_tree$AverageGrantham, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w), method = 'spearman')
Grantham = cor.test(data_tree$AverageGrantham,data_tree$GenerationLength_d, method = 'spearman')

#GenLen & DnDs
DnDsPIC = cor.test(pic(data_tree$DnDs, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w), method = 'spearman')
DnDs = cor.test(as.numeric(data_tree$DnDs),data_tree$GenerationLength_d, method = 'spearman')

#GenLen & KnKs
KnKsPIC = cor.test(pic(data_tree$KnKs, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w), method = 'spearman')
KnKs = cor.test(data_tree$KnKs,data_tree$GenerationLength_d, method = 'spearman')

estimate = data.frame(GranthamPIC$estimate,Grantham$estimate,DnDsPIC$estimate,
                      DnDs$estimate,KnKsPIC$estimate,KnKs$estimate);names(estimate) = c('GranthamPIC','Grantham','DnDsPIC','DnDs','KnKsPIC','KnKs')
pvalue = data.frame(GranthamPIC$p.value,Grantham$p.value,DnDsPIC$p.value,
                    DnDs$p.value,KnKsPIC$p.value,KnKs$p.value);names(pvalue) = c('GranthamPIC','Grantham','DnDsPIC','DnDs','KnKsPIC','KnKs')

fit_PICS = rbind (estimate,pvalue); rownames(fit)= c('estimate', 'p-value') ############ Табличка с результатами Пиков

######### МАНН-УИТНИ


wilcox.test(data$AverageGrantham,data$GenerationLength_d)
wilcox.test(data$DnDs,data$GenerationLength_d)
wilcox.test(data$KnKs,data$GenerationLength_d)
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

hist(data$GenerationLength_d,breaks=10)

######### ЛЯМБДА

GenLen <- data_tree$GenerationLength_d
names(GenLen) <- rownames(data_tree)

Grantham <- data_tree$AverageGrantham
names(Grantham) <- rownames(data_tree)

DnDs <- as.matrix(as.numeric(data_tree$DnDs))
rownames(DnDs) <- rownames(data_tree)

phylosig(tree_w, GenLen, method = "lambda", test = TRUE) 
phylosig(tree_w, Grantham, method = "lambda", test = TRUE)
phylosig(tree_w, DnDs, method = "lambda", test = TRUE)

################### PGLS

MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)

AverageGrantham_GenerationLength_d = pgls(scale(AverageGrantham) ~ log2(GenerationLength_d), MutComp, lambda="ML")
summary(pgls(scale(AverageGrantham) ~ log2(GenerationLength_d), MutComp, lambda="ML"))

KnKs_GenerationLength_d= pgls(scale(KnKs) ~ log2(GenerationLength_d), MutComp, lambda="ML")
summary(pgls(scale(KnKs) ~ log2(GenerationLength_d), MutComp, lambda="ML"))

DnDs_GenerationLength_d= pgls(scale(DnDs) ~ log2(GenerationLength_d), MutComp, lambda="ML")
summary(pgls(scale(DnDs) ~ log2(GenerationLength_d), MutComp, lambda="ML"))

Intercept = data.frame(AverageGrantham_GenerationLength_d$sterr[[1]],KnKs_GenerationLength_d$sterr[[1]],DnDs_GenerationLength_d$sterr[[1]]);names(Intercept) = c('AverageGrantham_GenerationLength_d','KnKs_GenerationLength_d','DnDs_GenerationLength_d')
GenerationLength_d = data.frame(AverageGrantham_GenerationLength_d$sterr[[2]],KnKs_GenerationLength_d$sterr[[2]],DnDs_GenerationLength_d$sterr[[2]]);names(GenerationLength_d) = c('AverageGrantham_GenerationLength_d','KnKs_GenerationLength_d','DnDs_GenerationLength_d')

Intercept = c(AverageGrantham_GenerationLength_d$sterr[1],KnKs_GenerationLength_d$sterr[1],DnDs_GenerationLength_d$sterr[1]);
names(Intercept) = c('AverageGrantham_GenerationLength_d','KnKs_GenerationLength_d','DnDs_GenerationLength_d')

PGLS= rbind (Intercept,GenerationLength_d); rownames(PGLS)= c('Intercept', 'log2(GenerationLength_d)') ############ Табличка с результатами PGLS

######################## ДЕЛЕНИЕ ПО СЕМЕЙСТВАМ

OrderFreq = data.frame(table(data$Order));
FrequentOrders = OrderFreq[OrderFreq$Freq >= 3,]$Var1; length(FrequentOrders) # 3!!!
OrderFreq = OrderFreq[OrderFreq$Var1 %in% FrequentOrders,]
names(OrderFreq)=c('Order','NumberOfSpecies')

Mammalia = data[data$Order %in% FrequentOrders,]
agg = aggregate(list(Mammalia$AverageGrantham,Mammalia$GenerationLength_d,Mammalia$DnDs,Mammalia$KnKs), by = list(Mammalia$Order), FUN = median)
names(agg) = c('Order','AverageGrantham','GenerationLength_d','DnDs','KnKs')

test1 = cor.test(agg$AverageGrantham,agg$GenerationLength_d,method = 'spearman') 
plot(agg$GenerationLength_d,agg$AverageGrantham)
test2 = cor.test(agg$DnDs,agg$GenerationLength_d,method = 'spearman') 
plot(agg$DnDs,agg$AverageGrantham)
test3 = cor.test(agg$KnKs,agg$GenerationLength_d,method = 'spearman')
plot(agg$KnKs,agg$AverageGrantham)

OrderGrantham_GenLen = c(test1$p.value,test1$estimate)
OrderDnDs_GenLen = c(test2$p.value,test2$estimate)
OrderKnKs_GenLen = c(test3$p.value,test3$estimate)

OrderFinal = data.frame(OrderGrantham_GenLen,OrderDnDs_GenLen,OrderKnKs_GenLen); rownames(OrderFinal) = c('p-value','rho')

FamFreq = data.frame(table(data$family));
FrequentFamilies = FamFreq[FamFreq$Freq >= 3,]$Var1; length(FrequentFamilies) # 3!!!
FamFreq = FamFreq[FamFreq$Var1 %in% FrequentFamilies,]
names(FamFreq)=c('Family','NumberOfSpecies')

# Family
Mammalia1 = data[data$family %in% FrequentFamilies,]
agg1 = aggregate(list(Mammalia1$AverageGrantham,Mammalia1$GenerationLength_d,Mammalia1$DnDs,Mammalia1$KnKs), by = list(Mammalia1$family), FUN = median)
names(agg1) = c('Family','AverageGrantham','GenerationLength_d','DnDs','KnKs')

test4 = cor.test(agg1$AverageGrantham,agg1$GenerationLength_d,method = 'spearman')
plot(agg1$GenerationLength_d,agg1$AverageGrantham)
test5 = cor.test(agg1$DnDs,agg1$GenerationLength_d,method = 'spearman') 
plot(agg1$DnDs,agg1$AverageGrantham)
test6 = cor.test(agg1$KnKs,agg1$GenerationLength_d,method = 'spearman') 
plot(agg1$KnKs,agg1$AverageGrantham)

FamilyGrantham_GenLen = c(test4$p.value,test4$estimate)
FamilyDnDs_GenLen = c(test5$p.value,test5$estimate)
FamilyKnKs_GenLen = c(test6$p.value,test6$estimate)

FamilyFinal = data.frame(FamilyGrantham_GenLen,FamilyDnDs_GenLen,FamilyKnKs_GenLen); rownames(FamilyFinal) = c('p-value','rho')

agg = merge(agg,OrderFreq)
agg = agg[order(agg$GenerationLength_d),]##"Eulipotyphla","Rodentia","Didelphimorphia","Lagomorpha","Chiroptera","Carnivora","Cetartiodactyla","Primates"
agg$AverageGrantham = round(agg$AverageGrantham,2)
agg$GenerationLength_d = round(agg$GenerationLength_d,0)

#################### Боксплоты по отрядам и семействам

#Order=c("Eulipotyphla","Rodentia","Didelphimorphia","Lagomorpha","Chiroptera","Carnivora","Cetartiodactyla","Cetartiodactyla")#    
#median_DnDs = c(median(Data[Data$Order == 'Eulipotyphla',]$DnDs),median(Data[Data$Order == 'Rodentia',]$DnDs),median(Data[Data$Order == 'Didelphimorphia',]$DnDs),
#                median(Data[Data$Order == 'Lagomorpha',]$DnDs),median(Data[Data$Order == 'Chiroptera',]$DnDs),median(Data[Data$Order == 'Carnivora',]$DnDs),
#                median(Data[Data$Order == 'Cetartiodactyla',]$DnDs),median(Data[Data$Order == 'Cetartiodactyla',]$DnDs))
#median_KnKs = c(median(Data[Data$Order == 'Eulipotyphla',]$KnKs),median(Data[Data$Order == 'Rodentia',]$KnKs),median(Data[Data$Order == 'Didelphimorphia',]$KnKs),
#                median(Data[Data$Order == 'Lagomorpha',]$KnKs),median(Data[Data$Order == 'Chiroptera',]$KnKs),median(Data[Data$Order == 'Carnivora',]$KnKs),
#                median(Data[Data$Order == 'Cetartiodactyla',]$KnKs),median(Data[Data$Order == 'Cetartiodactyla',]$KnKs))

#median_GenerationLength_d = c(median(Data[Data$Order == 'Eulipotyphla',]$GenerationLength_d),median(Data[Data$Order == 'Rodentia',]$GenerationLength_d),median(Data[Data$Order == 'Didelphimorphia',]$GenerationLength_d),
#                              median(Data[Data$Order == 'Lagomorpha',]$GenerationLength_d),median(Data[Data$Order == 'Chiroptera',]$GenerationLength_d),median(Data[Data$Order == 'Carnivora',]$GenerationLength_d),
 #                             median(Data[Data$Order == 'Cetartiodactyla',]$GenerationLength_d),median(Data[Data$Order == 'Cetartiodactyla',]$GenerationLength_d))
#median_AverageGrantham = c(median(Data[Data$Order == 'Eulipotyphla',]$AverageGrantham),median(Data[Data$Order == 'Rodentia',]$AverageGrantham),median(Data[Data$Order == 'Didelphimorphia',]$AverageGrantham),
#                           median(Data[Data$Order == 'Lagomorpha',]$AverageGrantham),median(Data[Data$Order == 'Chiroptera',]$AverageGrantham),median(Data[Data$Order == 'Carnivora',]$AverageGrantham),
#                           median(Data[Data$Order == 'Cetartiodactyla',]$AverageGrantham),median(Data[Data$Order == 'Cetartiodactyla',]$AverageGrantham))
#number_of_species = c(length(unique(Data[Data$Order == 'Eulipotyphla',]$Species)),length(unique(Data[Data$Order == 'Rodentia',]$Species)),length(unique(Data[Data$Order == 'Didelphimorphia',]$Species)),
#                      length(unique(Data[Data$Order == 'Lagomorpha',]$Species)),length(unique(Data[Data$Order == 'Chiroptera',]$Species)),length(unique(Data[Data$Order == 'Carnivora',]$Species)),
#                      length(unique(Data[Data$Order == 'Cetartiodactyla',]$Species)),length(unique(Data[Data$Order == 'Cetartiodactyla',]$Species)))

#Genetics_Ecology = data.frame(Order,median_GenerationLength_d,median_AverageGrantham,median_DnDs,median_KnKs,number_of_species) 
#Genetics_Ecology = Genetics_Ecology[order(Genetics_Ecology$median_GenerationLength_d),]# генетика ~ продолжительность жизни по семействам


#All_statistics = c(GenerationLength_FractionOfSyn_fit$estimate,GenerationLength_FractionOfNonsyn_fit$estimate,GenerationLength_SummOfAllGrantham_fit$estimate,
#                   GenerationLength_MedianOfAllNonsyn_fit$estimate,GenerationLength_AverageGrantham_fit$estimate,GenerationLength_KnKs_fit$estimate,GenerationLength_DnDs_fit$estimate)

#All_pvalue = c(GenerationLength_FractionOfSyn_fit$p.value,GenerationLength_FractionOfNonsyn_fit$p.value,GenerationLength_SummOfAllGrantham_fit$p.value,
#               GenerationLength_MedianOfAllNonsyn_fit$p.value,GenerationLength_AverageGrantham_fit$p.value,GenerationLength_KnKs_fit$p.value,GenerationLength_DnDs_fit$p.value)

#Variable = c('FractionOfSyn','FractionOfNonsyn','SummOfAllGrantham','MedianOfAllNonsyn','AverageGrantham','KnKs','DnDs')

#Corr_with_GenLen = data.frame(Variable, All_statistics,All_pvalue) 

boxplot(Data[Data$Order == "Eulipotyphla",]$KnKs,+
          Data[Data$Order == "Rodentia",]$KnKs, +
          Data[Data$Order == "Chiroptera",]$KnKs, +
          Data[Data$Order == "Carnivora",]$KnKs, +
          Data[Data$Order == "Primates",]$KnKs, +
          Data[Data$Order =="Cetartiodactyla",]$KnKs, names = c('Eulipotyphla',"Rodentia",'Chiroptera','Carnivora',"Primates",'Cetartiodactyla'), outline = FALSE, xlab = "GenLenght", ylab = "KnKs",notch = TRUE)

boxplot(Data[Data$Order == "Eulipotyphla",]$AverageGrantham,+
          Data[Data$Order == "Rodentia",]$AverageGrantham, +
          Data[Data$Order == "Chiroptera",]$AverageGrantham, +
          Data[Data$Order == "Carnivora",]$AverageGrantham, +
          Data[Data$Order == "Primates",]$AverageGrantham, +
          Data[Data$Order =="Cetartiodactyla",]$AverageGrantham, names = c('Eulipotyphla',"Rodentia",'Chiroptera','Carnivora',"Primates",'Cetartiodactyla'), outline = FALSE, xlab = "GenLenght", ylab = "AverageGrantham",notch = TRUE)

TempData = subset(Data, Data$Order == 'Eulipotyphla'|Data$Order =="Rodentia"|Data$Order =='Chiroptera'|Data$Order =='Carnivora'|Data$Order =="Primates"|Data$Order =='Cetartiodactyla')

ggplot(TempData, aes(x=Order, y=KnKs,fill=Order),notch=T) + 
  geom_violin(trim=FALSE)+ 
  stat_summary(fun=median, geom="crossbar", size=0.2, color="red")+
  labs(title="KnKs~GenLen",x="Order", y = "KnKs")+
  theme_classic()


ggplot(TempData, aes(x=Order, y=AverageGrantham,fill=Order),notch=T) + 
  geom_violin(trim=FALSE)+ 
  stat_summary(fun=median, geom="crossbar", size=0.2, color="red")+
  labs(title="AverageGrantham~GenLen",x="Order", y = "AverageGrantham")+
  theme_classic()


median(Data[Data$Order == "Primates",]$GenerationLength_d)
median(Data[Data$Order =="Cetartiodactyla",]$GenerationLength_d)# продолжительность жизни китопарнокопытных чуть больше, чем у приматов
# Carnivora Cetartiodactyla Chiroptera Eulipotyphla Primates Rodentia 

nrow(Mammalia)
for (i in 1:nrow(Mammalia))   {  Mammalia$FamilyShort[i] = paste(unlist(strsplit(Mammalia$family[i],''))[c(1:3)],collapse = '')  }
for (i in 1:nrow(Mammalia))   {  Mammalia$OrderShort[i] = paste(unlist(strsplit(Mammalia$Order[i],''))[c(1:3)],collapse = '')  }

library(boxplotdbl) # install.packages('boxplotdbl') #здесь надо разбиение по КК 0 и 1
X = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$AverageGrantham)
Y = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$GenerationLength_d)
par(mar = c(4, 4, 4, 4))
boxplotdou(Y,X,ylim = c(0,100), xlim = c(0,20000), name.on.axis = FALSE, cex = 1, pch = 0, cex.lab = 1, cex.axis = 1, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

X = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$DnDs)
Y = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$GenerationLength_d)
par(mar = c(4, 4, 4, 4))
boxplotdou(Y,X, xlim = c(0,15000), name.on.axis = FALSE, cex = 1, pch = 0, cex.lab = 1, cex.axis = 1, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

X = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$KnKs)
Y = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$GenerationLength_d)
par(mar = c(4, 4, 4, 4))
boxplotdou(Y,X,ylim = c(0,1), xlim = c(0,18500), name.on.axis = FALSE, cex = 1, pch = 0, cex.lab = 1, cex.axis = 1, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

X = data.frame(as.factor(Mammalia$OrderShort),Mammalia$AverageGrantham)
Y = data.frame(as.factor(Mammalia$OrderShort),Mammalia$GenerationLength_d)
par(mar = c(4, 4, 4, 4))
boxplotdou(Y,X, xlim = c(0,6500), name.on.axis = FALSE, cex = 1, pch = 0, cex.lab = 1, cex.axis = 1, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

X = data.frame(as.factor(Mammalia$OrderShort),Mammalia$DnDs)
Y = data.frame(as.factor(Mammalia$OrderShort),Mammalia$GenerationLength_d)
par(mar = c(4, 4, 4, 4))
boxplotdou(Y,X, xlim = c(0,6500), name.on.axis = FALSE, cex = 1, pch = 0, cex.lab = 1, cex.axis = 1, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

X = data.frame(as.factor(Mammalia$OrderShort),Mammalia$KnKs)
Y = data.frame(as.factor(Mammalia$OrderShort),Mammalia$GenerationLength_d)
par(mar = c(4, 4, 4, 4))
boxplotdou(Y,X, xlim = c(0,6500), name.on.axis = FALSE, cex = 1, pch = 0, cex.lab = 1, cex.axis = 1, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

##########################

