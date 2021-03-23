rm(list=ls(all=TRUE)) 

wd = getwd()
wd = paste(wd, '/mtdna-mammalian-evolution/Body/2Derived',sep='')
setwd(wd)

library(ape)
library(geiger)
library(caper)

tree = read.tree("../../Body/1Raw/FcC_supermatrix.part.treefile.txt") 
data = read.csv ("../../Body/2Derived/Distances_KnKs_Ecology_RG.csv", sep = "\t")

data = data[is.na(data$AverageGrantham)==F,]

data$Species = sub(" ", "_", data$Species, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)


data = data [-251,]# удаляем повторяющийся вид
row.names(data) = data$Species

tree_w = treedata(tree, data[, c('Species', 'AverageGrantham', 'KnKs', 'GenerationLength_d')], 
                  sort=T, warnings=T)$phy
data_tree = as.data.frame(treedata(tree_w, data[, c('Species', 'AverageGrantham', 'KnKs', 'GenerationLength_d')], 
                                   sort=T, warnings=T)$data)

setdiff(tree_w$tip.label, data_tree$Species)

data_tree$AverageGrantham = as.numeric(as.character(data_tree$AverageGrantham))
data_tree$GenerationLength_d = as.numeric(as.character(data_tree$GenerationLength_d))
data_tree$KnKs = as.numeric(as.character(data_tree$KnKs))


cor.test(pic(data_tree$AverageGrantham, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w), method = 'spearman')
#rho = -0.003478547 p-value = 0.9566

plot(pic(data_tree$AverageGrantham, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w))

cor.test(pic(data_tree$KnKs, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w), method = 'spearman')
#rho = 0.00416989 p-value = 0.948

plot(pic(data_tree$KnKs, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w))

MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)

model = pgls(scale(GenerationLength_d) ~ scale(AverageGrantham), MutComp, lambda="ML")
summary(model)

plot(pgls(scale(GenerationLength_d) ~ scale(AverageGrantham), MutComp, lambda="ML"))

# lambda [ ML]: 0.930
#(Intercept)            -0.093601   0.227325 -0.4117  0.68088  
#scale(AverageGrantham) -0.052622   0.029922 -1.7586  0.07988 

crunch(scale(GenerationLength_d) ~ scale(AverageGrantham), MutComp)
plot(crunch(scale(GenerationLength_d) ~ scale(AverageGrantham), MutComp))
# p-value = 0.1891
#scale(AverageGrantham)  0.03644    0.02768   1.317    0.189

######################Домашние виды
#Monodelphis_domestica - серый опоссум, экзотическое дом. животное
#Mus_musculus - мышь, а Mus musculus domestica - лабораторная мышь, надо подумать