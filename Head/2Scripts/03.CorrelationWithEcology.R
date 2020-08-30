rm(list=ls(all=TRUE))

file = read.table("/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Grantham.txt", sep='\t', header = TRUE) 

    string.counter<-function(strings, pattern){  
  counts<-NULL
  for(i in 1:length(strings)){
    counts[i]<-length(attr(gregexpr(pattern,strings[i])[[1]], "match.length")[attr(gregexpr(pattern,strings[i])[[1]], "match.length")>0])
  }
  return(counts)
}

################################### 
  for(i in 1:nrow(file)){
  vec = unlist(strsplit(file$NewDistance[i],';'))
  N <- c()
  Output <- c()
  
  if (length(vec)==0){N=N+1;next
  }else{
    
  for (j in 1:length(vec)){
    
    if (string.counter(vec[j], pattern=">")) {Output = c(Output, vec[j])
    } 
  }
  file[i, "only_distances"] <- paste(Output, collapse=";")
  } 
}
#################################### 

for(i in 1:nrow(file)){
  vec = unlist(strsplit(file$only_distances[i],';'))
    if (length(vec)==0){N=N+1;next
  }else{
    
Distances <- c()
  for (j in 1:length(vec)){
AA_with_distance = vec[j]
Distance = gsub("(.*)\\:",'',AA_with_distance)# оставляет все после знака :
Distances = c(Distances, Distance)
 
  }
  file[i, "distances"] <- paste(Distances, collapse=";") 
}
}
################################################
for(i in 1:nrow(file)){
  
  vec = unlist(strsplit(file$NewDistance[i],';'))
  Number <- c()
  if (length(vec)==0){N=N+1;next
  }else{
    
  for (j in 1:length(vec)){

    if (string.counter(vec[j], pattern="Non_Synon:")) {
      Number = gsub("(.*)\\:",'',vec[j])
      }
    file[i, "All_Non_Synon"] <- paste(Number, collapse=";")
  }
 }
}
###################################################

for(i in 1:nrow(file)){
  vec = unlist(strsplit(file$distances[i],';'))
  if (length(vec)==0){N=N+1;next
  }else{
    
  Distances_mean <- 0
  for (j in 1:length(vec)){
    num = as.numeric(vec[j])
    Distances_mean = (Distances_mean + num)
    
  }
  file[i, "distances_mean"] <- paste((as.numeric(Distances_mean/length(vec))), collapse=";") 
}
}

For_each <- aggregate(as.numeric(file[,11]),list(file$Species),mean,na.rm=TRUE)


write.table(For_each,file = "/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Distance_for_each_species.csv",quote = F, row.names = FALSE,
            sep = '\t')



library(ape)
library(gdata)
library(ggplot2)

data = read.table("/home/anastasia/mtdna-mammalian-evolution/Body/2Derived/Distance_for_each_species.csv", sep='\t', header = TRUE) 
GenLength <- read.xls("/home/anastasia/mtdna-mammalian-evolution/Body/1Raw/GenerationLengthForMammals.xlsx")# табличка с продолжительностью жизни от Алины

merge(data, GenLength, by='Species')

data$Group.1 = sub("_", " ", data$Group.1, ignore.case = FALSE, perl = FALSE,
    fixed = FALSE, useBytes = FALSE)

colnames(data) <- c("Species","Distance") 

NewData = merge(data,GenLength,
      by.x = "Species", by.y = "Scientific_name",all = FALSE,no.dups = TRUE,) 


cor.test(x = NewData$GenerationLength_d, y = NewData$Distance, method = "spearm")

fit  <- cor.test(x = NewData$GenerationLength_d, y = NewData$Distance)

str(fit)

fit$p.value

plot(x = NewData$GenerationLength_d, y = NewData$Distance)

ggplot(NewData, aes(x = log(GenerationLength_d), y = Distance,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(NewData, aes(x = log(AdultBodyMass_g), y = Distance,col = factor(Order) ))+
  geom_point(size = 2)

ggplot(NewData, aes(x = log(GenerationLength_d), y = Distance)) + 
  geom_point() +
  geom_smooth()

ggplot(NewData, aes(x = log(AdultBodyMass_g), y = Distance)) + 
  geom_point() +
  geom_smooth()








































