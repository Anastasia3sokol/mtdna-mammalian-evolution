library(Biostrings)
library(seqinr)
setwd('C:/USers/Professional/Desktop')
rm(list=ls(all=TRUE))


wd = getwd()
wd = paste(wd, '/mtdna-mammalian-evolution/Body/1Raw',sep='')
setwd(wd)

##### read tables for families and generation length ########
families = read.table('MitGenomics.txt',sep='\t', header = TRUE)
families = families[,c(1,37)]

GenerationL = read.table('GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
GenerationL = GenerationL[ ,c(2,11)] # Take information about name and GenLength
GenerationL$Scientific_name = gsub(' ','_',GenerationL$Scientific_name)


GenFam = merge(GenerationL, families, by.x = 'Scientific_name', by.y = 'Species' ) 
GenFam = GenFam[complete.cases(GenFam), ] #### remove rows with NA
GenFam = GenFam[-which(duplicated(GenFam$Scientific_name)),] ### delete repeated species

##### remove families with sppecies < 3 ######
VecOfFamilies = unique(GenFam$ECO.Family)
for (i in VecOfFamilies) {
  b = nrow(GenFam[GenFam$ECO.Family == i,])
  if ( b < 3){
    GenFam = GenFam[GenFam$ECO.Family != i,]
  }
}

VecOfFamilies = unique(GenFam$ECO.Family)### new vector for families with species >= 3
VecOfFamilies = as.character(VecOfFamilies)


wd = paste(wd, '/fasta_codons', sep='')
setwd(wd) ### Change direcotry for Codons


###### read file with sequence CytB #######
Testdf = readDNAStringSet('CytB.fasta')
Testdf = data.frame(Testdf); names(Testdf) = c('Sequence') 
Testdf$Scientific_name  = as.factor(rownames(Testdf))
row.names(Testdf) = c(1:nrow(Testdf))

#### merge dataframe GenLen with dataframe Sequence
MergeDf = merge(GenFam, Testdf, by = 'Scientific_name' )

wd = gsub('1Raw/fasta_codons','2Derived/DifferenceFamiliesCytB',wd) # change directory for our output 
setwd(wd)


####### Create a new data frame with 3 species: High and Low Generation time, and Outgroup####
for (fam in VecOfFamilies) {
  #fam = 'Felidae'
  
  numofloop = 1
  DifferenceInOneFamily = data.frame()
  
  
  OneFamilyDf = MergeDf[MergeDf$ECO.Family == fam,]
  NumberOfSp = nrow(OneFamilyDf)
  
  for (sp1 in 1:(NumberOfSp-2)){
    #sp1 = 1
    Species1 = OneFamilyDf[sp1,]
    
    for (sp2 in (sp1+1):(NumberOfSp-1)){
      #sp = 2
      Species2 = OneFamilyDf[sp2,]
      
      for(sp3 in (sp2+1):NumberOfSp){
        #sp3 = 3
        Species3 = OneFamilyDf[sp3,]
        
        
        ######### Make Sequences readable #########################
        
        VecOfNuc1 = unlist(strsplit(Species1$Sequence,''))
        VecOfNuc2 = unlist(strsplit(Species2$Sequence,''))
        VecOfNuc3 = unlist(strsplit(Species3$Sequence,''))
        
        StartNuc = 1
        
        Codons1 = c()
        Codons2 = c()
        Codons3 = c()
        
        ######## Create vectors for differences ##################
        
        DifBtwSp1AndSp2 = ''
        DifBtwSp1AndSp3 = ''
        DifBtwSp2AndSp3 = ''
        
        ############# Create a new data frame with information about Species and  difference between every species for One Gene ##############################
        
        DifferenceInOneTrio = data.frame(Species1$Scientific_name, Species2$Scientific_name, Species3$Scientific_name, Species1$ECO.Family, Species1$GenerationLength_d,
                                         Species2$GenerationLength_d, Species3$GenerationLength_d, '', '', '', 'CytB') 
        names(DifferenceInOneTrio) = c('Species1', 'Species2', 'Species3', 'Family', 'GenLenSp1', 'GenLenSp2','GenLenSp3', 'AllCodonSubstBetweenSp1AndSp2' ,
                                       'AllCodonSubstBetweenSp1AndSp3', 'AllCodonSubstBetweenSp2AndSp3', 'Gene')
        
        
        ############# Loop for differences ############ 
        
        for (i in 1:length(VecOfNuc1)/3){
          ##### i = 1
          Codons1 = c(VecOfNuc1[StartNuc:(StartNuc+2)])
          Codons2 = c(VecOfNuc2[StartNuc:(StartNuc+2)])
          Codons1 = paste(Codons1, collapse = '')
          Codons2 = paste(Codons2, collapse = '')
          if (Codons1 != Codons2){
            DifBtwSp1AndSp2 = paste(DifBtwSp1AndSp2,StartNuc, ':', Codons1, '>', Codons2, ';', sep = "")
          }
          StartNuc = StartNuc + 3
        }
        
        DifferenceInOneTrio$AllCodonSubstBetweenSp1AndSp2 = DifBtwSp1AndSp2
        
        StartNuc = 1
        
        for (i in 1:length(VecOfNuc1)/3){
          ##### i = 1
          Codons1 = c(VecOfNuc1[StartNuc:(StartNuc+2)])
          Codons3 = c(VecOfNuc3[StartNuc:(StartNuc+2)])
          Codons1 = paste(Codons1, collapse = '')
          Codons3 = paste(Codons3, collapse = '')
          if (Codons1 != Codons3){
            DifBtwSp1AndSp3 = paste(DifBtwSp1AndSp3, StartNuc, ':', Codons1, '>', Codons3, ';', sep = "")
          }
          StartNuc = StartNuc + 3
        }
        
        DifferenceInOneTrio$AllCodonSubstBetweenSp1AndSp3 = DifBtwSp1AndSp3
        
        StartNuc = 1
        
        for (i in 1:length(VecOfNuc2)/3){
          ##### i = 1
          Codons2 = c(VecOfNuc2[StartNuc:(StartNuc+2)])
          Codons3 = c(VecOfNuc3[StartNuc:(StartNuc+2)])
          Codons2 = paste(Codons2, collapse = '')
          Codons3 = paste(Codons3, collapse = '')
          if (Codons2 != Codons3){ DifBtwSp2AndSp3 = paste(DifBtwSp2AndSp3, StartNuc, ':', Codons2, '>', Codons3, ';', sep = "")}
          StartNuc = StartNuc + 3
        }
        
        
        DifferenceInOneTrio$AllCodonSubstBetweenSp2AndSp3 = DifBtwSp2AndSp3
      
        if (numofloop == 1) {
          DifferenceInOneFamily = DifferenceInOneTrio
        } else {
          DifferenceInOneFamily = rbind(DifferenceInOneFamily,DifferenceInOneTrio)
        }
        
        numofloop = numofloop + 1
      }
    }
  }
  
  #### write file for a family in CytB #####
  filename = paste('DifferenceBetween', fam ,'.txt', sep = '')
  write.table(DifferenceInOneFamily,file = filename)

}



