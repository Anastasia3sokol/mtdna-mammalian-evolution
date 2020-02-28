library(Biostrings)
library(seqinr)

rm(list=ls(all=TRUE))

wd <- getwd()
wd <- paste(wd, '/Mammalian mtDNA Evolution/Body/1Raw/', sep='')
setwd(wd)


##### Open file with Information About Generation Length in Animals #######
GenerationL = read.table('GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
GenerationL = GenerationL[ ,c(2,11)] # Take information about name and GenLength
GenerationL$Scientific_name = gsub(' ','_',GenerationL$Scientific_name)


SizeOfRandom = 50
RandomNumbers = sample(1:650,SizeOfRandom, replace = FALSE)


wd = paste(wd, '/Codons', sep='')
setwd(wd) ### Change direcotry for Codons
DifferenceBetweenSpecies = data.frame()

files = list.files(pattern=".*\\.fasta")

for ( f in files){
  # f = 'ATP6_NT.fasta'
  Testdf = readDNAStringSet(f)
  Testdf = data.frame(Testdf); names(Testdf) = c('Sequence') 
  Testdf$Scientific_name  = as.factor(rownames(Testdf))
  row.names(Testdf) = c(1:nrow(Testdf))
  
  #### merge dataframe GenLen with dataframe Sequence
  MergeDf = merge(GenerationL, Testdf, by = 'Scientific_name' )
  
  
  ####### Create a new data frame with 3 species: High and Low Generation time, and Outgroup####
  
  for ( v in 1:(SizeOfRandom-1)){
    #v = 1
    NumberofSp1 = RandomNumbers[v]
    Sp1 = MergeDf[NumberofSp1,]
    for (s in (v+1):SizeOfRandom){
      #s=2
      NumberofSp2 = RandomNumbers[s]
      Sp2 = MergeDf[NumberofSp2,]
      if (Sp1$GenerationLength_d > Sp2$GenerationLength_d){
        HighGenLength = MergeDf[NumberofSp1, ]
        LowGenLength = MergeDf[NumberofSp2, ]
      } else { 
        HighGenLength = MergeDf[NumberofSp2,]
        LowGenLength = MergeDf[NumberofSp1,] 
      }
      OutGroup = MergeDf[MergeDf$Scientific_name == 'Ornithorhynchus_anatinus', ]
      
      
      ######### Make Sequences readable #########################
      
      VecOfNucH = unlist(strsplit(HighGenLength$Sequence,''))
      VecOfNucL = unlist(strsplit(LowGenLength$Sequence,''))
      VecOfNucO = unlist(strsplit(OutGroup$Sequence,''))
      
      StartNuc = 1
      CodonsH = c()
      CodonsL = c()
      CodonsO = c()
    
      ######## Create vectors for differences ##################
      
      DifBtwSpLAndSpH = c()
      DifBtwOutAndSpL = c()
      DifBtwOutAndSpH = c()
      
      ############# Create a new data frame with information about Species and  difference between every species for One Gene ##############################
      
      DifferenceInOneGene = data.frame(OutGroup$Scientific_name, LowGenLength$Scientific_name, HighGenLength$Scientific_name, '', '', '', f) 
      names(DifferenceInOneGene) = c('Outgroup', 'LowGenerationLength', 'HighGenerationLength', 'AllCodonSubstBetweenOutgroupAndSpL' ,
                                     'AllCodonSubstBetweenOutgroupAndSpH', 'AllCodonSubstBetweenSpLAndSpH','Gene')
      
      
      ############# Loop for differences ############ 
      
      for (i in 1:length(VecOfNucO)/3){
        ##### i = 1
        CodonsO = c(VecOfNucO[StartNuc:(StartNuc+2)])
        CodonsL = c(VecOfNucL[StartNuc:(StartNuc+2)])
        CodonsO = paste(CodonsO, collapse = '')
        CodonsL = paste(CodonsL, collapse = '')
        if (CodonsO != CodonsL){
          DifBtwOutAndSpL = paste(DifBtwOutAndSpL, CodonsO, '>', CodonsL, ';', sep = "")
        }
        StartNuc = StartNuc + 3
      }
      
      DifferenceInOneGene$AllCodonSubstBetweenOutgroupAndSpL = DifBtwOutAndSpL
      
      StartNuc = 1
      
      for (i in 1:length(VecOfNucO)/3){
        ##### i = 1
        CodonsO = c(VecOfNucO[StartNuc:(StartNuc+2)])
        CodonsH = c(VecOfNucH[StartNuc:(StartNuc+2)])
        CodonsO = paste(CodonsO, collapse = '')
        CodonsH = paste(CodonsH, collapse = '')
        if (CodonsO != CodonsH){
          DifBtwOutAndSpH = paste(DifBtwOutAndSpH, CodonsO, '>', CodonsH, ';', sep = "")
        }
        StartNuc = StartNuc + 3
      }
      
      DifferenceInOneGene$AllCodonSubstBetweenOutgroupAndSpH = DifBtwOutAndSpH
      
      StartNuc = 1
      for (i in 1:length(VecOfNucL)/3){
        ##### i = 1
        CodonsL = c(VecOfNucL[StartNuc:(StartNuc+2)])
        CodonsH = c(VecOfNucH[StartNuc:(StartNuc+2)])
        CodonsL = paste(CodonsL, collapse = '')
        CodonsH = paste(CodonsH, collapse = '')
        if (CodonsL != CodonsH){
          DifBtwSpLAndSpH = paste(DifBtwSpLAndSpH, CodonsL, '>', CodonsH, ';', sep = "")
        }
        StartNuc = StartNuc + 3
      }
      
      DifferenceInOneGene$AllCodonSubstBetweenSpLAndSpH = DifBtwSpLAndSpH
      if ( s == 2 & f == 'ATP6_NT.fasta'){
        DifferenceBetweenSpecies = DifferenceInOneGene
      } else {
        DifferenceBetweenSpecies = rbind(DifferenceBetweenSpecies, DifferenceInOneGene)
      }
    }
    }
}
  
  
wd = gsub('1Raw/Codons', '2Derived', wd)
setwd(wd)
write.table(DifferenceBetweenSpecies, file = "DifferenceBetweenSpecies50")
