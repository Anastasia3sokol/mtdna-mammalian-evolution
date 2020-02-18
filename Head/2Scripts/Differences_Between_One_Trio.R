library(Biostrings)
library(seqinr)

rm(list=ls(all=TRUE))
wd <- getwd()
wd <- paste(wd, '/Mammalian mtDNA Evolution/Body/1Raw/Codons/', sep='')
setwd(wd)
DifferenceBetweenSpecies = data.frame()
files = list.files(pattern=".*\\.fasta")
for ( f in files){
  Testdf = readDNAStringSet(f)
  Testdf = data.frame(Testdf); names(Testdf) = c('Sequence') 
  Testdf$Name = rownames(Testdf)
  
  ####### Create a new data frame with 3 species: High and Low Generation time, and Outgroup####
  
  HighGenLenght = Testdf['Elephas_maximus', ]
  LowGenLenght = Testdf['Mus_musculus', ]
  OutGroup = Testdf['Ornithorhynchus_anatinus', ]
  FinalDataFrame = rbind(HighGenLenght,LowGenLenght, OutGroup)
  str(FinalDataFrame)
  
  
  ######### Split Sequences #########################
  
  VecOfNucH = unlist(strsplit(FinalDataFrame[1,1],''))
  VecOfNucL = unlist(strsplit(FinalDataFrame[2,1],''))
  VecOfNucO = unlist(strsplit(FinalDataFrame[3,1],''))
  
  StartNuc = 1
  CodonsH = c()
  CodonsL = c()
  CodonsO = c()
  
  ######## Create vectors for differences ##################
  
  DifBtwSpLAndSpH = c()
  DifBtwOutAndSpL = c()
  DifBtwOutAndSpH = c()
  
  ############# Create a new data frame with information about Species and  difference between every species for One Gene ##############################
  
  DifferenceInOneGene = data.frame( 'Ornithorhynchus_anatinus','Mus_musculus','Elephas_maximus','','','',f) 
  names(DifferenceInOneGene) = c('Outgroup', 'LowGenerationLenght', 'HighGenerationLenght', 'AllCodonSubstBetweenOutgroupAndSpL' ,
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
  if ( f == 'ATP6_NT.fasta'){
    DifferenceBetweenSpecies = DifferenceInOneGene
  } else {
    DifferenceBetweenSpecies = rbind(DifferenceBetweenSpecies, DifferenceInOneGene)
  }
  
}
wd = gsub('1Raw/Codons', '2Derived', wd)
setwd(wd)
write.table(DifferenceBetweenSpecies, file = "DifferenceBetweenSpeciesOneTrio")
