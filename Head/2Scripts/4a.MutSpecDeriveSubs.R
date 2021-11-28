rm(list=ls(all=TRUE))

library(seqinr)
library(Biostrings)

#### read triples ####

trios = read.table('../../Body/1Raw/trios_d5.txt', header =TRUE)


### read our genetic code for cytB ###

seq = readDNAStringSet('../../Body/1Raw/cytb.fa')
seq = data.frame(seq); names(seq) = c('Sequence') 
seq$Scientific_name  = as.factor(rownames(seq))
row.names(seq) = c(1:nrow(seq))

### merge our data sets #####

trios_seq = merge(trios, seq, by.x = 'Ingroup0', by.y = 'Scientific_name')
colnames(trios_seq)[colnames(trios_seq) == 'Sequence'] <- 'SeqIn0'

trios_seq = merge(trios_seq, seq, by.x = 'Ingroup1', by.y = 'Scientific_name')
colnames(trios_seq)[colnames(trios_seq) == 'Sequence'] <- 'SeqIn1'

trios_seq = merge(trios_seq, seq, by.x = 'Outgroup', by.y = 'Scientific_name')
colnames(trios_seq)[colnames(trios_seq) == 'Sequence'] <- 'SeqOut'


VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')


### function to take subs from 3 species sequence that have rule: BAA or ABA, where A and B two dif codons. 
makemut <- function(x)
{
  StartNuc = 1
  numofcodon = 1
  Output=''
  
  VecOfNucOut = unlist(strsplit(x[3],''))
  VecOfNucIn0 = unlist(strsplit(x[1],''))
  VecOfNucIn1 = unlist(strsplit(x[2],''))
  
  for (cc in 1:(length(VecOfNucIn0)/3))
  {
    CodonsOut = c(VecOfNucOut[StartNuc:(StartNuc+2)])
    CodonsIn0 = c(VecOfNucIn0[StartNuc:(StartNuc+2)])
    CodonsIn1 = c(VecOfNucIn1[StartNuc:(StartNuc+2)])
    
    CodonsOut = paste(CodonsOut,collapse='')
    CodonsIn0 = paste(CodonsIn0, collapse = '')
    CodonsIn1 = paste(CodonsIn1, collapse = '')
    
    if  (grepl('-',CodonsIn0) | grepl('-',CodonsIn1) | grepl('-',CodonsOut)) {StartNuc = StartNuc + 3; numofcodon = numofcodon + 1}
    else  if (!grepl('-',CodonsIn0) & !grepl('-',CodonsIn1) & !grepl('-',CodonsOut)) 
    {
      A1 = seqinr::translate(unlist(strsplit(CodonsOut,'')), numcode = 2) #Amino Acid for Outgroup
      A2 = seqinr::translate(unlist(strsplit(CodonsIn0,'')), numcode = 2) #Amino Acid for Ingroup0
      A3 = seqinr::translate(unlist(strsplit(CodonsIn1,'')), numcode = 2) #Amino Acid for Ingroup1
      
      if (A1 == A2 & A1 == A3 & A2 == A3)
      {
        if (CodonsIn0 %in% VecOfSynFourFoldDegenerateSites & CodonsIn1 %in% VecOfSynFourFoldDegenerateSites 
            & CodonsOut %in% VecOfSynFourFoldDegenerateSites) 
        {
          ### take just 1 difference btw codons
          
          DifIn0In1 = 0
          DifIn0Out = 0
          DifIn1Out = 0
          total_sub = 0
          
          NucOut = unlist(strsplit(CodonsOut,''))
          NucIn0 = unlist(strsplit(CodonsIn0,''))
          NucIn1 = unlist(strsplit(CodonsIn1,''))
          
          for (n in 1:3)  ### count difference between  the same codons
          {#n = 1
            
            if(NucIn0[n] != NucOut[n]) {DifIn0Out = DifIn0Out+1; total_sub = total_sub + 1}
            if(NucIn1[n] != NucOut[n]) {DifIn1Out = DifIn1Out+1 ; total_sub = total_sub + 1}
            if(NucIn0[n] != NucIn1[n]) {DifIn0In1 = DifIn0In1+1 ; pl_sub = n ; total_sub = total_sub + 1}
          }
          if (total_sub == 2 & DifIn0In1 == 1) 
          {
            if (DifIn0In1 == 1 & DifIn0Out == 0)
            {
              Subs = paste(NucIn0[pl_sub], '>', NucIn1[pl_sub], sep='')
              Dif = paste(numofcodon,'.',Subs,'.', 'In0>In1', sep='')
              Output = paste(Output, ';', Dif, sep = '')
              numofcodon = numofcodon +1 ; StartNuc = StartNuc + 3
              
            }
            else if (DifIn0In1 == 1 & DifIn1Out == 0)
            {
              Subs = paste(NucIn1[pl_sub], '>', NucIn0[pl_sub], sep='')
              Dif = paste(numofcodon,'.',Subs,'.', 'In1>In0',sep='')
              Output = paste(Output, ';', Dif, sep = '')
              numofcodon = numofcodon + 1 ; StartNuc = StartNuc + 3
            }
          }
        else {numofcodon = numofcodon +1 ; StartNuc = StartNuc + 3}
        } else {numofcodon = numofcodon +1 ; StartNuc = StartNuc + 3}
      } else {numofcodon = numofcodon +1 ; StartNuc = StartNuc + 3} 
    }
  }
  return(Output)
}

trios_seq$MutSpec = apply(trios_seq[,5:7],1,makemut)

trios_seq = trios_seq[trios_seq$MutSpec != '' ,] ### delete empty rows

### delete unnecessary columns
trios_seq$SeqIn0 <- NULL
trios_seq$SeqIn1 <- NULL
trios_seq$SeqOut <- NULL

#### save table
write.csv(trios_seq,file = '../../Mut_Spec_4Fold.csv')
