rm(list=ls(all=TRUE))

library(seqinr)

subs = read.table('C:/Users/MitoClubHouse/Documents/data/DIliushchenko/SubBtwTrio_d5.txt', header = T, skip = 1)

subs = subs[,-1]

VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')


### make function, that takes only 4-fold degenerative sites
makecodon <- function(a)
{
    Indel = 0
    VecofSubs = ''
    
    SubstitutionsSp1Sp2 = unlist(strsplit(a,';'))
    
    for (i in 1:length(SubstitutionsSp1Sp2)) 
    {
      #i=1
      TwoCodons = SubstitutionsSp1Sp2[i] # TwoCodons = '430:GCC>GCT'
      TwoCodons = gsub(".*:",'',TwoCodons) # GCC>GCT
      TwoCodons = unlist(strsplit(TwoCodons,'>'))
      Codon1 <- unlist(strsplit(TwoCodons[1],'>'))
      Codon2  <- unlist(strsplit(TwoCodons[2],'>'))
      if  (grepl('-',Codon1) | grepl('-',Codon2)) {Indel=Indel+1; break}
      if  (!grepl('-',Codon1) & !grepl('-',Codon2)) 
      {
        A1 = translate(unlist(strsplit(Codon1,'')), numcode = 2) #### use translate from seqinr except Biostrings(problems with 'I')
        A2 = translate(unlist(strsplit(Codon2,'')),numcode = 2 ) ### numcode mitochondrial for vertebrates look at: https://rdrr.io/rforge/seqinr/man/translate.html
        if (A1 == A2) 
        {
          if (Codon1 %in% VecOfSynFourFoldDegenerateSites) 
            {
              Dif = 0
              Codon1 = unlist(strsplit(Codon1,''))
              Codon2 = unlist(strsplit(Codon2,''))
              for (j in 1:3) 
              {
                if(Codon1[j] != Codon2[j]) {Dif = Dif+1}
              }
              if (Dif == 1) 
                {
                Subs = paste(Codon1[3],Codon2[3], sep = '>')
                VecofSubs = paste(VecofSubs,Subs,sep=';')  
                }
            }
        }
    }
  }
  return(VecofSubs)
}
  


# apply function for each row
subs$subsIn0In1 = apply(as.matrix(subs$DifIn0In1),1,makecodon)
subs$subsIn0Out = apply(as.matrix(subs$DifIn0Out),1,makecodon)
subs$subsIn1Out = apply(as.matrix(subs$DifIn1Out),1,makecodon)


subs$DifIn0In1 <- NULL
subs$DifIn0Out <- NULL
subs$DifIn1Out <- NULL


subs = subs[!subs$subsIn0In1 == "",] # delete bad trios  with no subs btw In0 and In1 


write.csv(x = subs, file = '../../Body/2Derived/Mut_Spec_4Fold.csv')
