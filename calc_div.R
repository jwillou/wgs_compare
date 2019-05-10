setwd("/scratch/snyder/j//jwillou/wgs/compare/dxy/")
setwd("/Volumes/jwillou/wgs/compare/dxy/")

library(ape)
directories = read.table("../blast_21/compare.txt", sep="\t")
colnames(directories) = c("sp1", "sp2", "class")
directories$combined = paste(directories$sp1, directories$sp2, sep="_")

#iterate over directories list and estimate divergence
divOUT = NULL
for(f in 1:1){ #length(directories)){
  #read in each data file
  data = read.table(paste("../blast_21/", directories$combined[f], "/outblast21.txtT", sep=""), sep="\t", header=F)
  colnames(data) = c("qid", "sid", "percid", "alength", "mismatches", "gapopen", "qstart", "qend", "sstart", 
                     "send", "e", "bitscore", "nident", "mismatch", "qseq", "sseq")
  #prep aligned sequences
  q = strsplit(paste(data$qseq, "", collapse=""), "")
  s = strsplit(paste(data$sseq, "", collapse=""), "")
  temp = as.DNAbin(c(q, s))
  
  #estimate divergences
  raw = dist.dna(temp, model=c("raw"))
  K80 = dist.dna(temp, model=c("K80"))
  K81 = dist.dna(temp, model=c("K81"))
  F84 = dist.dna(temp, model=c("F84"))
  TN93 = dist.dna(temp, model=c("TN93"))
  
  #save output
  divOUT = rbind(divOUT, c(directories$combined[f], summary(temp)[1,1], raw, K80, K81, F84, TN93))
}

divOUT = as.data.frame(divOUT)
colnames(divOUT) = c("spp", "length", "raw", "K80", "K81", "F84", "TN93")

