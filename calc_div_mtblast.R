setwd("/scratch/snyder/j//jwillou/wgs/compare/divergence/")
#setwd("/Volumes/jwillou/wgs/compare/divergence/")

library(ape)
directories = read.table("../mtblast_21/mtDNAcompare.txt", sep="\t")
colnames(directories) = c("class", "sp1", "sp2")
directories$combined = paste(directories$sp1, directories$sp2, sep="_")

#iterate over directories list and estimate divergence
divOUT = NULL
for(f in 1:nrow(directories)){ 
  #read in each data file
  data = read.table(paste("../mtblast_21/", directories$combined[f], "/outblast21.txt", sep=""), sep="\t", header=F)
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
write.table(divOUT, "mtDNAdiv_estimates.csv", sep=",", row.names=F, col.names=T)
