setwd("/scratch/snyder/j//jwillou/wgs/compare/divergence/")
#setwd("/Volumes/jwillou/wgs/compare/divergence/")

#install.packages("ape", lib="/scratch/snyder/j/jwillou/wgs/compare/divergence/rlibs/", repos='http://cran.us.r-project.org')
.libPaths("/scratch/snyder/j/jwillou/wgs/compare/divergence/rlibs") 

library(ape)
directories = read.table("../blast_21/todivR.txt", sep="\t")
colnames(directories) = c("sp1_sp2")

#iterate over directories list and estimate divergence
write.table(t(c("spp", "length", "raw", "K80", "K81", "F84", "TN93")), "DNAdiv_estimates.csv", sep=",", row.names=F, col.names=F, append=F)
write.table(t(c("starting loop")), "DNAdiv_watch.csv", sep=",", row.names=F, col.names=F, append=F)
for(f in 1:nrow(directories)){
  #read in each data file
  data = read.table(paste("../blast_21/", directories$sp1_sp2[f], "/outblast21.txt", sep=""), sep="\t", header=F)
  colnames(data) = c("qid", "sid", "percid", "alength", "mismatches", "gapopen", "qstart", "qend", "sstart", 
                     "send", "e", "bitscore", "nident", "mismatch", "qseq", "sseq")
  
  #sub-sample if number of sites greater than upper limit
  limit = (2^31-1)/2  #this could be run again with smallest length value across all species pairs
  if(sum(data$alength)>limit){
    tokeep = sample(seq(1,nrow(data)), limit/1000, replace=F)
    tdata = data[tokeep,]
    while(sum(tdata$alength)>limit){
      tokeep = sample(seq(1,nrow(tdata)), (nrow(tdata)-1), replace=F)
      tdata = tdata[tokeep,]
    }
    data = tdata
  }
  
  #prep aligned sequences
  q = strsplit(paste(data$qseq, "", collapse=""), "")
  s = strsplit(paste(data$sseq, "", collapse=""), "")
  
  temp = as.DNAbin(c(q, s))
  write.table(paste("created DNAbin for", as.character(directories$sp1_sp2[f])), "DNAdiv_watch.csv", sep=",", row.names=F, col.names=F, append=T)

  #estimate divergences
  raw = K80 = K81 = F84 = TN93 = NULL
  raw = dist.dna(temp, model=c("raw"))
  K80 = dist.dna(temp, model=c("K80"))
  K81 = dist.dna(temp, model=c("K81"))
  F84 = dist.dna(temp, model=c("F84"))
  TN93 = dist.dna(temp, model=c("TN93"))
  
  #save output
  write.table(t(c(as.character(directories$sp1_sp2[f]), summary(temp)[1,1], raw, K80, K81, F84, TN93)), "DNAdiv_estimates.csv", sep=",", row.names=F, col.names=F, append=T)
}
