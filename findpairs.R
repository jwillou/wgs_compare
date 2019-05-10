setwd("~/Desktop/")
data = read.table("mtDNAspp.txt", header=F, sep="\t")
colnames(data) = c("G", "s", "G_s")

genus = unique(data$G)

clist = NULL
for(g in 1:length(genus)){
  t = data[data$G==as.character(genus[g]),,drop=F]
  if(nrow(t)>1){
    glist = t(combn(x=c(as.character(t$G_s)), m=2))
    clist = rbind(clist, glist)
  }
}
colnames(clist) = c("sp1", "sp2")
write.table(clist, "mtDNAsppPairs.txt", sep="\t", row.names=F, col.names = T)
