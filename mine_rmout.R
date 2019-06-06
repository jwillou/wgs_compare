setwd("/Volumes/jwillou/wgs/te/")

spp = read.table("specieslist.txt")
colnames(spp) = "species"

OUT = NULL
for(s in 1:nrow(spp)){
  rmout = scan("~/Desktop/genome_format_unmask.fa.tbl", what="character", sep=NULL)
  tlen  = as.numeric(strsplit(rmout[11], '\\(')[[1]][2])
  sines = as.numeric(rmout[38])
  lines = as.numeric(rmout[56])
  ltrs  = as.numeric(rmout[81])
  dnae  = as.numeric(rmout[112])
  ucls  = as.numeric(rmout[130])
  totl  = sines + lines + ltrs + dnae + ucls
  prop  = totl/tlen
  OUT   = rbind(OUT, c(as.character(spp[s,1]), tlen, sines, lines, ltrs, dnae, ucls, prop))
}
colnames(OUT) = c("species", "tlen", "sines", "lines", "ltrs", "dnae", "ucls", "prop")
write.table(OUT, "te_genomes.csv", row.names=F, col.names=T, sep=",", quote=F)
