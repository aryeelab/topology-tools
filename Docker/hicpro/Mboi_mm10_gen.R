require(HiTC)
require(rtracklayer)
require(BSgenome.Mmusculus.UCSC.mm10)

## Generate restriction fragment file
human_chr <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)[1:21]
resFrag <- getRestrictionFragmentsPerChromosome(resSite="GATC", chromosomes=human_chr, overhangs5=0, genomePack="BSgenome.Mmusculus.UCSC.mm10")
allRF <- do.call("c",resFrag)
names(allRF) <- unlist(sapply(resFrag, function(x){paste0("HIC_", seqlevels(x), "_", 1:length(x))}))
export(allRF, format="bed", con="Mboi_resfrag_mm10.bed")

## Generate chromosome size file
chr <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)[1:21]
chrom.size <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)[chr]
write.table(chrom_size, file="chrom_mm10.sizes", quote=FALSE, col.names=FALSE, sep="\t")
