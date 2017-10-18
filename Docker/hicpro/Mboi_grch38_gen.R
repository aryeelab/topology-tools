# Uses 1, 2 naming

require(HiTC)
require(rtracklayer)
require(BSgenome.Hsapiens.NCBI.GRCh38)

## Generate restriction fragment file
human_chr <- seqlevels(Hsapiens)[1:24]
resFrag <- getRestrictionFragmentsPerChromosome(resSite="GATC", chromosomes=human_chr, overhangs5=1, genomePack="BSgenome.Hsapiens.NCBI.GRCh38")
allRF <- do.call("c",resFrag)
names(allRF) <- unlist(sapply(resFrag, function(x){paste0("HIC_", seqlevels(x), "_", 1:length(x))}))
export(allRF, format="bed", con="Mboi_resfrag_grch38.bed")


## Generate chromosome size file
human_chr <- seqlevels(Hsapiens)[1:24]
chrom_size <- seqlengths(Hsapiens)[human_chr]
write.table(chrom_size, file="chrom_grch38.sizes", quote=FALSE, col.names=FALSE, sep="\t")

