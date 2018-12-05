# Uses chr1, chr2 naming

require(HiTC)
require(rtracklayer)
require(BSgenome.Hsapiens.UCSC.hg38)

## Generate restriction fragment file
human_chr <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)[1:24]
resFrag <- getRestrictionFragmentsPerChromosome(resSite="AAGCTT", chromosomes=human_chr, overhangs5=1, genomePack="BSgenome.Hsapiens.UCSC.hg38")
allRF <- do.call("c",resFrag)
names(allRF) <- unlist(sapply(resFrag, function(x){paste0("HIC_", seqlevels(x), "_", 1:length(x))}))
export(allRF, format="bed", con="Hindiii_resfrag_hg38.bed")


## Generate chromosome size file
human_chr <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)[1:24]
chrom_size <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[human_chr]
write.table(chrom_size, file="chrom_hg38.sizes", quote=FALSE, col.names=FALSE, sep="\t")

