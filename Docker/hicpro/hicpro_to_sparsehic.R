# Usage:
# Rscript hicpro_to_sparsehic.R --sample_id=imr90-rep1  --genome_size=chrom_hg19.sizes --genome_id=hg19 --cores=2

library(sparseHiC)
library(BiocParallel)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("--sample_id"), type="character",
                     help="Sample name")
parser <- add_option(parser, c("--genome_size"), type="character",
                     default=".", help="Chromosome length table. One line per chr. Two columns: chr name, length in bp")
parser <- add_option(parser, c("--genome_id"), type="character",
                     help="Genome name for metadata slot annotation (e.g. hg19)")
parser <- add_option(parser, c("-p", "--cores"), type="integer",
                     default=2, help="Number of cores to use")


# args <- parse_args(parser, args = c("--sample_id=imr90-rep1",
#                             "--genome_size=chrom_hg19.sizes",
#                             "--genome_id=hg19",
#                             "--cores=8"))

args <- parse_args(parser)

interactionsFile <- list.files(".", full.names=TRUE, pattern=paste0(args$sample_id, ".*_iced.matrix$"))
bedFile <- list.files(".", full.names=TRUE, pattern="abs.bed")
stopifnot( gsub("\\S+_(\\S+)_\\S+", "\\1", basename( bedFile )) ==
             gsub("\\S+_(\\S+)_\\S+", "\\1", basename( interactionsFile )) )
rs <- gsub("\\S+_(\\S+)_\\S+", "\\1", basename( interactionsFile ))

gsize <- read.table(args$genome_size, stringsAsFactors = FALSE)
colnames(gsize) <- c("chr", "length")

hic <- import.HiCPro(
  matrix.files=interactionsFile, bed.files=bedFile,
  resolutions=rs, sampleName=args$sample_id,
  manual.chr = gsize$chr, manual.dist = gsize$length,
  BPPARAM=MulticoreParam(args$cores), n=0, temp=FALSE )

hic@metaData <- data.frame(genomeBuild = args$genome_id)

saveRDS(hic, file = paste0(args$sample_id, ".sparsehic.rds"))
