# Usage:
# Rscript filter_matrix.R --input_matrix=MGH-1904_100000.matrix --bins=MGH-1904_100000_abs.bed --blacklist=blacklist.bed  --output_matrix=MGH-1904_100000.blacklist.matrix


library(readr)
library(optparse)
suppressPackageStartupMessages(library(GenomicRanges))

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input_matrix"), type="character",
                     help="Input matrix file in TSV sparse triplet format (i j count)")
parser <- add_option(parser, c("-o", "--output_matrix"), type="character",
                     help="Output matrix file in TSV sparse triplet format (i j count)")
parser <- add_option(parser, c("-b", "--bins"), type="character",
                     help="Bins (chr start end)")
parser <- add_option(parser, c("--blacklist"), type="character",
                     default=NULL, help="Regions to filter out")


args <- parse_args(parser, args = c(
                            "--input_matrix=MGH-1904_100000.matrix",
                            "--output_matrix=MGH-1904_100000.blacklist.matrix",
                            "--bins=MGH-1904_100000_abs.bed",
                            "--blacklist=blacklist.bed"))

args <- parse_args(parser)

args

# Read bins
message ("Reading bins: ", args$bins)
bins <- read_tsv(args$bins, 
                 col_names = c("chrom", "start", "end", "id"),
                 col_types = cols_only(chrom = col_character(),
                                       start = col_integer(),
                                       end = col_integer()))
bin_gr <- GRanges(seqnames=bins$chrom, IRanges(start=bins$start, end=bins$end))

# Read matrix
message ("Reading input matrix: ", args$bins)
tab <- read_tsv(args$input_matrix, 
                col_names = c("i", "j", "count"),
                col_types = cols_only(i = col_integer(),
                                      j = col_integer(),
                                      count = col_double()))
message (nrow(bins), " non-zero cells")

# Read blacklist
if (!is.null(args$blacklist)) {
    message ("Reading blacklist: ", args$blacklist)
    blacklist <- read_tsv(args$blacklist, 
                        col_names = c("chrom", "start", "end"),
                        col_types = cols_only(chrom = col_character(),
                                              start = col_integer(),
                                              end = col_integer()))
  blacklist_gr <- GRanges(seqnames=blacklist$chrom, IRanges(start=blacklist$start+1, end=blacklist$end))
  bin_drop_idx <- which(countOverlaps(bin_gr, blacklist_gr) > 0)
  message ("Excluding ", length(bin_drop_idx), " bins out of ", nrow(bins), " because of overlap with the blacklist")
}

# Drop matrix cells
message ("Dropping cells")
drop_idx <- tab$i %in% bin_drop_idx | tab$j %in% bin_drop_idx
tab <- tab[!drop_idx,]

# Sort matrix
message ("Sorting by i,j")
o <- order(tab$i, tab$j)
tab <- tab[o,]

# Save matrix in sparse triplet format
message ("Writing output matrix: ", args$output_matrix)
write.table(tab, file=args$output_matrix, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
