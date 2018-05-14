# Usage:
# Rscript bnbc.R --in_dir=hicpro_matrices --sample_names=imr90-rep1,imr90-rep2 --resolution=10000 --chromosome=19 --cores=2

library(readr)
library(Matrix)
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(bnbc))
library(optparse)
suppressPackageStartupMessages(library(foreach))

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--in_dir"), type="character",
                                        default=".", help="Input matrix folder")
parser <- add_option(parser, c("-o", "--out_dir"), type="character",
                                            default=".", help="Output folder for normalized matrices")
parser <- add_option(parser, c("-s", "--sample_names"), type="character",
                                        help="Sample names")
parser <- add_option(parser, c("--matrix_suffix"), type="character",
                                        default="_iced.matrix", help="Matrix file name suffix")
parser <- add_option(parser, c("-r", "--resolution"), type="integer",
                     help="Bin size / resolution in bp")
parser <- add_option(parser, c("-c", "--chromosome"), type="character",
                     help="Chromosome to normalize. Use 'inter_chromosomal' for inter-chromosomal contacts")
parser <- add_option(parser, c("--blacklist"), type="character",
                     default=NULL, help="Regions to exclude from normalization and from output matrices.")
parser <- add_option(parser, c("-p", "--cores"), type="integer",
                     help="Number of cores to use")


# args <- parse_args(parser, args = c("--in_dir=hicpro_matrices",
#                             "--sample_names=imr90-rep1,imr90-rep2",
#                             "--resolution=10000",
#                             "--chromosome=19",
#                             "--cores=2"))

# args <- parse_args(parser, args = c("--in_dir=colon_matrix_1000000",
#                               "--sample_names=BRD3179N,BRD3187N",
#                               "--resolution=1000000",
#                               "--chromosome=chr19",
#                               "--cores=2"))

args <- parse_args(parser)

sample_names <- strsplit(args$sample_names, ",")[[1]]
args

cpm <- function(x) {
  library_size <- sum(x)
  (x + 0.5) / ((library_size+1)/10^6)
}

bin_file <- file.path(args$in_dir, paste0(sample_names[1], "_", args$resolution, "_abs.bed"))
message("Input bin file:")
message(bin_file)

matrix_files <- file.path(args$in_dir, paste0(sample_names, "_", args$resolution, args$matrix_suffix))
message("\nInput matrix files:")
message(paste(matrix_files, collapse="\n"))

bins <- read_tsv(bin_file, 
                 col_names = c("chrom", "start", "end", "id"),
                 col_types = cols_only(chrom = col_character(),
                                      start = col_integer(),
                                      end = col_double())
)
genome_bin_gr <- GRanges(seqnames=args$chromosome, IRanges(start=bins$start, end=bins$end))

if (!is.null(args$blacklist)) {
  blacklist_tab <- read.delim(args$blacklist, header=FALSE, stringsAsFactors = FALSE)
  blacklist_gr <- GRanges(seqnames=blacklist_tab[,1], IRanges(start=blacklist_tab[,2]+1, end=blacklist_tab[,3]))
  bin_drop_idx <- which(countOverlaps(genome_bin_gr, blacklist_gr) > 0)
  message ("Excluding ", length(bin_drop_idx), " bins out of ", nrow(bins), " because of overlap with the blacklist")
} else {
  message("No blacklist regions provided")
}


if (args$chromosome=="inter_chromosomal") {
  # Only perform CPM normalization
  message("\nNormalized output matrix files:")
  foreach (i = 1:length(sample_names)) %do% {
    matrix_file <- matrix_files[i]
    tab <- read_tsv(matrix_file, 
                    col_names = c("i", "j", "count"),
                    col_types = cols_only(i = col_integer(),
                                          j = col_integer(),
                                          count = col_double())
    )
    if (!is.null(args$blacklist)) {
      drop_idx <- tab$i %in% bin_drop_idx | tab$j %in% bin_drop_idx
      tab <- tab[!drop_idx,]
    }
    tab$cpm <- cpm(tab$count)
    
    # Write out CPM normalized sparse interchromosomal matrix
    idx <- bins$chrom[tab$i] != bins$chrom[tab$j] 
    filename <- file.path(args$out_dir, paste0(sample_names[i], "_", args$resolution, "_interchromosomal.cpm.matrix"))
    message(filename)
    write.table(tab[idx, c("i", "j", "cpm")], file=filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
} else {
  # Do bnbc normalization
  chr_idx <- bins$chrom==args$chromosome
  message("Of ", nrow(bins), " bins, ", sum(chr_idx), " correspond to chromosome ", args$chromosome)
  
  if (sum(chr_idx)==0) {
    message ("No output files will be created since no bins correspond to this chromosome")
  } else {
    # Create a list of matrices corresponding to the region
    # specified by 'args$chromosome' for each sample
    mat_list <- lapply(matrix_files, function(matrix_file) {
      message ("Extracting chromosome ", args$chromosome, " from ", matrix_file)
      tab <- read_tsv(matrix_file, 
                      col_names = c("i", "j", "count"),
                      col_types = cols_only(i = col_integer(),
                                            j = col_integer(),
                                            count = col_double())
      )
      if (!is.null(args$blacklist)) {
        drop_idx <- tab$i %in% bin_drop_idx | tab$j %in% bin_drop_idx
        tab <- tab[!drop_idx,]
      }
      
      # Use 0-based indexing for dgTMatrix
      mat <- new("dgTMatrix", i = as.integer(tab$i-1), j = as.integer(tab$j-1), x = tab$count, Dim=c(nrow(bins), nrow(bins)))
      as.matrix(mat[chr_idx, chr_idx])
    })
    names(mat_list) <- sample_names
    
    # Create ContactGroup
    sample_data <- DataFrame(rownames=sample_names)
    bin_gr <- GRanges(seqnames=args$chromosome, IRanges(start=bins$start, end=bins$end)[chr_idx])
    cg <- ContactGroup(rowData = bin_gr, contacts = mat_list, colData = sample_data)
    
    # Within-sample normalization
    cg_cpm <- logCPM(cg)
    cg_smooth <- boxSmoother(cg_cpm, h=3, mc.cores=args$cores)
    
    # Between-sample normalization
    cg_bnbc <- suppressWarnings(bnbc(cg_smooth, batch=1, bstart=2, nbands=nrow(cg)-1, threshold=NULL, step=as.numeric(args$resolution), mean.only=TRUE, verbose=FALSE))
    names(contacts(cg_bnbc)) <- names(contacts(cg))
    
    message("\nNormalized output matrix files:")
    
    # Write out normalized sparse matrices
    for (sample in sample_names) {
      mat <- contacts(cg_bnbc)[[sample]]
      mat[lower.tri(mat)] <- 0
      
      # Convert back from chromosome-level to genome-wide coordinates
      genome_mat <- new("dgTMatrix", Dim=c(nrow(bins), nrow(bins)))
      genome_mat[chr_idx, chr_idx] <- mat
      triplet_mat <- as(genome_mat, "dgTMatrix")
      
      # Save sparse triplet matrix
      triplet_df <- data.frame(bin1_id=triplet_mat@i, bin2_id=triplet_mat@j, count=triplet_mat@x)
      o <- order(triplet_df$bin1_id, triplet_df$bin2_id)
      triplet_df <- triplet_df[o,]
      filename <- file.path(args$out_dir, paste0(sample, "_", args$resolution, "_", args$chromosome, ".bnbc.matrix"))
      message(filename)
      write.table(triplet_df, file=filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
  }
}

