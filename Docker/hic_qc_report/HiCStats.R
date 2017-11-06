#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(optparse))
 
option_list = list(
    make_option(c("-s", "--sampleNames"), type="character", default=NULL, 
                help="Either a sample name or a set of samples separated by a comma. Example: -s sample1,sample2",
                metavar="character"),
    make_option(c("-d", "--inputDir"), type="character", default=NULL, 
                help="Path to HiCPro outputs, separated by a comma if more than one. Example: -d path1,path2",
                metavar="character"),
    make_option(c("-o", "--outputPrefix", type="character", default="muestra",
                help="prefix for plot names [default: muestra]. Example -o muestra",
                metavar="character") ) )  
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


samples <- opt[["sampleNames"]]
hicProDir <- opt[["inputDir"]]
prefix <- opt[["outputPrefix"]]

if( length(opt_parser) > 3 ){
    print_help( opt_parser )
    stop("No samples or paths were specified, nothing to be done\n")
}

#hicProDir <- "/data/aryee/bernstein/hic/hicpro_out_par"
#prefix <- "muestras"
#samples <- "HiC_1904,HiC_2834,HiC_5328,HiC_8416"

mapFiles <- list.files( hicProDir, pattern=".mpairstat", full.names=TRUE, recursive=TRUE )
names( mapFiles ) <- gsub( ".mpairstat", "", basename(mapFiles) )
statFiles <- list.files( hicProDir, pattern="mRSstat", full.names=TRUE, recursive=TRUE )
names( statFiles ) <- gsub( ".mRSstat", "", basename(statFiles) )
mergeFiles <- list.files( hicProDir, pattern="_allValidPairs.mergestat", full.names=TRUE, recursive=TRUE )
names( mergeFiles ) <- gsub("_allValidPairs.mergestat", "", basename( mergeFiles ) )

if (!is.null(samples)) {
    samples <- strsplit(samples, ",")[[1]]
    mapFiles <- mapFiles[samples]
    statFiles <- statFiles[samples]
    mergeFiles <- mergeFiles[samples]
}


if( !all( names( statFiles ) == names( mergeFiles ) ) ){
    stop("Missing files: each '.mRSstat' must have a corresponding '.mergestat' file")
}

if( !identical( samples, "all" ) ){
    allSamples <- samples %in% names( statFiles )
    if( !all( allSamples ) ){
        stop( sprintf("Files are missing for the following samples: %s", paste( samples[!allSamples], collapse=", ")))
    }
}

if( !all( file.exists( mapFiles ) ) ){
  missingSample <- samples[!file.exists( mapFiles )]
  stop(sprintf("The '*.mpairstat' files could not be found for the following sample(s):\n%s",
               paste(missingSample, collapse="\n")))
}
if( !all( file.exists( statFiles ) ) ){
    missingSample <- samples[!file.exists( statFiles )]
    stop(sprintf("The '*.mRSstat' files could not be found for the following sample(s):\n%s",
                 paste(missingSample, collapse="\n")))
}
if( !all( file.exists( mergeFiles ) ) ){
    missingSample <- samples[!file.exists( mergeFiles )]
    stop(sprintf("The '*.mergestat' files could not be found for the following sample(s):\n%s",
                 paste(missingSample, collapse="\n")))
}

allMapData <- melt( lapply( mapFiles, read.table, fill=TRUE),
                    id.vars="V1",
                    value.name="numberOfReads")
colnames( allMapData ) <- c( "mapStat", "dummy", "numberOfReads", "sample" )
tmp <- allMapData %>% filter(mapStat=="Total_pairs_processed" & dummy=="V2") 
totalReadNum <- tmp$numberOfReads
names(totalReadNum) <- tmp$sample

allMergeData <- melt( lapply( mergeFiles, read.table, fill=TRUE),
                     id.vars="V1",
                     value.name="numberOfReads")
colnames( allMergeData ) <- c( "interactionType", "dummy", "numberOfReads", "sample" )
allMergeData <- allMergeData[,colnames(allMergeData) !="dummy"]

allStatData <- melt( lapply( statFiles, read.table, fill=TRUE ),
                    id.vars="V1", value.name="numberOfReads" )
colnames( allStatData ) <- c( "interactionType", "dummy", "numberOfReads", "sample" )
allStatData <- allStatData[,colnames(allStatData) !="dummy"]

validInteractionNumbers <- filter(allStatData, grepl( "pairs$", allStatData$interactionType))
validInteractionNumbers <- mutate( validInteractionNumbers,
       fractionOfReads=validInteractionNumbers$numberOfReads / totalReadNum[validInteractionNumbers$sample] )
p1 <- ggplot(data.frame( totalReadPairs=totalReadNum, sample=names( totalReadNum ) ),
       aes(sample, totalReadPairs)) +
    geom_bar( stat="identity" ) + xlab("") +
    ylab("Number of sequenced fragments") +
    theme(axis.text.x=element_text(angle=90, vjust=.5))
p2 <- ggplot( validInteractionNumbers,
       aes( sample, fractionOfReads, fill=interactionType ) ) +
    geom_bar(stat="identity") +
    ylab("Fraction of sequenced fragments") + xlab("") +
    scale_y_continuous(breaks=seq(0, 1, 0.1), limits=c(0,1)) +
    guides(fill=guide_legend(title="Interactions")) +
    theme(axis.text.x=element_text(angle=90, vjust=.5))
shortLongInteractions <-
    filter( allMergeData,
           interactionType %in% c("cis_shortRange", "cis_longRange", "trans_interaction") ) %>%
        mutate( fractionOfReads=numberOfReads/totalReadNum[.$sample] )
duplicatedReadInteractions <- merge(
    filter( allMergeData, interactionType == "valid_interaction_rmdup" ),
    filter( validInteractionNumbers, interactionType == "Valid_interaction_pairs" ), by="sample") %>%
    mutate( numberOfReads=numberOfReads.y - numberOfReads.x ) %>%
    mutate( fractionOfReads=numberOfReads/totalReadNum[.$sample] )
duplicatedReadInteractions$fractionOfReads <-
    duplicatedReadInteractions$numberOfReads / totalReadNum[duplicatedReadInteractions$sample]
duplicatedReadInteractions$interactionType <- "Duplicates"
shortLongInteractions <- rbind( shortLongInteractions,
      duplicatedReadInteractions[,c("interactionType", "numberOfReads", "fractionOfReads", "sample")] )
shortLongInteractions$interactionType <-
    factor( shortLongInteractions$interactionType,
           levels=c("Duplicates", "trans_interaction", "cis_shortRange", "cis_longRange") )
levels( shortLongInteractions$interactionType ) <-
    c( "duplicates", "trans", "short range cis", "long range cis" )
colnames(shortLongInteractions)[1] <- "Valid interactions"
p3 <- ggplot( shortLongInteractions, aes( sample, fractionOfReads, fill=`Valid interactions` ) )  +
    geom_bar( stat="identity" ) + xlab("") +
    scale_y_continuous(breaks=seq(0, 1, 0.1), limits=c(0,1)) +
    ylab("Fraction of sequenced fragments") +
    theme(axis.text.x=element_text(angle=90, vjust=.5))    

bar <- .8 + .12*length(samples)
save_plot( p1, file=sprintf("%s-numberOfPairs.png", prefix), base_height=4, base_aspect_ratio=bar )
save_plot( p2, file=sprintf("%s-pairValidity.png", prefix), base_height=4, base_aspect_ratio=bar+.4 )
save_plot( p3, file=sprintf("%s-cisTrans.png", prefix), base_height=4, base_aspect_ratio=bar+.4 )
