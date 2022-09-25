# Make a small 100k read pairs test dataset
zcat /data/aryee/corri/christian_data/Gel-ex-rep1_Primer7_S3_L001_R1_001.fastq.gz | head -n 400000 | gzip > ~/small_rcmc_r1.fq.gz
zcat /data/aryee/corri/christian_data/Gel-ex-rep1_Primer7_S3_L001_R2_001.fastq.gz | head -n 400000 | gzip > ~/small_rcmc_r2.fq.gz

# Make a small chr10 subset genome
mkdir genome
cd genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr10.fa.gz
gunzip chr10.fa.gz
echo "chr10	87000000	91000000" > test.bed
echo "chr10:87000000-91000000 4000000" > test.chrom.sizes
bedtools getfasta -fi chr10.fa -bed test.bed -fo test.fa
rm chr10.fa*
rm test.bed
gzip test.fa
cd ..






