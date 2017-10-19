# Run on erisone

module load aryee/samtools-1.3.1
BWA=/apps/lab/aryee/bwa-0.7.12/bwa
BT2_HOME=/apps/lab/aryee/bowtie2-2.3.3
GENOME_FASTA=/data/aryee/pub/genomes/grch38/bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa



########################################################################################
########################################################################################

##  100kb region

# Create a fasta for a portion of chr19
mkdir -p test_data/genome_fasta
cd test_data/genome_fasta
samtools faidx $GENOME_FASTA 19:10000000-10100000 > hs_chr19_100kb.fa
sed -i 's/>19:1/>19 1/' hs_chr19_100kb.fa

# Make a bowtie2 index
mkdir -p test_data/bowtie2_index
cd test_data/bowtie2_index
$BT2_HOME/bowtie2-build ../genome_fasta/hs_chr19_100kb.fa hs_chr19_100kb
tar cfz hs_chr19_100kb_bt2.tar.gz *
rm *.bt2

# Get names of reads that fall in small region of chr19
for fq in /PHShome/ma695/work/projects/gbm_topology/input/rao/fastq/imr90/*.fastq.gz;
do
    echo $fq;
    bam=`basename ${fq/.fastq.gz/.bam}`
    readnames=`basename ${fq/.fastq.gz/.readnames}`
    bsub -sla ccr_sc -q big-multi -n 12  "$BWA mem -t 12 $GENOME_FASTA $fq | awk '\$3==19 && \$4 > 10000000 && \$4 < 10100000 { print \$1 }' > $readnames";
done

cat *.readnames | sort | uniq > imr90_chr19_100kb.readnames
wc -l imr90_chr19_100kb.readnames
rm SRR*.readnames


# Extract reads that fall in this small region
for fq in /PHShome/ma695/work/projects/gbm_topology/input/rao/fastq/imr90/*.fastq.gz;
do
    echo $fq;
    small_fq=`basename ${fq/.fastq.gz/.small.fastq.gz}`
    bsub -sla ccr_sc -q medium "zcat $fq | grep -F -A3 --no-group-separator -f imr90_chr19_100kb.readnames | gzip -c > $small_fq";
done

mkdir imr90_chr19_100kb
mv *.fastq.gz imr90_chr19_100kb




########################################################################################
########################################################################################

##  10Mb region

# Create a fasta for a portion of chr19
samtools faidx $GENOME_FASTA 19:1-11000000 > hs_chr19_11mb.fa
sed -i 's/>19:1/>19 1/' hs_chr19_11mb.fa


# Get names of reads that fall in small region of chr19
for fq in /PHShome/ma695/work/projects/gbm_topology/input/rao/fastq/imr90/*.fastq.gz;
do
    echo $fq;
    bam=`basename ${fq/.fastq.gz/.bam}`
    readnames=`basename ${fq/.fastq.gz/.readnames}`
    bsub -sla ccr_sc -q big-multi -n 12  "$BWA mem -t 12 $GENOME_FASTA $fq | awk '\$3==19 && \$4 < 10000000 { print \$1 }' > $readnames";
done

cat *.readnames | sort | uniq > imr90_chr19_1-10mb.readnames
wc -l imr90_chr19_1-10mb.readnames
rm SRR*.readnames


# Extract reads that fall in this small region
for fq in /PHShome/ma695/work/projects/gbm_topology/input/rao/fastq/imr90/*.fastq.gz;
do
    echo $fq;
    small_fq=`basename ${fq/.fastq.gz/.small.fastq.gz}`
    bsub -sla ccr_sc -q medium "zcat $fq | grep -F -A3 --no-group-separator -f imr90_chr19_1-10mb.readnames | gzip -c > $small_fq";
done

mkdir imr90_chr19_1-10mb
mv *.fastq.gz imr90_chr19_1-10mb



### Test hicpro image

docker run --rm -it -v /Users/maryee/Dropbox/projects/topology_tools/test_data:/test_data aryeelab/hicpro

mkdir /bowtie2_index
tar zxvf /test_data/bowtie2_index/hs_chr19_100kb_bt2.tar.gz -C /bowtie2_index

/HiC-Pro_2.9.0/bin/HiC-Pro  -i /test_data/hic/ -o /hicpro_out -c /PHShome/ma695/work/projects/glioma_topology/config-hicpro-mboi-extr1r3_mincisdist1000.txt -p
