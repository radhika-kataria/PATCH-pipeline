#!/bin/bash -l
#SBATCH --output=/scratch/users/%u/%j.out
#SBATCH --job-name=blastn-redo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20000
#SBATCH --output=blastn_2021.out.%j

#discordant-reads-pipeline

module load apps/bwa/0.7.17-singularity 
module load apps/bedtools2/2.29.0
module load apps/samtools/1.10.0-singularity
module load apps/openjdk 
module load apps/trimmomatic/0.39
module load apps/fastqc/0.11.8

#Convert bam to fastq
export bam2fq=/scratch/users/k1802884/tools/biobambam2/2.0.87-release-20180301132713/x86_64-etch-linux-gnu/bin/bamtofastq 
bam=$(echo *.bam)
$bam2fq \
    collate=1 \
    exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
    filename=${bam} \
    inputformat=bam \
    F=${bam%%.bam}_F1.fq.gz \
    F2=${bam%%.bam}_F2.fq.gz \
    S=${bam%%.bam}_s.fq.gz \
    0=${bam%%.bam}_0.fq.gz \
    02=${bam%%.bam}_02.fq.gz \
    tryoq=1 \
    gz=1 \
    exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
    level=5

#Fastq quality control & trimming 
mkdir fastqc
fastqc  *_F1.fq.gz --outdir fastqc
fastqc  *_F2.fq.gz --outdir fastqc

for fastq in *_F1.fq.gz 
do
        base=${fastq%%_F1.fq.gz}
        trimmomatic PE \
        ${base}_F1.fq.gz \
        ${base}_F2.fq.gz \
        ${base}_trimmed_F1.fq.gz \
        ${base}_F1_UP.fq.gz \
        ${base}_trimmed_F2.fq.gz \
        ${base}_F2_UP.fq.gz \
        LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:28 
done

fq1=*_trimmed_F1.fq.gz
fq2=*_trimmed_F2.fq.gz

#Align to combined reference genome
mkdir discordant
bwa mem /scratch/users/k1802884/ensembl/dna/index/combined.fa ${fq1} ${fq2} -t 10 | samtools sort -o discordant/fuso_human_sorted.bam - 

#Extract discordant reads
samtools view -F 14 -h -q 20 discordant/fuso_human_sorted.bam -@ 10 > discordant/discordant.bam

#remove duplicate reads
java -jar /scratch/users/k1802884/azure/radhika/tools/picard.jar MarkDuplicates \
      I=discordant/discordant.bam \
      O=discordant/discordant_remove_duplicates.bam \
      M=discordant/marked_dup_metrics.txt \
      REMOVE_DUPLICATES=true

#get pathogen-human discordant reads
samtools view discordant/discordant_remove_duplicates.bam| grep -f pathogen_fasta_headers.txt | awk '{for(i=7;i<=NF;i++){if($i~/^NC=/){a=$i}} print $0,a}' > discordant/discordant_host_pathogen_reads.txt
awk '{print $1}' discordant/discordant_host_pathogen_reads.txt > discordant/discordant_readname_extract.txt
#Extract co-ordinates
bedtools bamtobed -i discordant/discordant_remove_duplicates.bam  | grep -f discordant/discordant_readname_extract.txt > discordant/coordinates.txt

