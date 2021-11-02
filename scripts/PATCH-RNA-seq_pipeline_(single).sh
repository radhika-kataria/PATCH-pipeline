#!/bin/bash -l
#SBATCH --output=/scratch/users/%u/%j.out
#SBATCH --job-name=blastn-redo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20000
#SBATCH --output=blastn_2021.out.%j
sample=$(printf '%q\n' "${PWD##*/}")
module load apps/bwa/0.7.17-singularity 
module load apps/bedtools2/2.29.0
module load apps/samtools/1.10.0-singularity
module load apps/openjdk 
module load apps/trimmomatic/0.39
module load apps/fastqc/0.11.8

#PATCH RNA-seq script (single end)
#Convert bam file to fastq format using bedtools
bedtools bamtofastq -i *.bam -fq ${sample}.fq

echo ${sample} bam to fastq finished 

#Check quality of fastq file using fastqc
mkdir fastqc
fastqc ${sample}.fq --outdir /home/radhika/rna_seq/colorectal/${sample}/fastqc

#Trim the fastq file with trimmomatic
trimmomatic SE -phred33 \
	${sample}.fq \
	${sample}_trimmed.fq \
	LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:28

#Check quality of fastq file post trimming
fastqc ${sample}_trimmed.fq --outdir /home/radhika/rna_seq/colorectal/${sample}/fastqc

#Align trimmed fastq files to the human ref genome GRCh38 using hisat2
hisat2 -x /home/radhika/datadrive/GRCh38/index/GRCh38_index \
		--known-splicesite-infile /home/radhika/datadrive/GRCh38/index/GRCh38_splice_sites.txt   \
		-p 12 \
		-U ${sample}_28_trimmed.fq | samtools view -bS - > ${d}.hisat.bam 

echo ${sample} Hisat2 alignment finished 

#Extract "unmapped"/ non-human reads using samtools flags
samtools view -F 4 ${sample}.hisat.bam > ${sample}_unmapped.bam
echo ${sample} Unmapped extracted

#Sort the bam file using samtools
samtools sort ${sample}_unmapped.bam > ${sample}_unmapped_sorted.bam

#Convert bam file to fastq format using bedtools
bedtools bamtofastq -i ${sample}_unmapped_sorted.bam -fq ${sample}_unmapped.fq

echo ${sample} Unmapped bam to fastq

#Denovo assemlbly of unmapped reads using SPAdes
/scratch/users/k1802884/tools/SPAdes-3.14.1-Linux/bin/rnaspades.py -s ${sample}_unmapped_S.fq --only-assembler -o spades

#Run kraken for pathogen classification 

mkdir kraken

/scratch/users/k1802884/tools/kraken2-2.0.8-beta/kraken2 --use-names  --db /scratch/users/k1802884/azure/radhika/kraken/refseq/standard spades/transcripts.fasta \
--report kraken/kraken_report.txt --classified-out kraken/kraken_classifications.txt >> kraken/output_kraken.txt

grep -e 'pathogen_of_interest' kraken/output_kraken.txt | awk '{print $2}' > kraken/pathogen_nodes.txt

/scratch/users/k1802884/tools/seqtk/seqtk subseq spades/transcripts.fasta kraken/pathogen_nodes.txt > kraken/pathogen_sequences.fasta

/scratch/users/k1802884/tools/ncbi-blast-2.10.1+/bin/blastn -query kraken_2021/fuso_seq_21.fasta -db pathogens_of_interest_ref_database.fa \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
-max_target_seqs 1 -max_hsps 1 -out kraken/pathogen_gene_annotation.blastn


#Using 
mkdir centrifuge
/scratch/users/k1802884/azure/radhika/tools/centrifuge/centrifuge -x /scratch/users/k1802884/azure/radhika/centrifuge/p_compressed+h+v \
-f spades/transcripts.fasta \
--report-file centrifuge/centrifuge_report.txt -S centrifuge/centrifuge_output.txt

grep -e 'pathogen_of_interest' centrifuge/centrifuge_report.txt | awk '{print $3}' | sort | uniq > centrifuge/tax_id_list.txt
awk -F' ' 'NR==FNR{c[$1]++;next};c[$3]' tax_id_list.txt entrifuge/centrifuge_output.txt | awk '{print $1}' | sort -u | uniq > centrifuge/pathogen_nodes.txt
/scratch/users/k1802884/tools/seqtk/seqtk subseq spades/transcripts.fasta centrifuge/pathogen_nodes.txt > centrifuge/pathogen_sequences.fasta

/scratch/users/k1802884/tools/ncbi-blast-2.10.1+/bin/blastn -query centrifuge/pathogen_sequences.fasta -db /pathogens_of_interest_ref_database.fa \
	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
	-max_target_seqs 1 -max_hsps 1 -out centrifuge/pathogen_gene_annotation.blastn

#Run BLASTn for pathogen classification 

mkdir blastn

/scratch/users/k1802884/tools/ncbi-blast-2.10.1+/bin/blastn -db /scratch/users/k1802884/azure/radhika/blast_nt/nt -query spades/transcripts.fasta -num_threads 16 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
-max_target_seqs 1 -max_hsps 1 -out blastn/blastn_output.txt

grep -e 'pathogen_of_interest' blastn/blastn_output.txt | awk '{print $1}' | sort -u | uniq > blastn/pathogen_nodes.txt
/scratch/users/k1802884/tools/seqtk/seqtk subseq spades/transcripts.fasta blastn/pathogen_nodes.txt > blastn/pathogen_reads.fasta

/scratch/users/k1802884/tools/ncbi-blast-2.10.1+/bin/blastn -query blastn/pathogen_reads.fasta -db pathogens_of_interest_ref_database.fa \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
-max_target_seqs 1 -max_hsps 1 -out blastn/pathogen_gene_annotation.blastn

