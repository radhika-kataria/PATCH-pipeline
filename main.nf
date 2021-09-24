#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { FASTQC } from './modules/nf-core/modules/fastqc/main' addParams( options: [:] )
include { FASTQC_UMITOOLS_TRIMGALORE } from './modules/nf-core/modules/trimgalore/main' addParams( options: [:] )
include { BWA_MEM } from './modules/nf-core/modules/bwa/mem/main' addParams( options: [:] )
include { SAMTOOLS_VIEW } from './modules/nf-core/modules/samtools/view/main' addParams( options: [:] )
include { SAMTOOLS_SORT } from './modules/nf-core/modules/samtools/sort/main' addParams( options: [:] )
include { SAMTOOLS_FASTQ } from './modules/nf-core/modules/samtools/fastq/main' addParams( options: [:] )
include { SPADES } from './modules/nf-core/modules/spades/main' addParams( options: [:] )
include { KRAKEN2_KRAKEN2 } from './modules/nf-core/modules/kraken2/kraken2/main' addParams( options: [:] )
include { SEQTK_SUBSEQ } from './modules/nf-core/modules/seqtk/subseq/main.nf' addParams( options: [:] )
include { BLAST_BLASTN } from './modules/nf-core/modules/blast/blastn/main.nf' addParams( options: [:] )
include { CENTRIFUGE } from './modules/local/centrifuge/centrifuge.nf' addParams( options: [:] )

process BAM2FASTQ {
	memory '6GB'
	cpus '1'
	time '12h'
	container 'opengenomics/biobambam2'
	//scratch true
	
	input:
	path(name)

	output:
	path '*.{1,2}.fastq'
	val(name)

	shell:
	'''
	bamtofastq \
	    collate=1 \
	    exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
	    filename=!{name} \
	    inputformat=bam \
	    F=!{name}_R1.fastq.gz \
	    F2=!{name}_R2.fastq.gz \
	    S=!{name}_s.fq.gz \
	    0=!{name}_0.fq.gz \
	    02=!{name}_02.fq.gz \
	    tryoq=1 \
	    gz=1 \
	    exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
	    level=5
	'''
}
workflow {
	ch_bams = Channel.fromPath( params.input )
	BAM2FASTQ(ch_bams)
	FASTQC(BAM2FASTQ.out)
	FASTQC_UMITOOLS_TRIMGALORE(BAM2FASTQ.out)
	BWA_MEM(FASTQC_UMITOOLS_TRIMGALORE.out)
	SAMTOOLS_VIEW(BWA_MEM.out)
	SAMTOOLS_SORT(SAMTOOLS_VIEW.out)
	SAMTOOLS_FASTQ(SAMTOOLS_SORT.out)
	SPADES(SAMTOOLS_FASTQ.out)
	KRAKEN2_KRAKEN2(SPADES.out)
	SEQTK_SUBSEQ(KRAKEN2_KRAKEN2.out)
	BLAST_BLASTN(SEQTK_SUBSEQ.out)
	CENTRIFUGE(BLAST_BLASTN.out)
	SEQTK_SUBSEQ(CENTRIFUGE.out)
	BLAST_BLASTN(SEQTK_SUBSEQ.out)
	SEQTK_SUBSEQ(BLAST_BLASTN.out)
	BLAST_BLASTN(SEQTK_SUBSEQ.out)
}
