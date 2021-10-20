#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { INPUT_CHECK } from './subworkflows/local/input_check'
include { FASTQC } from './modules/nf-core/modules/fastqc/main' addParams( options: [:] )
include { TRIMGALORE } from './modules/nf-core/modules/trimgalore/main' addParams( options: [:] )
include { BWAMEM2_MEM } from './modules/nf-core/modules/bwamem2/mem/main' addParams( options: [:] )
include { SAMTOOLS_VIEW } from './modules/nf-core/modules/samtools/view/main' addParams( options: [:] )
include { SAMTOOLS_SORT } from './modules/nf-core/modules/samtools/sort/main' addParams( options: [:] )
//include { SAMTOOLS_FASTQ } from './modules/nf-core/modules/samtools/fastq/main' addParams( options: [:] )
//include { SPADES } from './modules/nf-core/modules/spades/main' addParams( options: [:] )
//include { KRAKEN2_KRAKEN2 } from './modules/nf-core/modules/kraken2/kraken2/main' addParams( options: [:] )
//include { SEQTK_SUBSEQ } from './modules/nf-core/modules/seqtk/subseq/main' addParams( options: [:] )
//include { BLAST_BLASTN } from './modules/nf-core/modules/blast/blastn/main' addParams( options: [:] )
//include { CENTRIFUGE } from './modules/local/centrifuge/centrifuge' addParams( options: [:] )

workflow {
	INPUT_CHECK()
	ch_raw_short_reads = INPUT_CHECK.out.raw_short_reads
	//FASTQC(ch_raw_short_reads)
	TRIMGALORE(ch_raw_short_reads)
	//BWAMEM2_MEM(TRIMGALORE.out.reads, '/sw/data/uppnex/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex') 
	//SAMTOOLS_SORT(BWAMEM2_MEM.out)
	//SAMTOOLS_VIEW(SAMTOOLS_SORT.out)
	//SAMTOOLS_SORT(SAMTOOLS_VIEW.out)
	//SAMTOOLS_FASTQ(SAMTOOLS_SORT.out)
	//SPADES(SAMTOOLS_FASTQ.out)
	//KRAKEN2_KRAKEN2(SPADES.out)
	//SEQTK_SUBSEQ(KRAKEN2_KRAKEN2.out)
	//BLAST_BLASTN(SEQTK_SUBSEQ.out)
	//CENTRIFUGE(BLAST_BLASTN.out)
	//SEQTK_SUBSEQ(CENTRIFUGE.out)
	//BLAST_BLASTN(SEQTK_SUBSEQ.out)
	//SEQTK_SUBSEQ(BLAST_BLASTN.out)
	//BLAST_BLASTN(SEQTK_SUBSEQ.out)
}
