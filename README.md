# PAThogen CHaracterisation (PATCH) pipeline

PATCH pipeline implementation [Nextflow][url_nextflow] pipeline for processing host transcriptomic and genomic sequencing data to identify pathogen-derived reads to a functional level in cancer datasets.

The pipeline was written by the [Cancer Bioinformatics][url_cb] and [Translational Systems Biology][url_sb] group at [King's College London][url_kcl], UK.

![patch_pipeline][pipeline]
## Pipeline summary - Pathogen characterisation
1. After QC steps ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc), [`trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic)), sequencing reads were aligned to the host reference genome ([`HISAT2`](https://github.com/DaehwanKimLab/hisat2) for transcriptomics and [`bwa`](http://bio-bwa.sourceforge.net/) for whole genome sequencing).
2. Extracting unaligned reads ([`SAMtools`](http://www.htslib.org/doc/samtools.html))
3. De novo assembley of host unmapped reads ([`SPAdes`](https://github.com/ablab/spades))
4. Pathogen classification using 3 tools: [`Kraken2`](https://ccb.jhu.edu/software/kraken2/), [`BLASTn`](https://www.ncbi.nlm.nih.gov/books/NBK279690/), [`Centrifuge`](https://ccb.jhu.edu/software/centrifuge/manual.shtml) where the consensus of two or more is taken forward.
5. Classified reads from the pathogen of interest are extracted and functionally annotated using [`BLASTn`](https://www.ncbi.nlm.nih.gov/books/NBK279690/) against indexed RefSeq for transcripts/genomes of the pathogen of interest. 

## Pipeline summary - Pathogen integration  
1. A custom combined reference genome is created using the host and pathogen of interest reference genoemes ([`bwa`](http://bio-bwa.sourceforge.net/))
2. Whole genome sequencing data is aligned to the combined reference genome ([`bwa`](http://bio-bwa.sourceforge.net/))
3. Discordant reads where one read maps to the pathogen of interest and it's mate to the host reference genome are extracted ([`SAMtools`](http://www.htslib.org/doc/samtools.html))
4. Filtering of duplicated reads and alignemnt quality (MAPQ scores), ([`Picard tools`](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates), [`SAMtools`](http://www.htslib.org/doc/samtools.html))
5. As before - classified reads from the pathogen of interest are extracted and functionally annotated using [`BLASTn`](https://www.ncbi.nlm.nih.gov/books/NBK279690/) against indexed RefSeq for transcripts/genomes of the pathogen of interest. 
6. Discordant read coordinates extracted ([`Bedtools`](https://bedtools.readthedocs.io/en/latest/))  


## Credits

The pipeline was written by the [Cancer Bioinformatics][url_cb] and [Translational Systems Biology][url_sb] group at [King's College London][url_kcl], UK.

Pipeline development and implementation by [Radhika Kataria](radhika.kataria@kcl.ac.uk), [Jelmar Quist](jelmar.quist@kcl.ac.uk), [Anargyros Megalios](argymeg@gmail.com), [Thomas Hardiman](thomas.hardiman@kcl.ac.uk), [Theo Portlock](theo.portlock@kcl.ac.uk) and [Sunjae Lee](sunjae.lee@kcl.ac.uk). 

Study concept and design [Radhika Kataria](radhika.kataria@kcl.ac.uk), [Anita Grigoriadis](anita.grigoriadis@kcl.ac.uk), [Saeed Shoaie](saeed.shoaie@kcl.ac.uk)

[url_cb]: http://cancerbioinformatics.co.uk/
[url_sb]: https://www.kcl.ac.uk/people/saeed-shoaie-1
[url_fastqc]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc
[url_kcl]: https://www.kcl.ac.uk/
[url_nextflow]: http://www.nextflow.io
[url_nextflow_tuto]: http://www.nextflow.io/docs/latest/getstarted.html#get-started

[pipeline]: https://github.com/radhika-kataria/PATCH-pipeline/blob/main/PATCH-github-image.png

