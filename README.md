# PcTF differential expression
Code for analyzing RNAseq data from breast cell lines expressing PcTF, and code for making figures. 

### Citation and data download:
Olney, K. C., Nyer, D. B., Wilson Sayres, M. A. & Haynes, K. Activation of tumor suppressor genes in breast cancer cells by a synthetic chromatin effector. (2017). bioRxiv 186056; doi: https://doi.org/10.1101/186056 

Fastq and CuffDiff output files are available at the National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) database (accession GSE103520)

### Methods:
RNA-seq reads were quality-checked using FastQC before and after trimming and filtering. Reads shorter than 50 bases, and base reads below the PHRED-scaled threshold quality of 10 at the 5’ end or 25 at the trailing 3’ end as well as below the average quality of 30 (within a sliding window of 4 bases) were filtered or clipped using TrimmomaticSE (33). A combined reference genome index and dictionary for GRCH38.p7 (1-22, X, MT, and non-chromosomal sequences) plus the full coding region of PcTF was generated using Spliced Transcripts Alignment to Reference (STARv2.5.2b)and the picard tools (version 1.1.19). Trimmed RNA-seq reads were mapped, and splice junctions extracted, using STARv2.5.2b read aligner. Bamtools2.4.0 was used to determine alignment quality (‘stats’ command), sort reads, mark duplicates, add read groups, and to index the BAM read files. CuffDiff (Cufflinks package), was used to run pairwise comparisons to identify significant (q ≤ 0.05) differential expression. CummeRbund was used to calculate distances between features (JSD plots). R ggplot2 and VennDiagrams were used to generate heat maps and Venn diagrams respectively

### Contents:
1. Download or obtain data 
2. FastQC to check the quality of the raw reads
3. Trim fastq files for quality and to remove adaptors. 
4. FastQC to check the quality of the trimmed reads
5. Obtain reference genome and gene annotation files
6. Generate genome indexes
7. STAR: map transcript reads to the reference genome and identify splice junctions 
8. check quality of raw BAM files
9. Sort BAM files, to be in the same order as the reference for downstream analysis 
10. check quality of sorted BAM files
11. Mark duplicates, duplicates will not be removed but will be marked for quality checks
12. check quality of mark duplicates BAM files
13. Add read groups to BAM files, this done to keep the sample ids organized when creating the merged vcf file
14. check quality of add read group BAM files
15. Index BAM files
16. Identify genes that differentially expressed 


### Publicly available packages:
fastqc		http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Trimmomatic	http://www.usadellab.org/cms/?page=trimmomatic

STAR		https://github.com/alexdobin/STAR

bamtools	https://github.com/pezmaster31/bamtools

cuffdiff	http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/  

cummeRbund	http://compbio.mit.edu/cummeRbund/manual_2_0.html 

ggplot2		http://ggplot2.tidyverse.org/reference/ 

VennDiagram	https://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf  

--------------------------------------
RNAseq processing for differential expression analysis 
--------------------------------------
#### 1. Download data
GEO accession number
Fastq and CuffDiff output files are available at the National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) database (accession GSE103520)

#### 2. Create and view fastqc reports
Fastqc reads raw sequence data from high throughput sequencers and runs a set of quality checks to produce a report. Best reports are those whose "per base sequence quality" are included in the green area of the graph & kmer content is good or average.

 	Example command: $fastqc sampleID.fastq
	
Reports were saved in fastq_files directory in Project /Project/fastq_files
This command will create two outputs: an .html file & an .zip file. Will output sampleID_fastqc.html and sampleID_fastqc.zip files

#### 3. Trim raw fastq files for quality and to remove adaptors 
The selection of trimming steps and their associated parameters are supplied on the command line.
The parameters selected were slidingwindow:4:30 leading10 trailing25 minlen50 phred33. They were chosen based on a better per sequence quality base and kmer content.
		
	Example command: $java -jar trimmomatic-0.36.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:/adapters/TruSeq3-SE.fa:2:30:10 LEADING:10 TRAILING:25 SLIDINGWINDOW:4:30 MINLEN:50

#### 4. Create and view fastqc reports for the trimmed fastq files

	Example command: $fastqc sampleID_minlen50_sliding430_leading30_trailing40.fastq

#### 5. Obtain reference genome and gene annotation file 
Obtain reference genome and gene annotation file to be used for mapping reads. 
Version GRCh38.p7 from gencode. http://www.gencodegenes.org/releases/current.html
Nucleotide sequence of the GRCh38.p7 genome assembly version on all regions, including reference chromosomes, scaffolds, assembly patches and haplotypes

 	Example command $wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.p7.genome.fa.gz

Comprehensive gene annotation gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz Contains the comprehensive gene annotation on the reference chromosomes, scaffolds, assembly patches and alternate loci (haplotypes)

 	Example command $wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz

#### 6. Generate reference genome index and dictionary 
A .dict dictionary of the contig names and sizes and a .fai fasta index file allow efficient random access to the reference bases for downstream analysis and mapping 

 Index reference genome using STAR
 
 	Example command: $STAR --runMode genomeGenerate --genomeDir /GRCh38.p7/ --genomeFastaFiles GRCh38.p7.genome.fa.gz --sjdbGTFfile gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz --runThreadN 8
 

Create reference dictionrary using Picard tools CreateSequenceDictionary 
 
 	Example command: $java -Xmx14g -jar picard.jar CreateSequenceDictionary R=GRCh38.p7.genome.fa O=GRCh38.p7.genome.fa.dict

#### 7. STAR: map transcript reads to the reference genome, output as a .bam
The user supplies the genome files generated in the pervious step (generate genome indexes), as well as the RNA-seq reads (sequences) in the form of FASTA or FASTQ files. STAR maps the reads to the genome, and writes several output files, such as alignments (SAM/BAM), mapping summary statistics, splice junctions, unmapped reads, signal (wiggle) tracks etc. 
Mapping is controlled by a variety of input parameters (options)

STAR highly recommends using --sjdbGTFfile which specifies the path to the file with annotated transcripts in the standard GTF format. STAR will extract splice junctions from this file and use them to greatly improve accuracy of the mapping. 
Compatibility with Cufflinks/Cuffdiff.For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute, which STAR will generate with --outSAMstrandField intronMotif option. As required, the XS strand attribute will be generated for all alignments that contain splice junctions. The spliced alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions) will be suppressed.
 
STAR: maps fastq files to the reference genome and identifies splice junctions.
 
 	Example comamnd: $STAR --genomeDir Project_genome --genomeLoad LoadAndKeep --readFilesIn sampleID.fastq --outSAMtype BAM Unsorted --outFileNamePrefix sampleID_pass1. --runThreadN 8

 #### 8. Check initial quality stats on .bam files 
 Will print basic statistics from input BAM file(s)

 	Example command: $bamtools stats -in sampleID_pass2.Aligned.out.bam > sampleID_pass2.txt

#### 9. Sort bam files
An appropriate @HD-SO sort order header tag will be added
 
 	Example command: $bamtools sort -in sampleID.bam -out sampleID.sorted.bam

 #### 10. Mark duplicates & Remove duplicates
 Mark duplicates: "Flags" where the duplicate reads are
 
 	Example command: $java -Xmx8g -jar picard.jar MarkDuplicates INPUT=sampleID.sorted.bam OUTPUT=sampleID.sorted.markdup.bam METRICS_FILE=sampleID.markdup.picardMetrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

#### 11. Add or replace read groups
For each sample, add a read group to the mark duplicate BAM files (a read group is a "tag" such as a sample ID)
 
	Example command: $java -Xmx8g -jar picard.jar AddOrReplaceReadGroups INPUT=sampleID.sorted.rmdup.bam OUTPUT=sampleID.sorted.rmdup.addReadGr.bam RGLB=sampleID RGPL=machineUsed RGPU=laneUsed RGSM=sampleName RGCN=location RGDS=species VALIDATION_STRINGENCY=LENIENT

#### 12. Generate statistics on final BAM sample files
For each sample get the read stats for the remove duplicates and add read groups BAM files.
Statistics on the BAM files should be the same as before the previous step when read groups were modified.
Compare stat results for each sample to the markdup.bam file (sanity check: is there the same number of reads as the original BAM file?) Compare stat results for each sample to the original bam file (sanity check: is there the same number of reads as the original BAM file?)If there is more than 15% of the reads being marked as duplicates may need to consider removing that sample
 
 	Example command: bamtools stats -in sampleID.sorted.markdup.addReadGr.bam

#### 13. Index BAM files 
For each sample index the processed BAM files that are sorted, have marked duplicates, and have read groups added. These will be used to identify callable loci. 
Indexing is used to "sort" by chromosome and region 
Output will be sampleID.sorted.markdup.addReadGr.bam.bai
 
 	Example command: $bamtools index -in sampleID.sorted.markdup.addReadGr.bam

#### 14. Differential expression using Cuffdiff
Identify differentially expressed genes and transcripts using cuffdiff. Cuffdiff is used to find significant changes in transcript expression, splicing, and promoter use. Cuffdiff uses a gene annotation file downloaded with a sbatch script from Gencode GRCh38
 
	Example command: $cuffdiff -use-sample-sheet -o diff_out -b reference.fa -p 8 --library-type fr-firststrand -L set1,set2 -u refence.gtf sampleSet.txt
