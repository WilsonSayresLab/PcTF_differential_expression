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

      fastqc            http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
      Trimmomatic       http://www.usadellab.org/cms/?page=trimmomatic
      STAR		https://github.com/alexdobin/STAR
      bamtools          https://github.com/pezmaster31/bamtools
      cuffdiff          http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/  
      cummeRbund	http://compbio.mit.edu/cummeRbund/manual_2_0.html 
      ggplot2		http://ggplot2.tidyverse.org/reference/ 
      VennDiagram	https://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf  
--------------------------------------
RNAseq processing for differential expression analysis 
--------------------------------------
### 1. Download data
GEO accession number
Fastq and CuffDiff output files are available at the National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) database (accession GSE103520)

### 2. Create and view fastqc reports
 Fastqc reads raw sequence data from high throughput sequencers and runs a set of quality checks to produce a report. 
 Best reports are those whose "per base sequence quality" are included in the green area of the graph & kmer content is good or average.

 	Example command: $fastqc sampleID.fastq
	
 Reports were saved in fastq_files directory in Project /Project/fastq_files
 This command will create two outputs: an .html file & an .zip file. Will output sampleID_fastqc.html and sampleID_fastqc.zip files

### 3. Trim raw fastq files for quality and to remove adaptors 
The selection of trimming steps and their associated parameters are supplied on the command line.
The parameters selected were slidingwindow:4:30 leading10 trailing25 minlen50 phred33. They were chosen based on a better per sequence quality base and kmer content.
		
	Example command: $java -jar trimmomatic-0.36.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:/adapters/TruSeq3-SE.fa:2:30:10 LEADING:10 TRAILING:25 SLIDINGWINDOW:4:30 MINLEN:50

--------------------------------------
 5. Create and view fastqc reports for the trimmed fastq files
--------------------------------------
 Fastqc reads sequence data from high throughput sequencers and runs a set of quality checks to produce a report. 
 Best reports are those whose "per base sequence quality" are included in the green area of the graph & kmer content is good or average.

 Example command: $fastqc sampleID_minlen50_sliding430_leading30_trailing40.fastq
 fastqc                                                        Babraham bioinformatics program that that checks for quality of reads 
 sampleID.fastq                                                path and name of sampleID in fastq format, may also be in fastq.gz format

 Reports were saved in fastq_files directory in Project /Project/fastq_files
 This command will create two outputs: an .html file & an .zip file. Will output sampleID_fastqc.html and sampleID_fastqc.zip files

 SBATCH script -
sbatch NameOfSBATCHscript.sh 

 Move fastqc reports to desktop to visualize them as you can't open html in a terminal.
 Open new terminal as this will not work if logged into a HPC (high performance computing) cluster
 Example command: $scp user@saguaro.a2c2.asu.edu:/Project/fastq_files/sampleID_raw_fastqc.html /Users/Desktop/
 scp                                                           secure copy  (linux command)                   
 /path/to/fastqc.html                                          path to where the files are located
 /path/where/you/want/fastqc.html                              path to where you would like to copy the files to 

 Reports were saved in Desktop in folder /trimmed_FASTQC
--------------------------------------
 6. Obtain reference genome and gene annotation file 
--------------------------------------
 Obtain reference genome and gene annotation file to be used for mapping reads. 
 Use the most relavent current version GRCh38.p7 from gencode. http://www.gencodegenes.org/releases/current.html
 Genome sequence (GRCh38.p7)	ALL	GRCh38.p7.genome.fa.gz
 Nucleotide sequence of the GRCh38.p7 genome assembly version on all regions, including reference chromosomes, scaffolds, assembly patches and haplotypes
 The sequence region names are the same as in the GTF/GFF3 files
 Example command $wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.p7.genome.fa.gz
 wget                                                          stands for "web get"
 ftp                                                           stands for File transfer program
 /GRCh38.p7.genome.fa.gz										human reference genome in fasta format

 Comprehensive gene annotation	ALL	gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz
 It contains the comprehensive gene annotation on the reference chromosomes, scaffolds, assembly patches and alternate loci (haplotypes)
 This is a superset of the main annotation file
 Example command $wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz
 wget                                                          stands for "web get"
 ftp                                                           stands for File transfer program
 gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz			human reference gene annotation file in .gtf format 

 SBATCH script -
sbatch NameOfSBATCHscript.sh 

--------------------------------------
 7. Generate reference genome index and dictionary 
--------------------------------------
 A .dict dictionary of the contig names and sizes and a .fai fasta index file allow efficient random access to the reference bases for downstream analysis and mapping 
 You have to generate these files in order to be able to use a Fasta file as reference.
 In this step user supplied the reference genome sequences (FASTA files) and annotations(GTF file), from which STAR generate genome indexes that are utilized in the 2nd step. 
 The genome indexes are saved to disk and need only be generated once for each genome/annotation combination. 
 Note - files can not be be zipped (i.e file.gz) use $gunzip file.gz 
 Example command: $gunzip reference.gz
 gunzip														linux command to unzip zipped files
 reference.gz													path and name to zipped file

 Index reference genome using STAR
 Example command: $STAR --runMode genomeGenerate --genomeDir /GRCh38.p7/ --genomeFastaFiles GRCh38.p7.genome.fa.gz --sjdbGTFfile gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz --runThreadN 8
 STAR															path to calling STAR read aligner package
 --runMode 
 genomeGenerate												option directs STAR to run genome indices generation job.
 --genomeDir 
 /path/to/genomeDir
 --genomeFastaFiles 
 /path/to/genome/fasta 
 --sjdbGTFfile /path/to/annotations.gtf
 --sjdbOverhang ReadLength-1
 --runThreadN 													NumberOfThreads. Option defines the number of threads to be used for genome generation, it has to be set to the number of available cores on the server node.

 Create reference dictionrary using Picard tools CreateSequenceDictionary 
 Example command: $java -Xmx14g -jar picard.jar CreateSequenceDictionary R=GRCh38.p7.genome.fa O=GRCh38.p7.genome.fa.dict
 java 															java program and will need java installed in order to run
 -Xmx14g 														declares memory 
 -jar 															tells the program that the following tool is in .jar format
 picard.jar 													picard tools from the Broad Institute. A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.
 CreateSequenceDictionary										creates a sequence dictionary for a reference sequence. The output file contains a header but no SAMRecords, and the header contains only sequence records.
 R=GRCh38.p7.genome.fa 										R stands for reference. Path and name of reference genome 
 O=GRCh38.p7.genome.fa.dict									O stand for output. Path and name of the output reference genome dictionary

 SBATCH script -
sbatch NameOfSBATCHscript.sh 

--------------------------------------
 8. STAR: map transcript reads to the reference genome, output as a .bam
--------------------------------------
 STAR read aligner is a 2 pass process. 
 The user supplies the genome files generated in the pervious step (generate genome indexes), as well as the RNA-seq reads (sequences) in the form of FASTA or FASTQ files. 
 STAR maps the reads to the genome, and writes several output files, such as alignments (SAM/BAM), mapping summary statistics, splice junctions, unmapped reads, signal (wiggle) tracks etc. 
 Mapping is controlled by a variety of input parameters (options)

 STAR highly recommends using --sjdbGTFfile which specifies the path to the file with annotated transcripts in the standard GTF format. 
 Where STAR will extract splice junctions from this file and use them to greatly improve accuracy of the mapping. 
 While this is optional, and STAR can be run without annotations, using annotations is highly recommended whenever they are available.
 However this option should not be included for projects that include hybrids, as this might cause a bias towards the reference. 

 Compatibility with Cufflinks/Cuffdiff.
 For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute, which STAR will generate with --outSAMstrandField intronMotif option. 
 As required, the XS strand attribute will be generated for all alignments that contain splice junctions. The spliced alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions) will be suppressed.
 If you have stranded RNA-seq data, you do not need to use any specific STAR options. Instead, you need to run Cufflinks with the library option --library-type options. 
 For example, cufflinks ... --library-type fr-firststrand should be used for the standard dUTP protocol, including Illumina’s stranded Tru-Seq. This option has to be used only for Cufflinks runs and not for STAR runs.
 In addition, it is recommended to remove the non-canonical junctions for Cufflinks runs using --outFilterIntronMotifs RemoveNoncanonical.

 Unstranded vs. stranded RNAseq
 stranded vs. unstranded RNAseq data is where the strand specificity of origin for each transcript is defined. 
 Without strand information it is difficult and sometimes impossible to accurately quantify gene expression levels 
 for genes with overlapping genomic loci that are transcribed from opposite strands (Zhao et al. 2015)
 Stranded RNA-seq provides a more accurate estimate of transcript expression compared with non-stranded RNA-seq, 
 and is therefore the recommended RNA-seq approach for future mRNA-seq studies.
 If your data is unstranded then you will want to include the STAR options --outSAMstrandField intronMotif option.
 STAR will generate the XS strand attribute for all alignments that contain splice junctions. 
 The spliced alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions) will be suppressed. 
 If you have stranded RNA-seq data, you do not need to use any specific STAR options. Instead, you need to run Cufflinks with the library option --library-type options (Dobin et al. 2012)

 First pass: maps fastq files to the reference genome and identifies splice junctions.
 Example comamnd: $STAR --genomeDir Project_genome --genomeLoad LoadAndKeep --readFilesIn sampleID.fastq --outSAMtype BAM Unsorted --outFileNamePrefix sampleID_pass1. --runThreadN 8
 STAR 									                        path to calling STAR read aligner package
 --genomeDir 							                        Specifies path to the directory (henceforth called ”genome directory” where the genome indices are stored. This directory has to be created (with mkdir) before STAR run and needs to have writing permissions. 
 Project_genome												path and directory to reference genome
 --genomeLoad													mode of shared memory usage for the genome files
 LoadAndKeep													load genome into shared and keep it in memory after run
 --sjdbGTFfile													path to the GTF file with annotations
 --readFilesIn 						                        if pair end reads, include path to both reads with a " " space inbetween _1 _2, /geuvadis_fastq/sampleID_1.fastq /home/kcolney/map_geuvadis/geuvadis_fastq/sampleID_2.fastq 
 sampleID.fastq												name and path to sample in fastq format. If paired end samples include both pairs and separate with a space i.e (sample_1.fastq sample_2.fastq)
 --outSAMtype 							                        indicate which output format, BAM unsorted
 BAM Unsorted													output unsorted Aligned.out.bam file. The paired ends of an alignment are always adjacent, and multiple alignments of a read are adjacent as well. This ”unsorted” file can be directly used with downstream software such as HTseq, without the need of name sorting. The order of the reads will match that of the input FASTQ(A) files only if one thread is used
 --outSAMstrandField 											* This option depends on the data * For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute, which STAR will generate with --outSAMstrandField intronMotif option
 intronMotif													* This option depends on the data * As required, the XS strand attribute will be generated for all alignments that contain splice junctions. The spliced alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions) will be suppressed.
 --outFilterIntronMotifs 										* This option depends on the data * it is recommended to remove the non-canonical junctions for Cufflinks runs using
 RemoveNoncanonical											* This option depends on the data * remove the non-canonical junctions
 --outFileNamePrefix 					                        define the sample id prefix, sampleID_pass1. (bam, will be added by the STAR program)
 --runThreadN 							                        for computing purpose allocate the number of threads, 8 

 SBATCH script - 
sbatch NameOfSBATCHscript.sh 

 For the most sensitive novel junction discovery, it is recommend to run STAR in the 2-pass mode. It does not increase the number of detected novel junctions, but allows to detect more splices reads mapping to novel junctions. 
 The basic idea is to run 1st pass of STAR mapping with the usual parameters, then collect the junctions detected in the first pass, and use them as ”annotated” junctions for the 2nd pass mapping.
 Second pass: maps fastq files between each other and identifies splice junctions.
 STAR will utilize annotations formatted as a list of splice junctions coordinates in a text file: --sjdbFileChrStartEnd /path/to/sjdbFile.txt. 
 This file should contains 4 columns separated by tabs:
 Chr \tab Start \tab End \tab Strand=+/-/.
 Example command: $STAR --genomeDir Project_genome --readFilesIn sampleID_1.fastq --outSAMtype BAM Unsorted --outFileNamePrefix sampleID_pass2. --sjdbFileChrStartEnd /sampleID1_pass1.SJ.out.tab /sampleID2_pass1.SJ.out.tab /sampleID3_pass1.SJ.out.tab --runThreadN 14
 STAR 									                        STAR read aligner package
 --genomeDir 							                        define where the genome is location, /star_genome 
 Project_genome												path and directory to reference genome
 --genomeLoad													mode of shared memory usage for the genome files
 LoadAndKeep													load genome into shared and keep it in memory after run
 --sjdbGTFfile													path to the GTF file with annotations
 --readFilesIn 						                        if pair end reads, include path to both reads with a " " space inbetween _1 _2, /geuvadis_fastq/sampleID_1.fastq /home/kcolney/map_geuvadis/geuvadis_fastq/sampleID_2.fastq 
 sampleID_1.fastq												name and path to sample in fastq format. If paired end samples include both pairs and separate with a space i.e (sample_1.fastq sample_2.fastq)
 --outSAMtype 							                        indicate which output format, BAM unsorted
 BAM Unsorted													output unsorted Aligned.out.bam file. The paired ends of an alignment are always adjacent, and multiple alignments of a read are adjacent as well. This ”unsorted” file can be directly used with downstream software such as HTseq, without the need of name sorting. The order of the reads will match that of the input FASTQ(A) files only if one thread is used
 --sjdbFileChrStartEnd					                        path to the pass_1.SJ.out.tab files made in the first pass
 sampleID1_pass1.SJ.out.tab									4 columns separated by tabs: Chr \tab Start \tab End \tab Strand=+/-/. Here Start and End are first and last bases of the introns (1-based chromosome coordinates). This file can be used in addition to the --sjdbGTFfile, in which case STAR will extract junctions from both files.
 sampleID2_pass1.SJ.out.tab									List all the samples from the first pass or all the samples in a group (i.e population, cases and controls, hybrids, males and females)
 --outFileNamePrefix 					                        define the sample id prefix, sampleID_pass1. (bam, will be added by the STAR program)
 --runThreadN 							                        for computing purpose allocate the number of threads, 14 

 SBATCH script - 
sbatch NameOfSBATCHscript.sh 

--------------------------------------
 9. Check initial quality stats on .bam files 
--------------------------------------
 Will print basic statistics from input BAM file(s)
 Total reads:       
 Mapped reads:      
 Forward strand:    
 Reverse strand:    
 Failed QC:         
 Duplicates:        
 Paired-end reads:  
 'Proper-pairs':    
 Both pairs mapped: 
 Read 1:            
 Read 2:            
 Singletons:   
     
 Example command: $bamtools stats -in sampleID_pass2.Aligned.out.bam > sampleID_pass2.txt
 bamtools											            package 
 stats												            command to get general alignment statistics
 -in												            indicates input file
 sampleID.bam										            path and name to bam file 
 >												                directs output     
 sampleID_pass2.txt                                            indicated output file name

 SBATCH script - 
sbatch NameOfSBATCHscript.sh 

--------------------------------------
 10. Sort bam files
--------------------------------------
 An appropriate @HD-SO sort order header tag will be added or an existing one updated if necessary.
 Example command: $bamtools sort -in sampleID.bam -out sampleID.sorted.bam
 bamtools								                        package 
 sort								                            command to add header tags
 -in									                        indicates input file
 sampleID.bam						                            path and name to bam file 
 -out 									                        indicated output file
 sampleID.sorted.bam					                        the sorted BAM file 

 SBATCH script - 
sbatch NameOfSBATCHscript.sh 

--------------------------------------
 11. Generate stats on sorted bam files
--------------------------------------
 Will print basic statistics from input BAM file(s)   
 Example command: $bamtools stats -in sampleID_pass2.Aligned.sorted.bam > sampleID_pass2.Aligned.sorted.txt
 bamtools											            package 
 stats												            command to get general alignment statistics
 -in												            indicates input file
 sampleID.bam										            path and name to bam file 
 >												                directs output     
 sampleID_pass2.txt                                            indicated output file name

 SBATCH script - 
sbatch NameOfSBATCHscript.sh 

 Compare sorted.bam stats to the original .bam stats, there should be no differences between them
 We do this to step (bamtools stats) every time we do anything to our bam files as a quality control check
--------------------------------------
 12. Mark duplicates & Remove duplicates
--------------------------------------
 Mark duplicates: "Flags" where the duplicate reads are
 Example command: $java -Xmx8g -jar picard.jar MarkDuplicates INPUT=sampleID.sorted.bam OUTPUT=sampleID.sorted.markdup.bam METRICS_FILE=sampleID.markdup.picardMetrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
 java												            program called 
 -Xmx8g 											            declares memory 
 picard.jar 										            path to picard jar file 
 MarkDuplicates 									            command to create sequence dictionary 
 INPUT=											            path to input file, sorted bam file per sample	
 OUTPUT= 											            path and name or output file .markdup to indicate this file will contain duplicates that have been marked
 METRICS_FILE=										            file to write duplication metrics to save as sampleID.markdup.picardMetrics.txt
 REMOVE_DUPLICATES=false 							            If true do not write duplicates to the output file instead of writing them with appropriate flags set. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
 ASSUME_SORTED=true 								            BAM files are sorted because we sorted them in step 6
 VALIDATION_STRINGENCY=LENIENT						            setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.

 SBATCH script - 
sbatch NameOfSBATCHscript.sh 

--------------------------------------
 13. Generate statistics for bam files with duplicates marked 
--------------------------------------
 For each sample get the read stats for the mark duplicates BAM files 
 Example command: $bamtools stats -in sampleID_pass2.Aligned.sorted.dupsMark.bam > sampleID_pass2.Aligned.sorted.dupsMark.txt

 SBATCH script - 
sbatch NameOfSBATCHscript.sh 

 Compare stat results for each sample to the original bam file (sanity check: is there the same number of reads as the original BAM file?)
 If there is more than 15% of the reads being marked as duplicates may need to consider removing that sample
--------------------------------------
 14. Add or replace read groups
--------------------------------------
 For each sample, add a read group to the mark duplicate BAM files (a read group is a "tag" such as a sample ID)
 Example command: $java -Xmx8g -jar picard.jar AddOrReplaceReadGroups INPUT=sampleID.sorted.rmdup.bam OUTPUT=sampleID.sorted.rmdup.addReadGr.bam RGLB=sampleID RGPL=machineUsed RGPU=laneUsed RGSM=sampleName RGCN=location RGDS=species VALIDATION_STRINGENCY=LENIENT
 java												            program called 
 -Xmx8g 											            declares memory 
 picard.jar 										            path to picard jar file 
 AddOrReplaceReadGroups							            Replaces all read groups in the INPUT file with a single new read group and assigns all reads to this read group in the OUTPUT BAM
 INPUT=											            path to input file, sorted bam file per sample	
 OUTPUT= 											            path and name or output file .markdup to indicate this file will contain duplicates that have been marked
 RGLB=												            Read Group Library Required (sampleID)
 RGPL=												            Read Group platform (e.g. illumina, solid) Required 
 RGPU=	 											            Read Group platform unit (eg. run barcode) Required	(laneUsed)
 RGSM=	 											            Read Group sample name Required (sampleName or sampleID)
 RGCN= 											            Read Group sequencing center name Default value: null (i.e. ASU)
 RGDS=												            Read Group description Default value: null (speciesName)
 VALIDATION_STRINGENCY=LENIENT						            setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.

 SBATCH script - 
sbatch NameOfSBATCHscript.sh 

--------------------------------------
 15. Generate statistics on read group BAM sample files
--------------------------------------
 For each sample get the read stats for the remove duplicates and add read groups BAM files.
 Statistics on the BAM files should be the same as before the previous step when read groups were modified.
 Compare stat results for each sample to the markdup.bam file (sanity check: is there the same number of reads as the original BAM file?)
 Example command: bamtools stats -in sampleID.sorted.markdup.addReadGr.bam

 SBATCH script - 
sbatch NameOfSBATCHscript.sh 

--------------------------------------
 16. Index BAM files 
--------------------------------------
 For each sample index the processed BAM files that are sorted, have marked duplicates, and have read groups added. These will be used to identify callable loci. 
 Indexing is used to "sort" by chromosome and region 
 Output will be sampleID.sorted.markdup.addReadGr.bam.bai
 Example command: $bamtools index -in sampleID.sorted.markdup.addReadGr.bam
 bamtools											            package 
 index												            Generates index for BAM file
 -in												            indicates input file
 sampleID.sorted.markdup.addReadGrbam				            path and name to bam file that has been sorted, duplicate reads were removed, and the read groups were added 

 SBATCH script - 
sbatch NameOfSBATCHscript.sh 

---------------------------------------
 18. Differential expression using Cuffdiff
---------------------------------------
 Identify differentially expressed genes and transcripts using cuffdiff. 
 Cuffdiff is used to find significant changes in transcript expression, splicing, and promoter use. 
 Cuffdiff uses a gene annotation file downloaded with a sbatch script from UCSC hg19
 Example command: $cuffdiff -use-sample-sheet -o diff_out -b reference.fa -p 8 --library-type fr-firststrand -L set1,set2 -u refence.gtf sampleSet.txt
 cuffdiff 											            program in cufflinks package that identifies if genes are differentially expressed 
 -use-sample-sheet 								            tells the program to use a text file containing the path to the samples and the group_id 
 -o 												            output
 diff_out 											            path and name of output directory where all of the output files from the cuffdiff program will be placed
 -b 												            Providing Cufflinks with the multifasta file your reads were mapped to via this option instructs it to run the bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates
 reference.fa 										            path and name of reference.fasta which the files were mapped to
 -p 												            Use this many threads to align reads. The default is 1
 8 												            declared 8 threads 
 --library-type												* This option depends on the data * unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute
 fr-firststrand												* This option depends on the data * should be used for the standard dUTP protocol, including Illumina’s stranded Tru-Seq.
 -L 												            Specify a label for each sample, which will be included in various output files produced by Cuffdiff.
 set1,set2 										            Must have at least 2 labels. (i.e male,female)	
 -u 												            Tells Cufflinks to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome
 refence.gtf 										            path and name to the gene annotation file. 
 sampleSet.txt										            text file containing the list of samples

 SBATCH script - 
sbatch NameOfSBATCHscript.sh 

