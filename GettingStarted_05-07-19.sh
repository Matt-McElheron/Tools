#!/bin/bash

# DATE: 05-07-19
# TASK: Alignment of sequenced M. bovis isolates' genomes to a reference genome with additional analysis
# AUTHOR: Matt mcElheron

# Files to be analysed were in .gz format
# .gz means the file has been compressed to reduce its size, and must be decompressed using gunzip to be used.
# -k allows retention of input files. 
# $1 means it parses the first commandline argument
gunzip -k *.gz

# Files are now stored in FASTQ format (.fastq)
# Fastq files contain the read information produced from sequencing
# Each read is represented by four lines:
# 1: A header line, beginning with an "@" character, followed by a read ID and other optional information
# 2: The read itsef, made up of characters representing the bases read.
# 3: An optional comment line, beginning with a '+' character
# 4: A line encoding the quality scores of each base read in line 2.
# Line 4 has to be the same length of line 2, meaning quality of corrosponding base is represented using (ASCII) 0x21 to 0x7e
# Hence, from lowest to highest quality, quality scores look like this: !"#$%&'()*+,-./0123456789:;<=>?



# Total number of reads in each file is then calculated
# A line count is performed to find the total number of lines, which is divided by 4 to find the total number of reads in a file
# This number can then be doubled, as pair-end sequencing provided both forward and reverse reads, giving the total number of reads.
# -l this option allows for counting of lines within a file, with other options allowing word/character count etc.
wc -l file.fastq

# The total number or reads can then be multiplied by the known read length (150) and divided by the expected genome size
# This gives a good rough estimate of the average coverage
# It wont be exact as some reads will be less than 150, and the genome of the isolate may differ in size.


# The fastQC program is then used to assess the quality of the reads in each file
# This produces a .HTML file allowing visualization of different aspects of quality
# This includes basic information such as sequence length, GC content and sequences flagged as poor quality
# The report indicates any need for adapter trimming, duplicate removal, potential contamination removal etc.
# The report can be used to help explain any potential downstream problems.
fastqc *.fastq

# While often uneeded due to modern sequence technology, a trim step is often used for quality assurance
# Adapters are short sequences ligated during PCR-based sequencing methods, not originating from the samples DNA.
# The cutadapt program is used to remove adapter sequences from high-throughput sequencing reads
# The program requires adapters to be removed, files in and file output name to be specified
# -b AGATCGGAAGAG	specifying removal this primer, which may be ligated to 5' or 3' end of the first read 
# -B AGATCGGAAGAG	5'/3 adapter to be removed from second read in a pair
# -b CTGTCTCTTATA	specifying removal this primer, which may be ligated to 5' or 3' end of the first read
# -B CTGTCTCTTATA	5'/3 adapter to be removed from second read in a pair
# -o 				output file for first read after trimming
# -p 				output file for second read after trimming
cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -b CTGTCTCTTATA -B CTGTCTCTTATA -o out1 -p out2 forward.fastq reverse.fastq

# Before aligning reads to a reference genome, the genome itself needs to be indexed using bwa index
# This constructs an FM-index - a compressed full-text substring index
# This allows the frequency and locations of a pattern within the compressed text to be found
# The indexing is often perform within a new directory, as many new files are created when indexing.
# Programs requiring an indexed genome need only the original genome as input, and know to search the dir for the accompanying files
# -p STR is optional and adds STR as a prefix for output database
bwa index reference_genome.fna

# The trimmed read files can then be aligned to the reference genome, using bwa mem
# mem aligns 70bp-1Mbp query sequences using the BWA-MEM algorithm
# The algorithm seeds alignments with "maximal exact matches" (mem), and extends seeds using the affine-gap SW algorithm
# Seeding in bioinformatics is essentially finding a match between a query sequence and a hit sequence
# bwa mem can be used on paired end or mate-end reads, depending on how may read files are inputted
# bwa mem is generally industry standard and more suited for illumina reads than bwa-backtrack or bwa-aln
# This outputs a .sam file
# -t specifies number of threads used
# -o specifies output file
bwa mem -t threadINT reference_genome.fna forward_reads.fastq reverse_reads.fastq -o alignment.sam 

# Sequence Alignment Map (SAM) is a text-base format
# The format contains a header section beginning with an '@' character and relevent information
# The alignment section contains a huge amount of information (with further optional information fields also)
# Most notable are the sequence itself, mapping quality, and single character base quality
# Different pieces of information are stored under flags, allowing combinations of flags to be called
# SAM files can be stored in binary as BAM files




# SAMtools is used to analyse SAM files
# Samtools view is used to view maps within the file, with parameters used to specify target sequences
# -f INT	is used to include reads with all of the flags in INT present
# -f 4		is used to view all unmapped reads
# -o 		output file
samtools view -f 4 -o unmapped.sam alignment.sam

# The number unmapped reads is then counted by extracted lines containing the correct header found once per read line counting
# grep retrieves all lines with the quoted term
# the '|' character sends the product into the next command, in this case a line count
# "NB501589" is here an example and specific to the data
grep "NB501589:" unmapped.sam | wc -l

# The above is repeated to count total number of maps using the alignment file
grep "NB501589:" alignment.sam | wc -l

#Proportion or percentage mapped/unmapped reads is then found through simple math


# The original SAM file is converted to a BAM file
# Binary alignment files are compressed versions of SAM files, but still allow analysis using SAMtools.
# This is performed using samtools view
# -b 	specifies the output file to be a BAM file
# -o 	output file
# -S 	ignored (input format is automatically detected)
samtools view -S -b -o alignment.bam alignment.sam

# This BAM file is then sorted using samtools sort
# BAM files can be sorted by different variables
# -n 		sort by read name
# -t TAG	sort by value of TAG. Uses position as secondary index or read name if -n is used
# -o 		output file
# -O BAM 	output format, (SAM, BAM, CRAM), tries to match suffix if not specified
samtools sort -O BAM -o alignment.sorted.bam alignment.bam

# The sorted BAM files then had statistics viewed using samtools flagstats
# This gives general statistics, such as reads mapped, and also duplicate number
# Duplicates may need to be removed in low quality sequencing data
# This tells us number of unmapped reads without the need for extracting using "view -f 4"
# This command usually prints stats to the screen, and does not have a specify output file option. so '>' is used
samtools flagstat alignment.sorted.bam > alignment.sorted.bam.stats

# Once sorted, the sorted BAM files must be indexed
# samtools index allows a coordinate-sorted BAM (or CRAM) file to be indexed for fast random access
# This index is needed when region arguments are used to limit samtools view and similar commands to regions of interest
samtools index alignment.sorted.bam alignment.index


# Duplicates are removed using samtools rmdup
# In NGS, PCR is used. Genomes are broken up and amplified, before being sequenced.
# However, sometimes there is a bias in which fragments are amplified the most. 
# Hence, potential duplicates can be highlighted and removed.
# Only the duplicate with the highest mapping quality will be kept
# Does not work for unpaired reads
# -s removes duplicates for single-end reads. By default, the command works for paired-end only
# -S treat paired-end and single-end reads
samtools rmdup -sS alignment.sorted.bam alignment.rmdups.sorted.bam


# Next, variant call analysis is performed. This uses a reference genome to inspect SNPs, INDELs and other structural variations
# The Variant Call Format (VCF) is used to store variant call information. BCF is the binary version of a VCF.
# It contains meta-info lines, a header line, and positional information.
# A BCF file is made using the samtools mpileup command, which is then be used to create a VCF file.
# -o	specifies output
# -u 	output is uncompressed, preferrable for piping
# -g 	specifies .bcf output
# -f 	reference genome
samtools mpileup -g -f ~/Desktop/WICKLOW/Reference_Mbovis/GCF_000195835.2_ASM19583v2_genomic.fna -o alignment.bcf alignment.rmdups.sorted.bam

# The BCF file is then used to call variants using BCFtools
# This produces a vcf file which includes information about variations found within the genome, such as SNPs
# This is done using the call command
# -o 		specifies output file
# -Ov 		specifies uncompressed VCF file type
# -Oz 		specifies compressed VCF file type
# -Ou 		specifies uncompressed BCF file type
# -Ob 		specifies compressed BCF file type
# -c 		tells command to use consensus caller
# -mv 		tells command to use multiallelic calling
# --ploidy	gives predefined ploidy
bcftools call -Ov --ploidy 1 -o alignment.vcf -mv alignment.bcf


