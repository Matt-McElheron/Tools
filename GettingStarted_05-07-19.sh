#!/bin/bash

# DATE: 05-07-19
# TASK: Alignment of sequenced M. bovis isolates' genomes to a reference genome with additional analysis
# AUTHOR: Matt mcElheron

# Files to be analysed were in .gz format
# .gz means the file has been compressed to reduce its size, and must be decompressed using gunzip to be used.
# -k allows retention of input files. 
gunzip *.gz

# Files are now stored in FASTQ format (.fastq)
# Fastq files contain the read information produced from sequencing
# Each read is represented by four lines:
# 1: A header line, beginning with an "@" character, followed by a read ID and other optional information
# 2: The read itsef, made up of characters representing the bases read.
# 3: An optional comment line, beginning with a '+' character
# 4: A line encoding the quality scores of each base read in line 2.
# Line 4 has to be the same length of line 2, meaning quality of corrosponding base is represented using (ASCII) 0x21 to 0x7e
# Hence, from lowest to highest quality, quality scores look like this: !"#$%&'()*+,-./0123456789:;<=>?



# Total number of reads in each file was then calculated
# A line count was performed to find the total number of lines, which was divided by 4 to find the total number of reads in a file
# This number can then be doubled, as pair-end sequencing provided both forward and reverse reads, giving the total number of reads.
# -l this option allows for counting of lines within a file, with other options allowing word/character count etc.
wc -l file.fastq

# The total number or reads can then be multiplied by the known read length (150) and divided by the expected genome size
# This gives a good rough estimate of the average coverage
# It wont be exact as some reads will be less than 150, and the genome of the isolate may differ in size.


# The fastQC program was then used to assess the quality of the reads in each file
# This produces a .HTML file allowing visualization of different aspects of quality
# This includes basic information such as sequence length, GC content and sequences flagged as poor quality
# The report indicates any need for adapter trimming, duplicate removal, potential contamination removal etc.
# The report can be used to help explain any potential downstream problems.
fastqc *.fastqc

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
bwa mem -t threadINT reference_genome.fna forward_reads reverse_reads -o alignment.sam 

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
# -b 		outputs a BAM file
# -o 		output file
samtools view -f 4 -b -o unmapped.bam align.sam






# Examine sorted reads to create BCF file
# BCF file is binary format of VCF file (Variant Calling Format)
# -f description
bcftools mpileup -f ~/Desktop/WICKLOW/Reference_Mbovis/GCF_000195835.2_ASM19583v2_genomic.fna align_23.sorted.bam -o 

# Calling variants from BCF file
# -m
# -v
# -O
# -o
bcftools call -mv -Ov -o align_23.vcf