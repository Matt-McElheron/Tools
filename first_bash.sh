#!/bin/bash

# DATE: 11-07-19
# TASK: beginning my first bash
# AUTHOR: Matt mcElheron





# Here we make a for loop to iterate over the input
# $@ takes the list of arguments inputted when calling the script
for file in $@; do
	echo $file


	# File to be analysed is in .gz format
	# .gz means the file has been compressed to reduce its size, and must be decompressed using gunzip to be used.
	# -k allows retention of input files. 
	# -c write on standard output, keep original files unchanged, removing .gz and keeping the .fastq
	# gunzip -c "$input_file" > "${file/.fastq*/.fastq}"
	#gunzip -k $file


	# Total number of reads in each file is then calculated
	# A line count is performed to find the total number of lines, which is divided by 4 to find the total number of reads in a file
	# This number can then be doubled, as pair-end sequencing provided both forward and reverse reads, giving the total number of reads.
	# -l this option allows for counting of lines within a file, with other options allowing word/character count etc.
#	wc -l $file
	

	# The fastQC program is then used to assess the quality of the reads in each file
	# This produces a .HTML file allowing visualization of different aspects of quality
	# This includes basic information such as sequence length, GC content and sequences flagged as poor quality
	# The report indicates any need for adapter trimming, duplicate removal, potential contamination removal etc.
	# The report can be used to help explain any potential downstream problems.
#	fastqc $file


# Here we end the for loop
#done

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
	#if [ -f trimmed_"${file//_*}"_R1.fastq ]; then
	#	continue
	#else:
	#cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -b CTGTCTCTTATA -B CTGTCTCTTATA -o trimmed_"${file//_*}"_R1.fastq -p trimmed_"${file//_*}"_R2.fastq "${file//_*}"_R1.fastq "${file//_*}"_R2.fastq

done
