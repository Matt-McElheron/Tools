#!/bin/bash

# DATE: 11-07-19
# TASK: beginning my first bash
# AUTHOR: Matt McElheron

# Command line structure
# REFERENCE_GENOME ThreadINT *.fastq.gz is general input



# Before aligning reads to a reference genome, the genome itself needs to be indexed using bwa index
# This constructs an FM-index - a compressed full-text substring index
# This allows the frequency and locations of a pattern within the compressed text to be found
# The indexing is often perform within a new directory, as many new files are created when indexing.
# Programs requiring an indexed genome need only the original genome as input, and know to search the dir for the accompanying files
# -p STR is optional and adds STR as a prefix for output database
bwa index $1


# Here we make a for loop to iterate over the input
# ${@:3} takes the list of arguments inputted when calling the script, skipping the initial inputs which aren't reads
for fastq_file in ${@:3}; do


	# File to be analysed is in .gz format
	# .gz means the file has been compressed to reduce its size, and must be decompressed using gunzip to be used.
	# -k allows retention of input files. 
	# -c write on standard output, keep original files unchanged, removing .gz and keeping the .fastq
	# gunzip -c "$input_file" > "${file/.fastq*/.fastq}"
	if [ -f ${fastq_file//.gz} ]; then
		continue
	else


		# Here we unzip the reverse read with the forward, so they can be used together downstream
		gunzip $fastq_file
		gunzip "${fastq_file//R1*}"R2"${fastq_file//*R1}"


		
	fi


	# This removes the characters stated, here ".gz"
	# This variable is set so that the reverse read can be more easily called using ${unzipped_file//*R1}" in one go
	# Otherwise, ${fastq_file//*R1}" would still have .gz suffix
	unzipped_file="${fastq_file//.gz}"


	# The fastQC program is then used to assess the quality of the reads in each file
	# This produces a .HTML file allowing visualization of different aspects of quality
	# This includes basic information such as sequence length, GC content and sequences flagged as poor quality
	# The report indicates any need for adapter trimming, duplicate removal, potential contamination removal etc.
	# The report can be used to help explain any potential downstream problems.
	fastqc $unzipped_file


	# Script then check to see if the argument is a forward read or a reverse read
	# These commands take the forward read, grab the associated reverse read
	# Hence, when looking at the reverse read argument, this if checkpoint prevents double commands
	# The statment checks to see if the output has already been created
	if [ -f trimmed_$unzipped_file ]; then
		continue
	else


	# While often uneeded due to modern sequence technology, a trim step is often used for quality assurance
	# Adapters are short sequences ligated during PCR-based sequencing methods, not originating from the samples DNA.
	# The cutadapt program is used to remove adapter sequences from high-throughput sequencing reads
	# The program requires adapters to be removed, files in and file output name to be specified
	# The command has to take in the forward and reverse reads together and not repeat the command for already trimmed files
	# -b AGATCGGAAGAG	specifying removal this primer, which may be ligated to 5' or 3' end of the first read 
	# -B AGATCGGAAGAG	5'/3 adapter to be removed from second read in a pair
	# -b CTGTCTCTTATA	specifying removal this primer, which may be ligated to 5' or 3' end of the first read
	# -B CTGTCTCTTATA	5'/3 adapter to be removed from second read in a pair
	# -o 				output file for first read after trimming
	# -p 				output file for second read after trimming
		cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -b CTGTCTCTTATA -B CTGTCTCTTATA -o trimmed_$unzipped_file -p trimmed_"${unzipped_file//R1*}"R2"${unzipped_file//*R1}" $unzipped_file "${unzipped_file//R1*}"R2"${unzipped_file//*R1}"


	# The trimmed read files can then be aligned to the reference genome, using bwa mem
	# mem aligns 70bp-1Mbp query sequences using the BWA-MEM algorithm
	# The algorithm seeds alignments with "maximal exact matches" (mem), and extends seeds using the affine-gap SW algorithm
	# Seeding in bioinformatics is essentially finding a match between a query sequence and a hit sequence
	# bwa mem can be used on paired end or mate-end reads, depending on how may read files are inputted
	# bwa mem is generally industry standard and more suited for illumina reads than bwa-backtrack or bwa-aln
	# This outputs a .sam file
	# -t specifies number of threads used
	# -o specifies output file

		bwa mem -t $2 $1 trimmed_$unzipped_file trimmed_"${unzipped_file//R1*}"R2"${unzipped_file//*R1}" -o "${unzipped_file//-*}".sam 


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
		samtools view -f 4 -o unmapped_"${unzipped_file//-*}".sam "${unzipped_file//-*}".sam


	# PERFORM UNMAPPED ANALYSIS SCRIPT HERE??









	# The original SAM file is converted to a BAM file
	# Binary alignment files are compressed versions of SAM files, but still allow analysis using SAMtools.
	# This is performed using samtools view
	# -b 	specifies the output file to be a BAM file
	# -o 	output file
	# -S 	ignored (input format is automatically detected)
		samtools view -S -b -o "${unzipped_file//-*}".bam "${unzipped_file//-*}".sam


	# This BAM file is then sorted using samtools sort
	# BAM files can be sorted by different variables
	# -n 		sort by read name
	# -t TAG	sort by value of TAG. Uses position as secondary index or read name if -n is used
	# -o 		output file
	# -O BAM 	output format, (SAM, BAM, CRAM), tries to match suffix if not specified
		samtools sort -O BAM -o "${unzipped_file//-*}".sorted.bam "${unzipped_file//-*}".bam


	# The sorted BAM files then had statistics viewed using samtools flagstats
	# This gives general statistics, such as reads mapped, and also duplicate number
	# Duplicates may need to be removed in low quality sequencing data
	# This tells us number of unmapped reads without the need for extracting using "view -f 4"
	# This command usually prints stats to the screen, and does not have a specify output file option. so '>' is used
		samtools flagstat "${unzipped_file//-*}".sorted.bam > "${unzipped_file//-*}".sorted.bam.stats


	# Once sorted, the sorted BAM files must be indexed
	# samtools index allows a coordinate-sorted BAM (or CRAM) file to be indexed for fast random access
	# This index is needed when region arguments are used to limit samtools view and similar commands to regions of interest
		samtools index "${unzipped_file//-*}".sorted.bam "${unzipped_file//-*}".index


	# Duplicates are removed using samtools rmdup
	# In NGS, PCR is used. Genomes are broken up and amplified, before being sequenced.
	# However, sometimes there is a bias in which fragments are amplified the most. 
	# Hence, potential duplicates can be highlighted and removed.
	# Only the duplicate with the highest mapping quality will be kept
	# Does not work for unpaired reads
	# -s removes duplicates for single-end reads. By default, the command works for paired-end only
	# -S treat paired-end and single-end reads
		samtools rmdup -sS "${unzipped_file//-*}".sorted.bam "${unzipped_file//-*}".rmdups.sorted.bam


	# Next, variant call analysis is performed. This uses a reference genome to inspect SNPs, INDELs and other structural variations
	# The Variant Call Format (VCF) is used to store variant call information. BCF is the binary version of a VCF.
	# It contains meta-info lines, a header line, and positional information.
	# A BCF file is made using the samtools mpileup command, which is then be used to create a VCF file.
	# -o	specifies output
	# -u 	output is uncompressed, preferrable for piping
	# -g 	specifies .bcf output
	# -f 	reference genome
		samtools mpileup -g -f $1 -o "${unzipped_file//-*}".bcf "${unzipped_file//-*}".rmdups.sorted.bam


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
		bcftools call -Ov --ploidy 1 -o "${unzipped_file//-*}".vcf -mv "${unzipped_file//-*}".bcf


	fi

############# Final contam compilation etc

done

