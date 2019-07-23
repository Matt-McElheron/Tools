#!/bin/bash

# DATE: 11-07-19
# TASK: beginning my first bash
# AUTHOR: Matt McElheron

# Command line structure
# REFERENCE_GENOME ThreadINT blastamountINT master_contam_file.csv *.fastq.gz is general input

# Here we take in the command line arguments, excluding the fastq files:
reference_genome=$1
thread_number=$2
num_blasts=$3
contam_file=$4
contam_file_summary="${contam_file//.csv}""_summary.csv"
contam_file_species="${contam_file//.csv}""_species.csv"


# These are universal adapter sequences removed from the reads
adap1=AGATCGGAAGAG
adap2=CTGTCTCTTATA

# Before aligning reads to a reference genome, the genome itself needs to be indexed using bwa index
# This constructs an FM-index - a compressed full-text substring index
# This allows the frequency and locations of a pattern within the compressed text to be found
# The indexing is often perform within a new directory, as many new files are created when indexing.
# Programs requiring an indexed genome need only the original genome as input, and know to search the dir for the accompanying files
# -p STR is optional and adds STR as a prefix for output database
# We skip this step if the referene genome has already been indexed, by checking for index products.
if [ -f $reference_genome".fai" ]; then
	echo "Reference genome $reference_genome already indexed, skipping indexing step."
else
	bwa index $reference_genome
fi

# Here we make a for loop to iterate over the input
# ${@:5} takes the list of arguments inputted when calling the script, skipping the initial inputs which aren't reads
for fastq_file in ${@:5}; do

	# File to be analysed is in .gz format
	# .gz means the file has been compressed to reduce its size, and must be decompressed using gunzip to be used.
	# -k allows retention of input files. 
	# -c write on standard output, keep original files unchanged, removing .gz and keeping the .fastq
	# gunzip -c "$input_file" > "${file/.fastq*/.fastq}"
	if [ -f "${fastq_file//-*}".vcf ]; then
		continue
	fi

	# Here we unzip the reverse read with the forward, so they can be used together downstream
	
	gunzip $fastq_file
	gunzip "${fastq_file//R1*}"R2"${fastq_file//*R1}"

	# This removes the characters stated, here ".gz"
	# This variable is set so that the reverse read can be more easily called using ${unzipped_forward//*R1}" in one go
	# Otherwise, ${fastq_file//*R1}" would still have .gz suffix
	unzipped_forward="${fastq_file//.gz}"
	unzipped_reverse="${unzipped_forward//R1*}"R2"${unzipped_forward//*R1}"

	# This prefix is going to be used to name files after both reads are combined into sam files, bam files, vcf files etc..
	prefix="${unzipped_forward//-*}"
	echo $prefix

	# The fastQC program is then used to assess the quality of the reads in each file
	# This produces a .HTML file allowing visualization of different aspects of quality
	# This includes basic information such as sequence length, GC content and sequences flagged as poor quality
	# The report indicates any need for adapter trimming, duplicate removal, potential contamination removal etc.
	# The report can be used to help explain any potential downstream problems.
	fastqc $unzipped_forward $unzipped_reverse


	# Script then check to see if the argument is a forward read or a reverse read
	# These commands take the forward read, grab the associated reverse read
	# Hence, when looking at the reverse read argument, this if checkpoint prevents double commands
	# The statment checks to see if the output has already been created
	if [ -f $prefix.vcf ]; then
		continue
	fi

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
	cutadapt -b $adap1 -B $adap1 -b $adap2 -B $adap2 -o trimmed_$unzipped_forward -p trimmed_$unzipped_reverse $unzipped_forward $unzipped_reverse


	# The trimmed read files can then be aligned to the reference genome, using bwa mem
	# mem aligns 70bp-1Mbp query sequences using the BWA-MEM algorithm
	# The algorithm seeds alignments with "maximal exact matches" (mem), and extends seeds using the affine-gap SW algorithm
	# Seeding in bioinformatics is essentially finding a match between a query sequence and a hit sequence
	# bwa mem can be used on paired end or mate-end reads, depending on how may read files are inputted
	# bwa mem is generally industry standard and more suited for illumina reads than bwa-backtrack or bwa-aln
	# This outputs a .sam file
	# -t specifies number of threads used
	# -o specifies output file
	bwa mem -t $thread_number $reference_genome trimmed_$unzipped_forward trimmed_$unzipped_reverse -o $prefix.sam 

	#  The trimmed files can be removed to save space
	rm trimmed_$unzipped_forward trimmed_$unzipped_reverse

	# The fastq files can be compressed to save space
	gzip $unzipped_forward $unzipped_reverse


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
	samtools view -f 4 -o $prefix"_unmapped.sam" $prefix.sam


	# PERFORM UNMAPPED ANALYSIS SCRIPT HERE??
	# Here a python script is called on each unmapped read file
	# The file has an amount of randomly selected reads Blasted to determine their source
	# requires unmapped read file (.sam), amount to be blasted, and the minimum read length to be blasted, and output file as parameters
	python3 ~/Desktop/Matt_py/pandas_blast.py $prefix"_unmapped.sam" $num_blasts 60 $contam_file



	# The original SAM file is converted to a BAM file
	# Binary alignment files are compressed versions of SAM files, but still allow analysis using SAMtools.
	# This is performed using samtools view
	# -b 	specifies the output file to be a BAM file
	# -o 	output file
	# -S 	ignored (input format is automatically detected)
	samtools view -S -b -o $prefix.bam $prefix.sam


	# This BAM file is then sorted using samtools sort
	# BAM files can be sorted by different variables
	# -n 		sort by read name
	# -t TAG	sort by value of TAG. Uses position as secondary index or read name if -n is used
	# -o 		output file
	# -O BAM 	output format, (SAM, BAM, CRAM), tries to match suffix if not specified
	samtools sort -O BAM -o $prefix.sorted.bam $prefix.bam


	# The sorted BAM files then had statistics viewed using samtools flagstats
	# This gives general statistics, such as reads mapped, and also duplicate number
	# Duplicates may need to be removed in low quality sequencing data
	# This tells us number of unmapped reads without the need for extracting using "view -f 4"
	# This command usually prints stats to the screen, and does not have a specify output file option. so '>' is used
	samtools flagstat $prefix.sorted.bam > $prefix.sorted.bam.stats


	# Once sorted, the sorted BAM files must be indexed
	# samtools index allows a coordinate-sorted BAM (or CRAM) file to be indexed for fast random access
	# This index is needed when region arguments are used to limit samtools view and similar commands to regions of interest
	samtools index $prefix.sorted.bam $prefix.index


	# Duplicates are removed using samtools rmdup
	# In NGS, PCR is used. Genomes are broken up and amplified, before being sequenced.
	# However, sometimes there is a bias in which fragments are amplified the most. 
	# Hence, potential duplicates can be highlighted and removed.
	# Only the duplicate with the highest mapping quality will be kept
	# Does not work for unpaired reads
	# -s removes duplicates for single-end reads. By default, the command works for paired-end only
	# -S treat paired-end and single-end reads
	samtools rmdup -sS $prefix.sorted.bam $prefix.rmdups.sorted.bam


	# Next, variant call analysis is performed. This uses a reference genome to inspect SNPs, INDELs and other structural variations
	# The Variant Call Format (VCF) is used to store variant call information. BCF is the binary version of a VCF.
	# It contains meta-info lines, a header line, and positional information.
	# A BCF file is made using the samtools mpileup command, which is then be used to create a VCF file.
	# -o	specifies output
	# -u 	output is uncompressed, preferrable for piping
	# -g 	specifies .bcf output
	# -f 	reference genome
	samtools mpileup -g -f $reference_genome -o $prefix.bcf $prefix.rmdups.sorted.bam


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
	bcftools call -Ov --ploidy 1 -o $prefix.vcf -mv $prefix.bcf

	# The accesory files can here be removed for storage reasons
	rm $prefix.sam $prefix.sorted.bam $prefix.bam $prefix.index $prefix.rmdups.sorted.bam $prefix.bcf

done


# The script for combining unmapped read data is then called
python3 ~/Desktop/Matt_py/Master_unmapped.py $contam_file $contam_file_summary

# This script finds out which contam is which
python3 ~/Desktop/Matt_py/taxid.py $contam_file_summary $contam_file_species
