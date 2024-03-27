#!/bin/bash
#SBATCH --ntasks=
#SBATCH --time=
#SBATCH --cpus-per-task=
#SBATCH --partition=
#SBATCH --output=
#SBATCH --error=
#SBATCH --job-name=
#SBATCH --mem=
#SBATCH --mail-type=
#SBATCH --mail-user=

---

###### CONTEXT

###### This (SLURM) script - Script1 - is written by Manon de Visser and serves as an example 
###### It is published with the Manuscript and it performs upstream filename + folder structure arrangements

###### With each command comes an explanation, as to facilitate any user to customize the script for ones own need
###### This script creates depth files from (deduplicated) BAM files and uses them to create an input file for the cross-correlation/PAV script (Script2)

---

### NOTES ON THE INPUT FILENAMES AND FOLDER STRUCTURE REQUIRED

### The working directory should have one or multiple subdirectories, which in turn have the BAM files of samples of two (or three) different phenotype/genotype classes (in the same folder)
### Depth files will be created, a sample ID column will be added, the files will be merged into one per spp/subdirectory, ...
### ... the format will be changed into CSV, the output will be split up based on the target names, and lastly some cleanup steps will get rid of intermediate files

### !!! Make sure that the filenames of the input BAMs have a sample identifier, as well as a phenotype/genotype class identifier embedded

---

# STEPS OF RUNNING THIS SCRIPT

# Load the modules you need for running on your HPC / in your Linux environment
# Load SAMtools / make sure you have this software installed

module load [MODULE]
module load SAMtools


# Provide your working directory

cd /path/to/working/directory


# Start up a for loop that will go through your subdirectory/subdirectories with the BAM files

for directory in */
do


# General checks

	echo "Starting loop for species/analysis '$directory'"
	cd $directory
	echo "Contents of working directory are as follows:"
  ls -l


# Run samtools' depth function using the 'all' argument to output read counts / depths for all the positions per target (including zeros)

	for file in *bam; do samtools depth -a $file >> ${file%.bam}.depth; done


# Add a column with the sample/file name, so that later on it is known which line belongs to which sample + phenotype/genotype class

	for file in *.depth; do awk 'NR == 1 {print $0 "\t" FILENAME; next;}{print $0 "\t" FILENAME;}' $file >> ${file%.depth}.samplename.depth; done


# Take all the samplename.depth files and merge them into one overall file. 
# [Be aware: depending on which (parts of) filenames you decide to use, figures created by PAV-spotter/Script 2 may or may not end up having a .jpg extension. In case they do not, you can simply add the extension afterwards to make the file visible, for instance by using this oneliner: 'for f in file; do mv "$f" "$f.jpg"; done' ]

	cat *.samplename.depth >> allgenes_allsamples.txt


# Sort the file by target (so the information will be available per target consecutively)

	sort -s -k 1,1 allgenes_allsamples.txt >> allgenes_allsamples_sorted.txt


# Change the format to CSV (replace tabs with commas)

  sed 's/\t/,/g' allgenes_allsamples_sorted.txt >> allgenes_allsamples_sorted.csv


# Create a new directory, copy the final/overall CSV to it, and then split it up in multiple files (one for each separate target, which will thus include the data for that target of all individual samples of a certain 'species')

	mkdir allgenes_split
	scp allgenes_allsamples_sorted.csv ./allgenes_split
	cd allgenes_split
	awk -F\, '{print>$1}' allgenes_allsamples_sorted.csv


# Create a file with the sample names of the input individuals (needed for downstream analysis)

	cat allgenes_allsamples_sorted.csv | awk -F\, '{print $4}' | sort | uniq >> individuals.txt


# Clean up (but make sure you don't lose critical data - to prevent that from happening, the command that deletes the BAM files is silenced with a '#' for now)

	rm allgenes_allsamples_sorted.csv
	cd ../
	#rm *dedup.*bam
	gzip *depth*
	gzip *.txt
	gzip *.csv


# direct two directories back to make sure the for loop can continue with the next 'species' 

	cd ../
 
  echo "Finished loop for species/analysis '$directory'. Moving on to (potential) next species/analysis loop"

done

echo "No next subdirectory detected, Script1 finished running. Thus, input files for PAV-spotter should have been created. Please check whether those files are in place and correct. These should at least be; 1) separate loose files per each target, with information of all samples/individuals for that particular target, in a subfolder called 'allgenes_split', and 2) the 'individuals.txt' file, which should be created in that same subfolder and which must contain the sample names with their phenotype/genotype 'class identifiers'. If all seems in order, you can proceed to PAV-spotter/Script2. Good luck!"
