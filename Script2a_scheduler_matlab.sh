#!/bin/bash

#SBATCH --ntasks=
#SBATCH --time=
#SBATCH --cpus-per-task=
#SBATCH --partition=
#SBATCH --output=
#SBATCH --error=
#SBATCH --mem=
#SBATCH --job-name=
#SBATCH --mail-type=
#SBATCH --mail-user=


### Author: Manon de Visser | Date: 2024-03-28 | Contact: devisser.manon@gmail.com

### Load any modules required, at least MATLAB:

module load [module]
module load MATLAB

### To run PAV-spotter, customize to your need the following variables in the command below (and see the manuscript and GitHub page for more details on what these variables entail):

# current_cd          --> Make sure your working directory corresponds with where you have locally downloaded Script 2 / the PAV-spotter repository
# save_file_name      --> Specify the desired name of your output file
# contig_width        --> Specify in the same format as below a number for minimum 'width' (could also be 0 in case irrelevant)
# reads_threshold     --> Specify in the same format as below a number for minimum 'depth' (could also be 0 in case irrelevant)
# categories          --> Specify in the same format as below the two, or three, identifiers for each phenotypic class by strings
# ctrl_category       --> Specify in the same format as below which out of the provided class identifiers is considered the control group
# common_identifier   --> Specicy in the same format as below a string that is a common identifier for all input files/targets
# plot_figures        --> Specify if you would like to output figures for each target by setting this argument to either 'true' or 'false'
# cd [PATH]           --> Adjust this path so that it corresponds with where you have locally downloaded Script2 / the PAV-spotter repository


matlab -nodisplay -r "current_cd = cd; save_file_name='OUTFILE'; contig_width = 0; reads_threshold = 0; categories = {'ST','FT','HC'}; ctrl_category = 'HC'; common_identifier = 'DN'; sum_for_each_subspecies=true; plot_figures = true; cd ~/path/to/PAV-spotter/; Script2_PAVSpotter.m ; exit"
