#!/bin/bash


---


###### CONTEXT

###### This (bash) script - Script3 - is written by Manon de Visser and serves as an example
###### It is published with the Manuscript and it performs downstream analyses

###### With each command comes an explanation, which makes it easier to customize the script for ones own need
###### The script takes as input a file that comes out of cross-correlation/PAV script - Script2

###### In this current script, the PAV-spotter output is further analyzed/manipulated in order to find signals of Presence Absence Variation...
###### ... using (dis)similarity thresholds that can be changed as desired


---


### NOTES ON THE INPUT FILENAMES AND FOLDER STRUCTURE REQUIRED

### Depending on your situation, you likely compare either two or three phenotype/genotype classes. 
### We work specifically with three here (important to note/repeat)


---

# STEPS OF RUNNING THIS SCRIPT:

# First, sort the output file, which is in csv format, in such a way that you have all the three comparisons lines per target together

sort -t ',' -k 2,2 -k 3,3 OUTFILE.csv > all_data_sorted



# Make sure that any 'NaN' characters, which originated through a zero by zero devision in the cross-correlation script, are replaced with what they should be: a zero (this command is silenced for now, but if you encounter NaN issues you can unsilence it)

#sed -i 's/NaN/0/g' all_data_sorted



# Create three different files that have only the content of one specific comparison of phenotypes/genotypes, so file1 is genotype1-2, file2 is genotype2-3, file3 is genotype 1-3

sed -n '1~3p' all_data_sorted > genotype1-2
sed -n '2~3p' all_data_sorted > genotype3-1
sed -n '3~3p' all_data_sorted > genotype3-2



# Paste those together into one matrix in two steps, so that eventually you will have a file that has one line per target with the cross-correlation values for each comparison in columns. 
# Cleanup after

cat genotype3-1 | cut -d ',' -f3,4,5 | paste -d , genotype1-2 - > intermediate_file
cat genotype3-2 | cut -d ',' -f3,4,5 | paste -d , intermediate_file - > matrix
rm intermediate_file
rm genotype*



# Now the interesting part begins: here we will decide a certain threshold of what we consider 'different from the average difference' between comparisons. 
# You can choose your own threshold value here. Below a threshold of 0.8 / 80% similarity is applied (as is elaborated on in the Manuscript)

cat matrix | awk 'BEGIN{ OFS=FS=","}{if($5<=0.8 || $5==0)$(NF+1)="DIFF";else if($5>0.8)$(NF+1)="NODIFF";print}' > diff1
cat diff1 | awk 'BEGIN{ OFS=FS=","}{if($9<=0.8 || $9==0)$(NF+1)="DIFF";else if($9>0.8)$(NF+1)="NODIFF";print}' > diff2
cat diff2 | awk 'BEGIN{ OFS=FS=","}{if($12<=0.8 || $12==0)$(NF+1)="DIFF";else if($12>0.8)$(NF+1)="NODIFF";print}' > matrix_diff
rm diff*
rm matrix



# We are interested to see which genotype differences show the expected patterns. 
# In our case the pattern [ Geno1!≈≈Geno2 and Geno3!≈≈Geno1, but Geno3≈≈Geno2 ] points towards absence from one chromosome (Geno1 deviant), 
# whereas the pattern [ Geno3!≈≈Geno1 and Geno3!≈≈Geno2, but Geno1≈≈Geno2 ] points towards absence from the other chromosome (Geno3 deviant)
# For more information on this, we refer again to the Manuscript

cat matrix_diff | awk -F "," 'BEGIN{ OFS=FS=","}{if ($13=="DIFF" && $14=="DIFF" && $15=="NODIFF")$(NF+1)="GENO1_diff";else $(NF+1)="FALSE";print}' > pattern1
cat pattern1 | awk -F "," 'BEGIN{ OFS=FS=","}{if ($13=="NODIFF" && $14=="DIFF" && $15=="DIFF")$(NF+1)="GENO3_diff";else $(NF+1)="FALSE";print}' > matrix_patterns
rm matrix_diff
rm pattern1



# If preferred, one can add a header to the new file (make sure those are correct and each column has a header - customize this to match your own needs)

sed -i -e '1iAnalysis,Target,Genotype1,Genotype2,Cross_1-2,PEAK_ALL,Genotype3,Genotype1,Cross_3-1,Genotype3,Genotype2,Cross_3-2,Cross1-2DIFF,Cross3-1DIFF,Cross3-2DIFF,Geno1_deviant,Geno3_deviant' matrix_patterns



# Pull out the targets for which a particular pattern was found (again, customize this to match your own needs)

awk -F ',' '$16 == "GENO1_diff" {print $2}' matrix_patterns > 1A_list
awk -F ',' '$17 == "GENO3_diff" {print $2}' matrix_patterns > 1B_list
awk -F ',' '$16 == "GENO1_diff" {print $2,$6}' matrix_patterns > 1A_list_maxpeak
awk -F ',' '$17 == "GENO3_diff" {print $2,$6}' matrix_patterns > 1B_list_maxpeak



echo "Script3 finished running. Good luck interpreting the results!"
