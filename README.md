# PAVspotter

Hello! Welcome to the PAV-spotter pipeline.

In order to run the scripts succesfully, you may need to manually adjust some identifiers, path- and filenames, and settings. 
So make sure you download the scripts and adjust them to your own need. We briefly explain how best to do this below and
refer to the Manuscript and the explanatory text in each of the Scripts for more detail.



Step 1: Running Script1

- Have your BAM files ready
- Fill in (or remove) the #SBATCH lines + manually adjust Script1 for running locally
- Run Script1 to generate input files/format required for PAV-spotter, Script2



Step 2: Running Script2

- Fill in (or remove) the #SBATCH lines + manually adjust Script2a
- Manually adjust Script2b if necessary (in general, this is likely unnecessary) 
- Run Script2a, which will call automatically Script2b, to generate output files/format required for Script3



Step 3: Running Script3

- 



