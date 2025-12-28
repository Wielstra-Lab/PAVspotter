# PAVspotter


Hello! Welcome to the PAV-spotter pipeline.

---

In order to run the scripts succesfully, you may need to manually adjust some;
- identifiers, 
- pathnames, 
- filenames, 
- settings, 
- and more
 
Thus, make sure you download the scripts and adjust them to your own need.

We briefly explain how best to do this below.
Also, we refer to the Manuscript and the explanatory text in each of the Scripts for more detail.

---

Running Script1 (batch/SLURM)

- Step 1: Have your BAM files ready
- Step 2: Fill in (or remove) the #SBATCH lines + customize Script1 manually for running remotely
- Step 3: Run Script1 to generate input files/format required for PAV-spotter (Script2)

Running Script2 (batch/SLURM)

- Step 4: Fill in (or remove) the #SBATCH lines + manually adjust Script2a for running remotely. Note there are two versions of Script2a: a version to initiate MATLAB, or a version to initiate Python - depending on your preferences 
- (Extra step: Customize Script2b manually if necessary, but this is less likely to be necessary - again note the '.m' version is the MATLAB version, and the '.py' version is the Python version) 
- Step 5: Run your preferred version of Script2a (MATLAB or Python), which will automatically kickstart Script2b, to generate output files/format required for Script3

Running Script3 (bash)

- Step 6: Manually adjust Script3 for running locally
- Step 7: Check if everything worked properly

Did it all run succesfully? Great job!
Are you still having trouble/questions? Feel free to contact the creators of PAV-spotter 

---

Good luck!

--

PS: An example batch script (to be used in a Windows environment) was added in this repository as well
It is meant for taking screenshots in IGV in an automated way, please see the Manuscript for more details on this


