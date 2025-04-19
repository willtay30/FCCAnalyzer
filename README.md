# Purpose
Analysis framework integrated with the FCC analysis software.  

This FCCAnalyzer framework relies on class definitions, functions and modules of the main FCC analysis framework, as described here: https://github.com/HEP-FCC/FCCAnalyses. This is necessary to read the official edm4hep Monte Carlo samples and to make use of the latest developments in terms of jet clustering and flavour tagging.  

This fork was originally meant to adapt the code at https://github.com/jeyserma/FCCAnalyzer to calculate the cross section and 
cross sectional uncertainty in the HZ decay spectrum, primarily H->bb, H->cc and H->ss. However, that code is now unsupported by 
the most recent updates to the FCC analysis software (mainly the jet clustering and tagging).  

**As such, only analysis/h_bb directory (and some various files) are being supported, other directories may or may not work!**  

To run the code here, one must correctly source your environment by running:  

```shell
source jsetup.sh
```
This allows one to use ROOT, the main frameworking for storing and processing data in this project.  

For more details on how this reposity works/used to work, consult the readme of https://github.com/jeyserma/FCCAnalyzer.  

For more details on the analysis/h_bb directory, condult that directoy's readme.  

One should also familiarize themselves with ROOT to understand the processing going on, found here: https://root.cern/manual/.  

# Combine environment

Some of the scripts at analysis/h_bb require the use of COMBINE, a CMS statistical analysis tool.  

Information about COMBINE can be found here: https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit  

**Please note that the sourcing for COMBINE and the sourcing for ROOT may conflict with each other!**  
**Only use the appropriate sourcing for whatever process you are doing!**  

To get COMBINE, one can follow one of two options:  

## Option 1 - Get the Standalone Version (Recommended)

The recommented method involves cloning your own version of COMBINE in your repository as that method has been tested and proven to work.  

The steps are outlined at the COMBINE website: https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#oustide-of-cmssw-recommended-for-non-cms-users  

This involves cloning the COMBINE directories and using make to setup the environment.
After which, one can use by using the appropriate sourcing:
```shell
cd HiggsAnalysis/CombinedLimit
source env_standalone.sh
```
This must be run in the specified directory in order to work.  


## Option 2 - Use Old Method

One may also be able to use COMBINE by following the procedure done previously, as described below. THis method is untested and therefore not recommended.  

To run fits with Combine, see https://github.com/jeyserma/FCCAnalyzer/tree/main/scripts/combine for instructions on how to install and use Combine.  

# Slurm (Work In Progress)

If resources on the local machine are limited, the batch system can be used to parallelize your analysis on multiple machines using the Slurm workload manager at SubMIT.  

Here is a sample SLURM script (.sbatch):
```shell
#!/bin/bash
#SBATCH -J h_ss_job                                       # job name
#SBATCH -o h_ss_output_%A_%a.txt              # output written to a text file with job number
#SBATCH -e h_ss_error_%A_%a.txt                # errors written to a text file with job number
#SBATCH -t 10:00:00                                        # time limit for the job (HH:MM:SS)
#SBATCH -n 1                                                   # number of tasks
#SBATCH --cpus-per-task=4                             # number of CPU cores per task
#SBATCH --mem=64GB                                   # memory allocated per node
#SBATCH -p submit                                          # partition name (pick between submit and submit-gpu)
#SBATCH --mail-type=BEGIN,END                  # send an email when done
#SBATCH --mail-user=isabellalynn622@gmail.com # email address to send confirmation to

# import the sourcing for using fccanalysis and location of the python folder
source /work/submit/jaeyserm/software/FCCAnalyses/setup.sh
export PYTHONPATH=$PYTHONPATH:/home/submit/isanford/FCCAnalyzer/

# running the script
fccanalysis run FCCAnalyzer/analyses/h_bb/h_ss.py
```

For you to run it, change the file paths for your desired script to run. This file should be stored in your home directory on submit (eg. /home/submit/<user>). Also edit the email address - it will send an email once the job is complete so you don't have to constantly log in and check.   

To run, login to submit from your terminal and run the command "sbatch [file name].sbatch". This will run the job on Slurm. Similar to what we've done before, there will be output text files directly viewable so you can see the progress of the script. Some helpful terminal commands:  

"sinfo" shows available partitions and their status  
"squeue" shows all jobs lined up in Slurm  
"squeue -u $USER" shows which jobes you have running on Slurm  
"scancel <job number>" will cancel said job  
"scontrol show job <job number>" will provide information on the current job  
"top -u $USER" will show live memory allocation (useful to determine if you will run out of memory)  

Please note that the job may take longer to run than normal execution and that the job may fail due to using more memory than specified.  

SLURM support is currently a work in progress.  

# Acknowledgements

Thank to the team behind the FCC analysis software, which is critical to our project.  

Special thank you to an Jan Eysermans, Luca Lavezzo and Christoph Paus from MIT. Their advice and support, along with their letting of us
using MIT computational resources, scripts, and data, made this project possible.  

Thank you to Christopher Palmer for his explanations, guidance and support through development.  

This code adapts much from code developed by Sarah Waldych, Jacob Lee, and Caitlin Kubina. Thank you to them for allowing us to use their 
code and aiding us in building on it. Their code can be found here: https://github.com/IOKnight/FCCAnalyzer .  

Contributers: Aniket Gullapalli, Isabella Sanford, William Taylor
