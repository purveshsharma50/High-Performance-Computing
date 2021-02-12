#!/bin/bash -l
#SBATCH --job-name=.         	    	 # Job name
#SBATCH --time=15-10:40:00                # Time limit days-hrs:min:sec
#SBATCH --partition=shared		 # Check partitions on the Wiki
#SBATCH --mem=30Gb			 # Memory size
#SBATCH -n 1				 # Number of processors in the node
#SBATCH --mail-user=p_s218@txstate.edu	 # Email me!
#SBATCH --mail-type=end,begin		 # Types of Emails (Start,End)

# Load Module
module load intel/18

# Go to directory
cd /gpfs/home/

# compiles program
icc -o HW2p_s218 HW2p_s218.cpp

# execute program
./HW2p_s218 > hw2final.txt

# Unload module
module unload intel/18

# End the Job
exit 0
