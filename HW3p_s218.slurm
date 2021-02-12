#!/bin/bash -l
#SBATCH --job-name=HW3p_s218	    	 # Job name
#SBATCH --time=15-10:40:00                # Time limit days-hrs:min:sec
#SBATCH --partition=shared		 # Check partitions on the Wiki
#SBATCH --mem=30Gb			 # Memory size
#SBATCH -n 28				 # Number of processors in the node
#SBATCH --mail-user=p_s218@txstate.edu	 # Email me!
#SBATCH --mail-type=end,begin		 # Types of Emails (Start,End)

# Load Module
module load gcc/7.2.0

# Go to directory
cd /home/p_s218/hpc/hw3

# compiles program
g++ HW3p_s218.cpp -fopenmp -o HW3p_s218 

# execute program
./HW3p_s218 > outputHW3.txt

# Unload module
module unload gcc/7.2.0

# End the Job
exit 0
