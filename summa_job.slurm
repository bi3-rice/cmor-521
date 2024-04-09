#!/bin/bash
#SBATCH --job-name=Summa
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00 
#SBATCH --mail-user=bi3@rice.edu 
#SBATCH --mail-type=ALL   
echo "My job ran on:" echo $SLURM_NODELIST
if [[ -d $SHARED_SCRATCH/$USER && -w $SHARED_SCRATCH/$USER ]]; then
    cd $SHARED_SCRATCH/$USER
    srun /home/bi3/summa.cpp
fi