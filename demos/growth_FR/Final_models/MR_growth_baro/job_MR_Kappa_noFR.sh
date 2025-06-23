#!/bin/bash
#SBATCH --time=13-00:00:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=gr_baro
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c).
#SBATCH --mem=180G                  # memory required per node - amount of memory (in bytes)
#SBATCH --account=col_jfwe223_uksr
#SBATCH --partition=CAL48M192_L# Partition/queue to run the job in. (REQUIRED)
#SBATCH --output=/mnt/gpfs2_4m/scratch/mme250/gr_paper/final_models/MR_growth_baro/output.%J.out # STDOUT
cd ../../../../python_codes
singularity exec --cleanenv /home/mme250/fenics.img  mpiexec -np $SLURM_NTASKS  python MyoFE.py LV_sim /home/mme250/MyoFE/demos/growth_FR/Final_models/MR_growth_baro/input_parameters_noFR.json
scontrol show job $SLURM_JOB_ID