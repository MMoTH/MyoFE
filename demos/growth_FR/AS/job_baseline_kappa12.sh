#!/bin/bash
#SBATCH --time=13-00:00:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=b_kp12
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c).
#SBATCH --mem=180G                  # memory required per node - amount of memory (in bytes)
#SBATCH --account=col_jfwe223_uksr
#SBATCH --partition=SKY32M192_L# Partition/queue to run the job in. (REQUIRED)
#SBATCH --output=/mnt/gpfs2_4m/scratch/mme250/AS_paper/output.%J.out # STDOUT
cd ../../../python_codes
singularity exec --cleanenv /home/mme250/fenics.img  mpiexec -np $SLURM_NTASKS  python MyoFE.py LV_sim /home/mme250/MyoFE/demos/growth_FR/AS/input_parameters_baseline_kappa12_cvein0.05.json
scontrol show job $SLURM_JOB_ID
