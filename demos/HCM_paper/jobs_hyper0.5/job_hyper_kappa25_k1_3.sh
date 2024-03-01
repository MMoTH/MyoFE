#!/bin/bash
#SBATCH --time=70:00:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=T25_K3
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c).
#SBATCH --mem=120G                  # memory required per node - amount of memory (in bytes)
#SBATCH --account=col_jfwe223_uksr
#SBATCH --partition=CAL48M192_L# Partition/queue to run the job in. (REQUIRED)
#SBATCH --output=../../../../../../mnt/gpfs2_4m/scratch/mme250/HCM_paper/hyper0.5/logs/output.%J.out # STDOUT
cd ../../../python_codes
singularity exec --cleanenv /home/mme250/fenics.img  mpiexec -np $SLURM_NTASKS  python MyoFE.py LV_sim /home/mme250/MYoFE/demos/HCM_paper/simulations_hyper0.5/hyper_kappa25_k1_3/sim_inputs/input_parameters.json
scontrol show job $SLURM_JOB_ID