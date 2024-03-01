#!/bin/bash
#SBATCH --time=70:00:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=kapa30_K1.8
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c).
#SBATCH --mem=120G                  # memory required per node - amount of memory (in bytes)
#SBATCH --account=col_jfwe223_uksr
#SBATCH --partition=CAL48M192_L# Partition/queue to run the job in. (REQUIRED)
#SBATCH -o output.%J.out # STDOUT
cd ../../../python_codes
singularity exec --cleanenv /home/mme250/fenics.img  mpiexec -np $SLURM_NTASKS  python MyoFE.py LV_sim /home/mme250/MYoFE/demos/HCM_paper/simulations_hyper0.3/hyper_kappa30_k1_1.8/sim_inputs/input_parameters.json
scontrol show job $SLURM_JOB_ID