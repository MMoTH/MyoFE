""" This file is aimed to generate a batch of jobs to run on LCC cluster """


import os
import json 

temp_json_input_str = '../../../demos/base/sim_inputs/base_instruction_lcc_apex.json'
job_folder_str = 'cores_study'
job_ext_str = '_core_sensetivity_'
if not os.path.isdir(job_folder_str):
    print('Making output dir')
    os.makedirs(job_folder_str)

with open(temp_json_input_str, 'r') as f:
    json_file = json.load(f)
    
num_of_iter = [1,2,3]
for c in num_of_iter:
    num_of_cores = int(2**c)
    for iter in num_of_iter:
        #change the json name
        folder_name = '../../../demos/' + str(c) + '_cores_'+str(iter)

        if not os.path.isdir(folder_name):
            print('Making output dir')
            os.makedirs(folder_name)
            os.makedirs(folder_name+'/sim_inputs')
            
        new_json_str = 'instruction_' + str(c) + '_cores_' + str(iter) + '.json'
        temp_instruction = json_file
        temp_instruction['output_handler']['output_data_path'][0] = \
            folder_name + '/sim_output/data.csv'
        temp_instruction['output_handler']['mesh_output_path'][0] = \
            folder_name + '/sim_output/mesh_output'
        new_json_str = folder_name + '/sim_inputs/' + new_json_str

        json_object = json.dumps(temp_instruction, indent = 4)
        with open(new_json_str, 'w') as fo:
            fo.write(json_object)
        #change the outputs
        #write the job file 
        job_str = job_folder_str + '/' +\
                 str(c)+job_ext_str+str(iter) + '.sh'
        with open(job_str , 'w') as jo:
            jo.write('#!/bin/bash\n')
            jo.write('#SBATCH --time=36:00:00             # Time limit for the job (REQUIRED).\n')
            jo.write('#SBATCH --job-name=Hossein_jobs     # Job name\n')
            jo.write(f'#SBATCH --ntasks={num_of_cores}                 # Number of cores to allocate. Same as SBATCH -n\n')
            jo.write('#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c).\n')
            jo.write('#SBATCH --mem=120G                  # memory required per node - amount of memory (in bytes)\n')
            jo.write('#SBATCH --account=col_jfwe223_uksr\n')
            jo.write('#SBATCH --partition=CAL48M192_L     # Partition/queue to run the job in. (REQUIRED)\n')
            jo.write('#SBATCH -o ./cores_study/output.%J.out # STDOUT\n')
            jo.write('cd ../FEniCS-Myosim/python_codes\n')
            jo.write("singularity exec --cleanenv /home/hsh245/fenics.img  mpiexec -np $SLURM_NTASKS  python MyoFE.py LV_sim /home/hsh245/FEniCS-Myosim/"+\
                folder_name + "/sim_inputs/" + new_json_str + '\n')
            jo.write('scontrol show job $SLURM_JOB_ID')


