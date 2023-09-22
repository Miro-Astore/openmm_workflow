#!/bin/sh
#SBATCH -p gpu
#SBATCH -C h100
#SBATCH -n 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=8
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-00:00:00

module purge
module load modules/2.1.1-20230405

#you can use newer cuda if you like but be careful it needs to be the one you installed openmm with.
#export MODULEPATH=/mnt/home/gkrawezik/modules/rocky8:$MODULEPATH
#module load cuda/12.1
#export OPENMM_CUDA_COMPILER=/mnt/home/gkrawezik/local/rocky8/cuda/12.1/bin/nvcc

#load in cuda and tell openmm where to find the cuda compiler 
module load cuda/11.8
export OPENMM_CUDA_COMPILER=/mnt/sw/nix/store/3xpm36w2kcri3j1m5j15hg025my1p4kx-cuda-11.8.0/bin/nvcc


conda activate openmm_workflow

#number of steps to run
n_steps=4
init="out"
equil_prefix="equil_pull"
prod_prefix="umbrella_production"

if [ ! -e prod_1.rst ] ; 
then 

    input_param=" -ff amber -t toppar.str -p ${init}.parm7 -c $init.rst7  -irst $equil_prefix.rst --restart-timer "
    python -m openmm_workflow.umbrella -i ${prod_prefix}.inp ${input_param} -orst prod_1.rst -odcd prod_1.dcd  --platform CUDA --restart-timer 

else 
    prev_steps=$(ls prod_*\.rst  | wc -l)
    curr_step=$(($prev_steps + 1 ))

    prod_name=prod_$curr_step
    input_param=" -ff amber -t toppar.str -p ${init}.parm7 -c ${init}.rst7  -irst prod_${prev_steps}.rst "
    python -m openmm_workflow.umbrella -i ${prod_prefix}.inp ${input_param} -orst ${prod_name}.rst -odcd ${prod_name}.dcd  --platform CUDA 

fi 

sbatch prod.slurm
