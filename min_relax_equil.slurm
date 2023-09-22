#!/bin/sh
#SBATCH -p gpu
#SBATCH -C a100
#SBATCH -n 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=8
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-00:00:00

module purge
module load modules/2.1.1-20230405


module load cuda/11.8
export OPENMM_CUDA_COMPILER=/mnt/sw/nix/store/3xpm36w2kcri3j1m5j15hg025my1p4kx-cuda-11.8.0/bin/nvcc

conda activate openmm_workflow

init="out"
relax_prefix="min_relax"
equil_prefix="equil"

# minimisation, relaxation
input_param=" -ff amber -t toppar.str -p ${init}.parm7 -c ${init}.rst7"
python prep_restraints.py -p ${init}.psf -ipdb ${init}.pdb
python -u openmm_workflow.openmm_run -i ${relax_prefix}.inp ${input_param} -orst ${relax_prefix}.rst -odcd ${relax_prefix}.dcd | tee ${equil_prefix}.out

# equilibration
input_param=" -ff amber -t toppar.str -p ${init}.parm7 -c ${init}.rst7 -irst ${relax_prefix}.rst --restart-timer" 
python -u openmm_workflow.openmm_run -i ${equil_prefix}.inp ${input_param} -orst ${equil_prefix}.rst -odcd ${equil_prefix}.dcd  | tee ${equil_prefix}.out

nbatch () {
    sbatch -J $( echo "$(pwd | awk -F/ '{print $(NF-2)$(NF-1)$NF}' | sed "s^_^^g" | sed "s^memb^^g" | sed "s^plumed^^g" | sed "s^CFTR^^g" | sed "s^simreadysystems^^" | sed "s^process_traj^^g"| sed "s^process_equil^^g")$(echo $1 | sed "s^\.pbs^^g"  | sed "s^\.slurm^^g" | sed "s^_^^g" | sed "s^memb^^g"  | sed "s^glob^^g" | sed "s^window^^g" | sed "s^long^^g" | sed "s^mart^^g" )" | tail -c 11 ) $1
    }
if [ -e equil.rst ] ;
then
    nbatch prod.slurm
fi
