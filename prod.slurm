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

export MODULEPATH=/mnt/home/gkrawezik/modules/rocky8:$MODULEPATH
#export OPENMM_CUDA_COMPILER=/mnt/home/gkrawezik/local/rocky8/cuda/12.1/bin/nvcc

module load cuda/12.1
#module load cuda/11.8
export OPENMM_CUDA_COMPILER=/mnt/sw/nix/store/3xpm36w2kcri3j1m5j15hg025my1p4kx-cuda-11.8.0/bin/nvcc
#
conda deactivate &> /dev/null

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/mnt/home/mastore/Software/cluster_conda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
eval "$__conda_setup"
else
if [ -f "/mnt/home/mastore/Software/cluster_conda/etc/profile.d/conda.sh" ]; then
. "/mnt/home/mastore/Software/cluster_conda/etc/profile.d/conda.sh"
else
export PATH="/mnt/home/mastore/Software/cluster_conda/bin:$PATH"
fi
fi
unset __conda_setup
# <<< conda initialize <<<
conda activate openmm_hpc2

init="step3_input"
equi_prefix="step4_equilibration"
prod_prefix="step5_production"
prod_step="step5"
istep="step5_1"

# Equilibration
#input_param=" -ff amber -t toppar.str -p ${init}.parm7 -c ${init}.rst7"
#python prep_restraints.py -p ${init}.psf -ipdb ${init}.pdb
#python -u openmm_run.py -i ${equi_prefix}.inp ${input_param} -orst ${equi_prefix}.rst -odcd ${equi_prefix}.dcd 

# Production
# The OpenMM check point file (.chk) cannot be used in a different machine environment.
# So please make sure if you are using the same GPU and CUDA version of machine while doing additional
# production steps with the check point file.

latest_rst=$(ls -v *.rst | tail -n 1 )

pstep=${equi_prefix}
#input_param="-t toppar.str -p ${init}.psf -ipdb ${init}.pdb -irst ${pstep}.rst --restart-timer "
input_param=" -ff amber -t toppar.str -p ${init}.parm7 -c ${init}.rst7 -irst ${pstep}.rst"
python -u openmm_run.py -i ${prod_prefix}.inp ${input_param} -orst ${istep}.rst -odcd ${istep}.dcd  --platform CUDA | tee ${istep}.out
