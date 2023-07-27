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
module load cuda/12.1
export OPENMM_CUDA_COMPILER=/mnt/home/gkrawezik/local/rocky8/cuda/12.1/bin/nvcc

#module load cuda/11.8
#export OPENMM_CUDA_COMPILER=/mnt/sw/nix/store/3xpm36w2kcri3j1m5j15hg025my1p4kx-cuda-11.8.0/bin/nvcc

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
conda activate openmm_hpc

init="out"
equil_prefix="equil"
prod_prefix="production"
prod_step="step5"

pstep=${equi_prefix}

if [ ! -e prod_1.rst ] ; 
then 
input_param=" -ff amber -t toppar.str -p ${init}.parm7 -c $init.rst7  -irst $equil_prefix.rst --restart-timer "
python -u openmm_run.py -i ${prod_prefix}.inp ${input_param} -orst prod_1.rst -odcd prod_1.dcd  --platform CUDA --restart-timer | tee prod_1.out 
else 
prev_steps=$(ls prod_*\.rst  | wc -l)
curr_step=$(($prev_steps + 1 ))
prod_name=prod_$curr_step
input_param=" -ff amber -t toppar.str -p ${init}.parm7 -c ${init}.rst7  -irst prod_${prev_steps}.rst "
python -u openmm_run.py -i ${prod_prefix}.inp ${input_param} -orst ${prod_name}.rst -odcd ${prod_name}.dcd  --platform CUDA | tee ${prod_name}.out 
fi 

nbatch () {
    sbatch -J $( echo "$(pwd | awk -F/ '{print $(NF-2)$(NF-1)$NF}' | sed "s^_^^g" | sed "s^memb^^g" | sed "s^plumed^^g" | sed "s^CFTR^^g" | sed "s^simreadysystems^^" | sed "s^process_traj^^g"| sed "s^process_equil^^g")$(echo $1 | sed "s^\.pbs^^g"  | sed "s^\.slurm^^g" | sed "s^_^^g" | sed "s^memb^^g"  | sed "s^glob^^g" | sed "s^window^^g" | sed "s^long^^g" | sed "s^mart^^g" )" | tail -c 11 ) $1
    }
 
line=$(squeue -h -j $SLURM_JOBID -O timeused) # Leading '-' aids parsing.
line=$(echo -$line)
parts=( 0 0 0 0 )
index=3
while [ ${#line} -gt 0 ]; do
  parts[${index}]=${line##*[-:]}
    line=${line%[-:]*}
      ((index--))
      done
root=$(echo ${parts[*]})
used_d=$(echo $root | awk  '{print $1}')
used_h=$(echo $root | awk  '{print $2}')
total_h=$(( $used_d * 24 + $used_h))
if $total_h
if [  $total_h -gt 0 ] ;
then
echo 'hi'
nbatch prod.slurm
fi
exit



