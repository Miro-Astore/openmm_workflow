# openmm_workflow
Here you will find a conda environment which forms a workflow for running common molecular dynamics tasks in openmm. To install perform the following:

1. Clone this repository and cd to it 

    `git clone git@github.com:miro-astore/openmm_workflow`

    `cd openmm_workflow`


* (optional) if you are on a HPC or you use modules it is recomended you make sure you load CUDA such as using 
`module load cuda`

2. Install dependencies and this environment using conda. 

    `conda create env -f environment.yml`

To perform an unbiased molecular dynamics in a slurm environment do the following. 
1. Copy the necessary `.inp` files.

``cp openmm_workflow/min_relax.inp .``

``cp openmm_workflow/unbiased_production.inp .``

* (optional) edit these files to reflect your desired minimisation, equilibration and unbiased molecular dynamics protocols.

2. Copy the necessary `.slurm` files if you are in a slurm environment. It is **recquired** that you edit these scripts to reflect your hpc and python environments. The scripts included here only serve as examples.   

``cp openmm_workflow/min_relax_equil.slurm . ``

``cp openmm_workflow/unbiased_production.slurm .``

 
2. Submit the minimisation and relaxation script which will start the workflow.

``sbatch min_relax_equil.slurm``

Openmm is extremely flexible, further customisation to the molecular dynamics protocols used in this package can be added by copying the base scripts like `unbiased.py` and `umbrella.py` in the `src` directory of this package.  
They can then be run in the following way after modification.

``` 
python customised_openmm_script.py -ff forcefield -i custom_inputs.inp 
```

To add further inputs to the molecular dynamics environment through the `.inp` interface, any extra parameters need to be reflected by a change in `src/openmm_workflow/utilities/omm_readinputs.py`.
