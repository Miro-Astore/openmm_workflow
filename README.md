# openmm_workflow
Here you will find a conda environment which forms a workflow for running common molecular dynamics tasks in openmm. To install perform the following:

1. Clone this repository and cd to it 
    `git clone git@github.com:miro-astore/openmm_workflow`
    `cd openmm_workflow`


(optional) if you are on a HPC or you use modules it is recomended you make sure you load CUDA such as using 
`module load cuda`

2. Install dependencies and this environment using conda. 
    `conda create env -f environment.yml`

To perform an unbiased molecular dynamics in a slurm environment do the following. 

