#!/bin/bash
#SBATCH -J runTCRdbScan_msCSF
#SBATCH -e /scratch/project_2005392/JoonaDawitAnalysisResults/error/runTCRdbScan%j_%A_%a.txt
#SBATCH -o /scratch/project_2005392/JoonaDawitAnalysisResults/output/runTCRdbScan%j_%A_%a.txt
#SBATCH -t 7-00:00:00
#SBATCH --account=project_2005392
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=200000
#SBATCH --partition=longrun
#SBATCH --array=1-8


# Load the right python environment
module load r-env/421
module load python-data

cd /scratch/project_2005392/JoonaDawitAnalysisResults/

# scan the VDJdb with tcrdist3
srun singularity_wrapper exec Rscript --no-save /scratch/project_2005392/JoonaDawitAnalysisResults/tcr_VDJDBScan.R "paired" ${SLURM_ARRAY_TASK_ID}

#srun singularity_wrapper exec Rscript --no-save /scratch/project_2005392/JoonaDawitAnalysisResults/tcr_VDJDBScan.R "beta"