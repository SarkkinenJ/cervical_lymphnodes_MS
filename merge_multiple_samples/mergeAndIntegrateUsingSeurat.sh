#!/bin/bash
#SBATCH -J mergeAndIntegrate2
#SBATCH -e error/mergeAndIntegrate2%j_%A_%a.txt
#SBATCH -o output/mergeAndIntegrate2%j_%A_%a.txt
#SBATCH -t 06:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100000
#SBATCH --partition=small
#SBATCH --array=1


# Load r-env-singularity
module load r-env-singularity


# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

echo "TMPDIR=/scratch/project_2005392" >> ~/.Renviron



cd combinedAnalysis

srun singularity_wrapper exec Rscript --no-save /scratch/project_2005392/JoonaDawitAnalysisResults/analyzeCombinedMSFNAdata.R

