#!/bin/bash
#SBATCH -J runSeuratSingleSample
#SBATCH -e error/runSeuratSingleSample%j_%A_%a.txt
#SBATCH -o output/runSeuratSingleSample%j_%A_%a.txt
#SBATCH -t 06:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=small
#SBATCH --array=1-9


# Load r-env-singularity
module load r-env-singularity


# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2005392" >> ~/.Renviron


# prepare the paths and sample names
scRNAseqSample=$(sed -n ${SLURM_ARRAY_TASK_ID}p msSamples)

# Path to where the gene expression matrices are for each sample
JoonaDataPath="/scratch/project_2005392/Joona/MS_FNA/scRNAseq/"
samp=${JoonaDataPath}${scRNAseqSample}

# prepare output directory
mkdir ${scRNAseqSample}_results
cd ${scRNAseqSample}_results


# Run the R script
srun singularity_wrapper exec Rscript --no-save /scratch/project_2005392/JoonaDawitAnalysisResults/analyzeSinglescRNAseqData.R ${samp}

