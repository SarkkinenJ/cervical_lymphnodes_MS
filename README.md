# Cervical lymph nodes in MS <br/>
The following summarizes the analysis steps used to perform the analysis and produce key results in the manuscript "Deep cervical lymph nodes of patients with multiple sclerosis show dysregulated B cells in the presence of Epstein-Barr virus" by Sarkkinen et al., (submitted).

Here, the analysis workflow is indicated along with a set of the main scripts used in the the study including for the scRNAseq, CITE-seq, and TCR- and BCR-analyses. To reproduce the key results, obtain the single cell gene expression, TCRseq and BCRseq data for each of the MS and Control sample, or/and the main Seurat objects from EGA (accession number to be provided). Put the obtained Seurat objects in the **data** directory. Then, clone this repository and load the script needed in R, install the necessary R packages if not installed, and follow the steps below, and the comments in the scripts. 

### scRNAseq analysis <br/>
1. Preprocess and analyze individual scRNAseq samples (in **single_sample_analysis**): use Slurm to run the shell script **runSeuratInBatches.sh** (modify as needed) in an High performance computing (HPC) environment or use the R script **analyzeSinglescRNAseqData.R** directly. 
   - Basic analysis using the Seurat clustering workflow (QC, dimensionality reduction, clustering)
   - Automatic cell annotation using singleR
   
2. Merge or integrate the individual scRNAseq samples (in **merge_multiple_samples**): use Slurm to run the shell script **mergeAndIntegrateUsingSeurat.sh** (modify as needed) in an HPC environment or use the R script **analyzeCombinedMSFNAdata.R** directly.
   - Merge samples into one Seurat object for analysis
   - Integrate samples into one Seurat object for analysis

3. Perform visual batch effect evaluation, group-wise comparison analyses and subclustering (in **multiSample_scRNAseq_analysis**): use the R scripts **AnalysesInMergedData_MSFNA.R** and **CD8_subclustering.R**. 
   - Evaluation of batch effect
   - Annotation of cells guided by CiteSeq and manual inspection
   - Cell type differential Abundance and differential expression analysis (Pseudobulk and single cell DE)
   - B cell subclustering
   - CD8 T cell subclustering
  
### scTCR and scBCR analysis <br/>
4. TCR and BCR analyses (in **scTCR_BCR_analyses**): use the R scripts **MSFNA_MergedData_TCRanalysis.R** and **MSFNA_MergedData_BCRanalysis.R**. 
   - Merge scRNAseq data with corresponding scTCR and BCR data
   - scTCR and scBCR comprehensive analysis (diversity, clonality, vgene usage, comparison between MS and CTRLs)
   - TCR scanning against the VDJdb database for specificity prediction using TCRdist3

### Intercellular communication analysis <br/>
5. Intercellular communication analysis using NicheNet (in **nichenet_analysis**): use the R scripts **nichenet_analysis.R** and **multinichenet_analysis.R**
    - Ligand identification and prioritizatin of ligand-receptor pairs 

<br/>
For comprehensive usage guidelines and further details, please reach out to the authors.

