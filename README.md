# Cervical lymph nodes in MS <br/>
Following summarizes the analysis steps used in the manuscript "Deep cervical lymph nodes of patients with multiple sclerosis show dysregulated B cells in the presence of Epstein-Barr virus" by Sarkkinen et al (to be updated).

The main data analysis scripts and the analysis workflow used in the study including for the scRNAseq, CITE-seq, and TCR- and BCR-analyses are listed here: 

## scRNAseq analysis
1. Preprocess and analyze individual scRNAseq samples (single_sample_analysis):
  - Basic analysis using the Seurat clustering workflow (QC, dimensionality reduction, clustering)
  - Automatic cell annotation using singleR

2. Merge or integrate the individual scRNAseq samples (merge_multiple_samples): 
  - Merge samples into one Seurat object for analysis
  - Integrate samples into one Seurat object for analysis

3. Perform visual batch effect evaluation, group-wise comparison analyses and subclustering:
  - Evaluation of batch effect
  - Annotation of cells guided by CiteSeq data and manual inspection
  - Cell type differential Abundance and differential expression analysis (Pseudobulk and single cell DE)
  - B cell subclustering
  - CD8 T cell subclustering
  - 
  
## scTCR and scBCR analysis
1. TCR and BCR analyses (TCR_BCR_analyses):

## sc Intercellular communication analysis
1. Intercellular communication analysis using NicheNet (nichenet_analysis):


