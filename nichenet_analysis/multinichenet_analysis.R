library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr) # v. 1.0.3

options(timeout = 600)

load("nichenet.RData")

colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()

colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()

sce = Seurat::as.SingleCellExperiment(seuratObj, assay = "RNA")

sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

sample_id = "sampleName"
group_id = "sampleGroup"
celltype_id = "rnaclustersvsAdtManAnn2MaxAnn2"
covariates = NA
batches = NA

# Define senders and receivers
senders_oi = c("GC_B_cells", "Memory_B_cells", "Plasmablasts", "Naive_B_cells")
receivers_oi = c("Treg_", "Naive_CD4", "Naive_CD8", "Memory_CD8", "Tfh.1", "Tfh.2", "Tfr")

# Subset sce to only sender and receiver cells
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]

min_cells = 10

abundance_expression_info = get_abundance_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, batches = batches)

abundance_expression_info$abund_plot_sample

contrasts_oi = c("'MS-CTRL','CTRL-MS'")
contrast_tbl = tibble(contrast = c("MS-CTRL","CTRL-MS"), group = c("MS","CTRL"))

# DE analysis

DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells)
DE_info$celltype_de$de_output_tidy %>% arrange(p_adj) %>% head()
celltype_de = DE_info$celltype_de$de_output_tidy

sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

sender_receiver_de %>% head(20)

logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.05

p_val_adj = FALSE # should be FALSE if <50 DE genes per group-celltype
empirical_pval = FALSE

top_n_target = 250

# Comment out when only one receiver celltype
verbose = TRUE
cores_system = 8
n.cores = min(cores_system, union(senders_oi, receivers_oi) %>% length())

# Ligand activity prediction (will take a while)

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose, 
  n.cores = n.cores
)))

ligand_activities_targets_DEgenes$de_genes_df %>% head(20)

ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

# Define prioritization weights for ranking communication patterns

prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)

prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
))

prioritization_tables$group_prioritization_tbl %>% head(20)

lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)


group_oi = "MS"

# zoom to specific receivers or senders with e.g. receivers_oi = c("GC_B_cells", "Memory_B_cells", "Naive_B_cells")
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, groups_oi = group_oi)

# Fig. S1H
plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_M_50)
plot_oi