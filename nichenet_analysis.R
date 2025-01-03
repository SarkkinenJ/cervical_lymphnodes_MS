library(nichenetr) # v. 2.0.1
library(Seurat)
library(SeuratObject)
library(tidyverse)

load("msFNAmerged.goodSamples.RData")

# Load in NicheNet networks
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
lr_network <- lr_network %>% distinct(from, to)

ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))

weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

# Load in data
seuratObj <- msFNAmerged.goodSamples
seuratObj <- alias_to_symbol_seurat(seuratObj, "human")
Idents(seuratObj) <- "rnaclustersvsAdtManAnn2MaxAnn2"

receivers <- c("Naive B cells", "GC B cells")

for (receiver in receivers) {
  
  # Define set of potential ligands for receivers, i.e., ligands in our dataset that also appear in NicheNet networks

  expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.1)
  
  sender_celltypes <- c("Intermediate CD4", "Naive CD4", "Tfh.1", "Tfh.2", "Tfr")
  list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.1)
  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
  
  all_ligands <- unique(lr_network$from)
  all_receptors <- unique(lr_network$to)
  
  expressed_ligands <- intersect(all_ligands, expressed_genes_sender)
  expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
  
  potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
    pull(from) %>% unique()
  
  # Define gene set of interest
  
  condition_oi <-  "MS"
  condition_reference <- "CTRL"
  
  seurat_obj_receiver <- subset(seuratObj, idents = receiver)
  
  DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                    ident.1 = condition_oi, ident.2 = condition_reference,
                                    group.by = "sampleGroup", 
                                    min.pct = 0.01) %>% rownames_to_column("gene")
  
  geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  
  # Define background genes
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  # Perform NicheNet ligand activity analysis
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                 background_expressed_genes = background_expressed_genes,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 potential_ligands = potential_ligands)
  
  ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>%
    mutate(rank = rank(desc(aupr_corrected)))
  
  # Perform prioritization of ligand-receptor pairs based on condition-specificities and expression values in cell types of interest
  
  lr_network_filtered <-  lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
  
  info_tables <- generate_info_tables(seuratObj,
                                      celltype_colname = "rnaclustersvsAdtManAnn2MaxAnn2",
                                      senders_oi = sender_celltypes,
                                      receivers_oi = receiver,
                                      lr_network = lr_network_filtered,
                                      condition_colname = "sampleGroup",
                                      condition_oi = condition_oi,
                                      condition_reference = condition_reference,
                                      scenario = "case_control")
  
  prioritizing_weights = c("de_ligand" = 1,
                           "de_receptor" = 1,
                           "activity_scaled" = 1,
                           "exprs_ligand" = 1,
                           "exprs_receptor" = 1,
                           "ligand_condition_specificity" = 2,
                           "receptor_condition_specificity" = 2)
  
  prior_table <- generate_prioritization_tables(info_tables$sender_receiver_info,
                                                info_tables$sender_receiver_de,
                                                ligand_activities,
                                                info_tables$lr_condition_de,
                                                prioritizing_weights,
                                                scenario = "case_control")
  

  # Create mushroom plots (Fig. 1F)
  
  legend_adjust <- c(0.7, 0.7)
  
  make_mushroom_plot(prior_table %>% filter(receiver == receiver),
                     top_n = 30,
                     size = "pct_expressed",
                     color = "scaled_avg_exprs",
                     true_color_range = TRUE,
                     show_rankings = FALSE,
                     show_all_datapoints = FALSE) +
    theme(legend.justification = legend_adjust,
          axis.title.x = element_text(hjust = 0.25))
  
  ggsave(filename = paste0(receiver, "_nichenet_mushroom.pdf"), device = "pdf", height = 6, width = 8, units = "in")
}

  
save.image(file = "nichenet.RData")
