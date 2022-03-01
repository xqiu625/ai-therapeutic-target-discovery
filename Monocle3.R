library(monocle3)
library(dplyr)
library(stringi)
library(reshape2)
library(ggplot2)
library(Seurat)
library(patchwork)
library(magrittr)

dat0 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/harmony_int/302005_Platelets_cluster_Anno.rds")
DefaultAssay(dat0) <- "RNA"

dat0@meta.data$Disease2 <- ifelse(dat0@meta.data$Disease %in% c('Flu', 'Lung_infect', 'Non_covid'), 
                                  'Non_covid_sepsis', dat0@meta.data$Disease)
dat0@meta.data$Category <- ifelse(dat0@meta.data$Category == "", "SLE", dat0@meta.data$Category)
dat0@meta.data$seurat_clusters <- paste0("C", dat0@meta.data$seurat_clusters)

exp_mtx <- as.matrix(GetAssayData(dat0[["RNA"]], slot = "counts"))
cell_meta <- dat0@meta.data
exp_mtx2 <- exp_mtx[, cell_meta$Cell_ID2]
gene_ano <- data.frame(gene_short_name=row.names(exp_mtx))
rownames(gene_ano) <- gene_ano$gene_short_name
cds <- new_cell_data_set(exp_mtx2,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_ano)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "Data_NO")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, resolution = 0.5)
cds <- learn_graph(cds)

saveRDS(cds, '/bigdata/godziklab/shared/Xinru/302005/v5_output/Monocle3/302005_Plateletsv5_monocle3.rds')

setwd("/bigdata/godziklab/shared/Xinru/302005/v5_output/Monocle3")
dpi = 300
png(file = "Platelet Traject by cluster.png", width = dpi * 6,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster= FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)+ theme(legend.position = "right")
dev.off()

