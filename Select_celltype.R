.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )

library(SeuratDisk)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(celldex)
library(SingleR)
library(HGNChelper)
library(tidyr)
library(stringr)

# datasets ----------------------------------------------------------------


df302004data15 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets/302004data15/302004data15_covid19.rds")
df302004data16 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets/302004data16/302004data16_covid19.rds")
df302004data17 <- LoadH5Seurat("/bigdata/godziklab/shared/Xinru/302005/datasets/302004data17/covid_portal_210320_with_raw.h5seurat")
df302004data18 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets/302004data18/302004data18_sepsis.rds")
df302005data01 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets/302005data01/302005data01_covid19.rds")
df302005data02 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets/302005data02/302005data02_covid19.rds")
df302005data03 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets/302005data03/labeled_clustered_integrated_sepsis_0729.rds")
df302005data04 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets/302005data04/302005data04_sle.rds")
DefaultAssay(df302004data17) <- "raw"


# metadata ----------------------------------------------------------------

Anno <- read.delim2('/bigdata/godziklab/shared/Xinru/302005/datasets/per_cell_type_annotation.tsv')
Anno <- Anno %>% 
  filter(pbmc == "Y") %>% 
  filter(Category != "HC_LPS")


df302004data15@meta.data$Cell_ID <- paste0(row.names(df302004data15@meta.data), "_5")
df302004data16@meta.data$Cell_ID <- paste0(row.names(df302004data16@meta.data), "_6")
df302004data17@meta.data$Cell_ID <- paste0(row.names(df302004data17@meta.data), "_7")
df302004data18@meta.data$Cell_ID <- paste0(row.names(df302004data18@meta.data), "_8")
df302005data01@meta.data$Cell_ID <- paste0(row.names(df302005data01@meta.data), "_9")
df302005data02@meta.data$Cell_ID <- paste0(row.names(df302005data02@meta.data), "_10")
df302005data03@meta.data$Cell_ID <- paste0(row.names(df302005data03@meta.data), "_11")
df302005data04@meta.data$Cell_ID <- paste0(row.names(df302005data04@meta.data), "_12")
  
# Integration with corrected gene names ----------------------------------------------------
setwd("/bigdata/godziklab/shared/Xinru/302005/datasets_v2/processed_datasets_v2/correct_gene_name")
file <- list.files(pattern = "*.txt")
gene_name <- do.call("rbind",lapply(file,
                                    FUN=function(file){read.delim2(file)}))
gene_name <- gene_name %>% distinct(x, .keep_all = T)

# Process seurat objects ready for Seurat v4 map--------------------------------------------------
process0 <- function(obj){
  meta <- obj@meta.data
  mat <- obj@assays$RNA@counts
  dat <- as.data.frame(row.names(mat))
  names(dat) <- c("x")
  dat<- left_join(dat, gene_name, by = "x")
  dat <- dat %>% filter(!is.na(Suggested.Symbol))
  mat <- mat[dat$x, ]
  row.names(mat) <- dat$Suggested.Symbol
  obj <- CreateSeuratObject(counts = mat)
  meta <- merge(meta, Anno, by = "Cell_ID")
  row.names(meta) <- meta$Cell_ID
  obj@meta.data <- meta
  return(obj)
}

process7 <- function(obj){
  meta <- obj@meta.data
  mat <- obj@assays$raw@counts
  dat <- as.data.frame(row.names(mat))
  names(dat) <- c("x")
  dat<- left_join(dat, gene_name, by = "x")
  dat <- dat %>% filter(!is.na(Suggested.Symbol))
  mat <- mat[dat$x, ]
  row.names(mat) <- dat$Suggested.Symbol
  obj <- CreateSeuratObject(counts = mat)
  meta <- merge(meta, Anno, by = "Cell_ID")
  row.names(meta) <- meta$Cell_ID
  obj@meta.data <- meta
  return(obj)
}



dat5 <- process0(df302004data15)
dat6 <- process0(df302004data16)
dat7 <- process7(df302004data17)
dat8 <- process0(df302004data18)
dat9 <- process0(df302005data01)
dat10 <- process0(df302005data02)
dat11 <- process0(df302005data03)
dat12 <- process0(df302005data04)

objs <- c(dat5, dat6, dat7, dat8, dat9, dat10, dat11, dat12)



for (i in 1: length(objs)){
  Platelets <- subset(objs[[i]], label %in% c("Platelets", "Cell_precursors", "Erythroblast"))
  saveRDS(Platelets, paste0("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/Platelets/",
                            "dat",i+4,"_Platelets_Cell_precursors_Erythroblast.rds"))
  
  B_cell <- subset(objs[[i]], label %in% c("Platelets", "B_cell"))
  saveRDS(B_cell, paste0("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/B_cell/",
                         "dat",i+4,"_Platelets_B_cell.rds"))
  
  T_cell <- subset(objs[[i]], label %in% c("Platelets", "T_cell"))
  saveRDS(T_cell, paste0("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/T_cell/",
                         "dat",i+4,"_Platelets_T_cell.rds"))
  
  Monocyte <- subset(objs[[i]], label %in% c("Platelets", "Monocyte"))
  saveRDS(Monocyte, paste0("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/Monocyte/",
                           "dat",i+4,"_Platelets_Monocyte.rds"))
  
  NK_cell <- subset(objs[[i]], label %in% c("Platelets", "NK_cell"))
  saveRDS(NK_cell, paste0("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/NK_cell/",
                          "dat",i+4,"_Platelets_NK_cell.rds"))
  
  Neutrophils <- subset(objs[[i]], label %in% c("Platelets", "Neutrophils"))
  saveRDS(Neutrophils, paste0("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/Neutrophils/",
                              "dat",i+4,"_Platelets_Neutrophils.rds"))
  
  DC <- subset(objs[[i]], label %in% c("Platelets", "DC"))
  saveRDS(DC, paste0("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/DC/",
                     "dat",i+4,"_Platelets_DC.rds"))
  
}




