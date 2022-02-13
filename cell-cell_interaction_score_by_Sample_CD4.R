.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )
library(Seurat)
library(dplyr)
library(cowplot)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(stringr)
library(iTALK)

data0 <- readRDS('/bigdata/godziklab/shared/Xinru/302005/datasets_v2/302005_platelet_CD4T_seurat_integrated.rds')


# add metadata ------------------------------------------------------------
Anno <- read.table("/bigdata/godziklab/shared/Xinru/302005/datasets/cell_annotation.tsv",
                   sep = "\t",
                   header = T,
                   stringsAsFactors = F
)

cell_meta <- data0@meta.data
cell_meta <- cell_meta %>% select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt,
                                  percent.rps, percent.rpl, percent.rrna, nCount_SCT,
                                  nFeature_SCT, pruned.labels)
cell_meta$Cell_ID <- row.names(cell_meta)

cell_meta <-  merge(cell_meta, Anno, by = "Cell_ID")
row.names(cell_meta)<- cell_meta$Cell_ID
dim(cell_meta)
data0@meta.data <- cell_meta
data1 <- subset(data0, pbmc == "Y")
DefaultAssay(data1) <- "RNA"
data1 <- subset(data1, Data_NO != "302005data01")
data1 <- subset(data1, pruned.labels %in%  c("Platelets", "T_cells"))

# Matrix Function ---------------------------------------------------------

as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

# Extract matrix ----------------------------------------------------------

df <- as.data.frame(table(data1@meta.data$Sample_ID))
samples <- factor(as.character(unique(df$Var1)))

for (i in 1:length(samples)){
  pt0 <- subset(data1, Sample_ID == samples[i])
  cell_type <- pt0@meta.data %>% dplyr::select(pruned.labels)
  names(cell_type) <- c("cell_type")
  New_matrix <- t(as.data.frame(as_matrix(pt0@assays[["RNA"]]@data)))
  expression_matrix <- merge(cell_type,New_matrix, by= 0)
  ntop= 100
  data <- as.data.frame(expression_matrix)
  data$cell_type <- as.character(data$cell_type)
  highly_exprs_genes <-rawParse(data,top_genes=ntop,stats='mean')
  
  res_cytokine <-FindLR(highly_exprs_genes,datatype='mean count',comm_type= "cytokine")
  res_checkpoint <-FindLR(highly_exprs_genes,datatype='mean count',comm_type= "checkpoint")
  res_growth_factor <-FindLR(highly_exprs_genes,datatype='mean count',comm_type= "growth factor")
  res_other <-FindLR(highly_exprs_genes,datatype='mean count',comm_type= "other")
  res_cat <- rbind(res_cytokine, res_checkpoint, res_growth_factor, res_other)
  
  platelet_res <- res_cat %>%
    dplyr::filter(cell_from == "Platelets") %>%
    dplyr::filter(cell_from_mean_exprs != 0) %>%
    dplyr::filter(cell_to_mean_exprs != 0)
  res_platelet <- res_cat %>%
    dplyr::filter(cell_to == "Platelets") %>%
    dplyr::filter(cell_from != "Platelets") %>%
    dplyr::filter(cell_from_mean_exprs != 0) %>%
    dplyr::filter(cell_to_mean_exprs != 0)
  
  platelet_res$LR_Pair <- paste(platelet_res$ligand, platelet_res$receptor, sep = " - ")
  platelet_res$CellType_Pair <- paste(platelet_res$cell_from, platelet_res$cell_to, sep = " - ")
  platelet_res$Score <- platelet_res$cell_from_mean_exprs * platelet_res$cell_to_mean_exprs
  platelet_res <- platelet_res %>% dplyr::select(CellType_Pair, LR_Pair, Score)
  
  res_platelet$LR_Pair <- paste(res_platelet$ligand, res_platelet$receptor, sep = " - ")
  res_platelet$CellType_Pair <- paste(res_platelet$cell_to, res_platelet$cell_from, sep = " - ")
  res_platelet$Score <- res_platelet$cell_from_mean_exprs * res_platelet$cell_to_mean_exprs
  res_platelet <- res_platelet %>% dplyr::select(CellType_Pair, LR_Pair, Score)
  
  platelet_inter <- rbind(platelet_res, res_platelet)
  if (dim(platelet_inter)[1] > 0){
    platelet_inter$Sample_ID <- samples[i]
    write.table(
      platelet_inter,
      file = paste0("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score_sample_v2/CD4T/", 
                    samples[i],'_platelet_CD4T_interaction_score.txt'),
      sep='\t',
      quote = FALSE,
      row.names = F
    )
    data2 <- platelet_inter
    data2$Score <- as.numeric(data2$Score)
    
    data2 <- data2 %>%
      group_by(CellType_Pair) %>%
      summarize(Score = mean(Score, na.rm = TRUE))
    data2$Sample_ID <- samples[i]
    data2$type <- c("AVG")
    write.table(
      data2,
      file = paste0("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score_sample_v2/CD4T/", 
                    samples[i],'_platelet_CD4T_interaction_AVG.txt'),
      sep='\t',
      quote = FALSE,
      row.names = F
    )
    data3 <- platelet_inter
    data3$Score <- as.numeric(data3$Score)
    data3 <- data3 %>%
      group_by(CellType_Pair) %>%
      summarize(Score = sum(Score, na.rm = TRUE))
    data3$Sample_ID <- samples[i]
    data3$type <- c("SUM")
    write.table(
      data3,
      file = paste0("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score_sample_v2/CD4T/", 
                    samples[i],'_platelet_CD4T_interaction_SUM.txt'),
      sep='\t',
      quote = FALSE,
      row.names = F
    )
  }
}


