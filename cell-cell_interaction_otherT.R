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

data0 <- readRDS('/bigdata/godziklab/shared/Xinru/302005/datasets_v2/302005_platelet_otherT_seurat_integrated.rds')
head(data0@meta.data)

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


# Extract matrix FT----------------------------------------------------------
#        CV     FT     HC HC_LPS     MD     ML     SV
# 597   4831  14042   4563    656   4940  11503  43502

Category_select <- c("FT")
pt <- subset(data1, Category %in% Category_select)
cell_type <- pt@meta.data %>% dplyr::select(pruned.labels)
head(cell_type)
names(cell_type) <- c("cell_type")
table(cell_type$cell_type)
New_matrix <- t(as.data.frame(pt@assays[["RNA"]]@data))
New_matrix[1:15, 1:15]
expression_matrix <- merge(cell_type,New_matrix, by= 0)
expression_matrix[1:15, 1:15]

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

dim(platelet_res)
platelet_res$LR_Pair <- paste(platelet_res$ligand, platelet_res$receptor, sep = " - ")
platelet_res$CellType_Pair <- paste(platelet_res$cell_from, platelet_res$cell_to, sep = " - ")
platelet_res$Score <- platelet_res$cell_from_mean_exprs * platelet_res$cell_to_mean_exprs
max(platelet_res$Score)
platelet_res <- platelet_res %>% dplyr::select(CellType_Pair, LR_Pair, Score)

dim(res_platelet)
res_platelet$LR_Pair <- paste(res_platelet$ligand, res_platelet$receptor, sep = " - ")
res_platelet$CellType_Pair <- paste(res_platelet$cell_to, res_platelet$cell_from, sep = " - ")
res_platelet$Score <- res_platelet$cell_from_mean_exprs * res_platelet$cell_to_mean_exprs
max(res_platelet$Score)
res_platelet <- res_platelet %>% dplyr::select(CellType_Pair, LR_Pair, Score)

platelet_inter <- rbind(platelet_res, res_platelet)
data2 <- platelet_inter
data2$Score <- as.numeric(data2$Score)
data2 <- data2 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = mean(Score, na.rm = TRUE))
names(data2)[2] <- Category_select
platelet_inter$Category <- Category_select

platelet_inter <- rbind(platelet_res, res_platelet)
data3 <- platelet_inter
data3$Score <- as.numeric(data3$Score)
data3 <- data3 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = sum(Score, na.rm = TRUE))
names(data3)[2] <- Category_select

data2$type <- c("AVG")
data3$type <- c("SUM")

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data2,
  file = paste0(Category_select,'_platelet_otherT_interaction_AVG.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data3,
  file = paste0(Category_select,'_platelet_otherT_interaction_SUM.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)
# Extract matrix CV----------------------------------------------------------
#        CV     FT     HC HC_LPS     MD     ML     SV
# 597   4831  14042   4563    656   4940  11503  43502

Category_select <- c("CV")
pt <- subset(data1, Category %in% Category_select)
cell_type <- pt@meta.data %>% dplyr::select(pruned.labels)
head(cell_type)
names(cell_type) <- c("cell_type")
table(cell_type$cell_type)
New_matrix <- t(as.data.frame(pt@assays[["RNA"]]@data))
New_matrix[1:15, 1:15]
expression_matrix <- merge(cell_type,New_matrix, by= 0)
expression_matrix[1:15, 1:15]

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

dim(platelet_res)
platelet_res$LR_Pair <- paste(platelet_res$ligand, platelet_res$receptor, sep = " - ")
platelet_res$CellType_Pair <- paste(platelet_res$cell_from, platelet_res$cell_to, sep = " - ")
platelet_res$Score <- platelet_res$cell_from_mean_exprs * platelet_res$cell_to_mean_exprs
max(platelet_res$Score)
platelet_res <- platelet_res %>% dplyr::select(CellType_Pair, LR_Pair, Score)

dim(res_platelet)
res_platelet$LR_Pair <- paste(res_platelet$ligand, res_platelet$receptor, sep = " - ")
res_platelet$CellType_Pair <- paste(res_platelet$cell_to, res_platelet$cell_from, sep = " - ")
res_platelet$Score <- res_platelet$cell_from_mean_exprs * res_platelet$cell_to_mean_exprs
max(res_platelet$Score)
res_platelet <- res_platelet %>% dplyr::select(CellType_Pair, LR_Pair, Score)

platelet_inter <- rbind(platelet_res, res_platelet)
data2 <- platelet_inter
data2$Score <- as.numeric(data2$Score)
data2 <- data2 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = mean(Score, na.rm = TRUE))
names(data2)[2] <- Category_select
platelet_inter$Category <- Category_select

platelet_inter <- rbind(platelet_res, res_platelet)
data3 <- platelet_inter
data3$Score <- as.numeric(data3$Score)
data3 <- data3 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = sum(Score, na.rm = TRUE))
names(data3)[2] <- Category_select

data2$type <- c("AVG")
data3$type <- c("SUM")

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data2,
  file = paste0(Category_select,'_platelet_otherT_interaction_AVG.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data3,
  file = paste0(Category_select,'_platelet_otherT_interaction_SUM.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)
# Extract matrix HC----------------------------------------------------------
#        CV     FT     HC HC_LPS     MD     ML     SV
# 597   4831  14042   4563    656   4940  11503  43502

Category_select <- c("HC")
pt <- subset(data1, Category %in% Category_select)
cell_type <- pt@meta.data %>% dplyr::select(pruned.labels)
head(cell_type)
names(cell_type) <- c("cell_type")
table(cell_type$cell_type)
New_matrix <- t(as.data.frame(pt@assays[["RNA"]]@data))
New_matrix[1:15, 1:15]
expression_matrix <- merge(cell_type,New_matrix, by= 0)
expression_matrix[1:15, 1:15]

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

dim(platelet_res)
platelet_res$LR_Pair <- paste(platelet_res$ligand, platelet_res$receptor, sep = " - ")
platelet_res$CellType_Pair <- paste(platelet_res$cell_from, platelet_res$cell_to, sep = " - ")
platelet_res$Score <- platelet_res$cell_from_mean_exprs * platelet_res$cell_to_mean_exprs
max(platelet_res$Score)
platelet_res <- platelet_res %>% dplyr::select(CellType_Pair, LR_Pair, Score)

dim(res_platelet)
res_platelet$LR_Pair <- paste(res_platelet$ligand, res_platelet$receptor, sep = " - ")
res_platelet$CellType_Pair <- paste(res_platelet$cell_to, res_platelet$cell_from, sep = " - ")
res_platelet$Score <- res_platelet$cell_from_mean_exprs * res_platelet$cell_to_mean_exprs
max(res_platelet$Score)
res_platelet <- res_platelet %>% dplyr::select(CellType_Pair, LR_Pair, Score)

platelet_inter <- rbind(platelet_res, res_platelet)
data2 <- platelet_inter
data2$Score <- as.numeric(data2$Score)
data2 <- data2 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = mean(Score, na.rm = TRUE))
names(data2)[2] <- Category_select
platelet_inter$Category <- Category_select

platelet_inter <- rbind(platelet_res, res_platelet)
data3 <- platelet_inter
data3$Score <- as.numeric(data3$Score)
data3 <- data3 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = sum(Score, na.rm = TRUE))
names(data3)[2] <- Category_select

data2$type <- c("AVG")
data3$type <- c("SUM")

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data2,
  file = paste0(Category_select,'_platelet_otherT_interaction_AVG.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data3,
  file = paste0(Category_select,'_platelet_otherT_interaction_SUM.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)
# Extract matrix HC_LPS----------------------------------------------------------
#        CV     FT     HC HC_LPS     MD     ML     SV
# 597   4831  14042   4563    656   4940  11503  43502

Category_select <- c("HC_LPS")
pt <- subset(data1, Category %in% Category_select)
cell_type <- pt@meta.data %>% dplyr::select(pruned.labels)
head(cell_type)
names(cell_type) <- c("cell_type")
table(cell_type$cell_type)
New_matrix <- t(as.data.frame(pt@assays[["RNA"]]@data))
New_matrix[1:15, 1:15]
expression_matrix <- merge(cell_type,New_matrix, by= 0)
expression_matrix[1:15, 1:15]

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

dim(platelet_res)
platelet_res$LR_Pair <- paste(platelet_res$ligand, platelet_res$receptor, sep = " - ")
platelet_res$CellType_Pair <- paste(platelet_res$cell_from, platelet_res$cell_to, sep = " - ")
platelet_res$Score <- platelet_res$cell_from_mean_exprs * platelet_res$cell_to_mean_exprs
max(platelet_res$Score)
platelet_res <- platelet_res %>% dplyr::select(CellType_Pair, LR_Pair, Score)

dim(res_platelet)
res_platelet$LR_Pair <- paste(res_platelet$ligand, res_platelet$receptor, sep = " - ")
res_platelet$CellType_Pair <- paste(res_platelet$cell_to, res_platelet$cell_from, sep = " - ")
res_platelet$Score <- res_platelet$cell_from_mean_exprs * res_platelet$cell_to_mean_exprs
max(res_platelet$Score)
res_platelet <- res_platelet %>% dplyr::select(CellType_Pair, LR_Pair, Score)

platelet_inter <- rbind(platelet_res, res_platelet)
data2 <- platelet_inter
data2$Score <- as.numeric(data2$Score)
data2 <- data2 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = mean(Score, na.rm = TRUE))
names(data2)[2] <- Category_select
platelet_inter$Category <- Category_select

platelet_inter <- rbind(platelet_res, res_platelet)
data3 <- platelet_inter
data3$Score <- as.numeric(data3$Score)
data3 <- data3 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = sum(Score, na.rm = TRUE))
names(data3)[2] <- Category_select

data2$type <- c("AVG")
data3$type <- c("SUM")

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data2,
  file = paste0(Category_select,'_platelet_otherT_interaction_AVG.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data3,
  file = paste0(Category_select,'_platelet_otherT_interaction_SUM.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)
# Extract matrix MD----------------------------------------------------------
#        CV     FT     HC HC_LPS     MD     ML     SV
# 597   4831  14042   4563    656   4940  11503  43502

Category_select <- c("MD")
pt <- subset(data1, Category %in% Category_select)
cell_type <- pt@meta.data %>% dplyr::select(pruned.labels)
head(cell_type)
names(cell_type) <- c("cell_type")
table(cell_type$cell_type)
New_matrix <- t(as.data.frame(pt@assays[["RNA"]]@data))
New_matrix[1:15, 1:15]
expression_matrix <- merge(cell_type,New_matrix, by= 0)
expression_matrix[1:15, 1:15]

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

dim(platelet_res)
platelet_res$LR_Pair <- paste(platelet_res$ligand, platelet_res$receptor, sep = " - ")
platelet_res$CellType_Pair <- paste(platelet_res$cell_from, platelet_res$cell_to, sep = " - ")
platelet_res$Score <- platelet_res$cell_from_mean_exprs * platelet_res$cell_to_mean_exprs
max(platelet_res$Score)
platelet_res <- platelet_res %>% dplyr::select(CellType_Pair, LR_Pair, Score)

dim(res_platelet)
res_platelet$LR_Pair <- paste(res_platelet$ligand, res_platelet$receptor, sep = " - ")
res_platelet$CellType_Pair <- paste(res_platelet$cell_to, res_platelet$cell_from, sep = " - ")
res_platelet$Score <- res_platelet$cell_from_mean_exprs * res_platelet$cell_to_mean_exprs
max(res_platelet$Score)
res_platelet <- res_platelet %>% dplyr::select(CellType_Pair, LR_Pair, Score)

platelet_inter <- rbind(platelet_res, res_platelet)
data2 <- platelet_inter
data2$Score <- as.numeric(data2$Score)
data2 <- data2 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = mean(Score, na.rm = TRUE))
names(data2)[2] <- Category_select
platelet_inter$Category <- Category_select

platelet_inter <- rbind(platelet_res, res_platelet)
data3 <- platelet_inter
data3$Score <- as.numeric(data3$Score)
data3 <- data3 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = sum(Score, na.rm = TRUE))
names(data3)[2] <- Category_select

data2$type <- c("AVG")
data3$type <- c("SUM")

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data2,
  file = paste0(Category_select,'_platelet_otherT_interaction_AVG.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data3,
  file = paste0(Category_select,'_platelet_otherT_interaction_SUM.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)
# Extract matrix ML----------------------------------------------------------
#        CV     FT     HC HC_LPS     MD     ML     SV
# 597   4831  14042   4563    656   4940  11503  43502

Category_select <- c("ML")
pt <- subset(data1, Category %in% Category_select)
cell_type <- pt@meta.data %>% dplyr::select(pruned.labels)
head(cell_type)
names(cell_type) <- c("cell_type")
table(cell_type$cell_type)
New_matrix <- t(as.data.frame(pt@assays[["RNA"]]@data))
New_matrix[1:15, 1:15]
expression_matrix <- merge(cell_type,New_matrix, by= 0)
expression_matrix[1:15, 1:15]

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

dim(platelet_res)
platelet_res$LR_Pair <- paste(platelet_res$ligand, platelet_res$receptor, sep = " - ")
platelet_res$CellType_Pair <- paste(platelet_res$cell_from, platelet_res$cell_to, sep = " - ")
platelet_res$Score <- platelet_res$cell_from_mean_exprs * platelet_res$cell_to_mean_exprs
max(platelet_res$Score)
platelet_res <- platelet_res %>% dplyr::select(CellType_Pair, LR_Pair, Score)

dim(res_platelet)
res_platelet$LR_Pair <- paste(res_platelet$ligand, res_platelet$receptor, sep = " - ")
res_platelet$CellType_Pair <- paste(res_platelet$cell_to, res_platelet$cell_from, sep = " - ")
res_platelet$Score <- res_platelet$cell_from_mean_exprs * res_platelet$cell_to_mean_exprs
max(res_platelet$Score)
res_platelet <- res_platelet %>% dplyr::select(CellType_Pair, LR_Pair, Score)

platelet_inter <- rbind(platelet_res, res_platelet)
data2 <- platelet_inter
data2$Score <- as.numeric(data2$Score)
data2 <- data2 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = mean(Score, na.rm = TRUE))
names(data2)[2] <- Category_select
platelet_inter$Category <- Category_select

platelet_inter <- rbind(platelet_res, res_platelet)
data3 <- platelet_inter
data3$Score <- as.numeric(data3$Score)
data3 <- data3 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = sum(Score, na.rm = TRUE))
names(data3)[2] <- Category_select

data2$type <- c("AVG")
data3$type <- c("SUM")

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data2,
  file = paste0(Category_select,'_platelet_otherT_interaction_AVG.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data3,
  file = paste0(Category_select,'_platelet_otherT_interaction_SUM.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)
# Extract matrix SV----------------------------------------------------------
#        CV     FT     HC HC_LPS     MD     ML     SV
# 597   4831  14042   4563    656   4940  11503  43502

Category_select <- c("SV")
pt <- subset(data1, Category %in% Category_select)
cell_type <- pt@meta.data %>% dplyr::select(pruned.labels)
head(cell_type)
names(cell_type) <- c("cell_type")
table(cell_type$cell_type)
New_matrix <- t(as.data.frame(pt@assays[["RNA"]]@data))
New_matrix[1:15, 1:15]
expression_matrix <- merge(cell_type,New_matrix, by= 0)
expression_matrix[1:15, 1:15]

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

dim(platelet_res)
platelet_res$LR_Pair <- paste(platelet_res$ligand, platelet_res$receptor, sep = " - ")
platelet_res$CellType_Pair <- paste(platelet_res$cell_from, platelet_res$cell_to, sep = " - ")
platelet_res$Score <- platelet_res$cell_from_mean_exprs * platelet_res$cell_to_mean_exprs
max(platelet_res$Score)
platelet_res <- platelet_res %>% dplyr::select(CellType_Pair, LR_Pair, Score)

dim(res_platelet)
res_platelet$LR_Pair <- paste(res_platelet$ligand, res_platelet$receptor, sep = " - ")
res_platelet$CellType_Pair <- paste(res_platelet$cell_to, res_platelet$cell_from, sep = " - ")
res_platelet$Score <- res_platelet$cell_from_mean_exprs * res_platelet$cell_to_mean_exprs
max(res_platelet$Score)
res_platelet <- res_platelet %>% dplyr::select(CellType_Pair, LR_Pair, Score)

platelet_inter <- rbind(platelet_res, res_platelet)
data2 <- platelet_inter
data2$Score <- as.numeric(data2$Score)
data2 <- data2 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = mean(Score, na.rm = TRUE))
names(data2)[2] <- Category_select
platelet_inter$Category <- Category_select

platelet_inter <- rbind(platelet_res, res_platelet)
data3 <- platelet_inter
data3$Score <- as.numeric(data3$Score)
data3 <- data3 %>%
  group_by(CellType_Pair) %>%
  summarize(Score = sum(Score, na.rm = TRUE))
names(data3)[2] <- Category_select

data2$type <- c("AVG")
data3$type <- c("SUM")

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data2,
  file = paste0(Category_select,'_platelet_otherT_interaction_AVG.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)

setwd("/bigdata/godziklab/shared/Xinru/302005/22-02/Interaction_score")
write.table(
  data3,
  file = paste0(Category_select,'_platelet_otherT_interaction_SUM.txt'),
  sep='\t',
  quote = FALSE,
  row.names = F
)
