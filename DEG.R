.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )

library(dplyr)
library(Seurat)
library(SeuratObject)
library(MAST)

dat0 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/harmony_int/302005_Platelets_cluster_Anno.rds")
DefaultAssay(dat0) <- "RNA"

dat0@meta.data$Disease2 <- ifelse(dat0@meta.data$Disease %in% c('Flu', 'Lung_infect', 'Non_covid'), 
                                  'Non_covid_sepsis', dat0@meta.data$Disease)
dat0@meta.data$Category <- ifelse(dat0@meta.data$Category == "", "SLE", dat0@meta.data$Category)

dat1 <- subset(dat0, label == "Platelets")
dat2 <- subset(dat0, label == "Erythroblast")
dat3 <- subset(dat0, label == "Cell_precursors")

a1 <- c("COVID-19", "Sepsis")
b1 <- c("COVID-19", "HC")
c1 <- c("Non_covid_sepsis", "HC")
d1 <- c("Sepsis", "HC")
e1 <- c("SLE", "HC")
f1 <- c("SLE", "COVID-19")
g1 <- c("SLE", "Non_covid_sepsis")
h1 <- c("SLE", "Sepsis")
i1 <- c("Sepsis", "Non_covid_sepsis")
j1 <- c("COVID-19", "Non_covid_sepsis")
disease_comparison <- list(a1, b1, c1, d1, e1, f1, g1, h1, i1, j1)
a2 <- c("HC", "ML")
b2 <- c("HC", "MD")
c2 <- c("HC", "SV")
d2 <- c("HC", "FT")
e2 <- c("ML", "MD")
f2 <- c("ML", "SV")
g2 <- c("ML", "FT")
h2 <- c("MD", "SV")
i2 <- c("MD", "FT")
j2 <- c("SV", "FT")
Category_comparison <- list(a2, b2, c2, d2, e2, f2, g2, h2, i2, j2)
diseases <- as.factor(unique(dat1@meta.data$Disease2))
categories <- as.factor(unique(dat1@meta.data$Category))


# Cell precursors ---------------------------------------------------------

for (i in 1: length(disease_comparison)){
  df_sub = subset(dat3, Disease2 %in% disease_comparison[[i]])
  Idents(df_sub) = 'Disease2'
  markers <- FindAllMarkers(object = df_sub, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
  setwd("/bigdata/godziklab/shared/Xinru/302005/v5_output/DEG/Cell_precursors")
  write.table(markers,
              file=(paste0(disease_comparison[[i]][1],"_",disease_comparison[[i]][2],'.txt')),quote = FALSE,sep = '\t')
}

for (i in 1: length(Category_comparison)){
  df_sub = subset(dat3, Category %in% Category_comparison[[i]])
  Idents(df_sub) = 'Category'
  markers <- FindAllMarkers(object = df_sub, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
  setwd("/bigdata/godziklab/shared/Xinru/302005/v5_output/DEG/Cell_precursors")
  write.table(markers,
              file=(paste0(Category_comparison[[i]][1],"_",Category_comparison[[i]][2],'.txt')),quote = FALSE,sep = '\t')
}

for (i in 1: length(diseases)){
  df_sub = subset(dat3, Disease2 %in% diseases[[i]])
  for (j in 1: length(Category_comparison)){
    if (Category_comparison[[j]] %in% (as.factor(unique(df_sub@meta.data$Category)))){
      df_sub = subset(df_sub, Category %in% Category_comparison[[j]])
      Idents(df_sub) = 'Category'
      markers <- FindAllMarkers(object = df_sub, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
      setwd("/bigdata/godziklab/shared/Xinru/302005/v5_output/DEG/Cell_precursors")
      write.table(markers,
                  file=(paste0(diseases[[i]],Category_comparison[[j]][1],
                               "_",diseases[[i]],Category_comparison[[j]][2],'.txt')),quote = FALSE,sep = '\t')
    }
  }
}

for (i in 1: length(categories)){
  df_sub = subset(dat3, Category %in% categories[[i]])
  for (j in 1: length(disease_comparison)){
    if (disease_comparison[[j]] %in% (as.factor(unique(df_sub@meta.data$Disease2)))){
      df_sub = subset(df_sub, Disease2 %in% disease_comparison[[j]])
      Idents(df_sub) = 'Disease2'
      markers <- FindAllMarkers(object = df_sub, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
      setwd("/bigdata/godziklab/shared/Xinru/302005/v5_output/DEG/Cell_precursors")
      write.table(markers,
                  file=(paste0(disease_comparison[[j]][1],categories[[i]],
                               "_",disease_comparison[[j]][2],categories[[i]],'.txt')),quote = FALSE,sep = '\t')
    }
  }
}


