.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )

library(SeuratDisk)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(KEGGREST)
library(HGNChelper)
library(stringr)
library(dplyr)
library(org.Hs.eg.db)
library(GO.db)
library(tidyr)
library(ggpubr)
library(rstatix)

dat0 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/harmony_int/302005_Platelets_cluster_Anno.rds")
DefaultAssay(dat0) <- "RNA"

dat0@meta.data$Disease2 <- ifelse(dat0@meta.data$Disease %in% c('Flu', 'Lung_infect', 'Non_covid'), 
                                  'Non_covid_sepsis', dat0@meta.data$Disease)
dat0@meta.data$Category <- ifelse(dat0@meta.data$Category == "", "SLE", dat0@meta.data$Category)

dat1 <- subset(dat0, label == "Platelets")
dat2 <- subset(dat0, label == "Erythroblast")
dat3 <- subset(dat0, label == "Cell_precursors")

df <- dat1
# module scores -----------------------------------------------------------
#COVID-19 module
covid19 <- as.data.frame(keggGet("hsa05171")[[1]]$GENE)
names(covid19) <- c("gene")
covid19 <- covid19 %>% 
  dplyr::filter(grepl(";", covid19$gene))
covid19$gene <- str_split_fixed(covid19$gene, ";", 2)[,1]
covid19v2 <- checkGeneSymbols(covid19$gene, species="human")
covid19v2 <- covid19v2 %>% 
  dplyr::filter(!is.na(covid19v2$Suggested.Symbol)) %>% 
  dplyr::distinct(Suggested.Symbol, .keep_all = T)
covid19_list <- list(as.list(covid19v2$Suggested.Symbol))
df <- AddModuleScore(df, 
                     features = covid19_list,
                     name = 'COVID-19',
                     search = TRUE)

addmodule <- function(i){
  go_id = GOID(GOTERM[Term(GOTERM) == i])
  allegs = get(go_id, org.Hs.egGO2ALLEGS)
  genes = as.data.frame(unlist(mget(allegs,org.Hs.egSYMBOL)))
  names(genes) <- c("gene")
  genev2 <- checkGeneSymbols(genes$gene, species="human")
  genev2$Suggested.Symbol <- str_split_fixed(genev2$Suggested.Symbol, "///", 2)[,1]
  genev2$Suggested.Symbol <- trimws(genev2$Suggested.Symbol)
  genev2_list <- list(as.list(genev2$Suggested.Symbol))
  AddModuleScore(df,
                 features = genev2_list,
                 name = i,
                 search = TRUE)
  
}

df <- addmodule("platelet activation")
df <- addmodule("coagulation")
df <- addmodule("cytokine activity")
df <- addmodule("acute inflammatory response")
df <- addmodule("positive regulation of hemostasis")
df <- addmodule("response to interferon-beta")
df <- addmodule("response to interferon-gamma")
df <- addmodule("response to type I interferon")
df <- addmodule("translational initiation")
df <- addmodule("response to oxidative stress")

OXPHOS <- list(c(
  "NDUFA10","NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7",
  "NDUFA8","NDUFA9","NDUFAB1","NDUFB10","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5",
  "NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFC1","NDUFC2","NDUFS4","NDUFS5","NDUFS6","NDUFV3",
  "NDUFAF1","NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2","SDHC","UQCR10",
  "UQCRB","UQCRC1","UQCRC2","UQCRH","COX10","COX15","COX4I1","COX4I2","COX5A","COX8A","COX8C",
  "ATP5C1","ATP5D","ATP7A","C1orf177","CHCHD10","COQ7","COQ9","CYCS","DLD","DNAJC15","FXN","GADD45GIP1",
  "GBAS","MECP2","MLXIPL","MSH2","MYOG","PARK7","PINK1","PMPCB","PPIF","SDHAF2","SLC25A23","SLC25A33",
  "SNCA","SURF1","TAZ","UQCRHL","VCP","ACTN3","AK2","APOC3"
))

Glycolysis <- list(c(
  "ENO1","ENO2","ENO3","ENTPD5","GAPDHS","GAPDH","GCK","GPD1","GPI","HIF1A","HK1",
  "HK2","HK3","HTR2A","IGF1","INSR","MYC","P2RX7","PFKFB1","PFKFB2","PFKFB3","PFKFB4",
  "PFKL","PFKM","PFKP","PGAM1","PGAM2","PGK1","PKLR","PPP2R5D","PRKAA1","PRKAA2","TPI1",
  "INS","ALDOA","ALDOB","ALDOC","ARNT"
))

MHC.ClassII <- list(c("CD74","HLA-DRA", "HLA-DQA1","HLA-DQA2","HLA-DPA1","HLA-DRB1","HLA-DPB1",
                      "HLA-DQB2","HLA-DRB5","HLA-DQB1","HLA-DMA","HLA-DMB"))
MHC.ClassI <- list(c("HLA-A", "HLA-B","HLA-C","HLA-E","HLA-F", "HLA-G"))
liteCoagulation <- list(c("F10",  "F11",  "F12",  "F13A1",  "F13B",  "F2",  "F5",  "F7",  "F8", 
                          "F9",  "FGA",  "FGB",  "FGG",  "GGCX",  "GP1BA",  "KLKB1",  "KNG1",  "LMAN1",  
                          "MCFD2",  "PLG",  "SERPINE1",  "SERPINF2",  "VKORC1",  "VWF"))

REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION <- list(c("A1BG", "A2M", "AAMP", "ABCC4", "ABHD12", "ABHD6", "ACTN1", "ACTN2", "ACTN4", 
                                                                 "ADRA2A", "ADRA2B", "ADRA2C", "AHSG", "AKT1", "ALB", "ALDOA", "ANXA5", "APBB1IP", 
                                                                 "APLP2", "APOA1", "APOH", "APOOL", "APP", "ARRB1", "ARRB2", "BCAR1", "BRPF3", "CALM1", 
                                                                 "CALU", "CAP1", "CD109", "CD36", "CD63", "CD9", "CDC37L1", "CDC42", "CFD", "CFL1", "CHID1", 
                                                                 "CLEC1B", "CLEC3B", "CLU", "COL1A1", "COL1A2", "CRK", "CSK", "CTSW", "CYB5R1", "CYRIB", 
                                                                 "DAGLA", "DAGLB", "DGKA", "DGKB", "DGKD", "DGKE", "DGKG", "DGKH", "DGKI", "DGKK", "DGKQ", 
                                                                 "DGKZ", "ECM1", "EGF", "ENDOD1", "F13A1", "F2", "F2R", "F2RL2", "F2RL3", "F5", "F8", "FAM3C", 
                                                                 "FCER1G", "FERMT3", "FGA", "FGB", "FGG", "FLNA", "FN1", "FYN", "GAS6", "GNA11", "GNA12", 
                                                                 "GNA13", "GNA14", "GNA15", "GNAI1", "GNAI2", "GNAI3", "GNAQ", "GNAT3", "GNB1", "GNB2", 
                                                                 "GNB3", "GNB4", "GNB5", "GNG10", "GNG11", "GNG12", "GNG13", "GNG2", "GNG3", "GNG4", "GNG5", 
                                                                 "GNG7", "GNG8", "GNGT1", "GNGT2", "GP1BA", "GP1BB", "GP5", "GP6", "GP9", "GRB2", "GTPBP2", 
                                                                 "HABP4", "HGF", "HRG", "HSPA5", "IGF1", "IGF2", "ISLR", "ITGA2B", "ITGB3", "ITIH3", "ITIH4", 
                                                                 "ITPR1", "ITPR2", "ITPR3", "KNG1", "LAMP2", "LAT", "LCK", "LCP2", "LEFTY2", "LGALS3BP", 
                                                                 "LHFPL2", "LY6G6F", "LYN", "MAGED2", "MANF", "MAPK1", "MAPK14", "MAPK3", "MGLL", "MMRN1", 
                                                                 "MPIG6B", "MPL", "NHLRC2", "OLA1", "ORM1", "ORM2", "P2RY1", "P2RY12", "PCDH7", "PCYOX1L", 
                                                                 "PDGFA", "PDGFB", "PDPK1", "PDPN", "PECAM1", "PF4", "PFN1", "PHACTR2", "PIK3CA", "PIK3CB", 
                                                                 "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R5", "PIK3R6", "PLA2G4A", "PLCG2", "PLEK", "PLG", 
                                                                 "PPBP", "PPIA", "PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCG", "PRKCH", "PRKCQ", "PRKCZ", "PROS1", 
                                                                 "PSAP", "PTK2", "PTPN1", "PTPN11", "PTPN6", "QSOX1", "RAB27B", "RAC1", "RAC2", "RAF1", "RAP1A", 
                                                                 "RAP1B", "RAPGEF3", "RAPGEF4", "RARRES2", "RASGRP1", "RASGRP2", "RHOA", "RHOB", "RHOG", "SCCPDH", 
                                                                 "SCG3", "SELENOP", "SELP", "SERPINA1", "SERPINA3", "SERPINA4", "SERPINE1", "SERPINF2", 
                                                                 "SERPING1", "SHC1", "SOD1", "SOS1", "SPARC", "SPP2", "SRC", "SRGN", "STX4", "STXBP2", "STXBP3", 
                                                                 "SYK", "SYTL4", "TAGLN2", "TBXA2R", "TEX264", "TF", "TGFB1", "TGFB2", "TGFB3", "THBS1", "THPO", 
                                                                 "TIMP1", "TIMP3", "TLN1", "TMSB4X", "TMX3", "TOR4A", "TRPC3", "TRPC6", "TRPC7", "TTN", "TUBA4A", 
                                                                 "VAV1", "VAV2", "VAV3", "VCL", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "VTI1B", "VWF", "WDR1", "YWHAZ"))

df <- AddModuleScore(df,features = OXPHOS,name = "OXPHOS", search = TRUE)
df <- AddModuleScore(df,features = Glycolysis,name = "Glycolysis", search = TRUE)
df <- AddModuleScore(df,features = MHC.ClassII,name = "MHC.ClassII", search = TRUE)
df <- AddModuleScore(df,features = MHC.ClassI,name = "MHC.ClassI", search = TRUE)
df <- AddModuleScore(df,features = liteCoagulation,name = "liteCoagulation", search = TRUE)
df <- AddModuleScore(df,features = REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION,
                     name = "REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION", 
                     search = TRUE)


head(df@meta.data$COVID.191)
head(df@meta.data$platelet.activation1)
head(df@meta.data$coagulation1)
head(df@meta.data$cytokine.activity1)
head(df@meta.data$acute.inflammatory.response1)
head(df@meta.data$positive.regulation.of.hemostasis1)
head(df@meta.data$response.to.interferon.beta1)
head(df@meta.data$response.to.interferon.gamma1)
head(df@meta.data$response.to.type.I.interferon1)
head(df@meta.data$translational.initiation1)
head(df@meta.data$response.to.oxidative.stress1)
head(df@meta.data$OXPHOS1)
head(df@meta.data$Glycolysis1)
head(df@meta.data$MHC.ClassI1)
head(df@meta.data$MHC.ClassII1)
head(df@meta.data$liteCoagulation1)
head(df@meta.data$REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION1)


saveRDS(df@meta.data, "/bigdata/godziklab/shared/Xinru/302005/v5_output/Module_score/Module_meta.rds")
# ggplot figures ----------------------------------------------------------
df_viz <- readRDS('/bigdata/godziklab/shared/Xinru/302005/v5_output/Module_score/Module_meta.rds')
df_viz$Category <- ifelse(df_viz$Disease == "SLE", "SLE", df_viz$Category)
# CV    FT    HC    MD    ML   SLE    SV
# 4268 12969  4381  6355 10491   603 33055
df_viz$Category <- factor(df_viz$Category,
                          levels = c("HC","CV", "ML", "MD", "SV", "FT", "SLE"))
df_viz$Disease2 <- factor(df_viz$Disease2,
                          levels = c("HC","Non_covid_sepsis","COVID-19","Sepsis","SLE"))
df_viz$Death <- factor(df_viz$Death,levels = c("Unknown", "S", "NS"))

names(df_viz) <- gsub("1$", "", names(df_viz))
modules <- c("COVID.19", "platelet.activation", "coagulation", "cytokine.activity", "acute.inflammatory.response", 
             "positive.regulation.of.hemostasis", "response.to.interferon.beta", "response.to.interferon.gamma", 
             "response.to.type.I.interferon", "translational.initiation", "response.to.oxidative.stress", "OXPHOS", 
             "Glycolysis", "MHC.ClassII", "MHC.ClassI", "liteCoagulation", 
             "REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION")

setwd("/bigdata/godziklab/shared/Xinru/302005/v5_output/Module_score/")
for (i in 1: length(modules)){
  module_name = modules[i]
  print(module_name)
  df1 <-df_viz %>% 
    dplyr::select(Category,modules[i]) %>% 
    dplyr::filter(!is.na(Category))
  df1[,2]<- as.numeric(df1[,2])
  names(df1) <- c("Category", "Module")
  stat.test <- df1 %>% wilcox_test(Module ~ Category)
  stat.test <- stat.test %>% 
    add_xy_position(fun = "max",x = "Category") %>% 
    filter(p.adj.signif != "ns")
  dpi = 300
  png(file = paste0(module_name,"_Category_module.png"), 
      width = dpi * 4,height = dpi * 4,
      units = "px",res = dpi,type = 'cairo')
  print(ggbarplot(df1, 
                  x = "Category", 
                  y =  "Module", 
                  fill = "Category", 
                  palette = c('#fee5d9','#fcbba1','#fc9272','#fb6a4a','#de2d26','#a50f15', '#bdbdbd'),
                  add = c("mean_se")) +
          # stat_pvalue_manual(stat.test,
          #                    y.position = 0.5,
          #                    step.increase = 0.05,
          #                    tip.length = 0.01,
          #                    label = "p.adj.signif") +
          theme_minimal() + 
          ggtitle(module_name) +
          theme(legend.position = "none"))
  dev.off()
}


# Disease -----------------------------------------------------------------

# COVID-19         Flu          HC      HC_LPS Lung_infect   Non_covid      Sepsis         SLE
# 33311          33        2450         636        1736         705        3172         136

setwd("/bigdata/godziklab/shared/Xinru/302005/v5_output/Module_score/")
for (i in 1: length(modules)){
  module_name = modules[i]
  print(module_name)
  df1 <-df_viz %>% 
    dplyr::select(Disease2,modules[i]) %>% 
    dplyr::filter(!is.na(Disease2))
  df1[,2]<- as.numeric(df1[,2])
  names(df1) <- c("Disease2", "Module")
  stat.test <- df1 %>% wilcox_test(Module ~ Disease2)
  stat.test <- stat.test %>% 
    add_xy_position(fun = "max",x = "Disease2") %>% 
    filter(p.adj.signif != "ns")
  dpi = 300
  png(file = paste0(module_name,"_Disease2_module.png"), 
      width = dpi * 4,height = dpi * 4,
      units = "px",res = dpi,type = 'cairo')
  print(ggbarplot(df1, 
                  x = "Disease2", 
                  y =  "Module", 
                  fill = "Disease2", 
                  palette = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c'),
                  add = c("mean_se")) +
          # stat_pvalue_manual(stat.test,
          #                    y.position = 0.5,
          #                    step.increase = 0.05,
          #                    tip.length = 0.01,
          #                    label = "p.adj.signif") +
          theme_minimal() + 
          ggtitle(module_name) +
          theme(legend.position = "none")+
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  ) 
  dev.off()
}


# Death -------------------------------------------------------------------

setwd("/bigdata/godziklab/shared/Xinru/302005/v5_output/Module_score/")
for (i in 1: length(modules)){
  module_name = modules[i]
  print(module_name)
  df1 <-df_viz %>% 
    dplyr::select(Death,modules[i]) %>% 
    dplyr::filter(!is.na(Death))
  df1[,2]<- as.numeric(df1[,2])
  names(df1) <- c("Death", "Module")
  stat.test <- df1 %>% wilcox_test(Module ~ Death)
  stat.test <- stat.test %>% 
    add_xy_position(fun = "max",x = "Death") %>% 
    filter(p.adj.signif != "ns")
  dpi = 300
  png(file = paste0(module_name,"_Death_module.png"), 
      width = dpi * 4,height = dpi * 4,
      units = "px",res = dpi,type = 'cairo')
  print(ggbarplot(df1, 
                  x = "Death", 
                  y =  "Module", 
                  fill = "Death", 
                  palette = c('#ffffbf','#99d594','#fc8d59'),
                  add = c("mean_se")) +
          # stat_pvalue_manual(stat.test,
          #                    y.position = 0.5,
          #                    step.increase = 0.05,
          #                    tip.length = 0.01,
          #                    label = "p.adj.signif") +
          theme_minimal() + 
          ggtitle(module_name) + 
          theme(legend.position = "none")+
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  )
  dev.off()
}


# Data_NO -----------------------------------------------------------------
df_viz$Data_NO <- factor(df_viz$Data_NO,levels = c("302004data01","302004data02","302004data04","302004data05","302004data15","302004data16","302004data17","302004data18","302005data01","302005data02","302005data03","302005data04"))
# 302004data01 302004data02 302004data04 302004data05 302004data15 302004data16 302004data17 302004data18 302005data01 302005data02 302005data03 302005data04
# 188         1028          744         4401         7459         1873        11397          452         5330         6334         2837          136

setwd("/bigdata/godziklab/shared/Xinru/302005/v5_output/Module_score/")
for (i in 1: length(modules)){
  module_name = modules[i]
  print(module_name)
  df1 <-df_viz %>% 
    dplyr::select(Data_NO,modules[i]) %>% 
    dplyr::filter(!is.na(Data_NO))
  df1[,2]<- as.numeric(df1[,2])
  names(df1) <- c("Data_NO", "Module")
  stat.test <- df1 %>% wilcox_test(Module ~ Data_NO)
  stat.test <- stat.test %>% 
    add_xy_position(fun = "max",x = "Data_NO") %>% 
    filter(p.adj.signif != "ns")
  dpi = 300
  png(file = paste0(module_name,"_Data_NO_module.png"), 
      width = dpi * 4,height = dpi * 8,
      units = "px",res = dpi,type = 'cairo')
  print(ggbarplot(df1, 
                  x = "Data_NO", 
                  y =  "Module", 
                  fill = "Data_NO", 
                  palette = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'),
                  add = c("mean_se")) +
          # stat_pvalue_manual(stat.test,
          #                    y.position = 0.5,
          #                    step.increase = 0.05,
          #                    tip.length = 0.01,
          #                    label = "p.adj.signif") +
          theme_minimal() + 
          ggtitle(module_name) + 
          theme(legend.position = "none")+
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  )
  dev.off()
}


