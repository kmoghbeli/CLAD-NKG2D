library(tidyverse)
library(Seurat)
library(harmony)
library(ggpubr)
library(ggprism)

options(Seurat.object.assay.version = "v3")

## CUSTOM GGPLOT THEME
GG_KM_THEME <- 
  list(
    theme_prism(palette = "colorblind_safe", 
                base_family = "sans",
                base_size = 16, 
                base_line_size = 1.5),
    ggprism::scale_color_prism(palette = "colorblind_safe", name = ""), 
    ggprism::scale_fill_prism(palette = "colorblind_safe", name = ""), 
    theme(plot.title = element_text(vjust = 0))
  )

datasets_path <- "../../../_DATASETS/"

data_dir <- "../data/"
figures_dir <- "../figures/"


############################################################

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 4)
counts <- readr::read_delim(paste0(datasets_path, "GSE166386_Eickelberg_LTx_Mouse_scRNAseq/GSE166386_count_matrix.tsv.gz")) %>% 
  column_to_rownames("gene")

metadata <- readr::read_delim(paste0(datasets_path, "GSE166386_Eickelberg_LTx_Mouse_scRNAseq/GSE166386_meta-data.tsv")) %>% 
  column_to_rownames("cell")

obj <- CreateSeuratObject(counts = counts, meta.data = metadata, min.cells = 3, min.features = 200)

# Initial QC (mitochondrial genes, low features, doublets)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^MT-|^mt-")

#Used below plots to get approx. nCount_RNA cut-off of ~20k
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))

obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

cat("After QC filtering:", nrow(obj), "genes x", ncol(obj), "cells\n\n")

seurat_objs <- SplitObject(obj, split.by = "sample_id")

for (i in 1:length(seurat_objs)) {
  seurat_objs[[i]] <- SCTransform(seurat_objs[[i]], 
                                  vst.flavor = "v2",
                                  vars.to.regress = "percent.mt",
                                  return.only.var.genes = FALSE,
                                  verbose = TRUE)
}

integ.features <- SelectIntegrationFeatures(object.list = seurat_objs, nfeatures = 3000)

## Merge
merged_seurat <- merge(x = seurat_objs[[1]],
                       y = seurat_objs[2:length(seurat_objs)],
                       merge.data = TRUE)

VariableFeatures(merged_seurat) <- integ.features

merged_seurat <- RunPCA(object = merged_seurat, assay = "SCT", npcs = 50, verbose = TRUE)

merged_seurat <- RunHarmony(object = merged_seurat,
                            assay.use = "SCT",
                            reduction = "pca",
                            dims.use = 1:50, 
                            lambda = NULL,  # Ridge regression penalty - when set to NULL, harmony tries to estimate
                            kmeans_init_nstart=20, kmeans_init_iter_max=100, 
                            group.by.vars = "sample_id", 
                            reduction.save = "harmony", 
                            plot_convergence = TRUE)

merged_seurat <- RunUMAP(object = merged_seurat, assay = "SCT", reduction = "harmony", dims = 1:50)
merged_seurat <- FindNeighbors(object = merged_seurat, assay = "SCT", reduction = "harmony", dims = 1:50, verbose = TRUE)
merged_seurat <- FindClusters(object = merged_seurat, resolution = 1.0)

merged_seurat <- PrepSCTFindMarkers(object = merged_seurat, verbose = TRUE)

merged_seurat %>% readr::write_rds(paste0(data_dir, "oliver_mouse.rds"))

#merged_seurat %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "oliver_mouse.h5Seurat"), overwrite = TRUE, verbose = TRUE)

###########################################################################################

merged_seurat <- readr::read_rds(paste0(data_dir, "oliver_mouse.rds"))

marker_genes <- 
  c(#"Itgam", # Cd11b - Myeloid marker (mouse)
    #"Adgre1", # F4/80, 
    "Pax5",  # B 
    "Ighd", "Cd27", # Naive (mouse) B cell markers (IgD+, CD27-)
    "Cd3d", "Cd3e", "Cd8a", # T
    "Klrb1c", "Prf1", "Klrk1", "Gzma", "Gzmb",  # NK 
    "Itga2", "Ncam1",  #NK-T
    "Cd83",  # DCs
    "Cd14", "Cd68",  # Macs - note that Cd16 never comes up 
    #"Itgax", # DCs
    "Ly6c1", 
    #"Cd74", # MHC-II mouse marker (used by Renthal 2022 to identify immune cells in TG)
    "Ptgs2", "Irf5", "Nos2",  # Mouse M1 Mac Markers 
    # "Stat1", "Retnla",  # Mouse M1 Mac Markers (less helpful)
    #"Il12a", "Il23a", "Cd163",  # M1 vs M2 (M1: IL-12 and IL23 high with CD163 neg and M2 the opposite)
    "Cd163",  # M2
    #"Arg1", # M2a
    "Socs3", "Cd86", # M2b
    "Ccr2", "Slamf6",   #M2c
    # "Tlr1", "Tlr8", "Scarb1", #M2c (less helpful)
    "Vegfa",    # M2d, 
    "Cx3cr1"  # Tissue-res Mac
  )

DimPlot(merged_seurat, group.by = "all_cells_cell_type_cluster", 
        reduction = "umap", label = TRUE, label.box = TRUE, label.size = 3, repel = TRUE, 
        #cols = DiscretePalette(length(unique(combined_tg$seurat_clusters)), palette = "glasbey")
) + NoLegend()

VlnPlot(merged_seurat, features = marker_genes, assay = "SCT", 
        group.by = "seurat_clusters", 
        #cols = viridis::viridis(4),
        stack = TRUE, flip = TRUE
) 

VlnPlot(merged_seurat, features = marker_genes, assay = "SCT", 
        group.by = "seurat_clusters", split.by = "condition", 
        #cols = viridis::viridis(4),
        stack = TRUE, flip = TRUE
) + scale_fill_prism(palette = "colorblind_safe")



cell_averages <- AverageExpression(merged_seurat, assays = "SCT", layer = "scale.data", group.by = "seurat_clusters",
                                   features = c("Cd3d", "Cd3e", "Cd8a", "Klrb1c", "Pax5", "Cd68"))[[1]] %>% 
  as.matrix() %>% t() %>% 
  as_tibble(rownames = "cluster") %>% 
  mutate(designation = case_when(Cd8a > 0 ~ "CD8", 
                                 Cd3d > 0 ~ "CD4", 
                                 #Klrb1c > 1 ~ "NK", 
                                 .default = "other"
                                 #.default = cluster
  ), 
  .after = cluster)


merged_seurat$t_cell_type <- factor(cell_averages$designation[as.numeric(merged_seurat$seurat_clusters)], 
                              levels = c("CD4", "CD8", "other"))

Idents(merged_seurat) <- merged_seurat$t_cell_type
t_cells <- subset(merged_seurat, idents = c("CD4", "CD8"))
cd8_cells <- subset(merged_seurat, idents = c("CD8"))


Idents(cd8_cells) <- cd8_cells$condition
dge.mouse_t <- FindMarkers(cd8_cells, ident.1 = "HLAB6", ident.2 = "B6B6", recorrect_umi = FALSE)

# EnhancedVolcano::EnhancedVolcano(dge.mouse_t, 
#                                  lab = rownames(dge.mouse_t), 
#                                  x = "avg_log2FC", y = "p_val_adj", 
#                                  pCutoff = 1e-40,
#                                  FCcutoff = 1.5, 
#                                  title = "CLAD vs Control",
#                                  subtitle = "Lung",
#                                  legendLabels = c("NS", expression(Log[2] ~ FC), "P-val", expression(P - value ~ and
#                                                                                                      ~ Log[2] ~ FC)),
#                                  legendPosition = "bottom",
#                                  legendLabSize = 10,
#                                  legendIconSize = 3,
#                                  legendDropLevels = TRUE,
#                                  colAlpha = 1, pointSize = 5, 
#                                  drawConnectors = TRUE, widthConnectors = 0.5,
#                                  labFace = "plain", boxedLabels = TRUE)

klrk1_res <- dge.mouse_t %>% rownames_to_column("gene") %>% filter("Klrk1" == gene)
klrk1_logfc <- signif(klrk1_res$avg_log2FC, digits = 3)
klrk1_fdr <- signif(klrk1_res$p_val_adj, digits = 3)

Idents(cd8_cells) <- cd8_cells$condition

cd8_cells$condition <- factor(ifelse("B6B6" == cd8_cells$condition, "Syn", "Allo"))

VlnPlot(cd8_cells, features = c("Klrk1"), assay = "SCT", 
        pt.size = 0, alpha = 0.5, 
        group.by = "condition") + 
  ggprism::scale_fill_prism(palette = "shades_of_gray", name = "") + 
  NoLegend() + 
  labs(title = "CD8 T Cells", y = "Klrk1 Expression Level") + 
  theme(axis.title.x = element_blank(), 
        axis.text = element_text(size = 20, face = "bold"), 
        axis.text.x = element_text(angle = 0, hjust = 0.5), 
        title = element_text(size = 24, face = "bold"), 
        axis.line =  element_line(colour = 'black', linewidth = 2), 
        #axis.ticks = element_line(colour = 'black', linewidth = 2), 
        aspect.ratio = 1.7,
        ) + 
  #geom_text(label = paste0("Avg LogFC: ", klrk1_logfc, "\np (adj) < ", klrk1_fdr), 
  geom_text(label = "****", 
    x=1.5,
    y=2.2,
    #label.padding = unit(0.55, "lines"), # Rectangle size around label
    #label.size = 0.35, 
    fontface = "bold", 
    size = 10,
    color = "black"
  )

ggsave("mouse_klrk1.pdf", plot = last_plot(), path = figures_dir, width = 5, height = 7)
