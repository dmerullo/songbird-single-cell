library(tidyverse)
library(Seurat)
library(cowplot)
library(egg)

#LOAD SEURAT OBJECT
#LOAD DATATABLE
#CLUSTER SEURAT OBJECT
#DATATABLE
#CLUSTERMERGED COMPRISED OF CLUSTERORIG
#HVC CLUSTERORIG GOING INTO CLUSTERMERGED
#X CLUSTERORIG GOING INTO CLUSTERMERGED
#CLUSTERPLOT UNLABELED
#CLUSTERPLOT UNLABELED REGION
#PCT
#REORDER IDENTS
#DOTPLOT
#GABA
#GABA DOTPLOT
#GABA CLUSTERPLOT
#LGE
#MARKERS
#MSN VS MSN-LIKE
#HEATMAP
#GABA PRE
#GABA PRE DOTPLOT
#GABA PRE CLUSTERPLOT
#LGE PRE
#LGE CLUSTERPLOT
#LGE PRE VIOLIN
#LGE PRE BLEND
#ESPN

#LOAD SEURAT OBJECT
load(
  file = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "postdoc/projects/HVCX_COMBINED/revision/REGION/res1.1/",
    "HVCX_seurat_REGION_res1.1.RData"
  )
)

#LOAD DATATABLE
load(
  file = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "postdoc/projects/HVCX_COMBINED/revision/REGION/res1.1/",
    "HVCX_seurat_REGION_res1.1_datatable.RData"
  )
)

#CLUSTER SEURAT OBJECT

#LOAD REGION
load(
  file = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "postdoc/projects/HVCX_COMBINED/revision/REGION/",
    "HVCX_seurat_REGION.RData"
  )
)

#CLUSTER PLOT
HVCX_seurat_REGION_res1.1 <- 
  FindNeighbors(
    HVCX_seurat_REGION,
    dims = 1:16
  )

HVCX_seurat_REGION_res1.1 <-
  FindClusters(
    HVCX_seurat_REGION_res1.1,
    resolution = 1.1,
  )

HVCX_seurat_REGION_res1.1 <-
  RunUMAP(
    HVCX_seurat_REGION_res1.1,
    dims = 1:16,
    umap.method = "umap-learn",
    metric = "correlation"
  )

HVCX_seurat_REGION_res1.1 <-
  RenameIdents(
    HVCX_seurat_REGION_res1.1,
    "0" = "1",
    "1" = "2",
    "2" = "3",
    "3" = "4",
    "4" = "5",
    "5" = "6",
    "6" = "7",
    "7" = "8",
    "8" = "9",
    "9" = "10",
    "10" = "11",
    "11" = "12",
    "12" = "13",
    "13" = "14",
    "14" = "15",
    "15" = "16",
    "16" = "17",
    "17" = "18",
    "18" = "19",
    "19" = "20",
    "20" = "21",
    "21" = "22",
    "22" = "23",
    "23" = "24",
    "24" = "25",
    "25" = "26",
    "26" = "27",
    "27" = "28",
    "28" = "29",
    "29" = "30",
    "30" = "31",
    "31" = "32",
    "32" = "33",
    "33" = "34",
    "34" = "35",
    "35" = "36",
    "36" = "37",
    "37" = "38",
    "38" = "39",
    "39" = "40",
    "40" = "41"
  )

#REORDER IDENTS
Idents(HVCX_seurat_REGION_res1.1) <- 
  factor(
    Idents(HVCX_seurat_REGION_res1.1),
    levels = rev(c(
      1, 5, 7, 35, 25, 12, 33, 40, 
      2, 3, 9, 17, 34, 15, 38, 8, 14, 28, 37, 18,
      13, 11, 23, 19, 
      4, 6, 16, 20, 26, 39, 30, 22, 24, 21, 41, 29, 10, 27, 36, 32, 31
    ))
  )

save(
  HVCX_seurat_REGION_res1.1,
  file = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "postdoc/projects/HVCX_COMBINED/revision/REGION/res1.1/",
    "HVCX_seurat_REGION_res1.1.RData"
  )
)

HVCX_seurat_REGION_res1.1_clusterplot <-
  UMAPPlot(
    HVCX_seurat_REGION_res1.1,
    label = TRUE
  )

HVCX_seurat_REGION_res1.1_clusterplot_region <-
  UMAPPlot(
    HVCX_seurat_REGION_res1.1,
    label = TRUE,
    split.by = "region"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "HVCX_seurat_REGION_res1.1_clusterplot.pdf"
  ),
  plot = HVCX_seurat_REGION_res1.1_clusterplot,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "HVCX_seurat_REGION_res1.1_clusterplot_region.pdf"
  ),
  plot = HVCX_seurat_REGION_res1.1_clusterplot_region,
  width = 24,
  height = 8,
  units = "in",
  dpi = 300
)

#DATATABLE
HVCX_seurat_REGION_res1.1_data_tibble <-
  t(
    as.matrix(
      GetAssayData(
        HVCX_seurat_REGION_res1.1, 
        slot = "data",
        assay = "RNA",
        features = c(
          "SIX3",
          "FOXP1",
          "FOXP2",
          "ZFHX4",
          "MEIS2",
          "PBX3",
          "TSHZ1",
          "PVALB",
          "ETV1",
          "BCAN",
          "TRPS1",
          "KCNC1",
          "KCNC2",
          "SST",
          "SOX6",
          "LHX6",
          "ARX",
          "SATB1",
          "NPY",
          "EFNA5",
          "ELFN1",
          "CALB1",
          "NR2F2",
          "ADARB2",
          "CNR1",
          "ZBTB16",
          "PROX1",
          "RELN",
          "PENK"
        )
      )
    )
  ) %>%
  as_tibble(rownames = "Cell")

DefaultAssay(HVCX_seurat_REGION_res1.1) <- "RNA"

HVCX_seurat_REGION_res1.1_data_fetch <-
  FetchData(
  HVCX_seurat_REGION_res1.1,
  vars = c(
    "SIX3",
    "FOXP1",
    "FOXP2",
    "ZFHX4",
    "MEIS2",
    "PBX3",
    "TSHZ1",
    "PVALB",
    "ETV1",
    "BCAN",
    "TRPS1",
    "KCNC1",
    "KCNC2",
    "SST",
    "SOX6",
    "LHX6",
    "ARX",
    "SATB1",
    "NPY",
    "EFNA5",
    "ELFN1",
    "CALB1",
    "NR2F2",
    "ADARB2",
    "CNR1",
    "ZBTB16",
    "PROX1",
    "RELN",
    "PENK"
  )
)

HVCX_seurat_REGION_res1.1_cells <-
  enframe(Cells(HVCX_seurat_REGION_res1.1), name = NULL) %>%
  rename("Cell" = value)

HVCX_seurat_REGION_res1.1_clusters <-
  enframe(Idents(HVCX_seurat_REGION_res1.1), name = NULL) %>%
  rename("Cluster" = value) 

HVCX_seurat_REGION_res1.1_UMAP <-
  as_tibble(HVCX_seurat_REGION_res1.1@reductions[["umap"]]@cell.embeddings)

HVCX_seurat_REGION_res1.1_metadata <-
  FetchData(
    HVCX_seurat_REGION_res1.1,
    vars = c(
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mito",
      "dataset",
      "region",
      "species",
      "cluster_orig"
    )
  ) %>%
  as_tibble()

HVCX_seurat_REGION_res1.1_datatable <-
  bind_cols(
    HVCX_seurat_REGION_res1.1_cells,
    HVCX_seurat_REGION_res1.1_metadata,
    HVCX_seurat_REGION_res1.1_clusters,
    HVCX_seurat_REGION_res1.1_UMAP,
  ) %>%
  rename(
    UMIs = `nCount_RNA`,
    Genes = `nFeature_RNA`,
    PctMito = percent.mito,
    ClusterMerged = Cluster,
    ClusterOrig = `cluster_orig`
  ) %>%
  add_count(ClusterMerged, name = "ClusterMerged_n") %>%
  add_count(ClusterOrig, name = "ClusterOrig_n") %>%
  select(1:11, `ClusterMerged_n`, `ClusterOrig_n`, everything())

HVCX_seurat_REGION_res1.1_datatable2 <-
  bind_cols(
    HVCX_seurat_REGION_res1.1_cells,
    HVCX_seurat_REGION_res1.1_metadata,
    HVCX_seurat_REGION_res1.1_clusters,
    HVCX_seurat_REGION_res1.1_UMAP,
    HVCX_seurat_REGION_res1.1_data_fetch
  ) %>%
  rename(
    UMIs = `nCount_RNA`,
    Genes = `nFeature_RNA`,
    PctMito = percent.mito,
    ClusterMerged = Cluster,
    ClusterOrig = `cluster_orig`
  ) %>%
  add_count(ClusterMerged, name = "ClusterMerged_n") %>%
  add_count(ClusterOrig, name = "ClusterOrig_n")

HVCX_seurat_REGION_res1.1_datatable <-
  HVCX_seurat_REGION_res1.1_datatable %>%
  mutate(
    ClusterOrig = case_when(
      ClusterOrig == "0" ~ "1",
      ClusterOrig == "1" ~ "2",
      ClusterOrig == "2" ~ "3",
      ClusterOrig == "3" ~ "4",
      ClusterOrig == "4" ~ "5",
      ClusterOrig == "5" ~ "6",
      ClusterOrig == "6" ~ "7",
      ClusterOrig == "7" ~ "8",
      ClusterOrig == "8" ~ "9",
      ClusterOrig == "9" ~ "10",
      ClusterOrig == "10" ~ "11",
      ClusterOrig == "11" ~ "12",
      ClusterOrig == "12" ~ "13",
      ClusterOrig == "13" ~ "14",
      ClusterOrig == "14" ~ "15",
      ClusterOrig == "15" ~ "16",
      ClusterOrig == "16" ~ "17",
      ClusterOrig == "17" ~ "18",
      ClusterOrig == "18" ~ "19",
      ClusterOrig == "19" ~ "20",
      ClusterOrig == "20" ~ "21",
      ClusterOrig == "21" ~ "22",
      ClusterOrig == "22" ~ "23",
      TRUE ~ as.character(ClusterOrig)
    )
  )

HVCX_seurat_REGION_res1.1_datatable2 <-
  HVCX_seurat_REGION_res1.1_datatable2 %>%
  mutate(
    ClusterOrig = case_when(
      ClusterOrig == "0" ~ "1",
      ClusterOrig == "1" ~ "2",
      ClusterOrig == "2" ~ "3",
      ClusterOrig == "3" ~ "4",
      ClusterOrig == "4" ~ "5",
      ClusterOrig == "5" ~ "6",
      ClusterOrig == "6" ~ "7",
      ClusterOrig == "7" ~ "8",
      ClusterOrig == "8" ~ "9",
      ClusterOrig == "9" ~ "10",
      ClusterOrig == "10" ~ "11",
      ClusterOrig == "11" ~ "12",
      ClusterOrig == "12" ~ "13",
      ClusterOrig == "13" ~ "14",
      ClusterOrig == "14" ~ "15",
      ClusterOrig == "15" ~ "16",
      ClusterOrig == "16" ~ "17",
      ClusterOrig == "17" ~ "18",
      ClusterOrig == "18" ~ "19",
      ClusterOrig == "19" ~ "20",
      ClusterOrig == "20" ~ "21",
      ClusterOrig == "21" ~ "22",
      ClusterOrig == "22" ~ "23",
      TRUE ~ as.character(ClusterOrig)
    )
  )

save(
  HVCX_seurat_REGION_res1.1_datatable,
  file = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "postdoc/projects/HVCX_COMBINED/revision/REGION/res1.1/",
    "HVCX_seurat_REGION_res1.1_datatable.RData"
  )
)

#CLUSTERMERGED COMPRISED OF CLUSTERORIG
HVCX_seurat_REGION_res1.1_datatable %>%
  group_by(region, ClusterMerged, ClusterOrig) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(ClusterMerged == 9)

#HVC CLUSTERORIG GOING INTO CLUSTERMERGED
HVCX_seurat_REGION_res1.1_datatable %>%
  filter(region %in% c("HVC", "RA")) %>%
  group_by(ClusterOrig, ClusterMerged) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(ClusterOrig == "GABA-Pre")

#X CLUSTERORIG GOING INTO CLUSTERMERGED
HVCX_seurat_REGION_res1.1_datatable %>%
  filter(region == "X") %>%
  group_by(ClusterOrig, ClusterMerged) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(ClusterOrig == 18)

HVCX_seurat_REGION_res1.1_datatable %>%
  filter(ClusterMerged == 7, ClusterOrig %in% c("4", "5")) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = ClusterOrig)) +
  geom_point()

#CLUSTERPLOT UNLABELED
plot_HVCX_seurat_REGION_res1.1_clusterplot_unlabeled <-
  UMAPPlot(
    HVCX_seurat_REGION_res1.1,
    label = TRUE,
    label.size = 10,
    pt.size = 1.5
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_clusterplot_unlabeled.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_clusterplot_unlabeled,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

#CLUSTERPLOT UNLABELED REGION
plot_HVCX_seurat_REGION_res1.1_clusterplot_unlabeled_region <-
  UMAPPlot(
    HVCX_seurat_REGION_res1.1,
    label = TRUE,
    label.size = 10,
    pt.size = 1.5,
    split.by = "region"
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_clusterplot_unlabeled_region.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_clusterplot_unlabeled_region,
  width = 24,
  height = 8,
  units = "in",
  dpi = 300
)

#CLUSTERPLOT labeled
plot_HVCX_seurat_REGION_res1.1_clusterplot_labeled <-
  UMAPPlot(
    HVCX_seurat_REGION_res1.1,
    label = TRUE,
    label.size = 10,
    pt.size = 1.5
  ) +
  scale_color_manual(
    values = rev(c(
      "light blue",
      "light blue",
      "light blue",
      "light blue",
      "light blue",
      "light blue",
      "light blue",
      "light blue",
      "red",
      "red",
      "red",
      "red",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "blue",
      "blue",
      "blue",
      "blue",
      "pink",
      "pink",
      "pink",
      "light green",
      "light green",
      "light green",
      "light green",
      "dark green",
      "dark green",
      "brown",
      "brown",
      "brown",
      "turquoise",
      "turquoise",
      "purple",
      "purple",
      "gray"
      )
    )
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "postdoc/projects/HVCX_COMBINED/revision/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_clusterplot_labeled.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_clusterplot_labeled,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

#CLUSTERPLOT REGION
plot_HVCX_seurat_REGION_res1.1_clusterplot_labeled_region <-
  UMAPPlot(
    HVCX_seurat_REGION_res1.1,
    label = TRUE,
    label.size = 10,
    pt.size = 1.5,
    split.by = "region",
    cols = rev(c(
      "light blue",
      "light blue",
      "light blue",
      "light blue",
      "light blue",
      "light blue",
      "light blue",
      "light blue",
      "red",
      "red",
      "red",
      "red",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "blue",
      "blue",
      "blue",
      "blue",
      "blue",
      "pink",
      "pink",
      "pink",
      "light green",
      "light green",
      "light green",
      "light green",
      "dark green",
      "dark green",
      "brown",
      "brown",
      "brown",
      "turquoise",
      "turquoise",
      "purple",
      "purple",
      "gray"
    ))
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_clusterplot_labeled_region.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_clusterplot_labeled_region,
  width = 24,
  height = 8,
  units = "in",
  dpi = 300
)

#DOTPLOT
plot_HVCX_seurat_REGION_res1.1_dotplot_all <-
  DotPlot(
    HVCX_seurat_REGION_res1.1,
    assay = "RNA",
    features = rev(c(
      "SLC17A6",
      "GAD2",
      "MLPH",
      "TLL1",
      "NECTIN3",
      "TBB2",
      "LRIG1",
      "MBP",
      "PDGFRA",
      "CSF1R",
      "FLT1",
      "HBAA",
      "RGS5"
    )),
    col.min = 0,
    col.max = 1,
    cols = c("gray90", "darkslateblue"),
    dot.scale = 10
  ) +
  scale_x_discrete(
    labels = c(
      "Slc17a6",
      "Gad2",
      "Lrig1",
      "Mbp",
      "Pdgfra",
      "Csf1r",
      "Flt1"
    )
  ) +
  scale_y_discrete(
    labels = rev(c(
      "1 (HVC Glut-1/4/5/6)",
      "5 (HVC Glut-2)",
      "12 (HVC Glut-3)",
      "23 (HVC Glut-1/2)",
      "26 (HVC Glut-7)",
      "14 (RA Glut-1/2/3)",
      "32 (X Glut+ INT)",
      "35 (X Glut+ Unknown)",
      "2 (X MSN)", 
      "3 (X MSN)",
      "7 (X MSN)",
      "21 (X MSN)",
      "31 (X PN)",
      "16 (HVC GABA-1)",
      "34 (RA GABA-1)",
      "11 (GABA-2/3)",
      "15 (GABA-7/8 & X INT)",
      "17 (GABA-3/4 & X INT)",
      "28 (GABA-5)",
      "33 (GABA-6)",
      "8 (Pre-1/2)",
      "10 (Pre-3/4)",
      "24 (Pre-5)",
      "36 (GABA-Pre)",
      "4 (HVC Astrocyte)",
      "6 (RA Astrocyte)",
      "18 (X Astrocyte)",
      "19 (HVC Oligodendrocyte)",
      "20 (RA Oligodendrocyte)",
      "27 (X Oligodendrocyte)",
      "13 (HVC/RA Oligo. Precursor)",
      "38 (X Oligo. Precursor)",
      "22 (HVC/X Microglia)",
      "30 (RA Microglia)",
      "9 (HVC/X Endothelial)",
      "25 (RA Endothelial)",
      "29 (RBC)",
      "37 (HVC/RA Mural)"
    ))
  ) +
  theme(
    axis.text.x = element_text(size = 24, hjust = 1, angle = 45, face = "italic"),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.key.size = unit(2, "lines"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_dotplot_all.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_dotplot_all,
  width = 15,
  height = 15,
  units = "in",
  dpi = 300
)

#PCT
HVCX_seurat_REGION_res1.1_datatable_pct <-
  HVCX_seurat_REGION_res1.1_datatable %>%
  select(1:13) %>%
  add_count(ClusterMerged, region, name = "nClusterMergedRegion") %>%
  distinct(ClusterMerged, region, .keep_all = TRUE) %>% 
  mutate(
    pct = nClusterMergedRegion / ClusterMerged_n * 100,
    ClusterMerged = factor(
      ClusterMerged,
      levels = rev(c(
        1, 5, 7, 35, 25, 12, 33, 40, 
        2, 3, 9, 17, 34, 15, 38, 8, 14, 28, 37, 18,
        13, 11, 23, 19, 
        4, 6, 16, 20, 26, 39, 30, 22, 24, 21, 41, 29, 10, 27, 36, 32, 31
      )
    )),
    region = factor(region, levels = rev(c(
      "HVC", "RA", "X"
    )))
    )

plot_HVCX_seurat_REGION_res1.1_pct <-
  HVCX_seurat_REGION_res1.1_datatable_pct %>%
  ggplot(
    aes(
      x = ClusterMerged,
      y = pct,
      fill = region
    )
  ) +
  coord_flip() +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  ylab("Percentage Composition") +
  scale_fill_manual(
    values = rev(c("tomato2", "slategray2", "darkgoldenrod2")), 
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_x_discrete(
    labels = rev(c(
      "1 (HVC Glut-1/4)",
      "5 (HVC Glut-2)",
      "7 (HVC Glut-3)",
      "35 (HVC Glut-5)",
      "25 (HVC Glut-1/2)",
      "12 (RA Glut-1/2/3 & HVC Glut-1)",
      "33 (X Glut+ INT)",
      "40 (X Glut+ Unknown)",
      "2 (X MSN)", 
      "3 (X MSN)",
      "9 (X MSN)",
      "17 (X MSN)",
      "34 (X PN)",
      "15 (HVC GABA-1)",
      "38 (RA GABA-1)",
      "8 (GABA-2/4)",
      "14 (GABA-3/4 & X INT)",
      "28 (GABA-5 & X INT)",
      "37 (GABA-6 & X INT)",
      "18 (GABA-7/8 & X INT)",
      "13 (Pre-1/Ependymal)",
      "11 (Pre-2)",
      "23 (Pre-3/4)",
      "19 (GABA-Pre)",
      "4 (HVC Astrocyte)",
      "6 (RA Astrocyte)",
      "16 (X Astrocyte)",
      "20 (HVC Oligodendrocyte)",
      "26 (RA Oligodendrocyte)",
      "39 (RA Oligodendrocyte)",
      "30 (X Oligodendrocyte)",
      "22 (HVC/X Oligo. Precursor)",
      "24 (RA Oligo. Precursor)",
      "21 (HVC/X Microglia)",
      "41 (HVC Microglia)",
      "29 (RA Microglia)",
      "10 (HVC/X Endothelial)",
      "27 (RA Endothelial)",
      "32 (HVC/RA RBC)",
      "36 (X RBC)",
      "31 (HVC/RA Mural)"
    ))
  ) +
  theme(
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 24),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "postdoc/projects/HVCX_COMBINED/revision/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_pct.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_pct,
  width = 12,
  height = 14,
  units = "in",
  dpi = 300
)

#GABA
HVCX_seurat_REGION_res1.1_GABA <-
  subset(
    HVCX_seurat_REGION_res1.1,
    slot = "data",
    idents = c(
      2, 3, 7, 21, 31, 16, 34, 11, 15, 17, 28, 33
    )
  )

#GABA DOTPLOT
plot_HVCX_seurat_REGION_res1.1_dotplot_GABA <-
  DotPlot(
    HVCX_seurat_REGION_res1.1_GABA,
    assay = "RNA",
    features = rev(c(
      "SIX3",
      "FOXP1",
      "FOXP2",
      "ZFHX4",
      "MEIS2",
      "PBX3",
      "TSHZ1",
      "PVALB",
      "ETV1",
      "BCAN",
      "TRPS1",
      "KCNC1",
      "KCNC2",
      "SST",
      "SOX6",
      "LHX6",
      "ARX",
      "SATB1",
      "NPY",
      "EFNA5",
      "ELFN1",
      "CALB1",
      "NR2F2",
      "ADARB2",
      "CNR1",
      "ZBTB16",
      "PROX1",
      "RELN",
      "PENK"
    )),
    col.min = 0,
    col.max = 10,
    cols = c("gray90", "darkslateblue"),
    dot.scale = 10
  ) +
  scale_y_discrete(
    labels = rev(c(
      "2 (X MSN)", 
      "3 (X MSN)",
      "7 (X MSN)",
      "21 (X MSN)",
      "31 (X PN)",
      "16 (HVC GABA-1)",
      "34 (RA GABA-1)",
      "11 (HVC/RA GABA-2/3)",
      "15 (HVC/RA GABA-7/8 & X INT)",
      "17 (HVC/RA GABA-3/4 & X INT)",
      "28 (HVC/RA GABA-5)",
      "33 (HVC/RA GABA-6)"
    ))
  ) +
  scale_x_discrete(
    labels = c(
      "Six3",
      "FoxP1",
      "FoxP2",
      "Zfhx4",
      "Meis2",
      "Pbx3",
      "Tshz1",
      "Pvalb",
      "Etv1",
      "Bcan",
      "Trps1",
      "Kcnc1",
      "Kcnc2",
      "Sst",
      "Sox6",
      "Lhx6",
      "Arx",
      "Satb1",
      "Npy",
      "Efna5",
      "Elfn1",
      "Calb1",
      "Nr2f2",
      "Adarb2",
      "Cnr1",
      "Zbtb16",
      "Prox1",
      "Reln",
      "Penk"
    )
  ) +
  theme(
    axis.text.x = element_text(size = 24, hjust = 1, angle = 45, face = "italic"),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.key.size = unit(2, "lines")
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_dotplot_GABA.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_dotplot_GABA,
  width = 20,
  height = 8,
  units = "in",
  dpi = 300
)

#GABA CLUSTERPLOT
HVCX_seurat_REGION_res1.1_GABA <-
  RunUMAP(
    HVCX_seurat_REGION_res1.1_GABA,
    dims = 1:15
  )

plot_HVCX_seurat_REGION_res1.1_GABA_clusterplot <-
  UMAPPlot(
    HVCX_seurat_REGION_res1.1_GABA,
    label = TRUE,
    label.size = 10,
    pt.size = 1.5
  ) +
  scale_color_manual(
    values = rev(
      c(
        "tomato1",
        "tomato2",
        "tomato3",
        "tomato4",
        "light green",
        "slateblue",
        "slateblue4",
        "darkgoldenrod",
        "darkgoldenrod1",
        "darkgoldenrod2",
        "darkgoldenrod3",
        "darkgoldenrod4"
      )
    )
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_GABA_clusterplot.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_GABA_clusterplot,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

#LGE
HVCX_seurat_REGION_res1.1_LGE <-
  subset(
    HVCX_seurat_REGION_res1.1,
    slot = "data",
    idents = c(
      2, 3, 7, 21, 31, 16, 34
    )
  )

HVCX_seurat_REGION_res1.1_LGE <-
  RunUMAP(
    HVCX_seurat_REGION_res1.1_LGE,
    dims = 1:15
  )

plot_HVCX_seurat_REGION_res1.1_LGE_clusterplot <-
  UMAPPlot(
    HVCX_seurat_REGION_res1.1_LGE,
    label = TRUE,
    label.size = 10,
    pt.size = 1.5
  ) +
  scale_color_manual(
    values = rev(c(
      "tomato1",
      "tomato2",
      "tomato3",
      "tomato4",
      "light green",
      "slateblue",
      "slateblue4"
      )
    )
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_LGE_clusterplot.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_LGE_clusterplot,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

#MARKERS
DefaultAssay(HVCX_seurat_REGION_res1.1) <-
  "RNA"

HVCX_seurat_REGION_res1.1_markers <-
  as_tibble(
    FindAllMarkers(
      HVCX_seurat_REGION_res1.1,
      assay = "RNA",
      slot = "data"
    )
  )

write_csv(
  HVCX_seurat_REGION_res1.1_markers,
  file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "HVCX_seurat_REGION_res1.1_markers.csv"
  )
)

HVCX_seurat_REGION_res1.1_markers <-
  read_csv(
    file.path(
      "~/OneDrive - University of Texas Southwestern/",
      "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
      "HVCX_seurat_REGION_res1.1_markers.csv"
    )
  )


HVCX_seurat_REGION_res1.1_markers %>%
  filter(
    cluster == 29
  ) %>%
  arrange(-avg_logFC)

HVCX_seurat_REGION_res1.1_markers %>%
  filter(
    cluster %in% c(8, 10, 24, 36, 29, 37)
  ) %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC) %>%
  arrange(cluster, -avg_logFC)

#MSN VS MSN-LIKE
DefaultAssay(HVCX_seurat_REGION_res1.1) <-
  "RNA"

HVCX_seurat_REGION_res1.1_DEG1_markers <-
    FindMarkers(
      HVCX_seurat_REGION_res1.1,
      assay = "RNA",
      slot = "data",
      ident.1 = c(2, 3, 8, 11),
      ident.2 = c(16, 38)
    ) %>%
  as_tibble(rownames = "gene")

write_csv(
  HVCX_seurat_REGION_res1.1_DEG1_markers,
  file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "HVCX_seurat_REGION_res1.1_DEG1_markers.csv"
  )
)

HVCX_seurat_REGION_res1.1_DEG1_markers <-
  read_csv(
    file.path(
      "~/OneDrive - University of Texas Southwestern/",
      "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
      "HVCX_seurat_REGION_res1.1_DEG1_markers.csv"
    )
  )

HVCX_seurat_REGION_res1.1_DEG1 <-
  subset(
    HVCX_seurat_REGION_res1.1,
    idents = c(2, 3, 8, 11, 16, 38)
  )

HVCX_seurat_REGION_res1.1_DEG1_avexp <- 
  log1p(
    AverageExpression(
      HVCX_seurat_REGION_res1.1_DEG1
    )
    $RNA
  ) %>%
  as_tibble(rownames = "gene") 

HVCX_seurat_REGION_res1.1_DEG1_avexp_stats <-
  inner_join(
    HVCX_seurat_REGION_res1.1_DEG1_avexp,
    HVCX_seurat_REGION_res1.1_DEG1_markers, 
    by = "gene"
  )

HVCX_seurat_REGION_res1.1_DEG1_avexp_stats %>%
  filter(
    p_val_adj < 0.5
  ) %>%
  arrange(
    -avg_logFC
  )

write_csv(
  HVCX_seurat_REGION_res1.1_DEG1_avexp_stats,
  file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "HVCX_markers_MSNvsHVC.csv"
  )
)

HVCX_markers_MSNvsHVC <-
  read_csv(
    file.path(
      "~/OneDrive - University of Texas Southwestern/",
      "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
      "HVCX_markers_MSNvsHVC.csv"
    )
  )

#TOP IN MSN
HVCX_markers_MSNvsHVC %>%
  filter(
    p_val_adj < 0.5,
  ) %>%
  group_by(gene) %>%
  filter(
    avg_logFC == max(avg_logFC)
  ) %>%
  ungroup() %>%
  arrange(
    -pct.1, -avg_logFC
  ) %>%
  top_n(n = 5, wt = avg_logFC) 

#TOP IN HVC
HVCX_markers_MSNvsHVC %>%
  filter(
    p_val_adj < 0.5,
  ) %>%
  group_by(gene) %>%
  filter(
    avg_logFC == max(avg_logFC)
  ) %>%
  ungroup() %>%
  arrange(
    -avg_logFC
  ) %>%
  top_n(n = 5, wt = -avg_logFC)

#HEATMAP
HVCX_seurat_REGION_res1.1_heatmap <-
  HVCX_seurat_REGION_res1.1_datatable %>%
  select(
    ClusterMerged,
    SLC17A6,
    GAD2,
    MLPH,
    TLL1,
    NECTIN3,
    TBB2,
    LRIG1,
    MBP,
    PDGFRA,
    CSF1R,
    FLT1,
    HBAA,
    RGS5
  ) %>%
  mutate(
    ClusterMerged = factor(
      ClusterMerged,
      levels = rev(c(
        1, 5, 12, 23, 26, 14, 32, 35,
        2, 3, 7, 21, 31, 16, 34,
        11, 15, 17, 28, 33,
        8, 10, 24, 36,
        4, 6, 18, 19, 20, 27, 13, 38, 22, 30, 9, 25, 29, 37
      ))
    )
  ) %>%
  pivot_longer(
    cols = -1,
    names_to = "Gene", 
    values_to = "Expression"
  ) %>%
  group_by(Gene, ClusterMerged) %>%
  mutate(
    Expression = mean(Expression),
  ) %>%
  distinct(ClusterMerged, .keep_all = TRUE)

plot_HVCX_seurat_REGION_res1.1_heatmap_SLC17A6 <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "SLC17A6") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  ) +
  scale_y_discrete(labels = rev(c(
    "1 (HVC Glut-1/4/5/6)",
    "5 (HVC Glut-2)",
    "12 (HVC Glut-3)",
    "23 (HVC Glut-1/2)",
    "26 (HVC Glut-7)",
    "14 (RA Glut-1/2/3)",
    "32 (X Glut+ INT)",
    "35 (X Glut+ Unknown)",
    "2 (X MSN)", 
    "3 (X MSN)",
    "7 (X MSN)",
    "21 (X MSN)",
    "31 (X PN)",
    "16 (HVC GABA-1)",
    "34 (RA GABA-1)",
    "11 (GABA-2/3)",
    "15 (GABA-7/8 & X INT)",
    "17 (GABA-3/4 & X INT)",
    "28 (GABA-5)",
    "33 (GABA-6)",
    "8 (Pre-1/2)",
    "10 (Pre-3/4)",
    "24 (Pre-5)",
    "36 (GABA-Pre)",
    "4 (HVC Astrocyte)",
    "6 (RA Astrocyte)",
    "18 (X Astrocyte)",
    "19 (HVC Oligodendrocyte)",
    "20 (RA Oligodendrocyte)",
    "27 (X Oligodendrocyte)",
    "13 (HVC/RA Oligo. Precursor)",
    "38 (X Oligo. Precursor)",
    "22 (HVC/X Microglia)",
    "30 (RA Microglia)",
    "9 (HVC/X Endothelial)",
    "25 (RA Endothelial)",
    "29 (RBC)",
    "37 (HVC/RA Mural)"
  ))) +
  scale_x_discrete(labels = c("Slc17a6")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_text(size = 76.8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_GAD2 <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "GAD2") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  ) +
  scale_x_discrete(labels = c("Gad2")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_MLPH <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "MLPH") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  ) +
  scale_x_discrete(labels = c("Mlph")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_TLL1 <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "TLL1") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  ) +
  scale_x_discrete(labels = c("TLL1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_NECTIN3 <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "NECTIN3") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  ) +
  scale_x_discrete(labels = c("Nectin3")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_TBB2 <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "TBB2") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_x_discrete(labels = c("Tbb2")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_LRIG1 <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "LRIG1") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  ) +
  scale_x_discrete(labels = c("Lrig1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_MBP <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "MBP") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 5),
    breaks = c(0, 5)
  ) +
  scale_x_discrete(labels = c("Mbp")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_PDGFRA <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "PDGFRA") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  ) +
  scale_x_discrete(labels = c("Pdgfra")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_CSF1R <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "CSF1R") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  ) +
  scale_x_discrete(labels = c("Csf1r")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_FLT1 <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "FLT1") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  ) +
  scale_x_discrete(labels = c("Flt1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_HBAA <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "HBAA") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 8),
    breaks = c(0, 8)
  ) +
  scale_x_discrete(labels = c("Hbaa")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap_RGS5 <-
  HVCX_seurat_REGION_res1.1_heatmap %>%
  filter(Gene == "RGS5") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  ) +
  scale_x_discrete(labels = c("Rgs5")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_HVCX_seurat_REGION_res1.1_heatmap <-
  ggarrange(
    plot_HVCX_seurat_REGION_res1.1_heatmap_SLC17A6,
    plot_HVCX_seurat_REGION_res1.1_heatmap_GAD2,
    plot_HVCX_seurat_REGION_res1.1_heatmap_MLPH,
    plot_HVCX_seurat_REGION_res1.1_heatmap_TLL1,
    plot_HVCX_seurat_REGION_res1.1_heatmap_NECTIN3,
    plot_HVCX_seurat_REGION_res1.1_heatmap_TBB2,
    plot_HVCX_seurat_REGION_res1.1_heatmap_LRIG1,
    plot_HVCX_seurat_REGION_res1.1_heatmap_MBP,
    plot_HVCX_seurat_REGION_res1.1_heatmap_PDGFRA,
    plot_HVCX_seurat_REGION_res1.1_heatmap_CSF1R,
    plot_HVCX_seurat_REGION_res1.1_heatmap_FLT1,
    plot_HVCX_seurat_REGION_res1.1_heatmap_HBAA,
    plot_HVCX_seurat_REGION_res1.1_heatmap_RGS5,
    nrow = 1
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_heatmap.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_heatmap,
  width = 60,
  height = 45,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

#GABA PRE
HVCX_seurat_REGION_res1.1_GABA_PRE <-
  subset(
    HVCX_seurat_REGION_res1.1,
    slot = "data",
    idents = c(
      8, 14, 28, 37, 18, 19, 15, 38, 2, 3, 9, 17, 34
    )
  )

#REORDER IDENTS GABA PRE
Idents(HVCX_seurat_REGION_res1.1_GABA_PRE) <- 
  factor(
    Idents(HVCX_seurat_REGION_res1.1_GABA_PRE),
    levels = rev(c(
      2, 3, 9, 17, 34, 15, 38, 8, 14, 28, 37, 18, 19
    ))
  )

#GABA PRE DOTPLOT
plot_HVCX_seurat_REGION_res1.1_dotplot_GABA_PRE <-
  DotPlot(
    HVCX_seurat_REGION_res1.1_GABA_PRE,
    assay = "RNA",
    features = rev(c(
      "SIX3",
      "FOXP1",
      "FOXP2",
      "ZFHX4",
      "MEIS2",
      "PBX3",
      "TSHZ1",
      "PVALB",
      "ETV1",
      "BCAN",
      "TRPS1",
      "KCNC1",
      "KCNC2",
      "SST",
      "SOX6",
      "LHX6",
      "ARX",
      "SATB1",
      "NPY",
      "EFNA5",
      "ELFN1",
      "CALB1",
      "NR2F2",
      "ADARB2",
      "CNR1",
      "ZBTB16",
      "PROX1",
      "RELN",
      "PENK"
    )),
    scale.max = 100,
    cols = c("gray90", "darkslateblue"),
    dot.scale = 10
  ) +
  scale_y_discrete(
    labels = rev(c(
      "2 (X MSN)", 
      "3 (X MSN)",
      "9 (X MSN)",
      "17 (X MSN)",
      "34 (X PN)",
      "15 (HVC GABA-1)",
      "38 (RA GABA-1)",
      "8 (GABA-2)",
      "14 (GABA-3/4 & X INT)",
      "28 (GABA-6 & X INT)",
      "18 (GABA-7/8 & X INT)",
      "37 (GABA-5 & X INT)",
      "19 (GABA-Pre)"
    ))
  ) +
  theme(
    axis.text.x = element_text(size = 24, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.key.size = unit(2, "lines"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_dotplot_GABA_PRE2.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_dotplot_GABA_PRE,
  width = 20,
  height = 8,
  units = "in",
  dpi = 300
)

#GABA PRE CLUSTERPLOT
HVCX_seurat_REGION_res1.1_GABA_PRE <-
  RunUMAP(
    HVCX_seurat_REGION_res1.1_GABA_PRE,
    dims = 1:16,
    umap.method = "umap-learn",
    metric = "correlation"
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_clusterplot <-
  UMAPPlot(
    HVCX_seurat_REGION_res1.1_GABA_PRE,
    label = TRUE,
    label.size = 10,
    pt.size = 1.5
  ) +
  scale_color_manual(
    values = rev(
      c(
        "tomato1",
        "tomato2",
        "tomato3",
        "tomato4",
        "light green",
        "slateblue",
        "slateblue4",
        "darkgoldenrod",
        "darkgoldenrod1",
        "darkgoldenrod2",
        "darkgoldenrod3",
        "darkgoldenrod4",
        "turquoise"
      )
    )
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/revision/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_GABA_PRE_clusterplot.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_GABA_PRE_clusterplot,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_GABA_PRE_clusterplot.png"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_GABA_PRE_clusterplot,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

#LGE PRE
HVCX_seurat_REGION_res1.1_LGE_PRE <-
  subset(
    HVCX_seurat_REGION_res1.1,
    slot = "data",
    idents = c(
      2, 3, 8, 11, 35, 16, 38, 41
    )
  )

HVCX_seurat_REGION_res1.1_LGE_PRE <-
  RenameIdents(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    "2" = "MSN",
    "3" = "MSN",
    "8" = "MSN",
    "11" = "MSN",
    "35" = "PN",
    "16" = "GABA-1",
    "38" = "GABA-1",
    "41" = "GABA-Pre"
    )

Idents(HVCX_seurat_REGION_res1.1_LGE_PRE) <- 
  factor(
    Idents(HVCX_seurat_REGION_res1.1_LGE_PRE),
    levels = rev(c(
      "MSN", "GABA-1", "PN", "GABA-Pre"
    ))
  )

#LGE PRE VIOLIN
HVCX_seurat_REGION_res1.1_LGE_PRE_FOXP1 <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "FOXP1",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 5),
    breaks = c(0, 5)
  ) +
  scale_x_discrete(
    labels = rev(c(
      "MSN",
      "GABA-1",
      "PN",
      "GABA-Pre"
    ))
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  labs(
    y = "Normalized Expression"
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )

HVCX_seurat_REGION_res1.1_LGE_PRE_FOXP2 <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "FOXP2",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 6),
    breaks = c(0, 6)
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )

HVCX_seurat_REGION_res1.1_LGE_PRE_PDE7B <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "PDE7B",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 6),
    breaks = c(0, 6)
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )

HVCX_seurat_REGION_res1.1_LGE_PRE_PCP4 <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "PCP4",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 6),
    breaks = c(0, 6)
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )

HVCX_seurat_REGION_res1.1_LGE_PRE_CNTN5 <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "CNTN5",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 6),
    breaks = c(0, 6)
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )

HVCX_seurat_REGION_res1.1_LGE_PRE_TAC1 <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "TAC1",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 6),
    breaks = c(0, 6)
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )

HVCX_seurat_REGION_res1.1_LGE_PRE_KCNT2 <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "KCNT2",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 5),
    breaks = c(0, 5)
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )

HVCX_seurat_REGION_res1.1_LGE_PRE_LRP1B <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "LRP1B",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 6),
    breaks = c(0, 6)
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )
 
HVCX_seurat_REGION_res1.1_LGE_PRE_GALNTL6 <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "GALNTL6",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 5),
    breaks = c(0, 5)
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )

HVCX_seurat_REGION_res1.1_LGE_PRE_TSHZ1 <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "TSHZ1",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 5),
    breaks = c(0, 5)
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )

HVCX_seurat_REGION_res1.1_LGE_PRE_RALYL <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "RALYL",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 5),
    breaks = c(0, 5)
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )

HVCX_seurat_REGION_res1.1_LGE_PRE_TRPM3 <-
  VlnPlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = "TRPM3",
    pt.size = 0,
    assay = "RNA",
    slot = "data"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(NA, 6),
    breaks = c(0, 6)
  ) +
  scale_fill_manual(
    values = rev(c(
      "red",
      "light green", 
      "purple",
      "turquoise"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24, angle = 0),
    plot.title = element_text(size = 24, face = "italic"),
    legend.position = "none"
  )

plot_HVCX_seurat_REGION_res1.1_LGE_PRE <-
  ggarrange(
    HVCX_seurat_REGION_res1.1_LGE_PRE_FOXP1,
    HVCX_seurat_REGION_res1.1_LGE_PRE_FOXP2,
    HVCX_seurat_REGION_res1.1_LGE_PRE_PCP4,
    HVCX_seurat_REGION_res1.1_LGE_PRE_PDE7B,
    HVCX_seurat_REGION_res1.1_LGE_PRE_CNTN5,
    HVCX_seurat_REGION_res1.1_LGE_PRE_KCNT2,
    HVCX_seurat_REGION_res1.1_LGE_PRE_TAC1,
    HVCX_seurat_REGION_res1.1_LGE_PRE_LRP1B,
    HVCX_seurat_REGION_res1.1_LGE_PRE_GALNTL6,
    HVCX_seurat_REGION_res1.1_LGE_PRE_TSHZ1,
    HVCX_seurat_REGION_res1.1_LGE_PRE_RALYL,
    HVCX_seurat_REGION_res1.1_LGE_PRE_TRPM3,
    nrow = 1
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_LGE_PRE.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_LGE_PRE,
  width = 18,
  height = 6,
  units = "in",
  dpi = 300
)

#LGE PRE BLEND
HVCX_seurat_REGION_res1.1_LGE_PRE_blend <- 
  FeaturePlot(
    HVCX_seurat_REGION_res1.1_LGE_PRE,
    features = c("DCX", "SOX4"),
    label.size = 10,
    min.cutoff = 0,
    pt.size = 1.5,
    blend = TRUE,
    combine = FALSE
  )

plot_HVCX_seurat_REGION_res1.1_LGE_PRE_blend <-
  HVCX_seurat_REGION_res1.1_LGE_PRE_blend[[3]] +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "none",
      plot.title = element_blank()
    )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_LGE_PRE_blend.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_LGE_PRE_blend,
  width = 8,
  height = 8,
  units = "in",
  limitsize = FALSE,
  dpi = 300
)

plot_HVCX_seurat_REGION_res1.1_LGE_PRE_blend_legend <-
  HVCX_seurat_REGION_res1.1_LGE_PRE_blend[[4]] +
  labs(
    x = "Relative Dcx Expression",
    y = "Relative Sox4 Expression"
  ) +
  scale_x_continuous(
    limits = c(0, 10),
    breaks = c(0, 10)
  ) +
  scale_y_continuous(
    limits = c(0, 10),
    breaks = c(0, 10)
  ) +
  theme(
    plot.title = element_blank(),
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24)
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_LGE_PRE_blend_legend.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_LGE_PRE_blend_legend,
  width = 4,
  height = 4,
  units = "in",
  limitsize = FALSE,
  dpi = 300
)

#ESPN
plot_HVCX_seurat_REGION_res1.1_dotplot_eSPN <-
  DotPlot(
    HVCX_seurat_REGION_res1.1,
    assay = "RNA",
    features = rev(c(
      "DRD1",
      "DRD2",
      "ADORA2A",
      "CASZ1",
      "OTOF",
      "CACNG5",
      "PCDH8",
      "ASIC2",
      "NTNG1",
      "VWCL2",
      "TSHZ1"
    )),
    col.min = 0,
    col.max = 1,
    cols = c("gray90", "darkslateblue"),
    dot.scale = 10
  ) +
  scale_x_discrete(
    labels = c(
      "Drd1",
      "Drd2",
      "Adora2a",
      "Casz1",
      "Otof",
      "Cacng5",
      "Pcdh8",
      "Ntng1",
      "Tshz1"
    )
  ) +
  scale_y_discrete(
    labels = rev(c(
      "1 (HVC-RA Glut-1/4/5/6)",
      "5 (Glut-2 Imm. HVC-RA)",
      "12 (Glut-3 HVC-X)",
      "23 (HVC Glut-1/2)",
      "26 (Glut-7 HVC-AV)",
      "14 (RA Glut-1/2/3)",
      "32 (X Glut+ INT)",
      "35 (X Glut+ Unknown)",
      "2 (X MSN)", 
      "3 (X MSN)",
      "7 (X MSN)",
      "21 (X MSN)",
      "31 (X PN)",
      "16 (HVC GABA-1)",
      "34 (RA GABA-1)",
      "11 (HVC/RA GABA-2/3)",
      "15 (HVC/RA GABA-7/8 & X INT)",
      "17 (HVC/RA GABA-3/4 & X INT)",
      "28 (HVC/RA GABA-5)",
      "33 (HVC/RA GABA-6)",
      "8 (HVC Glut-Pre-1/2)",
      "10 (HVC Glut-Pre-3/4)",
      "24 (HVC Glut-Pre-5)",
      "36 (HVC GABA-Pre)",
      "4 (HVC Astrocyte)",
      "6 (RA Astrocyte)",
      "18 (X Astrocyte)",
      "19 (HVC Oligodendrocyte)",
      "20 (RA Oligodendrocyte)",
      "27 (X Oligodendrocyte)",
      "13 (HVC/RA Oligo. Precursor)",
      "38 (X Oligo. Precursor)",
      "22 (HVC/X Microglia)",
      "30 (RA Microglia)",
      "9 (HVC/X Endothelial)",
      "25 (RA Endothelial)",
      "29 (RBC)",
      "37 (HVC/RA Mural)"
    ))
  ) +
  theme(
    axis.text.x = element_text(size = 24, hjust = 1, angle = 45, face = "italic"),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.key.size = unit(2, "lines"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/ZFBF/01/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_dotplot_eSPN.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_dotplot_eSPN,
  width = 15,
  height = 15,
  units = "in",
  dpi = 300
)















HVCX_seurat_REGION_res1.1_GABA_PRE <-
  subset(
    HVCX_seurat_REGION_res1.1,
    slot = "data",
    idents = c(
      2, 3, 7, 21, 31, 16, 34, 11, 15, 17, 28, 33, 36
    )
  )

#REORDER IDENTS GABA PRE
Idents(HVCX_seurat_REGION_res1.1_GABA_PRE) <- 
  factor(
    Idents(HVCX_seurat_REGION_res1.1_GABA_PRE),
    levels = rev(c(
      2, 3, 7, 21, 31, 16, 34, 36, 11, 15, 17, 28, 33
    ))
  )





#HEATMAP2
HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap <-
  HVCX_seurat_REGION_res1.1_datatable2 %>%
  filter(
    ClusterMerged %in% c(
      8, 14, 28, 37, 18, 19, 15, 38, 2, 3, 9, 17, 34
    )
  ) %>%
  select(
    ClusterMerged,
    SIX3,
    FOXP1,
    FOXP2,
    ZFHX4,
    MEIS2,
    PBX3,
    TSHZ1,
    PVALB,
    ETV1,
    BCAN,
    TRPS1,
    KCNC1,
    KCNC2,
    SST,
    SOX6,
    LHX6,
    ARX,
    SATB1,
    NPY,
    EFNA5,
    ELFN1,
    CALB1,
    NR2F2,
    ADARB2,
    CNR1,
    ZBTB16,
    PROX1,
    RELN,
    PENK
  ) %>%
  mutate(
    ClusterMerged = factor(
      ClusterMerged,
      levels = rev(c(
        2, 3, 9, 17, 34, 15, 38, 8, 14, 28, 37, 18, 19
      ))
    )
  ) %>%
  pivot_longer(
    cols = -1,
    names_to = "Gene", 
    values_to = "Expression"
  ) %>%
  group_by(Gene, ClusterMerged) %>%
  mutate(
    Expression = mean(Expression),
  ) %>%
  distinct(ClusterMerged, .keep_all = TRUE)

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_SIX3 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap %>%
  filter(Gene == "SIX3") %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = rev(c(
    "2 (X MSN)", 
    "3 (X MSN)",
    "9 (X MSN)",
    "17 (X MSN)",
    "34 (X PN)",
    "15 (HVC GABA-1)",
    "38 (RA GABA-1)",
    "8 (GABA-2/4)",
    "14 (GABA-3/4 & X INT)",
    "28 (GABA-5 & X INT)",
    "37 (GABA-6 & X INT)",
    "18 (GABA-7/8 & X INT)",
    "19 (GABA-Pre)"
  ))) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 76.8, face = "italic"),
    axis.text.y = element_text(size = 76.8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 76.8),
    legend.key.size = unit(2, "lines"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_genes <- 
  c(
    "FOXP1",
    "FOXP2",
    "ZFHX4",
    "MEIS2",
    "PBX3",
    "TSHZ1",
    "PVALB",
    "ETV1",
    "BCAN",
    "TRPS1",
    "KCNC1",
    "KCNC2",
    "SST",
    "SOX6",
    "LHX6",
    "ARX",
    "SATB1",
    "NPY",
    "EFNA5",
    "ELFN1",
    "CALB1",
    "NR2F2",
    "ADARB2",
    "CNR1",
    "ZBTB16",
    "PROX1",
    "RELN",
    "PENK"
  )

HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots <- 
  lapply(
    HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_genes, function(x) {
      HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap %>%
        filter(Gene == x) %>%
        ggplot(aes(Gene, ClusterMerged)) +
        geom_tile(aes(fill = Expression), color = "white")
    } +
      labs(fill = element_blank()) +
      theme(
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 76.8, face = "italic"),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 76.8),
        legend.key.size = unit(2, "lines"),
        legend.position = "bottom",
        legend.justification = "center",
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
  )

names(HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots) <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_genes

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_FOXP1 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$FOXP1 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 4),
    breaks = c(0, 4)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_FOXP2 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$FOXP2 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ZFHX4 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$ZFHX4 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_MEIS2 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$MEIS2 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 4),
    breaks = c(0, 4)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_PBX3 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$PBX3 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 4),
    breaks = c(0, 4)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_TSHZ1 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$TSHZ1 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_PVALB <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$PVALB +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ETV1 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$ETV1 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_BCAN <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$BCAN +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_TRPS1 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$TRPS1 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_KCNC1 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$KCNC1 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 4),
    breaks = c(0, 4)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_KCNC2 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$KCNC2 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_SST <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$SST +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_SOX6 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$SOX6 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_LHX6 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$LHX6 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ARX <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$ARX +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_SATB1 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$SATB1 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_NPY <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$NPY +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_EFNA5 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$EFNA5 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ELFN1 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$ELFN1 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_CALB1 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$CALB1 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_NR2F2 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$NR2F2 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ADARB2 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$ADARB2 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 5),
    breaks = c(0, 5)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_CNR1 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$CNR1 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ZBTB16 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$ZBTB16 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_PROX1 <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$PROX1 +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_RELN <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$RELN +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_PENK <-
  HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_plots$PENK +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  )

plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap <-
  ggarrange(
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_SIX3,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_FOXP1,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_FOXP2,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ZFHX4,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_MEIS2,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_PBX3,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_TSHZ1,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_PVALB,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ETV1,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_BCAN,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_TRPS1,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_KCNC1,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_KCNC2,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_SST,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_SOX6,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_LHX6,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ARX,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_SATB1,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_NPY,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_EFNA5,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ELFN1,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_CALB1,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_NR2F2,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ADARB2,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_CNR1,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_ZBTB16,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_PROX1,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_RELN,
    plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap_PENK,
    nrow = 1
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/revision/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_GABA_PRE_heatmap,
  width = 100,
  height = 30,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)
















HVCX_seurat_REGION_res1.1_heatmap2 <-
  HVCX_seurat_REGION_res1.1_datatable2 %>%
  select(
    ClusterMerged,
    SIX3,
    FOXP1,
    FOXP2,
    ZFHX4,
    MEIS2,
    PBX3,
    TSHZ1,
    PVALB,
    ETV1,
    BCAN,
    TRPS1,
    KCNC1,
    KCNC2,
    SST,
    SOX6,
    LHX6,
    ARX,
    SATB1,
    NPY,
    EFNA5,
    ELFN1,
    CALB1,
    NR2F2,
    ADARB2,
    CNR1,
    ZBTB16,
    PROX1,
    RELN,
    PENK
  ) %>%
  filter(
    ClusterMerged %in% c(
      2, 3, 9, 17, 34, 15, 38, 8, 14, 28, 37, 18, 19
    )
  ) %>%
  mutate(
    ClusterMerged = factor(
      ClusterMerged,
      levels = rev(c(
        2, 3, 9, 17, 34, 15, 38, 8, 14, 28, 37, 18, 19
      ))
    )
  ) %>%
  pivot_longer(
    cols = -1,
    names_to = "Gene", 
    values_to = "Expression"
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = c(
        "SIX3",
        "FOXP1",
        "FOXP2",
        "ZFHX4",
        "MEIS2",
        "PBX3",
        "TSHZ1",
        "PVALB",
        "ETV1",
        "BCAN",
        "TRPS1",
        "KCNC1",
        "KCNC2",
        "SST",
        "SOX6",
        "LHX6",
        "ARX",
        "SATB1",
        "NPY",
        "EFNA5",
        "ELFN1",
        "CALB1",
        "NR2F2",
        "ADARB2",
        "CNR1",
        "ZBTB16",
        "PROX1",
        "RELN",
        "PENK"
      )
    )
  ) %>%
  group_by(Gene, ClusterMerged) %>%
  mutate(
    Expression = mean(Expression),
  ) %>%
  distinct(ClusterMerged, .keep_all = TRUE)

plot_HVCX_seurat_REGION_res1.1_heatmap2 <-
  HVCX_seurat_REGION_res1.1_heatmap2 %>%
  ggplot(aes(Gene, ClusterMerged)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    name = "Normalized Expression",
    low = "white",  
    high = "steelblue",
    limits = c(0, 5),
    breaks = c(0, 5)
  ) +
  scale_y_discrete(
    labels = rev(c(
      "2 (X MSN)", 
      "3 (X MSN)",
      "9 (X MSN)",
      "17 (X MSN)",
      "34 (X PN)",
      "15 (HVC GABA-1)",
      "38 (RA GABA-1)",
      "8 (GABA-2/4)",
      "14 (GABA-3/4 & X INT)",
      "28 (GABA-5 & X INT)",
      "37 (GABA-6 & X INT)",
      "18 (GABA-7/8 & X INT)",
      "19 (GABA-Pre)"
    ))
  ) +
  theme(
    axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "top",
    legend.direction = "horizontal"
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "Postdoc/Projects/HVCX_COMBINED/revision/REGION/res1.1/",
    "plot_HVCX_seurat_REGION_res1.1_heatmap2.pdf"
  ),
  plot = plot_HVCX_seurat_REGION_res1.1_heatmap2,
  width = 20,
  height = 8,
  units = "in",
  dpi = 300
)





#HIERARCHICAL

#LOAD SEURAT OBJECT
load(
  file = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "postdoc/projects/HVCX_COMBINED/revision/REGION/res1.1/",
    "HVCX_seurat_REGION_res1.1.RData"
  )
)

HVCX_seurat_REGION_res1.1_markers <-
  read_csv(
    file.path(
      "~/OneDrive - University of Texas Southwestern/",
      "postdoc/projects/HVCX_COMBINED/revision/REGION/res1.1/",
      "HVCX_seurat_REGION_res1.1_markers.csv"
    )
  )

HVCX_seurat_REGION_res1.1_avexp <- 
  log1p(
    AverageExpression(
      HVCX_seurat_REGION_res1.1
    )
    $RNA
  ) %>%
  as_tibble(rownames = "gene") 

nsig <- 
  40

markers_int_top <-
  HVCX_seurat_REGION_res1.1_markers %>% 
  mutate(sign = avg_logFC > 0) %>%
  group_by(cluster, sign) %>%
  top_n(-1 * nsig, p_val_adj) %>%
  top_n(nsig, abs(avg_logFC)) %>%
  ungroup() %>% 
  distinct(gene, .keep_all = TRUE)

obj_int_avg_filt <-
  HVCX_seurat_REGION_res1.1_avexp %>%
  filter(gene %in% markers_int_top$gene)

test_df <-
  as.data.frame(obj_int_avg_filt)

test_matrix <-
  data.matrix(test_df)[,-1]

rownames(test_matrix) <-
  test_df$gene

library(colorspace)
library(dendextend)
library(pvclust)

pv <-
  pvclust(
    test_matrix, 
    method.dist = "correlation",
    method.hclust = "ward.D",
    nboot = 200, 
    parallel = 4L
  )

plot(pv)

text_colors <-
  qualitative_hcl(9)

dend <-
  as.dendrogram(pv)

dend <-
  dend %>%
  dendextend::set("branches_k_color", value = text_colors, k = 9) 

dend_gg <-
  as.ggdend(dend)

dend_seg <-
  dend_gg$segments

scale_factor = .05
buffer_size = .1
position_box_width = 1
dend_seg = dend_seg %>%
  mutate(y2 = max(y)-y,
         yend2 = max(yend) - yend) %>%
  mutate(y2 = y2 * scale_factor,
         yend2 = yend2 * scale_factor) %>%
  mutate(col = if_else(is.na(col), "grey50", col))

tree_width = max(dend_seg$yend2)

dend_leaves <- 
  dend_gg$labels %>%
  mutate(
    test = case_when(
      label == 1 ~ "1 (HVC Glut-1/4)",
      label == 5 ~ "5 (HVC Glut-2)",
      label == 7 ~ "7 (HVC Glut-3)",
      label == 35 ~ "35 (HVC Glut-5)",
      label == 25 ~ "25 (HVC Glut-1/2)",
      label == 12 ~ "12 (RA Glut-1/2/3 & HVC Glut-1)",
      label == 33 ~ "33 (X Glut+ INT)",
      label == 40 ~ "40 (X Glut+ Unknown)",
      label == 2 ~ "2 (X MSN)", 
      label == 3 ~ "3 (X MSN)",
      label == 9 ~ "9 (X MSN)",
      label == 17 ~ "17 (X MSN)",
      label == 34 ~ "34 (X PN)",
      label == 15 ~ "15 (HVC GABA-1)",
      label == 38 ~ "38  (RA GABA-1)",
      label == 8 ~ "8 (GABA-2/4)",
      label == 14 ~ "14 (GABA-3/4 & X INT)",
      label == 28 ~ "28 (GABA-5 & X INT)",
      label == 37 ~ "37 (GABA-6 & X INT)",
      label == 18 ~ "18 (GABA-7/8 & X INT)",
      label == 13 ~ "13 (Pre-1/Ependymal)",
      label == 11 ~ "11 (Pre-2)",
      label == 23 ~ "23 (Pre-3/4)",
      label == 19  ~ "19 (GABA-Pre)",
      label == 4 ~ "4 (HVC Astrocyte)",
      label == 6 ~ "6 (RA Astrocyte)",
      label == 16 ~ "16 (X Astrocyte)",
      label == 20 ~ "20 (HVC Oligodendrocyte)",
      label == 26 ~ "26 (RA Oligodendrocyte)",
      label == 39 ~ "39 (RA Oligodendrocyte)",
      label == 30 ~ "30 (X Oligodendrocyte)",
      label == 22 ~ "22 (HVC/X Oligo. Precursor)",
      label == 24 ~ "24 (RA Oligo. Precursor)",
      label == 21 ~ "21 (HVC/X Microglia)",
      label == 41 ~ "41 (RA Microglia)",
      label == 29 ~ "29 (HVC Microglia)",
      label == 10 ~ "10 (HVC/X Endothelial)",
      label == 27 ~ "27 (RA Endothelial)",
      label == 32 ~ "32 (X RBC)",
      label == 36 ~ "36 (HVC/RA RBC)",
      label == 31 ~ "32 (HVC/RA Mural)"
    )
  )

dend_leaves = dend_leaves %>% 
  mutate(pos_xmax = x+.4,
         pos_xmin = x-.4,
         perc_ymin1 = min(dend_seg$yend) - buffer_size,
         perc_ymax1 = perc_ymin1 - position_box_width,
         perc_ymin2 = perc_ymax1 - buffer_size,
         perc_ymax2 = perc_ymin2 - position_box_width,
         label_y = perc_ymax2 - buffer_size 
  )

flat_plot <-
  ggplot() +
  geom_segment(data = dend_seg,
               aes(x = -x,
                   xend = -xend,
                   y = -y,
                   yend = -yend,
                   color = col
               ),
               lineend = "square",
               size=1) + 
  coord_flip() +
  scale_color_identity() + 
  scale_fill_identity() +
  expand_limits(y=c(0, 15)) +
  theme_void() + 
  theme(legend.position = c(.1, .9)) +
  geom_text(data = dend_leaves,
            aes(x=-x,
                y=label_y,
                label=test),
            nudge_y = 2.5,
            hjust = 0
  )

ggsave(
  filename = file.path(
    "~/OneDrive - University of Texas Southwestern/",
    "postdoc/projects/HVCX_COMBINED/revision/REGION/res1.1/",
    "plot_tree.pdf"
  ),
  plot = flat_plot,
  width = 5,
  height = 8,
  units = "in",
  dpi = 300
)



