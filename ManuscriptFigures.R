## Explant Manuscript Figures

library(tidyverse)
library(Seurat)
library(sctransform)
library(scRepertoire)
library(ggpubr)
library(ggprism)
library(rstatix)
library(extrafont)

mean_sd_prop <- function(x) {
  df <- mean_sdl(x, mult = 1)
  
  if (df$ymin < 0) { df$ymin <- 1}
  if (df$ymax > 100) { df$ymax <- 99}
  
  return(df)
}

## CUSTOM GGPLOT THEME
GG_KM_THEME <- 
  list(
    ggprism::theme_prism(palette = "colorblind_safe", 
                         base_family = "sans",
                         base_size = 16, 
                         base_line_size = 1.5),
    ggprism::scale_color_prism(palette = "colorblind_safe", name = ""), 
    ggprism::scale_fill_prism(palette = "colorblind_safe", name = ""), 
    #theme(legend.position = "bottom"),
    theme(legend.position = "none", 
          plot.title = element_text(vjust = 0), 
          aspect.ratio = 2.5/1)
  )

data_dir <- "../data/"
figures_dir <- "../figures/"

human_gene_markers_yaml <- yaml::yaml.load("
t: 
  lineage: ['CD3', 'CD8A', 'CD4']
  naive: 
    surface: ['SELL', 'IL7R', 'CCR7', 'CXCR3']
    
  effector: 
    surface: ['IL2RA', 'TNFRSF8', 'CD44', 'CD69', 'IL2RB', 'TNFRSF4', 'LAG3', 'ICOS', 'KLRG1']
    intra: ['TBX21', 'PRDM1', 'ID2']
    cytokines: ['IFNG', 'IL2', 'PRF1', 'GZMA', 'GZMB', 'TNF', 'CCL3', 'CCL4', 'CCL5']
    
  effector_memory: 
    surface: ['CD44', 'KLRG1', 'B3GAT1']
    intra: ['EOMES', 'TBX21', 'PRDM1']
    cytokines: ['GZMB', 'IFNG', 'IL2', 'PRF1', 'TNF']
    low: ['SELL', 'CCR7']
  
  central_memory:
    surface: ['IL7R', 'CD27', 'CD28']
    intra: ['EOMES', 'TBX21']
    cytokines: ['IFNG', 'IL2', 'TNF']
    low: ['SELL', 'CCR7']
    
  resident_memory:
    surface: ['CD69', 'CXCR6', 'ITGAE', 'CTLA4']
  
  regulatory:
    surface: ['B3GAT1', 'CD28', 'KLRG1', 'LAG3', 'PDCD1', 'HLA-DRA']
    intra: ['FOXP3', 'IKZF1', 'EGR1', 'EGR2']
    cytokines: 
    low: ['IL2']
")


## Load in Single-cell Data
explantCLAD.with_trb <- readr::read_rds(paste0(data_dir, "explantCLAD.with_trb.rds"))

DefaultAssay(explantCLAD.with_trb) <- "SCT"
  
  
## READ IN FLOW DATA
flow_data_file <- "../CLAD Analysis V6.xls"

## Read in Flow Data
phenotype_flow_data <- readxl::read_xls(flow_data_file, sheet = "Phenotype") %>% 
  rename(sample_id = `...1`) %>%     # rename first column to "sample_id"
  drop_na(sample_id) %>% 
  filter(!grepl("WA", sample_id)) %>%   # Remove "WA" (warm autopsy) rows
  janitor::remove_empty("cols") %>%   # remove columns that only contain NA
  pivot_longer(!sample_id, names_to = "full_label", values_to = "value") %>% 
  mutate(value = ifelse(grepl("%", value),
                        readr::parse_number(value), 
                        readr::parse_number(value) * 100)) %>%
  mutate(panel = "phenotype")

function_pos_flow_data <- readxl::read_excel(flow_data_file, sheet = "Function") %>% 
  rename(sample_id = `...1`) %>%     # rename first column to "sample_id"
  drop_na(sample_id) %>% 
  filter(!grepl("WA", sample_id)) %>%   # Remove "WA" (warm autopsy) rows
  janitor::remove_empty("cols") %>%   # remove columns that only contain NA
  pivot_longer(!sample_id, names_to = "full_label", values_to = "value") %>% 
  mutate(value = ifelse(grepl("%", value),
                        readr::parse_number(value), 
                        readr::parse_number(value) * 100)) %>%
  mutate(panel = "function_pos")

function_diff_flow_data <- readxl::read_excel(flow_data_file, sheet = "Function w  neg samples") %>% 
  rename(sample_id = `...1`) %>%     # rename first column to "sample_id"
  drop_na(sample_id) %>% 
  filter(!grepl("WA", sample_id)) %>%   # Remove "WA" (warm autopsy) rows
  janitor::remove_empty("cols") %>%  # remove columns that only contain NA
  mutate(pos_neg = ifelse(grepl("neg", sample_id, ignore.case = TRUE), "neg", "pos")) %>% 
  mutate(sample_id = sub("[ _][PN][a-zA-Z]+\\.fcs", "", sample_id), 
         sample_id = gsub("[ -]", "_", sample_id)) %>% 
  pivot_longer(!c(sample_id, pos_neg), names_to = "full_label", values_to = "value") %>% 
  mutate(value = ifelse(grepl("%", value),
                        readr::parse_number(value), 
                        readr::parse_number(value) * 100)) %>%
  pivot_wider(names_from = pos_neg, values_from = value) %>% 
  drop_na(neg) %>% 
  mutate(value = pos - neg) %>% 
  select(-pos, -neg) %>% 
  mutate(panel = "function_diff") 

function_mfi_pos_flow_data <- readxl::read_excel(flow_data_file, sheet = "MFIs", na = c("", "n/a")) %>% 
  rename(sample_id = `...1`) %>%     # rename first column to "sample_id"
  drop_na(sample_id) %>% 
  filter(!grepl("WA", sample_id)) %>%   # Remove "WA" (warm autopsy) rows
  janitor::remove_empty("cols") %>%   # remove columns that only contain NA
  pivot_longer(!sample_id, names_to = "full_label", values_to = "value") %>% 
  filter(grepl("pos", sample_id, ignore.case = TRUE)) %>% 
  mutate(panel = "function_mfi_pos")

function_mfi_diff_flow_data <- readxl::read_excel(flow_data_file, sheet = "MFIs", na = c("", "n/a")) %>% 
  rename(sample_id = `...1`) %>%     # rename first column to "sample_id"
  drop_na(sample_id) %>% 
  filter(!grepl("WA", sample_id)) %>%   # Remove "WA" (warm autopsy) rows
  janitor::remove_empty("cols") %>%  # remove columns that only contain NA
  mutate(pos_neg = ifelse(grepl("neg", sample_id, ignore.case = TRUE), "neg", "pos")) %>% 
  mutate(sample_id = sub("[ _][PN][a-zA-Z]+\\.fcs", "", sample_id), 
         sample_id = gsub("[ -]", "_", sample_id)) %>% 
  pivot_longer(!c(sample_id, pos_neg), names_to = "full_label", values_to = "value") %>% 
  pivot_wider(names_from = pos_neg, values_from = value) %>% 
  drop_na(neg) %>% 
  mutate(value = pos - neg) %>% 
  select(-pos, -neg) %>% 
  mutate(panel = "function_mfi_diff") 


full_flow_data <- bind_rows(phenotype_flow_data, 
                            function_pos_flow_data, 
                            #function_diff_flow_data, 
                            function_mfi_pos_flow_data, 
                            #function_mfi_diff_flow_data
) %>% 
  mutate(disease = factor(case_when(
    grepl("IPF", sample_id, ignore.case = TRUE) ~ "IPF", 
    grepl("_C", sample_id, ignore.case = TRUE) | 
      grepl("C_", sample_id, ignore.case = TRUE) | 
      grepl("-C", sample_id, ignore.case = TRUE)~ "Control", 
    .default = "CLAD"
  ), levels = c("Control", "CLAD", "IPF")), 
  organ = factor(case_when(
    grepl("HLN", sample_id, ignore.case = TRUE) ~ "HLN", 
    .default = "Lung"
  ), levels = c("Lung", "HLN")), 
  full_label = sub("rla", "rHLA", full_label, ignore.case = TRUE),   # fix a typo in the file
  marker = str_squish(str_split_i(full_label, "of all|of", 1)),
  parent = str_squish(str_split_i(full_label, "of all|of", 2))
  ) %>% 
  select(panel, disease, organ, parent, marker, value, sample_id, full_label)

## Add CD4/CD8 Ratio
full_flow_data <- bind_rows(full_flow_data, 
                            full_flow_data %>% filter("phenotype" == panel, 
                                                      "CD3+" == parent,
                                                      "CD4+" == marker | "CD8+" == marker) %>% 
                              select(-full_label) %>% 
                              pivot_wider(names_from = marker, values_from = value) %>% 
                              mutate(marker = "CD4/CD8", 
                                     value = `CD4+` / `CD8+`, 
                                     full_label = "CD4/CD8 of CD3+") %>% 
                              select(-`CD4+`, -`CD8+`))

## Calculate the pairwise Wilcoxons and p-values
# https://www.datanovia.com/en/blog/how-to-add-p-values-onto-basic-ggplots/
cell_prop_wilcox <- 
  full_flow_data %>% 
  select(-sample_id, -full_label) %>% 
  group_by(panel, parent, marker, organ) %>% 
  drop_na(value) %>% 
  wilcox_test(value ~ disease) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p") %>% 
  add_significance("p.adj") %>% 
  #add_xy_position(x = "organ", dodge = 0.9, fun = "mean_sd", step.increase = 0.03) %>% 
  add_xy_position(x = "organ", dodge = 0.9, fun = "mean_sd", scales = "free", step.increase = 0.03) %>% 
  mutate(p.format = p_format(p, accuracy = 0.01, leading.zero = FALSE), 
         p.adj.format = p_format(p.adj, accuracy = 0.01, leading.zero = FALSE), 
         disease = group1)

## Helper function to ggplot flow data
flow_plot <- function(flow_data, 
                      wilcox_data, 
                      fig_panel, 
                      fig_parent, 
                      fig_marker, 
                      fig_y = NULL, 
                      fig_title = NULL, 
                      y_max = 100, y_steps = 25) { 
  
  w <- 0.9
  
  fig_y <- ifelse(is.null(fig_y), 
                  paste0("Mean % ", fig_marker), 
                  fig_y)
  
  fig_title <- ifelse(is.null(fig_title), 
                      paste0(fig_marker, " of ", fig_parent), 
                      fig_title)
  
  y_comparison_max <- max(y_max,
                          wilcox_data %>% 
                            filter(fig_panel == panel, 
                                   fig_parent == parent,
                                   fig_marker == marker) %>% 
                            pull(y.position))
  
  p <- full_flow_data %>% filter(fig_panel == panel, 
                                 fig_parent == parent,
                                 fig_marker == marker) %>% 
    drop_na(value) %>% 
    ggplot(aes(x = organ, y = value, fill = disease)) +
    stat_summary(geom = "bar", fun = mean, position = "dodge", width = w, alpha = 0.5) +
    stat_summary(geom = "errorbar", fun.data = mean_sdl, fun.args = list(mult = 1),  
                 position = position_dodge(width = w), 
                 width = 0.2, color = "darkslategray") + 
    geom_point(aes(color = disease), position = position_jitterdodge(jitter.width = 0.2, dodge.width = w)) + 
    stat_pvalue_manual(data = wilcox_data %>% 
                         filter(fig_panel == panel, 
                                fig_parent == parent,
                                fig_marker == marker), 
                       hide.ns = "p.adj", 
                       tip.length = 0, 
                       step.increase = 0.05,
                       step.group.by = "organ",
                       label = "{p.adj.signif}", 
                       label.size = 4, 
                       bracket.size = 1) + 
    xlab("") + ylab(fig_y) + 
    ggtitle(fig_title) + 
    #scale_y_continuous(breaks = seq(0, y_max, y_steps), limits = c(0, y_comparison_max)) + 
    GG_KM_THEME
  
  return(p)
}

## Old full_label parsing
# mini_label = sub("\nLymphocytes/Single Cells/Live", "", full_label), 
# mini_label = sub(" \\| Freq.*", "", mini_label), 
# mini_label = sub("^[^/]+/", "", mini_label), 
# mini_label = sub("T cells \\(Live ", "", mini_label),
# mini_label = sub("(All)", "", mini_label), 
# mini_label = gsub("TRM[ ]*", "", mini_label), 
# mini_label = sub("Treg[ ]*", "", mini_label),
# mini_label = gsub("[()]", "", mini_label), 
# mini_label = gsub(" /", "/", mini_label), 
# mini_label = sub("[ ]+$", "", mini_label), 
# mini_label = sub("Non-recipient.*", "rHLA-", mini_label), 
# mini_label = sub("Recipient.*", "rHLA+", mini_label), 
# marker = stringr::str_split_i(mini_label, "/", -1), 
# parent = stringr::str_split_i(mini_label, "/", -2)


# Aggregate to mean and median proportions
# flow_data_mean_props <- full_flow_data %>%
#   group_by(panel, disease, organ, parent, marker) %>%
#   summarise(mean = mean(value),
#             median = median(value),
#             sd = sd(value, na.rm = TRUE))

#### LEGENDPLEX FIGURES

## READ IN LEGENDPLEX DATA
lgplex_data_file <- "~/Library/CloudStorage/Dropbox/Project - 2022 - CLAD Explant/Data/CLAD legendplex data.xlsx"

lgplex_data <- readxl::read_xlsx(lgplex_data_file) %>% 
  rename(sample_id = well) %>%     # rename first column to "sample_id"
  drop_na(sample_id) %>%     # rename first column to "sample_id"
  janitor::remove_empty("cols") %>%   # remove columns that only contain NA
  pivot_longer(!sample_id, names_to = "molecule", values_to = "value") %>% 
  mutate(value = as.numeric(value), 
         sample_id = sub("\\.fcs", "", sample_id), 
         molecule = sub(" \\([A-Za-z0-9]+\\)", "", molecule), 
         molecule = ifelse(grepl("tnf", molecule, ignore.case = TRUE), "TNFa", molecule), 
         molecule = ifelse(grepl("ifn", molecule, ignore.case = TRUE), "IFNg", molecule), 
         treatment = str_to_title(str_extract(sample_id, "_\\w+$")), 
         treatment = sub("_", "", treatment), 
         sample_id = sub("_\\w+$", "", sample_id), 
         treatment = ifelse("Mica" == treatment, "MicA", treatment), 
         disease = ifelse(grepl("-LT$", sample_id, ignore.case = TRUE), "CLAD", "Control")) %>% 
  filter("Stim" != treatment) %>%   # No need to include the "Stim" data
  mutate(disease = factor(disease, levels = c("Control", "CLAD")), 
         treatment = factor(treatment, levels = c("Unstim", "MicA", "Block")), 
         molecule = factor(molecule))


lgplex_wilcox <- 
  lgplex_data %>% 
  select(-sample_id) %>% 
  group_by(disease, molecule) %>% 
  drop_na(value) %>% 
  wilcox_test(value ~ treatment) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p") %>% 
  add_significance("p.adj") %>% 
  add_xy_position(x = "disease", dodge = 0.9, fun = "mean_sd", step.increase = 0.03) %>% 
  #add_xy_position(x = "organ", dodge = 0.9, fun = "mean_sd", scales = "free") %>% 
  mutate(p.format = p_format(p, accuracy = 0.01, leading.zero = FALSE), 
         p.adj.format = p_format(p.adj, accuracy = 0.01, leading.zero = FALSE), 
         treatment = group1)


lgplex_main <- c("Granzyme B", "IL-2", "IFNg")
lgplex_supplement <- as.character(unique(lgplex_data$molecule)[!(unique(lgplex_data$molecule) %in% lgplex_main)])

w <- 0.9

lgplex_plots_main <- list()
lgplex_plots_supp <- list()

for(curr_molecule in c(lgplex_main, lgplex_supplement)) {
  
  plot <- 
  lgplex_data %>% 
    filter(molecule == curr_molecule) %>% 
    #mutate(molecule = factor(molecule, levels = lgplex_main)) %>% 
    drop_na(value) %>% 
    ggplot(aes(x = treatment, y = value, fill = disease)) +
    stat_summary(geom = "bar", fun = mean, position = "dodge", width = w, alpha = 0.5) +
    stat_summary(geom = "errorbar", fun.data = mean_sdl, fun.args = list(mult = 1),  
                 position = position_dodge(width = w), 
                 width = 0.2, color = "darkslategray") + 
    geom_point(aes(color = disease), position = position_jitterdodge(jitter.width = 0.2, dodge.width = w)) + 
    facet_grid(~ disease, scales = "free", switch = "x") + 
    stat_pvalue_manual(data = lgplex_wilcox %>% filter(molecule == curr_molecule),
                       hide.ns = "p.adj",
                       tip.length = 0,
                       #step.increase = 0.05,
                       #step.group.by = "organ",
                       label = "{p.signif}",
                       label.size = 4,
                       bracket.size = 1) +
    xlab("") + ylab("pg/mL") + 
    ggtitle(curr_molecule) + 
    GG_KM_THEME + theme(legend.position = "none", aspect.ratio = 1, 
                        axis.title.x = element_blank(), 
                        strip.placement = "outside", 
                        strip.background.x = element_rect(color = NA, fill = NA), 
                        strip.text = element_text(face = "bold", size = 18),
                        axis.text.x = element_text(size = 14, angle = 0)
                        #axis.text.y = element_text(face = "plain", size = 8, angle = 0), 
    )
  #scale_y_continuous(breaks = seq(0, y_max, y_steps), limits = c(0, y_comparison_max)) + 
  #GG_KM_THEME
  
  if (curr_molecule %in% lgplex_main) {
    lgplex_plots_main[[curr_molecule]] <- plot
  } else {
    lgplex_plots_supp[[curr_molecule]] <- plot
  }
}

## Main Figure Plot
#lgplex_plots_main[[2]] <- lgplex_plots_main[[2]] + theme(legend.position = "bottom")
cowplot::plot_grid(plotlist = lgplex_plots_main[1:3], ncol = 3, align = "h", axis = "b")
ggsave(paste0(figures_dir, paste0("lgplex_main.pdf")), plot = last_plot(), width = 16, height = 7)

## Supplement Figure Plot
#lgplex_plots_supp[[8]] <- lgplex_plots_supp[[8]] + theme(legend.position = "bottom")
cowplot::plot_grid(plotlist = lgplex_plots_supp, ncol = 3, align = "h", axis = "b")
ggsave(paste0(figures_dir, paste0("lgplex_supp.pdf")), plot = last_plot(), width = 16, height = 16)



#####---------##### START FIGURES CODE #####---------#####

#####---------##### PROGRESS MARKER #####---------#####

# Figure 1 - T cell phenotypes/memory subsets
## Donor vs Recipient Stromal Cells
flow_plot(full_flow_data, cell_prop_wilcox, 
          fig_panel = "phenotype", 
          fig_parent = "CD3- CD31+", 
          fig_marker = "rHLA+", 
          fig_y = "Mean % of Recipient Origin")

ggsave("flow_recip_stromal.pdf", plot = last_plot(), path = figures_dir, width = 5, height = 5, units = "in")

## Donor vs Recipient T cells
flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "phenotype", 
               fig_parent = "CD3+", 
               fig_marker = "rHLA+", 
               fig_y = "Mean % of Recipient Origin")

ggsave("flow_recip_t_cells.pdf", plot = last_plot(), path = figures_dir, width = 5, height = 5, units = "in")

## CD4, CD8
flow_plot(full_flow_data, cell_prop_wilcox, 
          fig_panel = "phenotype", 
          fig_parent = "CD3+", 
          fig_marker = "CD4+", 
          fig_y = "Mean % CD4+")

ggsave("flow_cd4_t_cells.pdf", plot = last_plot(), path = figures_dir, width = 5, height = 5, units = "in")

flow_plot(full_flow_data, cell_prop_wilcox,
          fig_panel = "phenotype", 
          fig_parent = "CD3+", 
          fig_marker = "CD8+", 
          fig_y = "Mean % CD8+")

ggsave("flow_cd8_t_cells.pdf", plot = last_plot(), path = figures_dir, width = 5, height = 5, units = "in")

## Memory Subsets
full_flow_data %>% filter(organ == "Lung",
                          panel == "phenotype", 
                          parent == "CD3+ CD4+",
                          marker %in% c("CCR7+ CD45RA+", "CCR7+ CD45RA-", "CCR7- CD45RA+", "CCR7- CD45RA-")) %>% 
  drop_na(value) %>% 
  ggplot(aes(x = disease, y = value, fill = marker)) +
  stat_summary(geom = "bar", fun = mean, position = "stack", width = 0.9, alpha = 0.8) +
  xlab("") + ylab("Mean Proportion of CD3+ CD4+") + 
  ggtitle("Lung CD4") + 
  GG_KM_THEME + 
  scale_fill_viridis_d(labels = c("CCR7+ CD45RA+" = "Naive", 
                                  "CCR7+ CD45RA-" = "TCM", 
                                  "CCR7- CD45RA+" = "TEMRA", 
                                  "CCR7- CD45RA-" = "TEM")) + 
  theme(legend.position = "right", axis.text.x = element_text(angle = 45), 
        aspect.ratio = 2.5,
        axis.title.x = element_blank(), 
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        legend.title = element_blank(),
        plot.margin = grid::unit(c(0,0,0,0),"mm")) 

ggsave("flow_lung_cd4_mem_subsets.pdf", plot = last_plot(), path = figures_dir, width = 4, height = 4, units = "in")

full_flow_data %>% filter(organ == "Lung",
                          panel == "phenotype", 
                          parent == "CD3+ CD8+",
                          marker %in% c("CCR7+ CD45RA+", "CCR7+ CD45RA-", "CCR7- CD45RA+", "CCR7- CD45RA-")) %>% 
  drop_na(value) %>% 
  ggplot(aes(x = disease, y = value, fill = marker)) +
  stat_summary(geom = "bar", fun = mean, position = "stack", width = 0.9, alpha = 0.8) +
  xlab("") + ylab("Mean Proportion of CD3+ CD8+") + 
  ggtitle("Lung CD8") + 
  GG_KM_THEME + 
  scale_fill_viridis_d(labels = c("CCR7+ CD45RA+" = "Naive", 
                                  "CCR7+ CD45RA-" = "TCM", 
                                  "CCR7- CD45RA+" = "TEMRA", 
                                  "CCR7- CD45RA-" = "TEM")) + 
  theme(legend.position = "right", axis.text.x = element_text(angle = 45), 
        aspect.ratio = 2.5,
        axis.title.x = element_blank(), 
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        legend.title = element_blank(),
        plot.margin = grid::unit(c(0,0,0,0),"mm")) 

ggsave("flow_lung_cd8_mem_subsets.pdf", plot = last_plot(), path = figures_dir, width = 4, height = 4, units = "in")

full_flow_data %>% filter(organ == "HLN",
                          panel == "phenotype", 
                          parent == "CD3+ CD4+",
                          marker %in% c("CCR7+ CD45RA+", "CCR7+ CD45RA-", "CCR7- CD45RA+", "CCR7- CD45RA-")) %>% 
  drop_na(value) %>% 
  ggplot(aes(x = disease, y = value, fill = marker)) +
  stat_summary(geom = "bar", fun = mean, position = "stack", width = 0.9, alpha = 0.8) +
  xlab("") + ylab("Mean Proportion of CD3+ CD4+") + 
  ggtitle("HLN CD4") + 
  GG_KM_THEME + 
  scale_fill_viridis_d(labels = c("CCR7+ CD45RA+" = "Naive", 
                                  "CCR7+ CD45RA-" = "TCM", 
                                  "CCR7- CD45RA+" = "TEMRA", 
                                  "CCR7- CD45RA-" = "TEM")) + 
  theme(legend.position = "right", axis.text.x = element_text(angle = 45), 
        aspect.ratio = 2.5,
        axis.title.x = element_blank(), 
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        legend.title = element_blank(),
        plot.margin = grid::unit(c(0,0,0,0),"mm")) 

ggsave("flow_hln_cd4_mem_subsets.pdf", plot = last_plot(), path = figures_dir, width = 4, height = 4, units = "in")

full_flow_data %>% filter(organ == "HLN",
                          panel == "phenotype", 
                          parent == "CD3+ CD8+",
                          marker %in% c("CCR7+ CD45RA+", "CCR7+ CD45RA-", "CCR7- CD45RA+", "CCR7- CD45RA-")) %>% 
  drop_na(value) %>% 
  ggplot(aes(x = disease, y = value, fill = marker)) +
  stat_summary(geom = "bar", fun = mean, position = "stack", width = 0.9, alpha = 0.8) +
  xlab("") + ylab("Mean Proportion of CD3+ CD8+") + 
  ggtitle("HLN CD8") + 
  GG_KM_THEME + 
  scale_fill_viridis_d(labels = c("CCR7+ CD45RA+" = "Naive", 
                                  "CCR7+ CD45RA-" = "TCM", 
                                  "CCR7- CD45RA+" = "TEMRA", 
                                  "CCR7- CD45RA-" = "TEM")) + 
  theme(legend.position = "right", axis.text.x = element_text(angle = 45), 
        aspect.ratio = 2.5,
        axis.title.x = element_blank(), 
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        legend.title = element_blank(),
        plot.margin = grid::unit(c(0,0,0,0),"mm")) 

ggsave("flow_hln_cd8_mem_subsets.pdf", plot = last_plot(), path = figures_dir, width = 4, height = 4, units = "in")


# Memory Subset statistical testing
cell_prop_wilcox %>% filter(panel == "phenotype", 
                            parent == "CD3+ CD4+", 
                            marker %in% c("CCR7+ CD45RA+", "CCR7+ CD45RA-", "CCR7- CD45RA+", "CCR7- CD45RA-"), 
                            p.adj <= 0.05) %>% 
  select(organ, parent, marker, group1, group2, p.adj)

cell_prop_wilcox %>% filter(panel == "phenotype", 
                            parent == "CD3+ CD8+", 
                            marker %in% c("CCR7+ CD45RA+", "CCR7+ CD45RA-", "CCR7- CD45RA+", "CCR7- CD45RA-"), 
                            p.adj <= 0.05) %>% 
  select(organ, parent, marker, group1, group2, p.adj)


## Clonal Overlap, Expansion, Etc (scRepertoire)

### Clonal Overlap
clonalOverlap(explantCLAD.with_trb, 
              cloneCall = "aa", 
              method = "overlap",
              chain = "TRB", 
              split.by = "organ.disease"
              ) + 
  scale_fill_viridis_c(option = "cividis", alpha = 0.75) + 
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.text.x = element_text(angle = 30, vjust = 0.5))

ggsave("fig_clonal_overlap.pdf", 
       plot = last_plot(),
       path = figures_dir, 
       width = 7, height = 5, units = "in", dpi = 300)


## Find CLAD HLN/Lung Overlapping Clones
clad_lung_TRBs <- explantCLAD.with_trb@meta.data %>% filter(organ.disease == "Lung_CLAD") %>% pull(TRBaa) %>% unique()
clad_hln_TRBs <- explantCLAD.with_trb@meta.data %>% filter(organ.disease == "HLN_CLAD") %>% pull(TRBaa) %>% unique()

clad_overlap_TRBs <- intersect(clad_lung_TRBs, clad_hln_TRBs)

clad_lung_overlap_barcodes <- explantCLAD.with_trb@meta.data %>% 
  filter(organ.disease == "Lung_CLAD", 
         TRBaa %in% clad_overlap_TRBs) %>% 
  rownames()

clad_hln_overlap_barcodes <- explantCLAD.with_trb@meta.data %>% 
  filter(organ.disease == "HLN_CLAD", 
         TRBaa %in% clad_overlap_TRBs) %>% 
  rownames()


Idents(explantCLAD.with_trb) <- explantCLAD.with_trb$organ.disease

# Lung
explantCLAD.with_trb.only_clad_lung = subset(explantCLAD.with_trb, idents = c("Lung_CLAD"))

explantCLAD.with_trb.only_clad_lung <- highlightClonotypes(explantCLAD.with_trb.only_clad_lung,
                                                           cloneCall = "TRBaa",
                                                           sequence = clad_overlap_TRBs)

explantCLAD.with_trb.only_clad_lung@meta.data$highlight <- ifelse(is.na(explantCLAD.with_trb.only_clad_lung@meta.data$highlight), 
                                                                  "Non-overlapped\nClones", "Overlapped\nClones")

Idents(explantCLAD.with_trb.only_clad_lung) <- explantCLAD.with_trb.only_clad_lung$highlight

overlapped_clones.dge.lung <- FindMarkers(explantCLAD.with_trb.only_clad_lung, 
                                          group.by = "highlight", 
                                          ident.1 = "Overlapped\nClones", 
                                          ident.2 = "Non-overlapped\nClones", 
                                          features = cytotox_genes_overlapped_clones, 
                                          recorrect_umi = FALSE
                                          )

# HLN
explantCLAD.with_trb.only_clad_hln = subset(explantCLAD.with_trb, idents = c("HLN_CLAD"))

explantCLAD.with_trb.only_clad_hln <- highlightClonotypes(explantCLAD.with_trb.only_clad_hln,
                                                           cloneCall = "TRBaa",
                                                           sequence = clad_overlap_TRBs)

explantCLAD.with_trb.only_clad_hln@meta.data$highlight <- ifelse(is.na(explantCLAD.with_trb.only_clad_hln@meta.data$highlight), 
                                                                 "Non-overlapped\nClones", "Overlapped\nClones")

Idents(explantCLAD.with_trb.only_clad_hln) <- explantCLAD.with_trb.only_clad_hln$highlight


clad_overlap_lung_vln <- 
  VlnPlot(explantCLAD.with_trb.only_clad_lung, 
        features = cytotox_genes_overlapped_clones,
        pt.size = 0, 
        stack = TRUE, flip = TRUE, 
        cols = viridis::viridis(n = 2),
        #group.by = "sample_type", 
        split.by = "highlight") + 
  NoLegend() + 
  ggtitle("CLAD Lung") + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 0, face = "bold", size = 14, hjust = 0.5))


clad_overlap_hln_vln <- 
  VlnPlot(explantCLAD.with_trb.only_clad_hln, 
          features = cytotox_genes_overlapped_clones,
          pt.size = 0, 
          stack = TRUE, flip = TRUE, 
          cols = viridis::viridis(n = 2),
          #group.by = "sample_type", 
          split.by = "highlight") + 
  NoLegend() + 
  ggtitle("Clad HLN") + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 0, face = "bold", size = 14, hjust = 0.5))

cowplot::plot_grid(clad_overlap_lung_vln, clad_overlap_hln_vln)

ggsave("overlap_cytotox.pdf", plot = last_plot(), path = figures_dir, width = 10, height = 7, units = "in")

# Idents(explantCLAD.with_trb) <- explantCLAD.with_trb$organ.disease
# 
# compareClonotypes(explantCLAD.with_trb, 
#                   numbers = 100, 
#                   samples = c("Lung_CLAD", "HLN_CLAD"), 
#                   cloneCall="aa", chain = "TRB",
#                   split.by = "ident",
#                   graph = "alluvial")


### Clonal Topo 
Idents(explantCLAD.with_trb) <- explantCLAD.with_trb$phenotype

clonalOverlay(explantCLAD.with_trb, 
                   reduction = "umap", 
                   freq.cutpoint = 5, 
                   bins = 10, 
                   facet = "disease") + 
  #ggprism::scale_fill_prism(palette = "colorblind_safe", name = "") + 
  #ggprism::scale_color_prism(palette = "colorblind_safe", name = "") + 
  scale_color_viridis_d() + 
  #ggpubr::theme_pubr() + 
  ggprism::theme_prism() + 
  ggtitle("T cell Clonal Density by Phenotype") + 
  theme(legend.position = "bottom", 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 21)) + 
  guides(color = guide_legend(override.aes = list(size=5)))


ggsave("fig_clonal_topo.pdf", 
       plot = last_plot(),
       path = figures_dir, 
       width = 7, height = 5, units = "in", dpi = 300)

### Expansion by Group
explantCLAD.with_trb$disease_organ_phenotype <- 
  factor(paste(explantCLAD.with_trb$disease, 
               explantCLAD.with_trb$sample_type, 
               explantCLAD.with_trb$phenotype,
               sep = "_"), 
         levels = c("Control_Lung_CD4 Naive", "Control_Lung_CD8 Effector", "Control_Lung_CD4 TRM", "Control_Lung_CD8 TRM",
                    "CLAD_Lung_CD4 Naive", "CLAD_Lung_CD8 Effector", "CLAD_Lung_CD4 TRM", "CLAD_Lung_CD8 TRM", 
                    "IPF_Lung_CD4 Naive", "IPF_Lung_CD8 Effector", "IPF_Lung_CD4 TRM", "IPF_Lung_CD8 TRM", 
                    "Control_HLN_CD4 Naive", "Control_HLN_CD8 Effector", "Control_HLN_CD4 TRM", "Control_HLN_CD8 TRM", 
                    "CLAD_HLN_CD4 Naive", "CLAD_HLN_CD8 Effector", "CLAD_HLN_CD4 TRM", "CLAD_HLN_CD8 TRM", 
                    "IPF_HLN_CD4 Naive", "IPF_HLN_CD8 Effector", "IPF_HLN_CD4 TRM", "IPF_HLN_CD8 TRM"))

Idents(explantCLAD.with_trb) <- explantCLAD.with_trb$sample_type
explantCLAD.with_trb.lung <- subset(explantCLAD.with_trb, idents = "Lung")

Idents(explantCLAD.with_trb.lung) <- explantCLAD.with_trb$disease

clonalHomeostasis(explantCLAD.with_trb.lung, 
                  cloneCall = "aa", chain = "TRB", 
                  cloneTypes = c(Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1), 
                  split.by = "disease") + 
  GG_KM_THEME + 
  scale_x_discrete(limits = c("Control", "CLAD", "IPF")) + 
  ggtitle("Clonal Expansion (Lung)") + 
  theme(axis.title.x = element_blank(), 
        legend.position = "right") + 
  scale_fill_viridis_d()

ggsave("fig_clonal_expansion_lung.pdf", 
       plot = last_plot(),
       path = figures_dir, 
       width = 7, height = 9, units = "in", dpi = 300)


Idents(explantCLAD.with_trb) <- explantCLAD.with_trb$sample_type
explantCLAD.with_trb.hln <- subset(explantCLAD.with_trb, idents = "HLN")

Idents(explantCLAD.with_trb.hln) <- explantCLAD.with_trb$disease

clonalHomeostasis(explantCLAD.with_trb.hln, 
                  cloneCall = "aa", chain = "TRB", 
                  cloneTypes = c(Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1), 
                  split.by = "disease") + 
  GG_KM_THEME + 
  scale_x_discrete(limits = c("Control", "CLAD", "IPF")) + 
  ggtitle("Clonal Expansion (HLN)") + 
  theme(axis.title.x = element_blank(), 
        legend.position = "right") + 
  scale_fill_viridis_d()

ggsave("fig_clonal_expansion_hln.pdf", 
       plot = last_plot(),
       path = figures_dir, 
       width = 7, height = 9, units = "in", dpi = 300)

# clonalProportion(explantCLAD.with_trb, 
#                   cloneCall = "aa", chain = "TRB", 
#                   split = c(10, 100, 1000, 10000, 30000, 1e+05), 
#                   split.by = "disease_organ_phenotype") + 
#   ggpubr::theme_pubr() + 
#   ggtitle("T cell Clonal Expansion") + 
#   theme(legend.position = "bottom", 
#         axis.text = element_text(face = "bold", size = 14), 
#         axis.text.x = element_text(angle = 90, vjust = 0.5), 
#         plot.title = element_text(face = "bold", hjust = 0.5, size = 21)) + 
#   guides(color = guide_legend(override.aes = list(size=5)))

# alluvialClonotypes(explantCLAD.with_trb, 
#                    cloneCall = "aa", chain = "TRB", 
#                    y.axes = c("sample_type", "disease", "phenotype"), 
#                    color = "disease") 

### Clonal Diversity
explantCLAD.with_trb$disease <- factor(explantCLAD.with_trb$disease, 
                                       levels = c("Control", "CLAD", "IPF"))

clonalDiversity(explantCLAD.with_trb, 
                cloneCall = "aa", 
                chain = "TRB", 
                split.by = "disease")

clonalDiversity(explantCLAD.with_trb, 
                cloneCall = "aa", 
                chain = "TRB", 
                exportTable = TRUE,
                split.by = "organ.disease") %>% 
  as_tibble() %>% 
  mutate(disease = case_match(Group, 
                              c("HLN_Control", "Lung_Control") ~ "Control", 
                              c("HLN_CLAD", "Lung_CLAD") ~ "CLAD", 
                              c("HLN_IPF", "Lung_IPF") ~ "IPF"), 
         organ = case_match(Group, 
                            c("Lung_Control", "Lung_CLAD", "Lung_IPF") ~ "Lung", 
                            c("HLN_Control", "HLN_CLAD", "HLN_IPF") ~ "HLN")) %>% 
  mutate(disease = factor(disease, levels = c("Control", "CLAD", "IPF")), 
         organ = factor(organ, levels = c("Lung", "HLN"))) %>% 
  {ggplot(., aes(x = organ, y = Shannon, fill = disease)) +
      geom_col(position = "dodge", width = 0.5, alpha = 0.5) +
      ggtitle("Clonal Diversity") + 
      coord_cartesian(ylim = c(min(.$Shannon), max(.$Shannon))) +
      ylab("Shannon Diversity Inddex") + 
      GG_KM_THEME + 
      theme(legend.position = "bottom", 
            axis.title.x = element_blank()) + 
      guides(color = guide_legend(override.aes = list(size=5)))}

ggsave("plot_clonal_diversity.pdf", plot = last_plot(), path = figures_dir, width = 5, height = 5, units = "in", dpi = 300)


clonalDiversity(explantCLAD.with_trb.lung, 
                cloneCall = "aa", 
                chain = "TRB", 
                exportTable = TRUE,
                split.by = "disease") %>% 
  as_tibble() %>% 
  mutate(disease = factor(Group,  levels = c("Control", "CLAD", "IPF"))) %>% 
  {ggplot(., aes(x = disease, y = Shannon, fill = disease)) +
      geom_col(position = "dodge", width = 0.5, alpha = 1) +
      ggtitle("Clonal Diversity") + 
      coord_cartesian(ylim = c(min(.$Shannon), max(.$Shannon))) +
      ylab("Shannon Diversity Inddex") + 
      GG_KM_THEME + 
      theme(legend.position = "none", 
            aspect.ratio = 1,
            axis.title.x = element_blank())}

ggsave("plot_clonal_diversity_lung.pdf", plot = last_plot(), path = figures_dir, width = 5, height = 5, units = "in", dpi = 300)


#### UMAPS / Feature Plots
Idents(explantCLAD.with_trb) <- explantCLAD.with_trb$seurat_cluster

FeaturePlot(explantCLAD.with_trb, 
            pt.size = 1.4, label = FALSE,
            features = c("CD4", "CD8A", "CD69", "ITGAE", "SELL", "S1PR1"))

ggsave("feature_plot.pdf", plot = last_plot(), path = figures_dir, width = 12, height = 12, units = "in", dpi = 300)

DimPlot(explantCLAD.with_trb, group.by = "seurat_cluster", label.box = TRUE, label = TRUE)

explantCLAD.with_trb$phenotype <- factor(
  case_match(as.numeric(explantCLAD.with_trb$seurat_cluster) - 1, 
             c(0, 3, 8) ~ "CD4 EM", 
             c(2) ~ "CD8 TRM",
             c(5, 6) ~ "CD8 EM", 
             .default = "CD4 TRM"
             ))

DimPlot(explantCLAD.with_trb, 
        group.by = "phenotype",
        split.by = "sample_type", 
        pt.size = 1.4, order = TRUE) + 
  ggtitle("") + 
  theme(legend.position = "bottom", legend.justification = "center") + 
  scale_color_viridis_d()

ggsave("plot_split_umap_phenotype.pdf", plot = last_plot(), path = figures_dir, width = 10, height = 7, units = "in", dpi = 300)

#VlnPlot(explantCLAD.with_trb, group.by = "seurat_cluster", features = c("CD4", "CD8A"))

## Figure 2 - TRMs

# TRMs
p <- flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "phenotype", 
               fig_parent = "CD3+", 
               fig_marker = "CD69+") + 
  scale_y_continuous(breaks = seq(0, 100, 25))

ggsave("flow_trm_t_cells.pdf", plot = p, path = figures_dir, width = 5, height = 5, units = "in")

# CD4 TRMs
flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "phenotype", 
               fig_parent = "CD3+ CD4+", 
               fig_marker = "CD69+") + 
  scale_y_continuous(breaks = seq(0, 100, 25))

ggsave("flow_trm_cd4.pdf", plot = last_plot(), path = figures_dir, width = 5, height = 5, units = "in")

# CD8 TRMs
flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "phenotype", 
               fig_parent = "CD3+ CD8+", 
               fig_marker = "CD69+") + 
  scale_y_continuous(breaks = seq(0, 100, 25))

ggsave("flow_trm_cd8.pdf", plot = last_plot(), path = figures_dir, width = 5, height = 5, units = "in")


## Figure 2 - Cytoxicity


cytotox_genes = c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "NKG7", "KLRK1", "IFNG", 
                  "KLRD1", "ZNF683", "HOPX", "CCL5", "CD74", "LAMP1", "TNF", "IL7R", 
                  "S1PR1", "CD27", "SELL", "CCR7")

cytotox_genes_overlapped_clones = c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "NKG7", "KLRK1", "CCL5", "CD74")


## Volcano Plots
Idents(explantCLAD.with_trb) <- explantCLAD.with_trb$organ.disease

dge.lungs.clad_v_control <- FindMarkers(explantCLAD.with_trb, ident.1 = "Lung_CLAD", ident.2 = "Lung_Control")
EnhancedVolcano::EnhancedVolcano(dge.lungs.clad_v_control, 
                                 lab = rownames(dge.lungs.clad_v_control), 
                                 selectLab = cytotox_genes_overlapped_clones,
                                 x = "avg_log2FC", y = "p_val_adj", 
                                 FCcutoff = 0.5, pCutoff =  1e-05, 
                                 title = "CLAD vs Control",
                                 subtitle = "Lung",
                                 legendLabels = c("NS", expression(Log[2] ~ FC), "P-val", expression(P - value ~ and
                                                                                                     ~ Log[2] ~ FC)),
                                 legendPosition = "bottom",
                                 legendLabSize = 10,
                                 legendIconSize = 3,
                                 legendDropLevels = TRUE,
                                 colAlpha = 1, pointSize = 3, 
                                 drawConnectors = TRUE, widthConnectors = 0.5,
                                 labFace = "plain", boxedLabels = TRUE)

ggsave("lung_volcano.pdf", 
       plot = last_plot(),
       path = figures_dir, 
       width = 7, height = 7, units = "in", dpi = 300)


dge.hln.clad_v_control <- FindMarkers(explantCLAD.with_trb, ident.1 = "HLN_CLAD", ident.2 = "HLN_Control")
EnhancedVolcano::EnhancedVolcano(dge.hln.clad_v_control, 
                                 lab = rownames(dge.hln.clad_v_control), 
                                 selectLab = cytotox_genes_overlapped_clones,
                                 x = "avg_log2FC", y = "p_val_adj", 
                                 FCcutoff = 0.5, 
                                 title = "CLAD vs Control",
                                 subtitle = "HLN",
                                 legendLabels = c("NS", expression(Log[2] ~ FC), "P-val", expression(P - value ~ and
                                                                                                     ~ Log[2] ~ FC)),
                                 legendPosition = "bottom",
                                 legendLabSize = 10,
                                 legendIconSize = 3,
                                 legendDropLevels = TRUE,
                                 colAlpha = 1, pointSize = 5, 
                                 drawConnectors = TRUE, widthConnectors = 0.5,
                                 labFace = "plain", boxedLabels = TRUE)

ggsave("hln_volcano.pdf", 
       plot = last_plot(),
       path = figures_dir, 
       width = 7, height = 7, units = "in", dpi = 300)


## Cytotoxic genes heatmap

### All T cells (average), Lung only
Idents(explantCLAD.with_trb) <- explantCLAD.with_trb$sample_type
explantCLAD.with_trb.lung <- subset(explantCLAD.with_trb, idents = c("Lung"))

Idents(explantCLAD.with_trb.lung) <- explantCLAD.with_trb.lung$disease
explantCLAD.lung.averages <- AverageExpression(explantCLAD.with_trb.lung, 
                                               return.seurat = TRUE, 
                                               group.by = "ident")

DoHeatmap(explantCLAD.lung.averages, 
          features = c("GZMA", "GZMB", "PRF1", "GNLY", "KLRK1"), 
          draw.lines = FALSE) + 
  theme(axis.text = element_text(face = "bold", size = 18)) + 
  guides(color = "none") + 
  ggtitle("Lung Explant CD4 and CD8 T Cells") + 
  guides(fill = guide_colourbar(title = "Mean Expression", 
                                title.theme = element_text(face = "bold", size = 14),
                                barheight = 1, 
                                barwidth = 20)) + 
  theme(legend.position = "bottom", 
        plot.title = element_text(face = "bold", hjust = 0, size = 14))

ggsave("avg_heatmap_cytotoxic_lung_all.png", 
       plot = last_plot(),
       path = figures_dir, 
       width = 9, height = 8, units = "in", dpi = 300)

# reorder nicely for display
explantCLAD.with_trb$organ.disease <- factor(explantCLAD.with_trb$organ.disease, 
                                             levels = c("Lung_Control", "Lung_CLAD", "Lung_IPF", 
                                                        "HLN_Control", "HLN_CLAD", "HLN_IPF"))

Idents(explantCLAD.with_trb) <- explantCLAD.with_trb$phenotype
explantCLAD.with_trb.CD4 <- subset(explantCLAD.with_trb, idents = c("CD4 TRM"))
explantCLAD.with_trb.CD8 <- subset(explantCLAD.with_trb, idents = c("CD8 TRM"))

# CD4
Idents(explantCLAD.with_trb.CD4) <- explantCLAD.with_trb.CD4$organ.disease

explantCLAD.CD4.averages <- AverageExpression(explantCLAD.with_trb.CD4, 
                                              return.seurat = TRUE, 
                                              group.by = "ident")

cd4heatmap <- DoHeatmap(explantCLAD.CD4.averages, 
                        features = cytotox_genes, 
                        draw.lines = FALSE) + 
  theme(axis.text = element_text(face = "bold", size = 18)) + 
  guides(color = "none") + 
  ggtitle("CD4") + 
  guides(fill = guide_colourbar(title = "Mean Expression", 
                                title.theme = element_text(face = "bold", size = 14),
                                barheight = 1, 
                                barwidth = 20)) + 
  theme(legend.position = "bottom", 
        plot.title = element_text(face = "bold", hjust = 0, size = 14))

ggsave("avg_heatmap_cytotoxic_cd4.png", 
       plot = last_plot(),
       path = figures_dir, 
       width = 9, height = 8, units = "in", dpi = 300)

# CD8
Idents(explantCLAD.with_trb.CD8) <- explantCLAD.with_trb.CD8$organ.disease

explantCLAD.CD8.averages <- AverageExpression(explantCLAD.with_trb.CD8, 
                                              return.seurat = TRUE, 
                                              group.by = "ident")

DoHeatmap(explantCLAD.CD8.averages, 
          group.colors = rep("white", 6),
          angle = 60,
          features = cytotox_genes, 
          raster = FALSE,
          draw.lines = FALSE) + 
  theme(axis.text = element_text(face = "bold", size = 18)) + 
  guides(color = "none") + 
  ggtitle("CD8") + 
  guides(fill = guide_colourbar(title = "Mean Expression", 
                                title.theme = element_text(face = "bold", size = 14),
                                barheight = 1, 
                                barwidth = 20)) + 
  theme(legend.position = "bottom", 
        plot.title = element_text(face = "bold", hjust = 0, size = 14))

ggsave("avg_heatmap_cytotoxic_cd8.pdf", 
       plot = last_plot(),
       path = figures_dir, 
       width = 5, height = 8, units = "in", dpi = 300)

# # Both CD4 and CD8
# Idents(explantCLAD.with_trb) <- paste0(explantCLAD.with_trb$organ.disease, "-", explantCLAD.with_trb$phenotype)
# 
# explantCLAD.averages <- AverageExpression(explantCLAD.with_trb, 
#                                           return.seurat = TRUE, 
#                                           group.by = "ident")
# 
# comboheatmap <- DoHeatmap(explantCLAD.averages, 
#                           features = cytotox_genes, 
#                           draw.lines = FALSE) + 
#   theme(axis.text = element_text(face = "bold", size = 18)) + 
#   guides(color = "none") + 
#   ggtitle("CD4") + 
#   guides(fill = guide_colourbar(title = "Mean Expression", 
#                                 title.theme = element_text(face = "bold", size = 14),
#                                 barheight = 1, 
#                                 barwidth = 20)) + 
#   theme(legend.position = "bottom", 
#         plot.title = element_text(face = "bold", hjust = 0, size = 14))
# 
# ggsave("avg_heatmap_cytotoxic_combo.png", 
#        plot = last_plot(),
#        path = figures_dir, 
#        width = 9, height = 8, units = "in", dpi = 300)


## Cytotoxic genes heatmap (CLAD only)
explantCLAD.with_trb.CLAD_only <- subset(explantCLAD.with_trb, subset = disease == "CLAD")

Idents(explantCLAD.with_trb.CLAD_only) <- explantCLAD.with_trb.CLAD_only$phenotype

explantCLAD.averages.clad_only <- AverageExpression(explantCLAD.with_trb.CLAD_only, 
                                                    return.seurat = TRUE, 
                                                    group.by = "ident")

DoHeatmap(explantCLAD.averages.clad_only, 
          features = cytotox_genes, 
          raster = FALSE,
          draw.lines = FALSE) + 
  theme(axis.text = element_text(face = "bold", size = 18)) + 
  guides(color = "none") + 
  ggtitle("CLAD") + 
  guides(fill = guide_colourbar(title = "Mean Expression", 
                                title.theme = element_text(face = "bold", size = 14),
                                barheight = 1, 
                                barwidth = 20)) + 
  theme(legend.position = "bottom", 
        legend.justification = c(0.7,0),
        plot.title = element_text(face = "bold", hjust = 0, size = 14))

ggsave("avg_heatmap_cytotoxic_clad_only.pdf", 
       plot = last_plot(),
       path = figures_dir, 
       width = 7, height = 7, units = "in", dpi = 300)

## Patient-level Violins
explantCLAD.with_trb$patient <- factor(explantCLAD.with_trb$patient, 
                                       levels = c("2020-10-LT", "2020-45-LT", "2020-47-LT", "2021-10-LT", "2021-47-LT",
                                                  "2020-5-C", "2020-56-C", "2021-38-C", 
                                                  "IPF1", "IPF2", "IPF3")) 

explantCLAD.with_trb$patient_label <- factor(paste0(explantCLAD.with_trb$disease, "_", 
                                                    "P", as.character(as.numeric(explantCLAD.with_trb$patient))))

VlnPlot(explantCLAD.with_trb, features = c("GZMA"), group.by = "patient_label")




#Now that we know we only have T cells/clones with TCR annotations, let's add to the metadata and calculate the 
#percentage that that individual clones makes up of all the T cells for that individual patient and sample


clone_summary_by_patient <- explantCLAD.with_trb@meta.data %>% 
  group_by(patient, TRBaa) %>% 
  mutate(clone_freq_for_patient = n()) %>% 
  ungroup() %>% 
  group_by(patient) %>% 
  mutate(total_cells_for_patient = n(), 
         clone_pct_tcells_for_patient = 100 * clone_freq_for_patient / total_cells_for_patient) %>% 
  select(disease, patient, TRBaa, clone_freq_for_patient, total_cells_for_patient, clone_pct_tcells_for_patient) %>% 
  distinct(patient, TRBaa, .keep_all = TRUE) %>% 
  arrange(patient, desc(clone_pct_tcells_for_patient)) %>% 
  mutate(cumul_pct = cumsum(clone_pct_tcells_for_patient))


# clone_summary_by_patient_organ <- explantCLAD.with_trb@meta.data %>% 
#   group_by(disease, patient, sample_type, TRBaa) %>% 
#   mutate(freq_for_patient_organ = n(), 
#          pct_clones_for_patient_organ = 100 * freq_for_patient_organ/sum((explantCLAD.with_trb$patient == patient) & 
#                                                                            (explantCLAD.with_trb$sample_type == sample_type))) %>% 
#   select(disease, patient, sample_type, seurat_cluster, phenotype, TRBaa, freq_for_patient_organ, pct_clones_for_patient_organ) %>% 
#   arrange(patient, sample_type, desc(pct_clones_for_patient_organ))


# Get top 10 percent of CLAD clones
clad_clones_top10pct <- clone_summary_by_patient %>% filter("CLAD" == disease, cumul_pct <= 10) %>% pull(TRBaa)

# Now let's do DGE: top CLAD clones vs all other CLAD T cells

cells_clad_top_clones <- names(which(explantCLAD.with_trb$disease == "CLAD" & 
                                       explantCLAD.with_trb@meta.data$TRBaa %in% clad_clones_top10pct))

cells_clad_NOT_top_clones <- names(which(explantCLAD.with_trb$disease == "CLAD" & 
                                           !(explantCLAD.with_trb@meta.data$TRBaa %in% clad_clones_top10pct)))

top_clone_degs_clad_vs_clad <- FindMarkers(explantCLAD.with_trb, 
                                           ident.1 = cells_clad_top_clones, 
                                           ident.2 = cells_clad_NOT_top_clones,
                                           only.pos = TRUE) %>% 
  filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))

# subset just the CLAD cells
explantCLAD.with_trb.only_clad = subset(explantCLAD.with_trb, subset = disease == "CLAD")

explantCLAD.with_trb.only_clad$highlight <- factor(ifelse(explantCLAD.with_trb.only_clad$TRBaa %in% clad_clones_top10pct, "Top Clone", "Other"), 
                                                   levels = c("Top Clone", "Other"))


Idents(explantCLAD.with_trb.only_clad) <- explantCLAD.with_trb.only_clad$highlight
explantCLAD.averages.only_clad <- AverageExpression(explantCLAD.with_trb.only_clad, 
                                                    return.seurat = TRUE, 
                                                    group.by = "ident")

DoHeatmap(explantCLAD.averages.only_clad, 
          group.bar.height = 0, 
          raster = FALSE,
          features = cytotox_genes, 
          draw.lines = FALSE) + 
  theme(axis.text = element_text(face = "bold", size = 18)) + 
  guides(color = "none") + 
  ggtitle("CLAD") + 
  guides(fill = guide_colourbar(title = "Mean Expression", 
                                title.position = "top",
                                title.theme = element_text(face = "bold", size = 14),
                                barheight = 1, 
                                barwidth = 10)) + 
  theme(legend.position = "bottom", 
        legend.title.align = 0.5,
        plot.title = element_text(face = "bold", hjust = 0, size = 14)) 

ggsave("avg_heatmap_cytotoxic_clad_top_clones.pdf", 
       plot = last_plot(),
       path = figures_dir, 
       width = 3.5, height = 7, units = "in", dpi = 300)

VlnPlot(explantCLAD.with_trb.only_clad, 
        features = cytotox_genes,
        pt.size = 0, 
        stack = TRUE, flip = TRUE, 
        cols = viridis::viridis(n = 2),
        group.by = "sample_type", 
        split.by = "highlight") + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 0, face = "bold", size = 14, hjust = 0.5))

ggsave("violin_top_clones_by_organ.png", 
       plot = last_plot(),
       path = figures_dir, 
       width = 8.5, height = 11, units = "in", dpi = 300)

VlnPlot(explantCLAD.with_trb.only_clad, 
        features = cytotox_genes,
        pt.size = 0, 
        stack = TRUE, flip = TRUE, 
        cols = viridis::viridis(n = 2),
        group.by = "phenotype", 
        split.by = "highlight") + 
  scale_x_discrete(labels = sub(" ", "\n", levels(explantCLAD.with_trb.only_clad$phenotype))) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 0, face = "bold", size = 14, hjust = 0.5))

ggsave("violin_top_clones_by_phenotype.png", 
       plot = last_plot(),
       path = figures_dir, 
       width = 8.5, height = 11, units = "in", dpi = 300)

## Flow
p <- flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "function_pos", 
               fig_parent = "CD3+ CD8+", 
               fig_marker = "GzmA+")

ggsave("flow_cd8_gzmA.pdf", plot = p, path = figures_dir, width = 5, height = 5, units = "in")

p <- flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "function_pos", 
               fig_parent = "CD3+ CD8+", 
               fig_marker = "GzmB+")

ggsave("flow_cd8_gzmB.pdf", plot = p, path = figures_dir, width = 5, height = 5, units = "in")

p <- flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "function_pos", 
               fig_parent = "CD3+ CD8+", 
               fig_marker = "GzmK+")

ggsave("flow_cd8_gzmK.pdf", plot = p, path = figures_dir, width = 5, height = 5, units = "in")

p <- flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "function_pos", 
               fig_parent = "CD3+ CD8+", 
               fig_marker = "Perf+")

ggsave("flow_cd8_perf.pdf", plot = p, path = figures_dir, width = 5, height = 5, units = "in")



## Figure - Functional
p <- flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "phenotype", 
               fig_parent = "CD3+ CD4+", 
               fig_marker = "FOXP3+ CD25+")

ggsave("flow_treg.pdf", plot = p, path = figures_dir, width = 5, height = 5, units = "in")

p <- flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "phenotype", 
               fig_parent = "CD3+ CD4+ CD69+ CD103+", 
               fig_marker = "FOXP3+ CD25+")

ggsave("flow_treg_trm.pdf", plot = p, path = figures_dir, width = 5, height = 5, units = "in")


## Figure 4F (GZMA and GZMB MFI)
p <- flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "function_mfi_pos", 
               fig_parent = "CD3+ CD8+", 
               fig_marker = "GZMA+", 
               fig_y = "Mean MFI") + 
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))

ggsave("mfi_gzma_cd8.pdf", plot = p, path = figures_dir, width = 5, height = 5, units = "in")

p <- flow_plot(full_flow_data, cell_prop_wilcox, 
               fig_panel = "function_mfi_pos", 
               fig_parent = "CD3+ CD8+", 
               fig_marker = "GZMB+", 
               fig_y = "Mean MFI") + 
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))

ggsave("mfi_gzmb_cd8.pdf", plot = p, path = figures_dir, width = 5, height = 5, units = "in")



full_flow_data %>% filter("function_mfi_pos" == panel) %>% 
  filter(marker %in% c("GZMA+", "GZMB+", "GZMK+", "Perf+")) %>% 
  mutate(parent = gsub(" ", "\n", parent), marker = gsub(" ", "\n", marker)) %>%
  drop_na(value) %>%
  ggplot(aes(x = organ, y = value, fill = disease)) +
  stat_summary(geom = "bar", fun = mean, position = "dodge", width = 0.9, alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.data = mean_sdl, fun.args = list(mult = 1),  
               position = position_dodge(width = 0.9), 
               width = 0.2, color = "darkslategray") + 
  geom_point(aes(color = disease), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) + 
  facet_grid(parent ~ marker) +
  xlab("") + ylab("Mean MFI") +
  #scale_y_continuous(breaks = seq(0,100,25), limits = c(0,100)) +
  #ggtitle("Function MFI (Mean Positive MFI)") + 
  ggprism::scale_color_prism(palette = "colorblind_safe", name = "") + 
  ggprism::scale_fill_prism(palette = "colorblind_safe", name = "") + 
  ggpubr::theme_pubr() + 
  theme(axis.title = element_text(face = "bold", size = 17),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14, angle = 90, face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 21),
        strip.text.y = element_text(angle = 0, size = 12, face = "bold"),
        strip.text.x = element_text(angle = 0, size = 12, face = "bold"),
        legend.position = "bottom")

ggsave("function_mfi_pos.eps", plot = last_plot(), path = figures_dir, width = 14, height = 10, units = "in")
ggsave("function_mfi_pos.pdf", plot = last_plot(), path = figures_dir, width = 14, height = 10, units = "in")


function_pos_flow <- 
  full_flow_data %>% 
  filter("function_pos" == panel, 
         parent == "CD3+ CD8+", 
         marker %in% c("GzmA+", "GzmB+", "GzmK+", "Perf+", "CD107a+", "IFNg+", "TNFa+"), 
         !(grepl("ccr7", parent, ignore.case = TRUE)))

function_pos_wilcox <- 
  function_pos_flow %>% 
  select(-sample_id, -full_label) %>% 
  group_by(panel, parent, marker, organ) %>% 
  drop_na(value) %>% 
  wilcox_test(value ~ disease) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p") %>% 
  add_significance("p.adj") %>% 
  add_xy_position(x = "organ", dodge = 0.9, fun = "mean_sd", step.increase = 0.03) %>% 
  #add_xy_position(x = "organ", dodge = 0.9, fun = "mean_sd", scales = "free") %>% 
  mutate(p.format = p_format(p, accuracy = 0.01, leading.zero = FALSE), 
         p.adj.format = p_format(p.adj, accuracy = 0.01, leading.zero = FALSE), 
         disease = group1)

function_pos_flow %>% 
  mutate(parent = gsub(" ", "\n", parent), marker = gsub(" ", "\n", marker)) %>%
  drop_na(value) %>%
  ggplot(aes(x = organ, y = value, fill = disease)) +
  stat_summary(geom = "bar", fun = mean, position = "dodge", width = 0.9, alpha = 0.5) +
  stat_summary(geom = "errorbar", 
               #fun.data = mean_sdl, fun.args = list(mult = 1),  
               fun.data = mean_sd_prop, 
               position = position_dodge(width = 0.9), 
               width = 0.2, color = "darkslategray") + 
  geom_point(aes(color = disease), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) + 
  stat_pvalue_manual(data = function_pos_wilcox, 
                     hide.ns = "p.adj", 
                     tip.length = 0, 
                     #step.increase = 0.05,
                     #step.group.by = "organ",
                     label = "{p.adj.signif}", 
                     label.size = 4, 
                     bracket.size = 1) + 
  facet_grid( ~ factor(marker, levels = c("GzmA+", "GzmB+", "GzmK+", "Perf+", "CD107a+", "IFNg+", "TNFa+"))) +
  xlab("") + ylab("Mean Proportion") +
  #scale_y_continuous(breaks = seq(0,100,25), limits = c(0,100)) +
  #ggtitle("Function MFI (Mean Positive MFI)") + 
  ggprism::scale_color_prism(palette = "colorblind_safe", name = "") + 
  ggprism::scale_fill_prism(palette = "colorblind_safe", name = "") + 
  ggpubr::theme_pubr() + 
  theme(axis.title = element_text(face = "bold", size = 17),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14, angle = 90, face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 21),
        strip.text.y = element_text(angle = 0, size = 12, face = "bold"),
        strip.text.x = element_text(angle = 0, size = 12, face = "bold"),
        legend.position = "bottom")

#ggsave("function_pos.eps", plot = last_plot(), path = figures_dir, width = 14, height = 10, units = "in")
ggsave("function_pos.pdf", plot = last_plot(), path = figures_dir, width = 10, height = 7, units = "in")



## Supplement



FeaturePlot(explantCLAD.with_trb, features = c("CD69", "ITGAE"), blend = TRUE, split.by = "organ.disease")

FeaturePlot(explantCLAD.with_trb, features = c("SELL", "S1PR1"))

ggsave("s1.png", 
       plot = last_plot(),
       path = figures_dir, 
       width = 11, height = 8.5, units = "in", dpi = 300)



