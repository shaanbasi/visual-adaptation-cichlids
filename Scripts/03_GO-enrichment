library(dplyr)

genes_Mz_up <- deg_tbl_Mz_vs_On_clean %>%
  filter(sig == "Mz higher", !is.na(gene_name)) %>%
  pull(gene_name) %>%
  unique()

genes_On_up <- deg_tbl_Mz_vs_On_clean %>%
  filter(sig == "On higher", !is.na(gene_name)) %>%
  pull(gene_name) %>%
  unique()

length(genes_Mz_up)
length(genes_On_up)

intersect(genes_Mz_up, genes_On_up)

background_genes <- deg_tbl_Mz_vs_On_clean %>%
  filter(!is.na(gene_name)) %>%
  pull(gene_name) %>%
  unique()


##GO profiler
library(gprofiler2)

gost_Mz <- gost(
  query = genes_Mz_up,
  organism = "drerio",
  custom_bg = background_genes,
  correction_method = "fdr",
  significant = FALSE
)

gost_On <- gost(
  query = genes_On_up,
  organism = "drerio",
  custom_bg = background_genes,
  correction_method = "fdr",
  significant = FALSE
)

head(gost_Mz$result$term_name)
head(gost_On$result$term_name)

head(genes_Mz_up)
head(genes_On_up)

nrow(gost_Mz$result)
nrow(gost_On$result)


##save
gost_Mz_final <- gost_Mz$result %>%
  filter(source %in% c("GO:BP", "GO:MF", "GO:CC")) %>%
  select(term_id, term_name, source, p_value, intersection_size, query_size)

gost_On_final <- gost_On$result %>%
  filter(source %in% c("GO:BP", "GO:MF", "GO:CC")) %>%
  select(term_id, term_name, source, p_value, intersection_size, query_size)

write.csv(gost_Mz_final, "GO_Mz_up_clean.csv", row.names = FALSE)
write.csv(gost_On_final, "GO_On_up_clean.csv", row.names = FALSE)


##
Mz_GO <- gost_Mz$result %>%
  filter(source %in% c("GO:BP", "GO:MF", "GO:CC"))

On_GO <- gost_On$result %>%
  filter(source %in% c("GO:BP", "GO:MF", "GO:CC"))

Mz_BP <- Mz_GO %>% filter(source == "GO:BP")
On_BP <- On_GO %>% filter(source == "GO:BP")
go_df <- bind_rows(Mz_GO, On_GO)

Mz_top <- Mz_BP %>%
  arrange(p_value) %>%
  select(term_name, p_value, intersection_size) %>%
  head(15)

On_top <- On_BP %>%
  arrange(p_value) %>%
  select(term_name, p_value, intersection_size) %>%
  head(15)

Mz_top
On_top


##vision related
Mz_vision <- Mz_BP %>%
  filter(grepl("phototransduction|opsin|retina|retinal|visual perception|eye development", term_name, ignore.case = TRUE)) %>%
  mutate(group = "Mz")

On_vision <- On_BP %>%
  filter(grepl("phototransduction|opsin|retina|retinal|visual perception|eye development", term_name, ignore.case = TRUE)) %>%
  mutate(group = "On")

Mz_vision
On_vision

vision_go <- bind_rows(Mz_vision, On_vision)
nrow(vision_go)

##celan up
On_BP_clean <- On_BP %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

Mz_BP_clean <- Mz_BP %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)



##VISUALISE
library(dplyr)
library(ggplot2)

Mz_plot <- Mz_BP %>%
  mutate(group = "Mz") %>%
  filter(p_value < 0.05)

On_plot <- On_BP %>%
  mutate(group = "On") %>%
  filter(p_value < 0.05)

go_plot <- bind_rows(Mz_plot, On_plot)

go_plot <- go_plot %>%
  mutate(gene_ratio = intersection_size / query_size,
         log_p = -log10(p_value))

top_go <- go_plot %>% ##not necessary for poster, maybe later?
  group_by(group) %>%
  arrange(p_value) %>%
  slice_head(n = 10) %>%
  ungroup()

##vision-specific separate plot
library(stringr)
library(tidytext)

top_n <- 10

vision_plot <- vision_go %>%
  filter(!is.na(p_value)) %>%
  
  group_by(group) %>%
  slice_min(order_by = p_value, n = top_n) %>%
  ungroup() %>%
  
  mutate(
    gene_ratio = intersection_size / query_size,
    neglogp = -log10(p_value),
    term_name = str_trunc(term_name, width = 35),
    
    term_f = reorder_within(term_name, gene_ratio, group)
  )

ggplot(vision_plot, aes(
  x = gene_ratio,
  y = term_f
)) +
  
  geom_point(aes(
    size = intersection_size,
    color = neglogp
  ), alpha = 0.9) +
  
  facet_wrap(~group, scales = "free_y") +
  
  scale_y_reordered() +
  
  scale_color_gradient(
    low = "#2c7bb6",
    high = "#d7191c",
    name = "Enrichment significance"
  ) +
  scale_size(range = c(2, 7), name = "Gene count") +
  
  theme_classic(base_size = 13) +
  
  labs(
    title = "Vision-related GO enrichment",
    x = "Gene ratio",
    y = NULL,
    size = "Gene count",
    color = "Enrichment significance"
  ) +
  
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    panel.spacing = unit(1.5, "lines")
  )


## POSTER specific - simple and readable
poster_GO <- ggplot(vision_plot, aes(
  x = gene_ratio,
  y = term_f
)) +
  
  geom_point(aes(
    size = intersection_size,
    color = neglogp
  ), alpha = 0.85) +
  
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  
  scale_y_reordered() +
  
  scale_color_gradient(
    low = "#2c7bb6",
    high = "#d7191c",
    name = "Enrichment\n(-log10 p)"
  ) +
  
  scale_size(range = c(1.5, 5), name = "Gene count") +
  
  theme_classic(base_size = 11) +
  
  labs(
    title = "Vision-related GO enrichment",
    x = "Gene ratio",
    y = NULL
  ) +
  
  theme(
    strip.text = element_text(face = "bold", size = 10),
    
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 9),
    
    axis.ticks = element_line(linewidth = 0.3),
    
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    
    legend.key.size = unit(0.4, "cm"),
    
    panel.spacing = unit(1, "lines"),
    
    plot.title = element_text(face = "bold", size = 13)
  )

poster_GO

ggsave(
  filename = "GO_vision_poster.png",
  plot = last_plot(),
  width = 10,
  height = 12,
  units = "in",
  dpi = 300
)
