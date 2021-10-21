# Exploring the data before inference in RevBayes

library(ape)
library(tidyverse)
library(ggtree)

tree <- read.tree("./data/itho_phylo.tre")
host_tree <- read.tree("./data/44tips.tre")
data <- read.nexus.data("./data/original_Ithomiini_mimicry_data_mtx.nex") %>% 
  as.data.frame() %>% 
  apply(2, as.numeric) %>% 
  t()
colnames(data) <- paste0("C",1:44)

write.csv(t(data),"./data/spp_col_net.csv", row.names = T)

ggt <- ggtree(tree, layout = 'circular') + geom_tiplab(size = 1) + theme_tree2()
gheatmap(ggt, data)

# histogram

common <- data %>%
  as_tibble() %>% 
  mutate(across(everything(), sum)) %>% 
  slice_head() %>% 
  pivot_longer(everything(),names_to = "color", values_to = "n_spp") %>% 
  arrange(desc(n_spp))

plot_spp_rich <- ggplot(common) + 
  geom_col(aes(reorder(color, desc(n_spp)),n_spp)) + 
  labs(x = "Color patterns", y = "Species richness") +
  #geom_hline(yintercept = 10, col = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270))

ggplot(common) + geom_density(aes(n_spp)) + theme_bw()


col_div <- t(data) %>%
  as_tibble() %>% 
  mutate(across(everything(), sum)) %>% 
  slice_head() %>% 
  pivot_longer(everything(),names_to = "spp", values_to = "n_col")

plot_col_div <- ggplot(col_div) + 
  geom_col(aes(reorder(spp, desc(n_col)),n_col)) + 
  labs(x = "Number of color patterns", y = "Species") +
  #geom_hline(yintercept = 10, col = "red") +
  theme_bw() +
  theme(axis.text.x = element_blank())


ggplot(col_div) + geom_density(aes(n_col)) + theme_bw()


# Only color patterns present in more than 10 species

topc <- common %>% 
  filter(n_spp >= 10) %>%   # we loose 17 spp 
  pull(color)               # the only way to keep all spp is to filter >= 2 (only 5 colors out)

data_topc <- data[,topc]
table(rowSums(data_topc))

