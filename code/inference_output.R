library(evolnets)
library(tidyverse)
library(MCMCpack)
library(coda)
library(treeio)
library(ggtree)
library(ape)

## Check convergence ----

chain1_full <- read.table("./output/log_and_txt_files/out.3itho_color.log", header = TRUE)[,c(1,5,8,9)]
chain2_full <- read.table("./output/log_and_txt_files/out.5itho_color.log", header = TRUE)[,c(1,5,8,9)]

colnames(chain1_full) <- colnames(chain2_full) <- c("iteration","clock", "lambda[02]", "lambda[20]")

max(chain1_full$iteration)
max(chain2_full$iteration)

chain1 <- filter(chain1_full, iteration >= 8000 & iteration <= 80000)
chain2 <- filter(chain2_full, iteration >= 8000 & iteration <= 80000)

gelman.diag(mcmc.list(as.mcmc(chain1), as.mcmc(chain2)))

# already converged at 28000 gen


## Character history ----

# Read history
# history1 <- read_history("./output/log_and_txt_files/out.3itho_color.history.txt", burnin = 0.1)
# history2 <- read_history("./output/log_and_txt_files/out.5itho_color.history.txt", burnin = 0.1)
# 
# (ne1 <- count_events(history1))
# (ne2 <- count_events(history2))
# effective_rate(history1,tree)
# (gl1 <- count_gl(history1))
# (gl2 <- count_gl(history2))
# 
# length(unique(history1$iteration))
# length(unique(history2$iteration))
# 
# history <- filter(history1, iteration %in% sample(unique(history1$iteration), 1000))
# length(unique(history$iteration))
# saveRDS(history, "./output/history_1000it.rds")

history <- readRDS("./output/history_1000it.rds")
n_events <- count_events(history)

# Read trees
tree_index <- treeio::read.beast.newick("./data/itho_color_Rev.tre")
tree <- tree_index@phylo

rate <- effective_rate(history, tree)

index_node <- tree_index@data
names(index_node$index) <- NULL

index_labels <- index_node %>% 
  mutate(node = as.numeric(node)) %>% 
  arrange(node) %>% 
  filter(node > Ntip(tree)) %>% 
  pull(index)

tree$node.label <- paste0("Index_",index_labels)  # index_labels
#saveRDS(tree, "./output/tree_node_index.rds")


#plot
ggtree(tree, ladderize = F, layout="circular", size = 0.2) +
  geom_nodelab(colour = "red", size = 1) +
  geom_tiplab(size = 0.8)

host_tree <- read.tree("./data/44tips.tre")
hosts <- paste0("C",1:Ntip(host_tree))
host_tree$tip.label <- hosts


## Ancestral states ----

at_nodes <- posterior_at_nodes(history, tree, host_tree)
pp_at_nodes <- at_nodes[[2]]

high_pp <- function(matrix, pt, weighted = TRUE){
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      if(matrix[i,j] < pt){
        matrix[i,j] = 0
      } else{
        if(weighted == FALSE){
          matrix[i,j] = 1
        }
      }
    }
  }
  
  matrix <- matrix[rowSums(matrix)!= 0, colSums(matrix)!= 0]
  
  return(matrix)
}

matrix_95 <- high_pp(pp_at_nodes, 0.95)
matrix_90 <- high_pp(pp_at_nodes, 0.90)


# # find modules to use for plotting colors
# net <- read.csv("./data/spp_col_net.csv", header = T, row.names = 1)
# mod_members <- mycomputeModules(net)
# 
# mod_list <- listModuleInformation(mod_members)[[2]]
# nmod <- length(mod_list)
# 
# # !!! fix species names ( add _ )
# all_mod <- tibble()
# for(m in 1:nmod){
#   members <- unlist(mod_list[[m]])
#   mtbl <- tibble(name = members,
#                  module = rep(paste0("M",m), length(members)))
# 
#   all_mod <- bind_rows(all_mod, mtbl)
# }
# 
# all_mod <- mutate(all_mod, type = case_when(name %in% hosts ~ "host",
#                                             TRUE ~ "symbiont"))
# 
# write.csv(all_mod, "./output/all_mod.csv", row.names = F)
all_mod <- read.csv("./output/all_mod.csv", header = T)
all_mod <- read.csv("./output/most_common_11.csv", header = T)
  
asr <- plot_ancestral_states(tree, at_nodes, all_mod,
                               #module_order = paste0("M",1:13),
                               layout = "circular",
                               threshold = 0.9, 
                               point_size = 2, 
                               dodge_width = 0.002)

plot_module_matrix2(t(as.matrix(net)), at_nodes, tree, host_tree,
                     all_mod, 
                     threshold = 0.95, 
                     point_size = 2, 
                     dodge_width = 0.015)

plotModuleWeb(mod_members)
