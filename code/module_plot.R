nodes <- nodes %>% arrange(type)

order1 <- (nodes[1:44,] %>% arrange(module))$name ##### adjust number
order2 <- (nodes[45:383,] %>% arrange(module))$name

edges <- get.data.frame(graph, what = "edges")

edge_list <- edges %>%
  inner_join(nodes, by = c("from" = "name")) %>%
  inner_join(nodes, by = c("to" = "name")) %>%
  mutate(Module = ifelse(module.x == module.y, module.x, NA))

plot_data <- edge_list %>% mutate(
  to = factor(to, levels = order1),
  from = factor(from, levels = order2))

ggplot(plot_data, aes(x = from, y = to, fill = Module)) +
  geom_tile() +
  theme_bw() +
  #scale_fill_discrete(na.value ="white") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())

