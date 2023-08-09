library(igraph)
library(data.table)

nodesTable <- fread("/cytoscape_sessions/clueGOResults_nodeTable_pairedGO_cov2.csv", sep=",")
network <- fread("/cytoscape_sessions/networkCluego_cov2.sif", sep="\t")
sharedGenes <- fread("sharedGenesAcrossStrains_goterms_clueGO_genelist.txt", header = F)
colnames(network) <- c("node1","interaction_type","node2")

nodeNames <- union(network$node1,network$node2)

nodeNames.df <- left_join(as.data.frame(nodeNames),nodesTable %>% dplyr::select(`shared name`,Term),by=c("nodeNames"="shared name"))

nodeNames.df <- nodeNames.df %>% mutate(Term=if_else(is.na(Term),nodeNames,Term))

network <- network %>% mutate(node1=gsub("%20","_",node1),node2=gsub("%20","_",node2)) %>% dplyr::select(-interaction_type)

nodes1 <- left_join(network %>% dplyr::select(node1),nodeNames.df, by=c("node1"="nodeNames")) %>% rename(nodeA=Term)
nodes2 <- left_join(network %>% dplyr::select(node2),nodeNames.df, by=c("node2"="nodeNames")) %>% rename(nodeB=Term)

networkNamed <- cbind(nodes1 %>% dplyr::select(nodeA),nodes2 %>% dplyr::select(nodeB))

networkNamed <- networkNamed %>% filter((nodeB %in% sharedGenes$V1) | str_detect(nodeB,"[[::lower::]]"))

net <- graph_from_data_frame(networkNamed, directed = F)

 

#Analysis
#Degree centrality:

geneList <- c("FTO","HNRNPC","YTHDF1","WTAP","RBM15","RBM15B","YTHDC1","HNRNPA2B1")


degree_scores <- degree(net,v=V(net),mode ="total", 
       loops = FALSE, normalized = FALSE)

hub_scores <- hub.score(net)
closeness_score <- closeness(net)
betweenness_score <- betweenness(net,v = V(net),directed = F)


df <-  betweenness_score %>% as.data.frame() %>% 
  rownames_to_column("name") %>% 
  rename(betweeness=".") %>% 
  filter(!str_detect(name,"[[::lower::]]")) %>%
  #filter(betweeness > 0) %>% 
  mutate(label=if_else(name %in% geneList,name,""))

#reorder(name,desc(betweeness))
df %>% ggplot(aes(desc(rank(betweeness)),betweeness, label=label)) + 
  geom_point(color = ifelse(df$name %in% geneList, "red", "grey50")) +
  #geom_boxplot(alpha=0.1, width=0.1) +
  geom_text_repel(nudge_x = 0.9, nudge_y = 0.1) + 
  theme_prism(base_size = 10) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(x="Rank", y="Betweeness centrality")

# Hub score plot
hubdf <- hub_scores$vector %>% as.data.frame() %>% 
  rownames_to_column("name") %>% 
  rename(hub_score=".") %>% 
  filter(!str_detect(name,"[[::lower::]]")) %>%
  #filter(hub_score > 0) %>% 
  mutate(label=if_else(name %in% geneList,name,""))

#reorder(name,desc(hub_score))
hubdf %>% ggplot(aes("hub",hub_score, label=label)) + 
  geom_point(color = ifelse(hubdf$name %in% geneList, "red", "grey50"), position = position_jitter(width = 0.01)) +
  geom_boxplot(alpha=0.4, width=0.1) +
  geom_text_repel(nudge_x = 0.1, nudge_y = 0.01) + theme_prism(base_size = 10) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 0.75) + 
  labs(x="Nodes", y="Hub centrality") -> pHubCentrality



hubdf %>% ggplot(aes(desc(rank(hub_score,na.last = T)),hub_score, label=label)) + 
  geom_point(color = ifelse(hubdf$name %in% geneList, "red", "grey50"), size=1, alpha=0.5) +
  #geom_boxplot(alpha=0.1, width=0.1) +
  geom_text_repel(nudge_x = 0.9, nudge_y = 0.05, max.overlaps = Inf) + 
  theme_prism(base_size = 10) + 
  scale_x_continuous(breaks = seq(-400,0,by=100),labels =  as.character(seq(0,400,by=100))) +
  scale_y_log10() +
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(x="Rank", y="Hub centrality") -> pHubRank
# Degree distribution

df %>% ggplot(aes(desc(rank(betweeness,na.last = T, ties.method = "first")),betweeness, label=label)) + 
  geom_point(color = ifelse(df$name %in% geneList, "red", "grey50"), size=1, alpha=0.5) +
  geom_text_repel(max.overlaps = Inf, nudge_x = 0.5, nudge_y = 0.5) +
  theme_prism(base_size = 10) + 
  scale_x_continuous(breaks = seq(-400,0,by=100),labels =  as.character(seq(0,400,by=100))) +
  labs(x="Rank", y="Betweeness centrality") -> pBetweenessRank


df.deg <- degree_scores %>% as.data.frame() %>% 
  rownames_to_column("name") %>% 
  rename(degree=".") %>% 
  filter(!str_detect(name,"[[::lower::]]")) %>%
  #filter(betweeness > 0) %>% 
  mutate(label=if_else(name %in% geneList,name,""))

df.deg %>% ggplot(aes(desc(rank(degree,na.last = T, ties.method = "first")),degree, label=label)) + 
  geom_point(color = ifelse(df$name %in% geneList, "red", "grey50"), size=1, alpha=0.5) +
  geom_text_repel(max.overlaps = Inf, nudge_x = 0.5, nudge_y = 0.5) +
  theme_prism(base_size = 10) + 
  scale_x_continuous(breaks = seq(-400,0,by=100),labels =  as.character(seq(0,400,by=100))) +
  labs(x="Rank", y="Connectivity degree") -> pDegreeRank

df.deg %>% mutate(type=if_else(label!="","m6A related","non-m6A related")) %>% 
  ggplot(aes(type,degree)) + 
  geom_boxplot(width=0.25, alpha=0.5,outlier.size = 0) + 
  #geom_point(position = position_jitter(width = 0.05), size=1, alpha=0.5) + 
  stat_compare_means(label.y = 9,method = "wilcox", paired = FALSE) + 
  labs(x="",y="Connectivity degree") + 
  theme_prism(base_size = 10,axis_text_angle = 45, border = T) +
  theme(panel.background = element_rect(fill="transparent"),
                                        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                                        ) -> boxplotDegree

hist(degree_scores)

# Generate random networks with same number of nodes and links:


n <- 1000

nodes <- length(V(net))
edges <- length(E(net))

set.seed(23)
#library(devtools)
#install_github("jfukuyama/phyloseqGraphTest")
#library(phyloseqGraphTest)

netDegreeMatrix <- NULL
for(i in 1:n){

  p <- runif(1,min = 0, max = 1)
  #newnet <- erdos.renyi.game(n=nodes, p.or.m = 0.5, type = "gnp",directed = F,loops = F)
  newnet <- sample_pa(n=nodes,directed = F)
  #newnet <- rewire(net,with = each_edge(prob = p,loops = F,multiple = T))
  #newnet <- sample_gnp(n =nodes ,p=p , directed = F,loops = F)
  degMat <- igraph::degree(newnet,mode = "total")
  netDegreeMatrix <- rbind(netDegreeMatrix,degMat)
}

point.estimate.mean.random <- mean(netDegreeMatrix)

point.estimate.mean.observed <- mean(degree_scores)

t.test(netDegreeMatrix, mu = point.estimate.mean.observed, alternative = "two.sided")
wilcox.test(netDegreeMatrix, mu = point.estimate.mean.observed, alternative = "two.sided")

netDegreeMatrix %>% as.data.frame() %>%
  pivot_longer(names_to = "network",values_to = "degree",dplyr::everything()) %>%
  ggplot(aes(degree)) + 
  geom_histogram(binwidth = 1,fill="white", color="black") +
  geom_vline(xintercept = unlist(point.estimate.mean.observed), color="red") +
  geom_vline(xintercept = unlist(point.estimate.mean.random), color="blue") +
  annotate(geom = "text",x=unlist(point.estimate.mean.observed)+0.15,y=250000,label="Mean degree \n observed",color="red", size=5, hjust=0) +
  annotate(geom = "text",x=unlist(point.estimate.mean.random)+0.15,y=250000,label="Mean degree \n random networks",color="blue", size=5, hjust=0) +
  annotate(geom = "text",x=8,y=25000,label="Wilcoxon \n p < 2.2e-16") +
  scale_x_continuous(limits = c(0,10), expand = expansion(mult = 0.001), breaks = seq(0,10, by=2), labels = seq(0,10,by=2)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  labs(x="Connectivity degree",y=str_wrap("Node degree distribution \n 1000 random networks",20)) +
  theme(aspect.ratio = 1) +
  theme_prism(base_size = 10) -> pCompareToRandomNetwork


#p <- pCompareToRandomNetwork + pHubCentrality + plot_annotation(tag_levels = 'A')

p <- pDegreeRank + inset_element(boxplotDegree, 0.4,0.3,1,1,align_to = 'full')

patchwork <- pCompareToRandomNetwork + plot_annotation(tag_levels = c('A'))

patchwork <- patchwork + p + plot_layout(tag_level = "new") + plot_annotation(tag_levels = c("A",""))

ggsave("randomNetworkAnalysis.pdf", plot = patchwork, width = 12, height = 5.5)
