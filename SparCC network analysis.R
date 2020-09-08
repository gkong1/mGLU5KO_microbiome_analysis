# Load the libraries
library(igraph)
library(qgraph)
library(vegan)
library(MCL)
library(SpiecEasi)

# Import in data
source("common_functions.R")
source("igraph_plot_functions.R")
asv_data <- read.csv("ASVs_taxa_as_rows.txt", sep = "\t", header = T, row.names = 1) # ASVs as rows, samples as columns

# Filtering low counts
asv_data.filt <- low.count.removal(asv_data, percent = 0.01)
asv_data.filt <- asv_data.filt$data.filter
asv_data.filt <- t(asv_data.filt)

# SparCC network
## Bootstrap sparCC
sparcc.boot.matrix <- sparccboot(asv_data.filt, R = 100, ncpus = 3)

## Bootstrap pval
sparcc.boot.pval <- pval.sparccboot(sparcc.boot.matrix)

### Number of significant correlations
length(which(sparcc.boot.pval$pvals < 0.05)) #604

### Replace "NaN" with 1
index <- which(sparcc.boot.pval$pvals == "NaN")
sparcc.boot.pval$pvals[index] <- 1

## Reshape pvalue output to matrix
sparcc.boot.pval.mat <- diag(0.5, nrow = dim(asv_data.filt)[2], ncol = dim(asv_data.filt)[2])
sparcc.boot.pval.mat[upper.tri(sparcc.boot.pval.mat, diag = FALSE)] <- sparcc.boot.pval$pvals
sparcc.boot.pval.mat <- sparcc.boot.pval.mat + t(sparcc.boot.pval.mat)

### Add taxa names to rows and columns
rownames(sparcc.boot.pval.mat) <- colnames(asv_data.filt)
colnames(sparcc.boot.pval.mat) <- colnames(asv_data.filt)

## Reshape bootstrapped correlations output to matrix
sparcc.boot.cor.mat <- diag(0.5, nrow = dim(asv_data.filt)[2], ncol = dim(asv_data.filt)[2])
sparcc.boot.cor.mat[upper.tri(sparcc.boot.cor.mat, diag = FALSE)] <- sparcc.boot.pval$cors
sparcc.boot.cor.mat <- sparcc.boot.cor.mat + t(sparcc.boot.cor.mat)

### Add taxa names to rows and columns
rownames(sparcc.boot.cor.mat) <- colnames(asv_data.filt)
colnames(sparcc.boot.cor.mat) <- colnames(asv_data.filt)

### Number of significant abs(correlations) >= 0.6
length(which(sparcc.boot.pval$cors >= 0.6 & sparcc.boot.pval$pvals < 0.05)) #46
length(which(sparcc.boot.pval$cors <= -0.6 & sparcc.boot.pval$pvals < 0.05)) #32

## Finally to have matrix of reads with p < 0.05 and abs(cor) >= 0.6
sparcc.final <- matrix(nrow = dim(asv_data.filt)[2], ncol = dim(asv_data.filt)[2])
rownames(sparcc.final) <- colnames(asv_data.filt)
colnames(sparcc.final) <- colnames(asv_data.filt)

### Set threshold
sparcc.cutoff <- 0.6
pval.cutoff <- 0.05

### Making the matrix
for (i in 1:nrow(sparcc.final)) {
  for (j in 1:ncol(sparcc.final)) {
    if (abs(sparcc.boot.cor.mat[i,j]) >= sparcc.cutoff) {
      if (sparcc.boot.pval.mat[i,j] < pval.cutoff) {
        sparcc.final[i,j] = sparcc.boot.cor.mat[i,j]
      }
    }
  }
}

sparcc.final.1 <- sparcc.final
sparcc.final.1[is.na(sparcc.final.1)] <- 0
diag(sparcc.final.1) <- 1

# Remove non-sig correlations
index <- which(colnames(sparcc.final.1) %in% zero_colnames)
sparcc.final.1 <- sparcc.final.1[,-index]
index <- which(rownames(sparcc.final.1) %in% zero_rownames)
sparcc.final.1 <- sparcc.final.1[-index,]

# Build network from adjacency
sparcc.net <- graph.adjacency(sparcc.final.1, weighted = TRUE,
                                 mode = "undirected",
                                 diag = FALSE)

## Hub detection 
# Use sparcc.net or whichever choosing for the rest of the method
net <- sparcc.net

# Hub detection
net.cn <- closeness(net)
net.bn <- betweenness(net)
net.pr <- page_rank(net)$vector
net.hs <- hub_score(net)$vector

# Sort the species based on hubbiness score
net.hs.sort <- sort(net.hs, decreasing = TRUE)
# Choose the top 5 keystone species
net.hs.top10 <- head(net.hs.sort, n = 10)

## Cluster detection
# Get clusters
set.seed(123) # For reproductibility
wt.com <- walktrap.community(net)

# Get membership of walktrap clusters
membership(wt.com)

# Calculate modularity
modularity(net, membership(wt.com))

### Find articulation points
net.ap <- articulation.points(net)

## qgraph method: Plot nodes' degrees, closeness and betweenness
centralityPlot(net, include = c("Closeness", "Betweenness","Degree"), decreasing = TRUE, )

# Execute this command after running Function 3
plot.net.cls(net, net.hs, wt.com, net.ap,
             outfile = "Network_hub_size_cluster", title = "ASV network")

# Export graph as graphml for Cytoscape
write.graph(net, file = "Final_ASV_network.graphml", format = "graphml")

# Export table of network and the edge weight
network_df <- as.data.frame(get.edgelist(net))
network_df$edgeweights = get.edge.attribute(net)$weight
colnames(network_df) <- c("Node1","Node2","weight")
write.table(network_df, file = "Network_edge_table.tsv", sep = "\t", quote = F, row.names = F)
