#-------------------------------------------
# ASV level
abund_file = "feature-table.tsv"
data <- read.csv(file = abund_file, sep = "\t", header = T, row.names = 1)
colnames(data) <- rownames(meta) # Set rownames to match
index <- seq(1,nrow(data), by=1)
index <- as.vector(sapply(index, pasting))
asv_id <- as.data.frame(index) # ASV_id and respective sequences
asv_id$seq <- rownames(data)

pasting <- function(x){
  paste("ASV", sep = "_", x)
}

rownames(data) <- asv_id$index

# Phyloseq alpha div analysis
library(phyloseq)
phy.obj <- phyloseq(otu_table(data, taxa_are_rows = T), sample_data(meta))

## Rarefy for alpha-div
set.seed(223)
phy.obj.rare = rarefy_even_depth(phy.obj, sample.size = min(sample_sums(phy.obj)))
alpha_div = plot_richness(phy.obj.rare, shape = "Genotype", color = "Genotype", x = "Genotype",
                          measures = c("Observed", "Shannon", "InvSimpson","Fisher"), 
                          title = "Alpha Diversity of mGlu5 mice microbiome") +  scale_color_manual(values = c('#388ECC','#F68B33'))
alpha_div 

# Store as a new data variable
alphadt = data.table::as.data.table(alpha_div$data)
# Subset to individual alpha diversity index
alphadt_obs <- alphadt[(variable == "Observed")]
kruskal.test(value ~ Genotype, data = alphadt_obs)
shapiro.test(alphadt_obs$value)
var.test(value ~ Genotype, data = alphadt_obs)

obs.p <- t.test(value ~ Genotype, data = alphadt_obs, var.equal = TRUE)
obs.p$p.value


alphadt_shan <- alphadt[(variable == "Shannon")]
kruskal.test(value ~ Genotype, data = alphadt_shan)

alphadt_InvSimpson <- alphadt[(variable == "InvSimpson")]
kruskal.test(value ~ Genotype, data = alphadt_InvSimpson)

alphadt_fisher <- alphadt[(variable == "Fisher")]
kruskal.test(value ~ Genotype, data = alphadt_fisher)


# Phyloseq beta div analysis
phy.obj.ra = transform_sample_counts(phy.obj, function(x) x / sum(x) )

### Jaccard distance
distjac = distance(phy.obj.ra, method = "jaccard")
ordjac = ordinate(phy.obj.ra, method = "PCoA", distance = distjac)
plot_ordination(phy.obj.ra, ordjac, color = "Genotype", shape = "Genotype") +
  ggtitle("Jaccard PCoA of the mGlu5 mice") +
  geom_jitter(aes(color = Genotype)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(name = "Component 1") +
  scale_y_continuous(name = "Component 2") +
  theme(legend.position = "right")

### Bray Curtis distance
distBC = distance(phy.obj.ra, method = "bray")
ordBC = ordinate(phy.obj.ra, method = "PCoA", distance = distBC)

#jpeg('images/Bray_curtis.jpg', units = "in", width = 9, height = 4.5, res = 300)
bray_pcoa <- plot_ordination(phy.obj.ra, ordBC, color = "Genotype", shape = "Genotype") +
  ggtitle("Bray Curtis PCoA of mGlu5 mice") +
  geom_jitter(aes(color = Genotype)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c('#388ECC','#F68B33'))

# Unweighted and weighted UniFrac analysis
## Make rooted tree
tree <- ape::read.tree("tree.nwk")
data1 <- data
rownames(data1) <- asv_id$seq
data1 <- t(data1)
data1.filter <- low.count.removal(data1, percent=0.01)
data1.filter = data1.filter$data.filter
data1.TSS <- t(apply(data1.filter, 1, TSS.divide))
phy.obj.ra1 = phyloseq(otu_table(data1.TSS, taxa_are_rows = F), sample_data(meta),phy_tree(tree))

### Unweighted unifrac distance
distUF = UniFrac(phy.obj.ra, weighted = F, normalized = T, parallel = F)
ordUF = ordinate(phy.obj.ra, method = "PCoA", distance = distUF)

#jpeg('UniFrac.jpg', units = "in", width = 9, height = 4.5, res = 300)
uf_pcoa <- plot_ordination(phy.obj.ra1, ordUF, color = "Genotype", shape = "Genotype") +
  ggtitle("Unweighted UniFrac PCoA of mGlu5 mice") +
  geom_jitter(aes(color = Genotype)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c('#388ECC','#F68B33'))
#scale_x_continuous(name = "Component 1") +
#scale_y_continuous(name = "Component 2") + 

### Weighted unifrac distance
distwUF = UniFrac(phy.obj.ra1, weighted = T, normalized = T, parallel = F)
ordwUF = ordinate(phy.obj.ra1, method = "PCoA", distance = distwUF)

#jpeg('wUniFrac.jpg', units = "in", width = 9, height = 4.5, res = 300)
plot_ordination(phy.obj.ra1, ordwUF, color = "Genotype", shape = "Genotype") +
  ggtitle("Weighted UniFrac PCoA of mGlu5 mice") +
  geom_jitter(aes(color = Genotype)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(name = "Component 1") +
  scale_y_continuous(name = "Component 2") 

vegan::adonis(distwUF ~ Genotype, data = meta)
vegan::adonis(distUF ~ Genotype, data = meta)
vegan::adonis(distBC ~ Genotype, data = meta)
vegan::adonis(distjac ~ Genotype, data = meta)