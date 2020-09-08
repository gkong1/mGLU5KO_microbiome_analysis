setwd("~/phd/carol - mGlu5_microbiome/R_analysis/")

source("~/phd/scripts/common_functions.R")
library(mixOmics)

abund_file = 
meta_file = "sample_metadata.txt"

data <- read.csv(file = "l7_matrix.csv", sep = ",", header = T, row.names = 1)
meta <- read.csv(file = meta_file, sep = "\t", header = T, row.names = 1)
meta$Genotype = factor(meta$Genotype, levels = c("WT","KO"))

index <- which(colnames(data) == "X")
l7data <- data[,1:(index - 1)]

l5data <- data[,1:27]
rownames(meta) == rownames(l7data) # Make sure sample names match up before proceeding

#----------------------------------------------------------------------
# Phylum level
# Data filter, TSS transf
l2data_clean <- data_cleanup(l2data, filter_percent = 0.01, offset = FALSE)

# Data prep for barplots
l2data_tidy <- tidy_data_barplot(l2data_clean, meta_file = meta, factor1 = "Genotype")
bacteroidetes <- plotIndiv_bar(l2data_tidy, x = Genotype, 
                               y = k__Bacteria.p__Bacteroidetes,
                               title = "Bacteroidetes",
                               ylabel = "Relative abundance (%)",
                               color_pick = c("darkgreen","darkmagenta"))

verrucomicrobia <- plotIndiv_bar(l2data_tidy, x = Genotype, 
                                 y = k__Bacteria.p__Verrucomicrobia,
                                 title = "verrucomicrobia",
                                 ylabel = "Relative abundance (%)",
                                 color_pick = c("darkgreen","darkmagenta"))

firmicutes <- plotIndiv_bar(l2data_tidy, x = Genotype, 
                            y = k__Bacteria.p__Firmicutes,
                            title = "Firmicutes",
                            ylabel = "Relative abundance (%)",
                            color_pick = c("darkgreen","darkmagenta"))

# One-way ANOVA test:
l2.good <- data.transf(l2data_clean)
l2.good$Genotype = meta$Genotype

### Firmicutes
shapiro.test(l2.good$k__Bacteria.p__Firmicutes)
var.test(k__Bacteria.p__Firmicutes ~ Genotype, data = l2.good)
res.l2.firm <- t.test(k__Bacteria.p__Firmicutes ~ Genotype, data = l2.good, var.equal = TRUE)
res.l2.firm

kruskal.test(k__Bacteria.p__Firmicutes ~ Genotype, data = l2.good)

#res.l2.firm <- aov(k__Bacteria.p__Firmicutes ~ Genotype, data = l2data_clr)
#summary(res.l2.firm)

### Bacteroidetes
shapiro.test(l2.good$k__Bacteria.p__Bacteroidetes)
var.test(k__Bacteria.p__Bacteroidetes ~ Genotype, data = l2.good)
res.l2.bac <- t.test(k__Bacteria.p__Bacteroidetes ~ Genotype, data = l2.good, var.equal = TRUE)
res.l2.bac

kruskal.test(k__Bacteria.p__Bacteroidetes ~ Genotype, data = l2.good)

### Verrucomicrobia
shapiro.test(l2.good$k__Bacteria.p__Verrucomicrobia)
var.test(k__Bacteria.p__Verrucomicrobia ~ Genotype, data = l2.good)
res.l2.veru <- t.test(k__Bacteria.p__Verrucomicrobia ~ Genotype, data = l2.good, var.equal = TRUE)
res.l2.veru

kruskal.test(k__Bacteria.p__Verrucomicrobia ~ Genotype, data = l2.good)

phyla_p <- c("0.1392","0.8306" ,"0.05923")
phyla_padj <- p.adjust(phyla_p, method="fdr")

mean(l2.data.TSS.1$k__Bacteria.p__Bacteroidetes)
mean(l2.data.TSS.1$k__Bacteria.p__Verrucomicrobia)
mean(l2.data.TSS.1$k__Bacteria.p__Firmicutes)
#----------------------------------------------------------------------------
# Family level
# Data filter, TSS transf
l5data_clean <- data_cleanup(l5data, filter_percent = 0.01, offset = FALSE)

# Data prep for barplots
l5data_tidy <- tidy_data_barplot(l5data_clean, meta_file = meta, factor1 = "Genotype")
s24.7 <- plotIndiv_bar(l5data_tidy, x = Genotype, 
                               y = k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__S24.7,
                               title = "S24.7",
                               ylabel = "Relative abundance (%)",
                               color_pick = c("darkgreen","darkmagenta"))

verrucomicrobiaceae <- plotIndiv_bar(l5data_tidy, x = Genotype, 
                       y = k__Bacteria.p__Verrucomicrobia.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae,
                       title = "Verrucomicrobiaceae",
                       ylabel = "Relative abundance (%)",
                       color_pick = c("darkgreen","darkmagenta"))

bacteroidaceae <- plotIndiv_bar(l5data_tidy, x = Genotype, 
                    y = k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae,
                    title = "Bacteroidaceae",
                    ylabel = "Relative abundance (%)",
                    color_pick = c("darkgreen","darkmagenta"))

prevotellaceae <- plotIndiv_bar(l5data_tidy, x = Genotype, 
                                y = k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae,
                                title = "Bacteroidaceae",
                                ylabel = "Relative abundance (%)",
                                color_pick = c("darkgreen","darkmagenta"))

erysipelotrichaceae <- plotIndiv_bar(l5data_tidy, x = Genotype, 
                                y = k__Bacteria.p__Firmicutes.c__Erysipelotrichi.o__Erysipelotrichales.f__Erysipelotrichaceae,
                                title = "Erysipelotrichaceae",
                                ylabel = "Relative abundance (%)",
                                color_pick = c("darkgreen","darkmagenta"))

# One-way ANOVA test:
l5.good <- data.transf(l5data_clean)
l5.good$Genotype = meta$Genotype

### S24.7
shapiro.test(l5.good$k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__S24.7)
var.test(k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__S24.7 ~ Genotype, data = l5.good)

res.l5.s24.7 <- t.test(k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__S24.7 ~ Genotype, 
                   data = l5.good, var.equal = TRUE)
res.l5.s24.7$p.value

### Verrucomicrobiaceae
shapiro.test(l5.good$k__Bacteria.p__Verrucomicrobia.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae)
var.test(k__Bacteria.p__Verrucomicrobia.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae ~ Genotype, data = l5.good)

res.l5.verrucom <- t.test(k__Bacteria.p__Verrucomicrobia.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae ~ Genotype, 
                  data = l5.good, var.equal = TRUE)
res.l5.verrucom

### bacteroidaceae
shapiro.test(l5.good$k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae)
var.test(k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae ~ Genotype, data = l5.good)

res.l5.bact <- t.test(k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae ~ Genotype, 
                   data = l5.good, var.equal = TRUE)
res.l5.bact$p.value

kruskal.test(k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae ~ Genotype, 
             data = l5.good)
wilcox.test(k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae ~ Genotype, 
             data = l5.good)

### prevotellaceae
shapiro.test(l5.good$k__Bacteria.p__Firmicutes.c__Erysipelotrichi.o__Erysipelotrichales.f__Erysipelotrichaceae)
var.test(k__Bacteria.p__Firmicutes.c__Erysipelotrichi.o__Erysipelotrichales.f__Erysipelotrichaceae ~ Genotype, data = l5.good)

res.l5.prev <- t.test(k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae ~ Genotype, 
                   data = l5.good, var.equal = TRUE)
res.l5.prev$p.value

### erysipelotrichaceae
shapiro.test(l5.good$k__Bacteria.p__Firmicutes.c__Erysipelotrichi.o__Erysipelotrichales.f__Erysipelotrichaceae)
var.test(k__Bacteria.p__Firmicutes.c__Erysipelotrichi.o__Erysipelotrichales.f__Erysipelotrichaceae ~ Genotype, data = l5.good)

res.l5.erys <- t.test(k__Bacteria.p__Firmicutes.c__Erysipelotrichi.o__Erysipelotrichales.f__Erysipelotrichaceae ~ Genotype, 
                   data = l5.good, var.equal = TRUE)
res.l5.erys$p.value

fam_p <- c(res.l5.s24.7$p.value, res.l5.verrucom$p.value, res.l5.bact$p.value, res.l5.prev$p.value, res.l5.erys$p.value)
fam_padj <- p.adjust(fam_p, method="BH")

mean(l5.data.TSS.1$k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__S24.7)
mean(l5.data.TSS.1$k__Bacteria.p__Verrucomicrobia.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae)
mean(l5.data.TSS.1$k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae)
mean(l5.data.TSS.1$k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae)
mean(l5.data.TSS.1$k__Bacteria.p__Firmicutes.c__Erysipelotrichi.o__Erysipelotrichales.f__Erysipelotrichaceae)

fam_pval_mat <- matrix(nrow = 2, ncol = 23)
colnames(fam_pval_mat) = colnames(l5.data.TSS.1)[1:23]
rownames(fam_pval_mat) = c("p.val","p.adj")

for (j in 1:ncol(fam_pval_mat)) {
  index <- which(l5.data.TSS.1$Genotype == "WT")
  wt.value <- l5.data.TSS.1[index,j]
  index <- which(l5.data.TSS.1$Genotype == "KO")
  ko.value <- l5.data.TSS.1[index,j]
  res <- t.test(wt.value, ko.value, var.equal = TRUE)
  fam_pval_mat[1,j] <- res$p.value
}

fam_padj_mat <- p.adjust(fam_pval_mat[1,], method="BH")
fam_pval_mat[2,] <- fam_padj_mat
which(fam_padj_mat < 0.05)
which(fam_pval_mat[1,] < 0.1)
#-------------------------------------------------------------
# Species level
l7data_clean <- data_cleanup(l7data, filter_percent = 0.01, offset = FALSE)

# Data prep for barplots
l7data_tidy <- tidy_data_barplot(l7data_clean, meta_file = meta, factor1 = "Genotype")
akker.mucin <- plotIndiv_bar(l7data_tidy, x = Genotype, 
                       y = k__Bacteria.p__Verrucomicrobia.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.g__Akkermansia.s__muciniphila,
                       title = "Akkermansia muciniphila",
                       ylabel = "Relative abundance (%)",
                       color_pick = c("darkgreen","darkmagenta"))

# One-way ANOVA test:
l7data_clr <- data.transf(l7data_clean)
l7data_clr$Genotype = meta$Genotype

l7.data.clr <- as.matrix(l7.good[,1:40])
### Akkermansia muciniphila
shapiro.test(l7.good$k__Bacteria.p__Verrucomicrobia.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.g__Akkermansia.s__muciniphila)
var.test(k__Bacteria.p__Verrucomicrobia.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.g__Akkermansia.s__muciniphila ~ Genotype, data = l7.good)

res.l7.akkm <- t.test(k__Bacteria.p__Verrucomicrobia.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.g__Akkermansia.s__muciniphila ~ Genotype, 
                    data = l7.good, var.equal = TRUE)
res.l7.akkm$p.value

#----------------------------------------------------------
# ASV level
asv_file = "feature-table.tsv"
data <- read.csv(file = asv_file, sep = "\t", header = T, row.names = 1)
asv_data <- data
## Check sample names to match
colnames(data) <- rownames(meta) # Set rownames to match

## Change ASV names
index <- seq(1,nrow(data), by=1)
pasting <- function(x){
  paste("ASV", sep = "_", x)
}
index <- as.vector(sapply(index, pasting))
asv_id <- as.data.frame(index) # ASV_id and respective sequences
asv_id$seq <- rownames(data)

data1 <- data
rownames(data1) <- asv_id$seq
rownames(data) <- asv_id$index

asv_data <- data

# Phyloseq alpha div analysis
library(phyloseq)
phy.obj <- phyloseq(otu_table(asv_data, taxa_are_rows = T), sample_data(meta))

## Rarefy for alpha-div
set.seed(223)
phy.obj.rare = rarefy_even_depth(phy.obj, sample.size = min(sample_sums(phy.obj)))
alpha_div = plot_richness(phy.obj.rare, shape = "Genotype", 
                          color = "Genotype", x = "Genotype",
                          measures = c("Observed", "Shannon", "InvSimpson","Fisher"), 
                          title = "Alpha Diversity of mGlu5 mice microbiome") +  
                          scale_color_manual(values = c('darkgreen','darkmagenta'))
alpha_div # use geom_jitter(size = 3, alpha = 0.5) to show all points

jpeg(filename = "../images/Rarefaction_curve.jpg", width = 5, height = 4, units = "cm", res = 600)
vegan::rarecurve(t(otu_table(phy.obj)), step = 50, cex = 0.5, 
                 xlab = "Sequencing depth", ylab = "Number of observed ASV",
                 main = "Rarefaction curve")
abline(v=1100, col = "blue")
dev.off()

## Extracting alpha div data for kruskal wallis test
alphadt = data.table::as.data.table(alpha_div$data)
alphadf <- alpha_div_extract(alphadt)
colnames(alphadf)
obs_plot <- ggplot(data = alphadf, aes(x = Genotype, y = Observed, fill = Genotype)) + 
  geom_boxplot() + 
  geom_jitter(position = position_jitter(0.3)) +  
  scale_fill_manual(values = c('darkgreen','darkmagenta')) +
  ggtitle("Observed") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
shannon_plot <- ggplot(data = alphadf, aes(x = Genotype, y = Shannon, fill = Genotype)) + 
  geom_boxplot() + 
  geom_jitter(position = position_jitter(0.3)) +  
  scale_fill_manual(values = c('darkgreen','darkmagenta')) +
  ggtitle("Shannon") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
invs_plot <- ggplot(data = alphadf, aes(x = Genotype, y = InvSimpson, fill = Genotype)) + 
  geom_boxplot() + 
  geom_jitter(position = position_jitter(0.3)) +  
  scale_fill_manual(values = c('darkgreen','darkmagenta')) +
  ggtitle("Inverse Simpson") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
fisher_plot <- ggplot(data = alphadf, aes(x = Genotype, y = Fisher, fill = Genotype)) + 
  geom_boxplot() + 
  geom_jitter(position = position_jitter(0.3)) +  
  scale_fill_manual(values = c('darkgreen','darkmagenta')) +
  ggtitle("Fisher") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

jpeg(filename = "../images/Alpha_div_plot.jpg", height = 10, width = 18, units = "cm", res = 600)
gridExtra::grid.arrange(obs_plot, shannon_plot, invs_plot, fisher_plot, nrow = 1)
dev.off()

kruskal.test(Observed ~ Genotype, data = alphadf)
kruskal.test(Shannon ~ Genotype, data = alphadf)
kruskal.test(InvSimpson ~ Genotype, data = alphadf)
kruskal.test(Fisher ~ Genotype, data = alphadf)

# Phyloseq beta div analysis
## Make rooted tree
tree <- ape::read.tree("tree.nwk")

## Cleanup data and import to phyloseq
data_clean1 <- data_cleanup(data1, samples_as_rows = FALSE,filter_percent = 0.01, offset = FALSE)
phy.obj.ra = phyloseq(otu_table(data_clean1, taxa_are_rows = F), sample_data(meta),phy_tree(tree))

### Bray Curtis distance
distBC = distance(phy.obj.ra, method = "bray")
ordBC = ordinate(phy.obj.ra, method = "PCoA", distance = distBC)

bray_pcoa <- plot_ordination(phy.obj.ra, ordBC, color = "Genotype", shape = "Genotype") +
  ggtitle("Bray Curtis PCoA") +
  geom_jitter(aes(color = Genotype), size = 3.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c('darkgreen','darkmagenta'))

jpeg('../images/Bray_curtis.jpg', units = "in", width = 5.5, height = 4.5, res = 300)
bray_pcoa
dev.off()
vegan::adonis(distBC ~ Genotype, data = meta)

### Unweighted unifrac distance
distUF = UniFrac(phy.obj.ra, weighted = F, normalized = T, parallel = F)
ordUF = ordinate(phy.obj.ra, method = "PCoA", distance = distUF)

uf_pcoa <- plot_ordination(phy.obj.ra, ordUF, color = "Genotype", shape = "Genotype") +
  ggtitle("Unweighted UniFrac PCoA") +
  geom_jitter(aes(color = Genotype), size = 3.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c('darkgreen','darkmagenta'))

jpeg('../images/Unweighted_UniFrac.jpg', units = "in", width = 5.5, height = 4.5, res = 300)
uf_pcoa
dev.off()
vegan::adonis(distUF ~ Genotype, data = meta)

# Filter low count samples
asv.data.filter = low.count.removal(asv_data, percent = 0.01)
asv.data.filter = asv.data.filter$data.filter



#----------------------------------------------------------------------------------
# Differential abundance with ALDEx2
library(ALDEx2)

## Differential abundance in ASVs
### Create Aldex object
ASVdata.aldex <- aldex(asv.data.filter, meta$Genotype, mc.samples = 128, test = "t", effect = TRUE, 
                    include.sample.summary = FALSE, verbose = FALSE, denom = "all")
par(mfrow=c(1,2))
aldex.plot(ASVdata.aldex, type = "MA", test = "welch", 
           xlab = "Log-ratio abundance", ylab = "Difference")
aldex.plot(ASVdata.aldex, type = "MW", test = "welch", 
           xlab = "Dispersion", ylab = "Difference")

## Differential abundance in Phylum
### Create Aldex object
l2data.aldex <- aldex(l2data, meta$Genotype, mc.samples = 128, test = "t", effect = TRUE, 
                    include.sample.summary = FALSE, verbose = FALSE, denom = "all")
par(mfrow=c(1,2))
aldex.plot(l2data.aldex, type = "MA", test = "welch", 
           xlab = "Log-ratio abundance", ylab = "Difference")
aldex.plot(l2data.aldex, type = "MW", test = "welch", 
           xlab = "Dispersion", ylab = "Difference")

## Differential abundance in Family
l5data.aldex <- aldex(l5data, meta$Genotype, mc.samples = 128, test = "t", effect = TRUE, 
                    include.sample.summary = FALSE, verbose = FALSE, denom = "all")
par(mfrow=c(1,2))
aldex.plot(l5data.aldex, type = "MA", test = "welch", 
           xlab = "Log-ratio abundance", ylab = "Difference")
aldex.plot(l5data.aldex, type = "MW", test = "welch", 
           xlab = "Dispersion", ylab = "Difference")

## Differential abundance in Species
l7.data.filter <- low.count.removal(t(l7data), percent = 0.01)
l7.data.filter = l7.data.filter$data.filter 

l7data.aldex <- aldex(l7.data.filter, meta$Genotype, mc.samples = 128, test = "t", effect = TRUE, 
                      include.sample.summary = FALSE, verbose = FALSE, denom = "all")
par(mfrow=c(1,2))
aldex.plot(l7data.aldex, type = "MA", test = "welch", 
                              xlab = "Log-ratio abundance", ylab = "Difference")
aldex.plot(l7data.aldex, type = "MW", test = "welch", 
                              xlab = "Dispersion", ylab = "Difference")
