library(ggplot2)
library(dplyr)

tidy_data_barplot <- function(data, meta_file, factor1){
  # Sort by abundance, Samples as rows, OTU as columns
  index = order(colSums(data),decreasing = T)
  data <- data[,order(colSums(data),decreasing = T)]
  data_tidy = cbind(as.data.frame(data), meta_file[factor1])
  colnames(data_tidy)[ncol(data_tidy)] = factor1
  #data_tidy[factor1] = factor(metadata[factor1], levels = f1_levels)
  #data_tidy = cbind(as.data.frame(data),as.character(meta$Genotype),as.character(meta$Age),as.character(meta$MouseID))
  #colnames(data_tidy)[(ncol(data_tidy)-2):ncol(data_tidy)] = c(factor1)
  #data_tidy$Genotype = factor(meta$Genotype, levels = c("WT","HD"))
  #data_tidy$Age = factor(meta$Age, levels = c("4","6","8","10",'12'))
  return(data_tidy)
  }

plotIndiv_bar <- function(data, x, y, title, ylabel) {
  ggplot(data, aes(x=!!enquo(x), y=!!enquo(y), fill=!!enquo(x))) + 
    geom_boxplot(color="black") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(enquo(title)) +
    ylab(enquo(ylabel)) +
    scale_x_discrete(element_blank()) +
    theme(legend.position = "none") + 
    #scale_fill_manual(values = color_pick ) + #= c('#388ECC','#F68B33')
    geom_jitter()
}

low.count.removal = function(
  data, # OTU count data frame of OTU as rows x samples as columns
  percent=0.01 # cutoff chosen
){
  keep.otu = which(rowSums(data)*100/(sum(rowSums(data))) > percent)
  data.filter = data[keep.otu,]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# Convert to relative abundance
TSS.divide = function(x){
  (x/sum(x))*100
}

# Cleanup raw data to return matrix of samples as rows, OTU as columns in relative abundance + filtered + offset
data_cleanup <- function(counts, samples_as_rows = TRUE, filter_percent, offset = TRUE, offset_number = 0.01){
  if (samples_as_rows == TRUE) {
    counts = t(counts)
  } else {
    counts = counts
  }
  data.filter = low.count.removal(counts, percent = filter_percent)
  data.filter = data.filter$data.filter
  if (offset == TRUE) {
    data.filter = data.filter + offset_number
  } else {
    data.filter = data.filter
  }
 
  data.TSS = t(apply(data.filter, 2, TSS.divide))  # function is applied to each row (i.e. each sample)
  return(data.TSS)
}

# Data transformation with CLR
data.transf <- function(data.TSS){ # take sample as rows, OTUs as columns
  data.clr = logratio.transfo(data.TSS, logratio = "CLR", offset = 0.1)
  data.good <- as.data.frame(matrix(ncol = ncol(data.clr), 
                                nrow = nrow( data.clr)))
  rownames(data.good) <- rownames(data.clr)
  colnames(data.good) <- colnames(data.clr)
  for( i in c(1:nrow(data.clr))){
     for( j in c(1:ncol(data.clr))){
       data.good[i,j] <- data.clr[i,j]
     }
  }
  return(data.good)
}

# alpha div data cleanup for kruskal wallis test
alpha_div_extract <- function(alphadt){
  # Subset to individual alpha diversity index
  obs <- alphadt[(variable == "Observed")]
  shan <- alphadt[(variable == "Shannon")]
  invsimpson <- alphadt[(variable == "InvSimpson")]
  fisher <- alphadt[(variable == "Fisher")]
  
  alpha_df <- data.frame(obs$value, shan$value, invsimpson$value, fisher$value, obs$Genotype)
  rownames(alpha_df) = obs$samples
  colnames(alpha_df) = c("Observed","Shannon","InvSimpson","Fisher","Genotype")
  return(alpha_df)
}
