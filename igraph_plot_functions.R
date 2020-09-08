plot.net <- function(net, scores, outfile, title) {
  # Convert node label from names to numerical IDs.
  features <- V(net)$name
  col_ids <- seq(1, length(features))
  V(net)$name <- col_ids
  node.names <- features[V(net)]
  # Nodes' color.
  V(net)$color <- "white"
  # Define output image file.
  outfile <- paste(outfile, "jpg", sep=".")
  # Image properties.
  jpeg(outfile, width = 4800, height = 9200, res = 300, quality = 100)
  par(oma = c(4, 1, 1, 1))
  # Main plot function.
  plot(net, vertex.size = (scores*5)+4, vertex.label.cex = 1)
  title(title, cex.main = 4)
  # Plot legend containing OTU names.
  labels = paste(as.character(V(net)), node.names, sep = ") ")
  legend("bottom", legend = labels, xpd = TRUE, ncol = 5, cex = 1.2)
  dev.off()
}

# Function 3: Plot network with clusters and node size scaled to hubbiness
plot.net.cls <- function(net, scores, cls, AP, outfile, title) {
  # Get size of clusters to find isolated nodes.
  cls_sizes <- sapply(igraph::groups(cls), length)
  # Randomly choosing node colors. Users can provide their own vector of colors.
  #colors <- sample(colours(), length(cls))
  colours = c("darkslategray3","chocolate3","darkolivegreen3","darkorchid3",
              "antiquewhite3","deeppink3","goldenrod3", "royalblue3","springgreen3")
  # Nodes in clusters will be color coded. Isolated nodes will be white.
  V(net)$color <- sapply(membership(cls),
                         function(x) {ifelse(cls_sizes[x]>1,
                                             colours[x], "white")})
  # Convert node label from names to numerical IDs.
  node.names <- V(net)$name
  col_ids <- seq(1, length(node.names))
  V(net)$name <- col_ids
  # To draw a halo around articulation points.
  AP <- lapply(names(AP), function(x) x)
  marks <- lapply(1:length(AP), function(x) which(node.names == AP[[x]]))
  # Define output image file.
  outfile <- paste(outfile, "jpg", sep=".")
  # Image properties.
  jpeg(outfile, width = 4800, height = 4800, res = 300, quality = 100)
  par(oma = c(4, 1, 1, 1))
  # Customized layout to avoid nodes overlapping.
  e <- get.edgelist(net)
  class(e) <- "numeric"
  ew <- get.edge.attribute(net)$weight
  l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(net),
                                         area=2*(vcount(net)^2),
                                         repulse.rad=(vcount(net)^3.1))
  # Main plot function.
  plot(net, vertex.size = (scores*15)+4, vertex.label.cex=2.5,
       width = ew,
       labels.font = 2,
       labels.dist = 1,
       vertex.label = node.names,
       vertex.label.color = "black",
       #mark.border="black",
       #mark.groups = marks,
       #mark.col = "white",
       #mark.expand = 10,
       #mark.shape = 1,
       layout=l)
  title(title, cex.main=4)
  # Plot legend containing OTU names.
  #labels = paste(as.character(V(net)), node.names, sep = ") ")
  #legend("bottom", legend = labels, xpd = TRUE, ncol = 5, cex = 1.2)
  dev.off()
}
