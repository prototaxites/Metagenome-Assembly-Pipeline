#!/usr/bin/env Rscript
if (!require("pacman", quietly = TRUE)) 
  install.packages("pacman")
pacman::p_load(TreeTools, ggnewscale, RColorBrewer, optparse, 
               phytools, svglite)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require(ggtree))
  BiocManager::install("ggtree")

library(TreeTools,ggnewscale, RColorBrewer, optparse, 
        phytools, svglite)

library(ggtree)
###################
# Functions #
getDescendantTips <- function(tree, node){
  all_descendants <- getDescendants(tree,node)
  return(tree$tip.label[all_descendants[all_descendants <= length(tree$tip.label)]])
}
###################
# tree <- "gtdbtk.backbone.bac120.classify.tree"
# metadata <- "bin_data.full.taxcols.csv"
# outfile <- "metadata_cladogram.png"

option_list = list(
  make_option(c("-t", "--tree"), type="character", default=NA, 
              help="Input tree (required)", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NA, 
              help="Metadata associated with tips (csv) (required)", 
              metavar="character"),
  make_option(c("-o", "--outfile"), type="character", 
              default="./metadata_cladogram.png", 
              help="Output file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (!is.na(opt$tree) & file.exists(opt$tree)) {
  tree <-  ape::read.tree(opt$tree)
} else {
  stop("Input tree not provided or does not exist. See script usage (--help)")
}
if (!is.na(opt$metadata) & file.exists(opt$metadata)) {
  data <- read.csv(opt$metadata, header = TRUE)
} else {
  stop("Metadata not provided. See script usage (--help)")
}

# tree <- ape::read.tree(tree)
# data <- read.csv(metadata, header = TRUE)


row.names(data) <- data$bin_id
trimmed_tree <- KeepTip(tree, data$bin_id)
data <- data[data$bin_id %in% TipLabels(trimmed_tree),]

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


n <- length(unique(data$phylum)) + 
  length(unique(data$class))
colours <- sample(col_vector, n)
i <- 0
taxon_nodes <- data.frame(level = character(), Taxon = character(), 
                           node = numeric(), colour = character())
for(phylum in unique(data$phylum)){
  i <- i + 1
  color <- colours[i]
  node <- MRCA(trimmed_tree, data$bin_id[data$phylum == phylum])
  if(length(getDescendantTips(trimmed_tree, node)) == length(data$bin_id[data$phylum == phylum])){
    new_line <- data.frame(level = "phylum", Taxon = phylum, node = node, colour = color)
    taxon_nodes <- rbind(new_line, taxon_nodes)
  } else {
    for(tax in data$bin_id[data$phylum == phylum]){
      node <- which(trimmed_tree$tip.label == tax)
      new_line <- data.frame(level = "phylum", Taxon = phylum, node = node, colour = color)
      taxon_nodes <- rbind(new_line, taxon_nodes)
    }
  }
}
for(class in unique(data$class)){
  i <- i + 1
  color <- colours[i]
  node <- MRCA(trimmed_tree, data$bin_id[data$class == class])
  if(length(getDescendantTips(trimmed_tree, node)) == length(data$bin_id[data$class == class])){
    new_line <- data.frame(level = "class", Taxon = class, node = node, colour = color)
    if(node %in% taxon_nodes$node){
      taxon_nodes <- taxon_nodes[which(taxon_nodes$node != node),]
    }
    taxon_nodes <- rbind(new_line, taxon_nodes)
  } else {
    for(tax in data$bin_id[data$class == class]){
      node <- which(trimmed_tree$tip.label == tax)
      new_line <- data.frame(level = "class", Taxon = class, node = node, colour = color)
      if(node %in% taxon_nodes$node){
        taxon_nodes <- taxon_nodes[which(taxon_nodes$node != node),]
      }
      taxon_nodes <- rbind(new_line, taxon_nodes)
    }
  }
}
# for(order in unique(data$order)){
#   i <- i + 1
#   color <- colours[i]
#   node <- MRCA(trimmed_tree, data$bin_id[data$order == order])
#   new_line <- data.frame(level = "order", Taxon = order, node = node, colour = color)
#   if(node %in% taxon_nodes$node){
#     taxon_nodes <- taxon_nodes[which(taxon_nodes$node != node),]
#   }
#   taxon_nodes <- rbind(new_line, taxon_nodes)
# }

data$bin_type_binary <- 0
data$bin_type_binary[data$bin_type == "MAG"] <- 1
data$fully_circ <- 0
data$fully_circ[data$contigs == data$circular] <- 1
data$size_mbp <- data$size / 10^6
data$log_coverage <- log10(data$mean_coverage)
data$mag_circ <- "Binned Assembly"
data$mag_circ[data$bin_type == "MAG" & data$fully_circ == 0 ] <- "MAG"
data$mag_circ[data$bin_type == "MAG" & data$fully_circ == 1 ] <- "Fully-Circularized MAG"


circ <- ggtree(trimmed_tree, layout = 'circular', branch.length = 'none')
if(length(TipLabels(trimmed_tree)) < 30){
  circ <- ggtree(trimmed_tree, branch.length = 'none')
}


width <- .04
p1 <- gheatmap(circ, data[,which(colnames(data) == "size_mbp"),drop=FALSE], 
               offset=0, width=width,
               colnames = FALSE) +
  scale_fill_gradient(low = "#cce7e8", high = "#042f66", name="Size (Mbp)",
                      na.value = "white")
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, data[,which(colnames(data) == "log_coverage"),drop=FALSE], 
               offset=width * 10, width=width,
               colnames = FALSE) +
  scale_fill_gradient(low = "#ECF1D4", high = "#C70E03", name="log10(Coverage)",
                      na.value = "white")

p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, data[,which(colnames(data) == "mag_circ"),drop=FALSE], 
               offset=2 * width * 10, width=.04,
               colnames = FALSE) +
  scale_fill_manual(name = "", 
                    breaks = c("MAG", "Fully-Circularized MAG"),
                    values = c("Binned Assembly" = "white", 
                               "MAG" = "grey", 
                               "Fully-Circularized MAG" = "black"),
                    na.value = "white")

p4 <- p3 + new_scale_fill()
p4 <- p4 + geom_hilight(data = taxon_nodes, aes(node = node, fill = Taxon),  
                            to.bottom = TRUE) +
  scale_fill_manual(values=taxon_nodes$colour)
ggsave(opt$outfile, p4, width = 10, height = 10)

