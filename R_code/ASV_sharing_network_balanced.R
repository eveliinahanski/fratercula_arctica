# Load libraries
library(phyloseq)
library(igraph)
library(dplyr)
library(tidyr)

set.seed(295)

# Load phyloseq object
phy <- readRDS("Puffin_phyloseq.RDS")

# Only include chicks
phy <- subset_samples(phy, Age_clean == "Chick")

# Define colony map for plotting
colony_map <- c(
  "Skomer" = "1",
  "Shiants" = "2",
  "Sule_Skerry" = "3",
  "Orkney" = "4",
  "Fair_Isle" = "5",
  "Isle_of_May" = "6"
)

# Function to compute average richness and Jaccard index between populations
compute_asv_metrics <- function(phy, pops, n_iter = 100, n_samples = 4) {
  otu <- as.data.frame(otu_table(phy))
  if (taxa_are_rows(phy)) otu <- t(otu)
  sample_info <- data.frame(sample_data(phy))
  
  all_pops <- intersect(unique(sample_info$Colony), pops)
  richness_list <- list()
  jaccard_list <- list()
  
  for (iter in 1:n_iter) {
    subsamples <- lapply(all_pops, function(pop) {
      candidates <- rownames(sample_info[sample_info$Colony == pop, ])
      if (length(candidates) >= n_samples) {
        sample(candidates, n_samples)
      } else {
        NA
      }
    })
    names(subsamples) <- all_pops
    valid_pops <- names(subsamples)[sapply(subsamples, function(x) !any(is.na(x)))]
    
    if (length(valid_pops) < 2) next
    
    otu_bin <- as.matrix(otu)
    otu_bin[otu_bin > 0] <- 1
    
    richness_iter <- sapply(valid_pops, function(pop) {
      samples <- subsamples[[pop]]
      sum(colSums(otu_bin[samples, , drop = FALSE]) > 0)
    })
    
    jaccard_iter <- matrix(NA, length(valid_pops), length(valid_pops),
                           dimnames = list(valid_pops, valid_pops))
    
    for (i in seq_along(valid_pops)) {
      for (j in seq_along(valid_pops)) {
        if (i < j) {
          s1 <- subsamples[[valid_pops[i]]]
          s2 <- subsamples[[valid_pops[j]]]
          asvs1 <- colSums(otu_bin[s1, , drop = FALSE]) > 0
          asvs2 <- colSums(otu_bin[s2, , drop = FALSE]) > 0
          intersect_count <- sum(asvs1 & asvs2)
          union_count <- sum(asvs1 | asvs2)
          jaccard <- ifelse(union_count > 0, intersect_count / union_count, 0)
          jaccard_iter[i, j] <- jaccard
          jaccard_iter[j, i] <- jaccard
        }
      }
    }
    
    richness_list[[iter]] <- richness_iter
    jaccard_list[[iter]] <- jaccard_iter
  }
  
  avg_richness <- Reduce("+", richness_list) / length(richness_list)
  avg_jaccard <- Reduce("+", jaccard_list) / length(jaccard_list)
  
  list(richness = avg_richness, jaccard = avg_jaccard)
}

# Function to generate network plot
generate_network_plot <- function(richness, jaccard, color, title, output_path) {
  diag(jaccard) <- 0
  valid_nodes <- rowSums(jaccard, na.rm = TRUE) > 0 | colSums(jaccard, na.rm = TRUE) > 0
  jaccard <- jaccard[valid_nodes, valid_nodes]
  richness <- richness[names(richness) %in% rownames(jaccard)]
  
  if (any(!is.finite(jaccard))) stop("Jaccard matrix has non-finite values.")
  if (nrow(jaccard) < 2) stop("Too few connected populations to plot.")
  
  # Node size scaling
  min_size <- 20
  max_size <- 50
  norm_size <- min_size + (richness - min(richness)) / (max(richness) - min(richness)) * (max_size - min_size)
  
  g <- graph_from_adjacency_matrix(jaccard, mode = "undirected", weighted = TRUE)
  
  # Edge thickness based on Jaccard index
  min_edge_width <- 1.5
  max_edge_width <- 10
  edge_weights <- E(g)$weight
  E(g)$width <- min_edge_width + (edge_weights - min(edge_weights)) / (max(edge_weights) - min(edge_weights)) * (max_edge_width - min_edge_width)
  
  V(g)$name <- colony_map[V(g)$name]
  layout <- layout_with_fr(g, niter = 1000)
  layout <- layout.norm(layout, -1, 1)
  
  png(output_path, width = 1200, height = 1200, res = 150)
  plot(
    g,
    layout = layout,
    vertex.label = V(g)$name,
    vertex.label.color = "white",
    vertex.label.cex = 2,
    vertex.size = norm_size,
    vertex.color = color,
    vertex.frame.color = color,
    edge.width = E(g)$width,
    edge.color = adjustcolor(color, alpha.f = 0.5),
    edge.curved = 0.2,
    main = title
  )
  dev.off()
}

# Subset for aerobes and anaerobes
aero <- read.csv("Puffin_taxa_aerotolerance.csv")
obligate_anaerobe_genera <- unique(aero$Genus[aero$Aerotolerance_cat == "obligate anaerobe"])
obligate_aerobe_genera <- unique(aero$Genus[aero$Aerotolerance_cat == "aerotolerant"])

anaerobes <- subset_taxa(phy, Genus %in% obligate_anaerobe_genera)
aerobes <- subset_taxa(phy, Genus %in% obligate_aerobe_genera)

# Common populations
pops <- intersect(unique(sample_data(phy)$Colony), names(colony_map))

# Plot 1: Full microbiota
metrics_all <- compute_asv_metrics(phy, pops)
generate_network_plot(
  richness = metrics_all$richness,
  jaccard = metrics_all$jaccard,
  color = "#655f12",
  title = "ASV Sharing Network Between Colonies - Full Microbiota",
  output_path = "asv_network_all.png"
)

# Plot 2: Anaerobes
metrics_anaerobes <- compute_asv_metrics(anaerobes, pops)
generate_network_plot(
  richness = metrics_anaerobes$richness,
  jaccard = metrics_anaerobes$jaccard,
  color = "#f45751",
  title = "ASV Sharing Network Between Colonies - Anaerobes",
  output_path = "asv_network_anaerobes.png"
)

# Plot 3: Aerobes
metrics_aerobes <- compute_asv_metrics(aerobes, pops)
generate_network_plot(
  richness = metrics_aerobes$richness,
  jaccard = metrics_aerobes$jaccard,
  color = "#4d8190",
  title = "ASV Sharing Network Between Colonies - Aerobes",
  output_path = "asv_network_aerobes.png"
)
