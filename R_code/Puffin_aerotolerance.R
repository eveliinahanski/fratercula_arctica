
aero_data <- read.csv("Puffin_taxa_aerotolerance.csv") 
aero_data<-aero_data[which(aero_data$Genus!=""),]
aero_data<-aero_data[which(aero_data$Genus!=" "),]
phy <- readRDS("Puffin_phyloseq.RDS")

# Load necessary libraries
library(phyloseq)
library(dplyr)
library(tidyr)
library(microbiome)  # For transformation

# Ensure OTU table sums to 1 per sample (before processing)
unique(rowSums(otu_table(phy)))  # Should be all 1s

# Aggregate at Genus Level (Ensure each Genus only appears once per Sample)
phy2 <- phy %>%
  aggregate_taxa(level = "Genus") %>%
  microbiome::transform(transform = "compositional")  # Ensure proportions sum to 1

# Convert to data frame
phy_melt <- psmelt(phy2)

# **Check if Genus is duplicated within a Sample**
dup_check <- phy_melt %>%
  group_by(Sample, Genus) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  warning("Duplicates detected: Some genera appear multiple times per sample!")
} else {
  print("No duplicate Genus per Sample. Proceeding...")
}

# **Aggregate abundance at Genus level before merging with aerotolerance**
phy_agg <- phy_melt %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

# **Merge with aerotolerance categories**
phyp <- merge(phy_agg, aero_data, by = "Genus", all.x = TRUE, all.y = FALSE)

# Select relevant columns for analysis
phe <- phyp[, c("Sample", "Abundance", "Aerotolerance_cat")]

# **Check if total abundance is still 1 per sample (before aerotolerance grouping)**
test_abundance <- phe %>%
  group_by(Sample) %>%
  summarise(total_abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance))

print(test_abundance)  # Should all be 1

# **Summarise abundance per Sample across Aerotolerance Categories**
aero_summary <- phe %>%
  group_by(Sample) %>%
  summarise(
    aerotolerant = sum(Abundance[Aerotolerance_cat == "aerotolerant"], na.rm = TRUE),
    anaerobic = sum(Abundance[Aerotolerance_cat == "obligate anaerobe"], na.rm = TRUE),
    unknown = sum(Abundance[Aerotolerance_cat == "Unknown"], na.rm = TRUE),
    .groups = "drop"
  )

# **Final Check: Ensure Total = 1 per Sample**
aero_summary$t <- round(aero_summary$aerotolerant + aero_summary$anaerobic + aero_summary$unknown, 6)
range(aero_summary$t)  # Should be exactly 1.000000 for all

# **Verify if any samples still have totals ≠ 1**
if (!all(abs(aero_summary$t - 1) < 1e-6)) {
  warning("Some samples do not sum to 1! Investigate further.")
} else {
  print("All samples correctly sum to 1. No issues detected.")
}

# Measure proportion of aerotolerant taxa out of taxa with known aerotolerance
aero_summary$aerotolerant_ofknowntaxa <- aero_summary$aerotolerant/(aero_summary$aerotolerant+aero_summary$anaerobic)

# Add other relevant info
sampledata<- data.frame(sample_data(phy))
sampledata <- sampledata[,c("sample_name", "EB_code", "read_count", "Colony", "Age_clean")]
sampledata$EB <- str_sub(sampledata$EB_code, start=1, end=2)
sampledata <- sampledata[,c("sample_name", "EB", "read_count", "Colony", "Age_clean")]

# Merge
data <- merge(sampledata, aero_summary, by.x = "sample_name", by.y="Sample")

# Make a plot with mean aerotolerance profiles for adult vs chicks
library(tidyverse)

# Reorder samples: adults first, ordered by aerotolerant abundance
data_indiv <- data %>%
  mutate(sample_id = row_number()) %>%
  arrange(Age_clean, desc(aerotolerant)) %>%
  mutate(sample_id = factor(sample_id, levels = sample_id)) %>%
  pivot_longer(cols = c(aerotolerant, anaerobic, unknown),
               names_to = "Oxygen_Tolerance", values_to = "Abundance")

# Set colours
tolerance_colors <- c(
  "aerotolerant" = "#4d8190",
  "anaerobic" = "#f45751",
  "unknown" = "#d3d3d3"
)

# Plot 1: Individuals only
p_indiv <- ggplot(data_indiv, aes(x = sample_id, y = Abundance, fill = Oxygen_Tolerance)) +
  geom_col(width = 1) +
  scale_fill_manual(values = tolerance_colors) +
  labs(x = NULL, y = "Relative abundance", fill = "Oxygen tolerance category") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )
ggsave("Puffin_aerotolerance_samples.png")

# Prep data
data_mean <- data %>%
  group_by(Age_clean) %>%
  summarise(
    aerotolerant = mean(aerotolerant),
    anaerobic = mean(anaerobic),
    unknown = mean(unknown),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(aerotolerant, anaerobic, unknown),
               names_to = "Oxygen_Tolerance", values_to = "Abundance")

# Plot 2: Mean bar per age group
p_mean <- ggplot(data_mean, aes(x = Age_clean, y = Abundance, fill = Oxygen_Tolerance)) +
  geom_col(width = 1) +
  scale_fill_manual(values = tolerance_colors) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_void() +
  theme(legend.position = "none")
ggsave("Puffin_aerotolerance_mean.png")

# Save a file with sample name and relative abundances
mini <- data[,c("sample_name", "aerotolerant", "anaerobic", "unknown")]
write.csv(mini, "sample_aerotolerance.csv")

# Plot

# Convert data to long format while keeping metadata
aero_long <- data %>%
  pivot_longer(cols = c(aerotolerant, anaerobic, unknown), 
               names_to = "Aerotolerance", 
               values_to = "Relative_abundance") %>%
  mutate(Age_Colony = paste(Age_clean, Colony, sep = ", "))  # Format x-axis labels correctly

# Standardize Aerotolerance category names for better readability
aero_long$Aerotolerance[aero_long$Aerotolerance == "aerotolerant"] <- "Aerotolerant"
aero_long$Aerotolerance[aero_long$Aerotolerance == "anaerobic"] <- "Anaerobic"
aero_long$Aerotolerance[aero_long$Aerotolerance == "unknown"] <- "Unknown"

# Standardize Colony names for better readability
aero_long$Colony[aero_long$Colony == "Isle_of_May"] <- "Isle of May"
aero_long$Colony[aero_long$Colony == "Sule_Skerry"] <- "Sule Skerry"
aero_long$Colony[aero_long$Colony == "Fair_Isle"] <- "Fair Isle"
aero_long$Colony[aero_long$Colony == "Orkney"] <- "Pentland Skerries"

# **Take the mean relative abundance for each Age_Colony category**
aero_avg <- aero_long %>%
  group_by(Age_Colony, Aerotolerance) %>%
  summarise(Relative_abundance = mean(Relative_abundance, na.rm = TRUE), .groups = "drop")

# Create the stacked bar plot with adjusted aesthetics
ggplot(aero_avg, aes(x = Age_Colony, y = Relative_abundance, fill = Aerotolerance)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Aerotolerant" = "#4d8190",  # Dark teal
                               "Anaerobic" = "#f45751",    # Red
                               "Unknown" = "grey62")) +    # Gray
  theme_minimal() +
  guides(fill = guide_legend(ncol = 1)) +  # Ensure legend stays in one column
  labs(x = "", y = "Relative abundance (%)",
       title = "", fill = "Aerotolerance") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        text=element_text(size=15),
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, vjust = .5, hjust=.5))  # Rotate x-axis labels

# Does age category predict the proportion of aerotolerant taxa out of taxa with known aerotolerance?

# scale read count
dataset <- data
dataset$read_count_scaled <- scale(dataset$read_count)

model <- brm(
  aerotolerant_ofknowntaxa ~ Age_clean + read_count_scaled + EB + unknown + (1|Colony), 
  data = dataset,
  family = zero_one_inflated_beta(),  
  warmup = 2000, iter = 10000, 
  chains = 1, cores = 1,  
  control = list(adapt_delta = 0.98, max_treedepth = 12),
  init = 0
)

pp_check(model)
ggsave("aerotolerance_ppcheck.png")
pp_check(model, type = "stat")
ggsave("aerotolerance_ppcheckStat.png")
