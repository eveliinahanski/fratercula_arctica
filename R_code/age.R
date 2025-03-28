library(brms)
library(rstan)

# Use rstan settings to avoid recompilation
rstan_options(auto_write = TRUE)  
options(mc.cores = parallel::detectCores())  

Sys.setenv(TMPDIR = "/scratch/project_2010022")  # Set temporary directory

library(dplyr)
library(stringr)

# Load data
dyadic <- read.csv("puffin_dyadic_data.csv")

# Scaling function
range.use <- function(x, min.use, max.use) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (max.use - min.use) + min.use
}

scalecols <- c("At_sea_distance", "Absolute_distance", "Jaccard_dissimilarity", "read_count_diff")

for (i in scalecols) {
  dyadic[[i]] <- range.use(dyadic[[i]], 0, 1)
}

# Convert dissimilarities to similarities
dyadic$Jaccard_similarity <- 1 - dyadic$Jaccard_dissimilarity

# Transform response variables for beta distribution
samplesize <- nrow(dyadic)
dyadic$Jaccard <- (dyadic$Jaccard_similarity * (samplesize - 1) + 0.5) / samplesize

# Select same age category pairs
dyadic <- dyadic[which(dyadic$age_pair=="Chick-Chick" | dyadic$age_pair=="Adult-Adult"),]

dyadic <- dyadic %>%
  rowwise() %>%
  mutate(pair_id = paste(sort(c(Sample1, Sample2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE)

# brm model using rstan backend
model <- brm(
  Jaccard ~ 1 + age_pair + colony_similarity + EB_pair + read_count_diff + (1 | mm(Sample1, Sample2)),  
  data = dyadic, 
  family = "beta",
  warmup = 20000, iter = 40000, 
  chains = 4, cores = 4, threads = threading(9), 
  control = list(adapt_delta = 0.95, max_treedepth = 13), 
  backend = "rstan",  # Use rstan instead of cmdstanr
  save_pars = save_pars(group = TRUE),
  init = 0
)

saveRDS(model, "model_age_similarity2.rds")
