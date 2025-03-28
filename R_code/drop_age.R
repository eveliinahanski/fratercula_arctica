library(phyloseq)
library(brms)
library(rstan)
library(dplyr)
library(stringr)
library(parallel)
library(doParallel)
library(data.table)

# Load data
micdata <- readRDS("micdata.RDS")
micdata <- filter_taxa(micdata, function(x) max(x) > 0, TRUE)

# Convert spaces in genus names with underscores
taxa_table <- data.frame(tax_table(micdata))
taxa_table$Genus <- gsub(' ', '_', taxa_table$Genus)

# Since taxa will be dropped at genus level, make sure each genus is present only once (even if within genus there would be species variation)
taxa_table$Species <- NA
taxa_table <- taxa_table[!duplicated(taxa_table),] 

# Update tax_table
taxa_table <- as.matrix(taxa_table)
tax_table(micdata) <- tax_table(taxa_table)

# Load dyadic data
data.dyad_REAL <- readRDS("puffin_dyadic_data.RDS")
data.dyad_REAL <- data.dyad_REAL[which(data.dyad_REAL$age_pair=="Adult-Adult" | data.dyad_REAL$age_pair=="Chick-Chick"),]

# Data wrangling step: Give all unclassified genera a name based on their family or order
tax <- as.data.frame(tax_table(micdata))

# Make a list of the bacterial genera present in the data and how many unique taxa (ASVs) each genus has
tax_G <- as.data.frame(table(tax[, 6]))
taxNone <- tax_G
taxNone$Var1 <- "none"
taxNone$Freq <- 0
taxNone <- taxNone[!duplicated(taxNone), ]
tax_GG <- rbind(tax_G, taxNone)

# Extract genera (select only a couple if optimising)
genuses <- tax_GG$Var1

# Make a key for the order of sample names and their associated individual IDs.
key <- data.frame(ID = sample_data(micdata)$sample_name, sample_name = sample_data(micdata)$sample_name)

# Optionally run in serial mode for debugging (set to TRUE to debug, set FALSE to run in parallel mode)
serial_mode <- TRUE

if (serial_mode) {
  # Serial mode
  for (i in 1:length(genuses)) {
    cat("Processing genus:", genuses[i], "\n")
    
    # Use tryCatch for error handling
    result <- tryCatch({
      gen.i <- genuses[i]
      cat("Genus:", gen.i, "- Starting...\n")
      
      taxa_tokeep <- rownames(tax[which(tax$Genus != gen.i), ])
      mic.i <- prune_taxa(taxa_tokeep, micdata)
      
      JACM.i <- as.matrix(phyloseq::distance(mic.i, method = "jaccard", binary = TRUE))
      BRAY.i <- as.matrix(phyloseq::distance(mic.i, method = "bray"))
      
      bray <- c(as.dist(BRAY.i))
      jac <- c(as.dist(JACM.i))
      
      data.dyad.i <- data.frame(Jaccard = jac, BrayCurtis = bray)
      
      list <- expand.grid(key$sample_name, key$sample_name)
      list <- list[which(list$Var1 != list$Var2), ]
      list$key <- apply(list, 1, function(x) paste(sort(x), collapse = ''))
      list <- subset(list, !duplicated(list$key))
      
      data.dyad.i$Sample_A <- list$Var2
      data.dyad.i$Sample_B <- list$Var1
      
      keyA <- key[, c("sample_name", "sample_name")]
      colnames(keyA) <- c("id1", "Sample_A")
      keyB <- key[, c("sample_name", "sample_name")]
      colnames(keyB) <- c("id2", "Sample_B")
      
      keyA <- keyA[match(data.dyad.i$Sample_A, keyA$Sample_A), ]
      keyB <- keyB[match(data.dyad.i$Sample_B, keyB$Sample_B), ]
      
      data.dyad.i$id1 <- keyA$id1
      data.dyad.i$id2 <- keyB$id2
      
      samplekey <- data.dyad_REAL[, c("sample_pair", "age_pair", "read_count_diff", "EB_pair", "colony_similarity")]
      data.dyad.i$sample_pair <- paste(data.dyad.i$id1, data.dyad.i$id2, sep = "-")
      data.dyad.i <- merge(data.dyad.i, samplekey, by = "sample_pair", all.x = TRUE, all.y = FALSE)
      data.dyad.i <- na.omit(data.dyad.i)
      data.dyad.i <- data.dyad.i[which(data.dyad.i$id1 != data.dyad.i$id2), ]
      
      data.dyad <- data.dyad.i
      
      data.dyad$id1 <- as.factor(data.dyad$id1)
      data.dyad$id2 <- as.factor(data.dyad$id2)
      
      scalecols <- c("Jaccard", "read_count_diff")
      range.use <- function(x, min.use, max.use) {
        (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (max.use - min.use) + min.use
      }
      
      for (col in scalecols) {
        data.dyad[, col] <- range.use(data.dyad[, col], 0, 1)
      }
      
      data.dyad$Jaccard <- 1 - data.dyad$Jaccard
      
      samplesize <- nrow(data.dyad)
      data.dyad$Jaccard <- (data.dyad$Jaccard * (samplesize - 1) + 0.5) / samplesize
      
      cat("Genus:", gen.i, "- Running model...\n")
      
      # Run model
      dropmodel <- brm(Jaccard ~ 1 + age_pair + colony_similarity + EB_pair + read_count_diff + (1 | mm(id1, id2)), data = data.dyad, family = "Beta",
                       control = list(adapt_delta = 0.9, max_treedepth = 12), warmup = 800, iter = 1600, cores = 1, chains = 1,
                       save_pars = save_pars(group = FALSE), init = 0)
      
      adapt_delta_save <- 0.9
      max_treedepth_save <- 12
      warmup_save <- 800
      iter_save <- 1600
      
      sampler_params <- rstan::get_sampler_params(dropmodel$fit, inc_warmup = FALSE)
      divergent_transitions <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
      total_divergent_transitions <- sum(divergent_transitions)
      rhat <- summary(dropmodel)$fixed[2, "Rhat"]
      bulk_ess <- summary(dropmodel)$fixed[2, "Bulk_ESS"]
      
      if (total_divergent_transitions > 10 || any(rhat > 1.05) || any(bulk_ess < 160)) {
        cat("Genus:", gen.i, "- Re-running model with updated parameters (adapt_delta = 0.95)...\n")
        
        dropmodel <- brm(Jaccard ~ 1 + age_pair + colony_similarity + EB_pair + read_count_diff + (1 | mm(id1, id2)), data = data.dyad, family = "Beta",
            control = list(adapt_delta = 0.95, max_treedepth = 13), warmup = 1500, iter = 3000, cores = 1, chains = 1,
            save_pars = save_pars(group = FALSE), init = 0)
        
        adapt_delta_save <- 0.95
        max_treedepth_save <- 13
        warmup_save <- 1500
        iter_save <- 3000
        
        sampler_params <- rstan::get_sampler_params(dropmodel$fit, inc_warmup = FALSE)
        divergent_transitions <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
        total_divergent_transitions <- sum(divergent_transitions)
        rhat <- summary(dropmodel)$fixed[2, "Rhat"]
        bulk_ess <- summary(dropmodel)$fixed[2, "Bulk_ESS"]
        
        if (total_divergent_transitions > 10 || any(rhat > 1.05) || any(bulk_ess < 300)) {
          cat("Genus:", gen.i, "- Re-running model with updated parameters (adapt_delta = 0.98)...\n")
          
          dropmodel <- brm(Jaccard ~ 1 + age_pair + colony_similarity + EB_pair + read_count_diff + (1 | mm(id1, id2)), data = data.dyad, family = "Beta",
              control = list(adapt_delta = 0.98, max_treedepth = 13), warmup = 3000, iter = 6500, cores = 1, chains = 1,
              save_pars = save_pars(group = FALSE), init = 0)
          
          adapt_delta_save <- 0.98
          max_treedepth_save <- 13
          warmup_save <- 3000
          iter_save <- 6500
        }}
      
      # Use the final model results
      
      sampler_params <- rstan::get_sampler_params(dropmodel$fit, inc_warmup = FALSE)
      divergent_transitions <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
      total_divergent_transitions <- sum(divergent_transitions)
      
      # Extract number of ASVs dropped when genus was excluded
      ASVs_dropped.i <- nrow(tax_table(micdata)) - nrow(tax_table(mic.i))
      
      # Extract model estimates and construct output dataframe
      fixed_summary <- summary(dropmodel)$fixed
      
      resdf.i <- data.frame(
        Genus_dropped = gen.i, 
        ASVs_dropped = ASVs_dropped.i,
        divergent_transitions = total_divergent_transitions,
        adapt_delta = adapt_delta_save,
        max_treedepth = max_treedepth_save,
        warmup = warmup_save,
        iterations = iter_save,
        Age_Estimate = fixed_summary[2, "Estimate"],
        Age_EstError = fixed_summary[2, "Est.Error"],
        Age_lCI = fixed_summary[2, "l-95% CI"],
        Age_uCI = fixed_summary[2, "u-95% CI"],
        Age_Rhat = fixed_summary[2, "Rhat"],
        Age_Bulk_ESS = fixed_summary[2, "Bulk_ESS"]
      )
      
    
      # Save the results for this genus to an RDS file
      saveRDS(resdf.i, paste0("dropResults_genus_", gen.i, ".rds"))
      
      cat("Genus:", gen.i, "- Completed!\n")
      
      return(resdf.i)
      
    }, error = function(e) {
      cat("Error encountered for genus:", genuses[i], "\n")
      cat("Error message:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(result)) {
      if (!exists("final_results")) {
        final_results <- result
      } else {
        final_results <- rbind(final_results, result)
      }
    }
  }
}

# Gather output files
directory <- "/scratch/project_2010022/drop"
rds_files <- list.files(path = directory, pattern = "\\.rds$", full.names = TRUE)
# Check if there are any .rds files in the directory
if (length(rds_files) == 0) {
  stop("No .rds files found in the specified directory.")
}

# Combine all .rds files using rbind
combined_data <- lapply(rds_files, readRDS)
combined_data <- rbindlist(combined_data, fill = TRUE)
combined_data <- combined_data[!duplicated(combined_data),]

saveRDS(combined_data, "puffin_Age_drop_model_output.rds")

# Calculate importance scores
dropResults <- combined_data
none <- dropResults[which(dropResults$Genus=="none"),] 
none <- abs(none$Age_uCI-none$Age_lCI) 
none <- none[1]
dropResults$Age_CIbr_baseline <- none 

dropResults$Age_CIbr <- abs(dropResults$Age_uCI-dropResults$Age_lCI)
dropResults$Age_CIbr_increase <- dropResults$Age_CIbr-dropResults$Age_CIbr_baseline

#Importance value is this increase in uncertainty divided by the square root of how many ASVs were dropped, and multiplied by 100 to increase the scale.
dropResults$IMPORTANCE_Age<-(dropResults$Age_CIbr_increase/dropResults$Age_CIbr_baseline)/sqrt(dropResults$ASVs_dropped)*100

# Plot distribution and save as PNG
hist(dropResults$IMPORTANCE_Age, main="Importance score distribution", xlab="Importance score", ylab="Frequency", col="blue")
png("importance_distribution.png")
dev.off()

# Scale importance scores to 0-1

dropResults$extra <- dropResults$IMPORTANCE_Age

dropResults <- data.frame(dropResults)

scalecols <- c("IMPORTANCE_Age", "extra")

range.use <- function(x, min.use, max.use) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (max.use - min.use) + min.use
}

for(i in 1:ncol(dropResults[,which(colnames(dropResults)%in%scalecols)])){
  dropResults[,which(colnames(dropResults)%in%scalecols)][,i] <- range.use(dropResults[,which(colnames(dropResults)%in%scalecols)][,i],0,1)
}

# Add taxonomy 
micdata <- readRDS("micdata.RDS")
taxa <- data.frame(tax_table(micdata))
taxa <- taxa[,c(2:6)]
taxa$Genus <- gsub(' ', '_', taxa$Genus)
taxa <- taxa[!duplicated(taxa),]

dropResults <- merge(dropResults, taxa, by.x="Genus_dropped", by.y="Genus", all.x=T, all.y=F)
dropResults$RHat <- "<1.05"
dropResults$RHat[dropResults$Age_Rhat>=1.05] <- ">=1.05"
dropResults$Bulk <- ">10%"
dropResults$Bulk[(dropResults$Age_Bulk_ESS/dropResults$iteratations)<=0.1] <- "<=10%"
dropResults$Fail <- FALSE
dropResults$Fail[dropResults$RHat==">=1.05" | dropResults$Bulk=="<=10%"] <- TRUE

saveRDS(dropResults, "puffin_Age_genus_importance_scaled.rds")