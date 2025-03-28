library(brms)
library(rstan)
library(phyloseq)

# Use rstan settings to avoid recompilation
rstan_options(auto_write = TRUE)  
options(mc.cores = parallel::detectCores())  

Sys.setenv(TMPDIR = "/scratch/project_2010022")  # Set temporary directory

library(dplyr)
library(stringr)

# Load data
phy <- readRDS("Puffin_phyloseq.RDS")
data <- data.frame(sample_data(phy))

# Scaling function
range.use <- function(x, min.use, max.use) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (max.use - min.use) + min.use
}

scalecols <- c("Shannon_asy", "Richness_asy", "read_count")

for (i in scalecols) {
  data[[i]] <- range.use(data[[i]], 0.000000001, 0.99999999)
}

# Extract extraction batch info
data$batch <- str_sub(data$EB_code, start = 1, end=2)

### Shannon diversity
data$Shannon_logit <- log(data$Shannon_asy / (1 - data$Shannon_asy))

# beta regression
shannon_beta <- brm(
  Shannon_asy ~ 1 + Age_clean + batch + read_count + (1 | Colony),  
  data = data, 
  family = "beta",
  warmup = 2000, iter = 8000, 
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

# logit transformed gaussian
shannon_logit <- brm(
  Shannon_logit ~ 1 + Age_clean + batch + read_count + (1 | Colony),  
  data = data, 
  family = gaussian(),
  warmup = 2000, iter = 8000, 
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

# student t (no nu prior)
shannon_t <- brm(
  Shannon_logit ~ 1 + Age_clean + batch + read_count + (1 | Colony),  
  data = data, 
  family = student,
  warmup = 2000, iter = 8000, 
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

# student t (with nu prior)
shannon_t_nu <- brm(
  Shannon_logit ~ 1 + Age_clean + batch + read_count + (1 | Colony),  
  data = data, 
  family = student,
  prior = c(
    prior(gamma(7, 1), class = "nu")  
  ),
  warmup = 2000, iter = 8000, 
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

# Compare models to find best fit
loo_beta <- loo(shannon_beta)
loo_logit <- loo(shannon_logit)
loo_t <- loo(shannon_t)
loo_t_nu <- loo(shannon_t_nu)
loo_compare(loo_beta, loo_logit, loo_t, loo_t_nu)
#               elpd_diff se_diff
# shannon_beta     0.0       0.0 
# shannon_t     -180.2      18.3 
# shannon_t_nu  -181.2      18.8 
# shannon_logit -224.0      23.2 
pp_check(shannon_beta) # still really bad fit. saving this plot and then trying a different approach
ggsave("Shannon_beta_pp_check.png")

# Going to try a few more models and also remove batch effect as it doesn't seem to have significant effects on diversity and is currently complicating models.

# Logit transformation for Gaussian and Student-T models
data$Shannon_logit <- log(data$Shannon_asy / (1 - data$Shannon_asy))

# 1.. Beta Regression Model
shannon_beta <- brm(
  Shannon_asy ~ 1 + Age_clean + read_count + (1 | Colony),  
  data = data, 
  family = "beta",
  warmup = 2000, iter = 8000, 
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

# 2. Zero-Inflated Beta Regression (Handles many near-zero values)
shannon_zib <- brm(
  bf(Shannon_asy ~ 1 + Age_clean + read_count + (1 | Colony), 
     phi ~ 1 + Age_clean, 
     zi ~ 1 + Age_clean),  # Zero-inflation term
  data = data,
  family = zero_inflated_beta(),
  warmup = 2000, iter = 8000,
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

# 3. Hurdle Lognormal Model (For zero-inflated right-skewed data)
shannon_hurdle <- brm(
  Shannon_asy ~ 1 + Age_clean + read_count + (1 | Colony),  
  data = data, 
  family = hurdle_lognormal(),
  warmup = 2000, iter = 8000, 
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

# 4. Student-T Model (With nu prior to control outliers)
shannon_t_nu <- brm(
  Shannon_logit ~ 1 + Age_clean + read_count + (1 | Colony),  
  data = data, 
  family = student,
  prior = c(prior(gamma(7, 1), class = "nu")),  # Prior to limit extreme values
  warmup = 2000, iter = 8000, 
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

# Compare Model Fits Using LOO Cross-Validation
loo_beta <- loo(shannon_beta)
loo_zib <- loo(shannon_zib)
loo_hurdle <- loo(shannon_hurdle)
loo_t_nu <- loo(shannon_t_nu)

# Display LOO comparison
loo_compare(loo_beta, loo_zib, loo_hurdle, loo_t_nu)
#                 elpd_diff se_diff
# shannon_zib       0.0       0.0 
# shannon_beta    -46.1      17.3 
# shannon_hurdle  -51.9      25.6 
# shannon_t_nu   -228.9      27.7

# Run Posterior Predictive Checks for Model Diagnostics
pp_check(shannon_beta)
pp_check(shannon_zib) 
pp_check(shannon_hurdle)
pp_check(shannon_t_nu)
# zero inflated beta gives best fit

summary(shannon_zib)
# standard deviation (sd(Intercept) = 0.50) suggests there’s colony-level variability.
# Going to check if allowing phi (precision) to vary by colony helps. This would let different colonies have different precision parameters,
# which might better account for colony-specific dispersion in Shannon diversity.
shannon_zib_phi_colony <- brm(
  bf(Shannon_asy ~ 1 + Age_clean + read_count + (1 | Colony), 
     phi ~ 1 + Age_clean + (1 | Colony),  # Allow precision to vary by colony
     zi ~ 1 + Age_clean),  # Zero-inflation component
  data = data,
  family = zero_inflated_beta(),
  warmup = 2000, iter = 8000,
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.995, max_treedepth = 15),  # Slightly higher adapt_delta
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

loo1 <- loo(shannon_zib)  
loo2 <- loo(shannon_zib_phi_colony)  

loo_compare(loo1, loo2)
#                        elpd_diff se_diff
# shannon_zib             0.0       0.0   
# shannon_zib_phi_colony -0.2       1.2 
# No meaningful difference -> stick to simpler shannon_zib model
pp_check(shannon_zib)
ggsave("Shannon_beta_zib_pp_check.png")
saveRDS(shannon_zib, "Shannon_beta_zib_model.RDS")
conditional_effects(shannon_zib, effects = "Age_clean", re_formula = NULL)
ggsave("Shannon_beta_zib_conditionaleffects.png")

### ASV richness

# Log transform for Gaussian models
data$Richness_log <- log(data$Richness_asy + 1)  # Avoid log(0) errors

# 1. Gaussian Regression on Raw Richness
richness_gaussian <- brm(
  Richness_asy ~ 1 + Age_clean + read_count + (1 | Colony),  
  data = data,
  family = gaussian(),
  warmup = 2000, iter = 8000,
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

# 2. Gaussian Regression on Log-Transformed Richness
richness_log_gaussian <- brm(
  Richness_log ~ 1 + Age_clean + read_count + (1 | Colony),  
  data = data,
  family = gaussian(),
  warmup = 2000, iter = 8000,
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

# 3. Hurdle Lognormal Model (For zero-inflated right-skewed data)
richness_hurdle <- brm(
  Richness_asy ~ 1 + Age_clean + read_count + (1 | Colony),  
  data = data, 
  family = hurdle_lognormal(),
  warmup = 2000, iter = 8000,
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "rstan",
  save_pars = save_pars(all = TRUE)
)

# 
loo_gaussian <- loo(richness_gaussian)
loo_log_gaussian <- loo(richness_log_gaussian)
loo_hurdle <- loo(richness_hurdle)

# Display LOO comparison
loo_compare(loo_gaussian, loo_log_gaussian, loo_hurdle)
#                       elpd_diff se_diff
# richness_log_gaussian   0.0       0.0  
# richness_gaussian     -19.7       1.0  
# richness_hurdle       -67.7      14.3 

pp_check(richness_gaussian) + ggtitle("Gaussian (Raw)")
pp_check(richness_log_gaussian) + ggtitle("Gaussian (Log)")
pp_check(richness_hurdle) + ggtitle("Hurdle Lognormal")

summary(richness_log_gaussian)
saveRDS(richness_log_gaussian, "richness_log_gaussian_model.RDS")
pp_check(richness_log_gaussian) + ggtitle("Gaussian (Log)")
ggsave("richness_log_gaussian_pp_check.png")
conditional_effects(richness_log_gaussian, effects = "Age_clean", re_formula = NULL)
ggsave("richness_log_gaussian_conditionaleffects.png")
