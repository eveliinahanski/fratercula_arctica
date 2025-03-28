library(brms)
library(rstan)
library(dplyr)

dyadic <- read.csv("puffin_dyadic_data.csv")

cmdstanr::set_cmdstan_path("/projappl/project_2008131/cmdstan-2.34.1")

#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }

scalecols <- c("At_sea_distance", "Jaccard_dissimilarity", "Aitchison_dissimilarity", "read_count_diff")

for(i in 1:ncol(dyadic[,which(colnames(dyadic)%in%scalecols)])){
  dyadic[,which(colnames(dyadic)%in%scalecols)][,i] <- range.use(dyadic[,which(colnames(dyadic)%in%scalecols)][,i],0,1)
}

# transform dissimilarities to similarities
dyadic$Jaccard <- 1-dyadic$Jaccard_dissimilarity
dyadic$Aitchison <- 1-dyadic$Aitchison_dissimilarity

samplesize <- nrow(dyadic)
dyadic$Jaccard <- (dyadic$Jaccard * (samplesize - 1) + 0.5) / samplesize
dyadic$Aitchison <- (dyadic$Aitchison * (samplesize - 1) + 0.5) / samplesize

# convert at-sea distance to at-sea proximity
dyadic$At_sea_distance <- 1-dyadic$At_sea_distance

# simplify age pairs
dyadic$age_pair[dyadic$age_pair=="Chick-Adult"]<- "Adult-Chick"

dyadic <- dyadic %>%
  rowwise() %>%
  mutate(pair_id = paste(sort(c(Sample1, Sample2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE)

#dyadic <- dyadic[which(dyadic$colony_similarity=="DIFF_COLONY"),]

# Chicks:
df <- dyadic[which(dyadic$age_pair=="Chick-Chick"),]

model <- brm(Jaccard ~ 1 + At_sea_distance + colony_similarity + read_count_diff + EB_pair +
               (1|mm(Sample1, Sample2)),  
             data = df, 
             family= "beta",
             warmup = 20000, iter = 40000, 
             chains = 4, cores = 4, threads = threading(9), 
             control = list(adapt_delta = 0.97, max_treedepth = 13), 
             backend = "cmdstanr",
             stan_model_args = list(stanc_options = list("O1")),
             save_pars = save_pars(group = TRUE),
             init = 0)
saveRDS(model, "at_sea_prox_chicks_Jaccard.rds")

# Adults:
df <- dyadic[which(dyadic$age_pair=="Adult-Adult"),]

model <- brm(Jaccard ~ 1 + At_sea_distance + colony_similarity + read_count_diff + EB_pair +
               (1|mm(Sample1, Sample2)),  
            data = df, 
            family= "beta",
           warmup = 20000, iter = 40000,
            cores=2, chains=2, 
           control = list(adapt_delta = 0.97), max_treedepth = 13,  
            backend = "cmdstanr",
           stan_model_args = list(stanc_options = list("O1")),
           save_pars = save_pars(group = TRUE),
            init = 0)
saveRDS(model, "at_sea_prox_adults_Jaccard.rds")
