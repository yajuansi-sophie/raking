rm(list = ls())

source("utils.R")

library(dplyr)
library(tidyr)
library(data.table)
library(rstan)
library(survey)
library(ggplot2)
library(knitr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

acs_data <- readRDS("Data/acs_simu.rds")

# ============= SETTING ================ # 
SEED = 1234
set.seed(SEED)
VARIABLES <- c("age_dc", "race_dc","educat", "sex", "opmres_x")
DEPENDENT <- "y" ## Used in simlation outcome name

acs_data = acs_data[, VARIABLES]

CALIBRATE_MARGIN = list( ~age_dc, ~race_dc, ~educat, ~sex, ~opmres_x)
SUBGROUP = c(CALIBRATE_MARGIN, list(~age_dc+opmres_x))
SREG = ~ age_dc + race_dc + educat + sex + opmres_x
YREG = ~ age_dc + race_dc + educat + sex + opmres_x

ITER = 1000
CHAINS = 4
STANFILE <- "stan/bayes_raking.stan"

SIMULATION_TIME = 300

# ====================================== # 
cat("The following code replicate simulation in section 3.1.", 
    "Simulation time: ", SIMULATION_TIME,
    "It may take a while. Please change SIMULATION_TIME in this script accordingly", sep = "\n")

invisible(readline(prompt="Press [enter] to continue. "))

cat("Use calibrate margins:", 
    VARIABLES, sep = "\n")

# ============== Outcome and inclusion  ================== #

YBETA = c(`(Intercept)` = 0.85, age_dc2 = 0.41, age_dc3 = 0.48, age_dc4 = 0, 
          age_dc5 = -0.63, race_dc2 = 1, race_dc3 = 0, race_dc4 = 1.14, 
          race_dc5 = 1.28, educat2 = 0, educat3 = 0, educat4 = -0.81, sex2 = 0.31, 
          opmres_x2 = -0.61, opmres_x3 = 0, opmres_x4 = -0.78, opmres_x5 = -1.38)

PBETA = c(`(Intercept)` = -4.31, age_dc2 = 0.26, age_dc3 = 0.46, age_dc4 = 0.57, 
          age_dc5 = 0.51, race_dc2 = 0.43, race_dc3 = -0.84, race_dc4 = 1.13, 
          race_dc5 = 0.68, educat2 = 0.47, educat3 = 0.64, educat4 = 1.19, 
          sex2 = 0.32, opmres_x2 = 0, opmres_x3 = -0.26, opmres_x4 = -0.46, 
          opmres_x5 = -0.48)

cat("Outcome model use logistic model with parameters:",
    "P(Y = 1 | X ) = logit(X'ybeta)",
    "ybeta: ", sep = "\n")
cat(knitr::kable(YBETA), sep="\n")

cat("Inclusion model use logistic model with parameters:",
    "P(I = 1 | X) = logit(X'pbeta)",
    "pbeta: ", sep = "\n")
cat(knitr::kable(PBETA), sep="\n")

pop = as.data.frame(sapply(acs_data, factor))

pop_contingency = pop %>%
  xtabs(~., data = .) %>%
  as.data.frame() %>%
  filter(Freq != 0) %>%
  mutate(id = 1:nrow(.)) %>%
  as.data.table()

template <- pop_contingency
setkeyv(pop_contingency, VARIABLES)
template$Freq <- NULL
setkeyv(template, VARIABLES)
pop_contingency = pop_contingency %>% arrange(id)

yloading = model.matrix(YREG, pop)
sloading = model.matrix(SREG, pop)
if (any (colnames(sloading) == names(PBETA))) {
   PBETA = PBETA[colnames(sloading)]
}
yprob = log_inv(yloading %*% YBETA)
sprob = log_inv(sloading %*% PBETA)

pop[[DEPENDENT]] = sapply(yprob, FUN = function(x) {rbinom(1, 1, x)})

cat("Finish preparing population data.", 
    "Simulation start.", sep = "\n")

# ============= Extract marginal information === #

cat("Extract marginals information.\n")

population.margins = lapply(CALIBRATE_MARGIN, FUN = function(x) {xtabs(x, pop)})
Nmargin = margin_vector(population.margins = population.margins)

print_marign(CALIBRATE_MARGIN, population.margins)
J = NROW(pop_contingency) # 986
D = NROW(Nmargin) # 21

cat("Contingency table has: ", J, " non-empty cells", "\n")

# ============= STAN Prepare =================== #

cat("Start complie Stan:", STANFILE, 
    "It may take a while", sep='\n')

bayes_raking = stan_model(model_name = 'bayes_raking', 
                          file = STANFILE)

# ============== Holding quantities ============ #

L_pop = loading_matrix(pop_contingency, SUBGROUP)
popy = template[pop] %>% group_by(id) %>% summarise(ysum = sum(!!sym(DEPENDENT)), Freq = n())
ymean_true = mean(pop[[DEPENDENT]])
ymarginal_true = (L_pop %*% popy$ysum) / (L_pop %*% popy$Freq)

ymean_bayes_summary = c()
ymean_raking_summary = c()
ymarginal_bayes_summary = c()
ymarginal_raking_summary = c()

bayes_raking_list = raking_list = list()
bayes_time = raking_time = 0
# ============================================== #
for (iter in seq_len(SIMULATION_TIME)) {
  set.seed(SEED + iter)
  cat("Seed for current simulation: ", SEED + iter, "\n")
  sample_data = sampling_from_pop(pop, sprob)
  
  sample_summary = bayes_sample_contingency(sample_data, VARIABLES,
                                            CALIBRATE_MARGIN, SUBGROUP, 
                                            SREG, template)
  
  non_empty_J = sum(sample_summary$ncell != 0)
  tmp_sam = template[sample_data] %>% 
    group_by(id) %>% 
    summarise(ysum = sum(!!sym(DEPENDENT)), Freq = n())
  y_success = tmp_sam$ysum
  y_total = tmp_sam$Freq
  
  ydesign_J = model.matrix(YREG, sample_summary$sample_contingency)[, YBETA != 0]
  ydesign_non_empty = ydesign_J[tmp_sam$id, ]
  py = NCOL(ydesign_J)
  
  sample_summary$pdesign_J = sample_summary$pdesign_J[, PBETA != 0]
  sample_summary$ps = NCOL(sample_summary$pdesign_J)
  
  data_list = c(sample_summary, list(J = J, Nmargin = as.vector(Nmargin),
                   non_empty_J = non_empty_J, y_success = as.vector(y_success), y_total = as.vector(y_total), 
                   py = py, ydesign_non_empty = ydesign_non_empty, ydesign_J = ydesign_J))
  
  cat("Start sampling with iteration: ", ITER, ", number of chains: ", CHAINS, 
      ", seed: ", SEED, "\n")
  
  ptm <- proc.time()
  braking_fit = sampling(bayes_raking, data = data_list, chains = CHAINS,
                         iter = ITER, seed = SEED, open_progress = FALSE,
                         show_messages = FALSE)
  bayes_time  = bayes_time + (proc.time() - ptm)
  
  ymean_bayes_raking = summary_overall(summary_wrap(braking_fit, "ymean"), ymean_true)
  ymarginal_bayes_raking = summary_marginal(summary_wrap(braking_fit, "ymarginals"), ymarginal_true)
  rownames(ymarginal_bayes_raking) = rownames(data_list$L_quant)
  
  stopifnot(all(rownames(ymarginal_true) == rownames(ymarginal_bayes_raking)))
  cat("Start using survey::rake\n")
  
  ptm = proc.time()
  design = svydesign(id = ~0, probs = NULL, data = sample_data)
  rclus = rake(design, sample.margins = CALIBRATE_MARGIN, population.margins = population.margins)
  raking_time = raking_time + (proc.time() - ptm)
  
  ymean_raking = summary_overall(raking_mean(rclus, DEPENDENT), ymean_true)
  ymarginal_raking = summary_marginal(
    raking_marginal(rclus, DEPENDENT, SUBGROUP), ymarginal_true)
  
  stopifnot(all(rownames(ymarginal_true) == rownames(ymarginal_raking)))
  
  
  bayes_raking_list[[iter]] = list(ymarginal = ymarginal_bayes_raking, 
                                 ymean = ymean_bayes_raking)
  raking_list[[iter]] = list(ymarginal = ymarginal_raking, 
                           ymean = ymean_raking)
  
  if (is.null(ymarginal_bayes_summary)) {
    ymean_bayes_summary = ymean_bayes_raking
    ymean_raking_summary = ymean_raking
    
    ymarginal_bayes_summary = ymarginal_bayes_raking
    ymarginal_raking_summary = ymarginal_raking
  } else {
    ymean_bayes_summary = ymean_bayes_summary + ymean_bayes_raking
    ymean_raking_summary = ymean_raking_summary + ymean_raking
    
    ymarginal_bayes_summary = ymarginal_bayes_summary + ymarginal_bayes_raking
    ymarginal_raking_summary = ymarginal_raking_summary + ymarginal_raking
  }
  rm(ymean_raking, ymean_bayes_raking, 
     ymarginal_raking, ymarginal_bayes_raking, data_list, sample_data, sample_summary,
     braking_fit, design, rclus)
  cat("Finish ", iter, "th simulation.\n")
}

ymean_bayes_summary = ymean_bayes_summary / iter
ymean_raking_summary = ymean_raking_summary / iter
ymar_bayes_plot = ymarginal_bayes_summary = ymarginal_bayes_summary / iter
ymar_raking_plot = ymarginal_raking_summary = ymarginal_raking_summary / iter

cat("Bayes raking average time: ", bayes_time[3] / iter, "\n")
cat("Raking average time: ", raking_time[3] / iter, "\n")

cat("True outcome mean: ", ymean_true, '\n')
cat("Bayes raking estimate: \n")
print_summary(ymean_bayes_summary)
cat("Raking estimate: \n")
print_summary(ymean_raking_summary)

ymar_bayes_plot = plot_prepare(ymar_bayes_plot)
ymar_raking_plot = plot_prepare(ymar_raking_plot)
ymar_bayes_plot['Method'] = "Bayes-raking"
ymar_raking_plot['Method'] = "Raking"


p1 =  rbind(ymar_bayes_plot, ymar_raking_plot) %>% 
  select(-Margin) %>%
  gather(Quantities, Value, -Method) %>%
  filter(Quantities %in% c("Coverage", "Abs.Bias", "RMSE", "StandardErr")) %>%
  ggplot(aes(x = Method, y = Value)) +
  geom_violin() + facet_wrap(~Quantities, scales = "free_y")
print(p1)
