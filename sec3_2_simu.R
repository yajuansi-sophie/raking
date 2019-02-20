rm(list = ls())
SEED = 1234
set.seed(SEED)

source("utils.R")

library(dplyr)
library(tidyr)
library(data.table) # MUST INCLUDE
library(rstan)
library(survey)
library(ggplot2)
library(knitr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

acs_data <- readRDS("Data/acs_simu.rds")

VARIABLES <- c("age_dc", "opmres_x", "childx_ca")
acs_data = acs_data[, VARIABLES]
DEPENDENT <- "y" ## Used in simlation outcome name

CALIBRATE_MARGIN = list(~age_dc, ~childx_ca, ~opmres_x, ~age_dc+opmres_x, ~childx_ca + opmres_x)
SUBGROUP = c(~age_dc, ~childx_ca, ~opmres_x, ~childx_ca + opmres_x)

SREG = ~ age_dc + opmres_x + childx_ca + age_dc:opmres_x + childx_ca:opmres_x
YREG = SREG

SIMULATION_TIME = 10
ITER = 2000
CHAINS = 4
PROJECTION_ON = FALSE
PROJECTION_STAN = "stan/projection_32.stan"
BASIS_ON = FALSE
BASIS_STAN = "stan/basis_32.stan"

INDEPENDENT = FALSE

BAYESRAKING_ON = TRUE
BAYESRAKING_STAN = "stan/bayes_raking_32.stan"
RAKING_ON = TRUE




# ====================================== # 
cat("The following code replicate simulation in section 3.2.", 
    "Simulation time: ", SIMULATION_TIME,
    "It may take a while. Please change SIMULATION_TIME in this script accordingly", sep = "\n")

invisible(readline(prompt="Press [enter] to continue. "))

cat("Use calibrate margins:", 
    VARIABLES, sep = "\n")

# ============== Outcome and inclusion  ================== #
PBETA = c(`(Intercept)` = -3.45, 
          age_dc2 = 0, age_dc3 = 0.53, age_dc4 = 0.32, age_dc5 = 0.76, 
          opmres_x2 = 0.12, opmres_x3 = 0, opmres_x4 = 0, opmres_x5 = -0.3, 
          childx_ca2 = 0.48, childx_ca3 = 0.17, childx_ca4 = -0.17, 
          `age_dc2:opmres_x2` = 0, `age_dc3:opmres_x2` = 0, 
          `age_dc4:opmres_x2` = 0, `age_dc5:opmres_x2` = -0.87, 
          `age_dc2:opmres_x3` = 0, `age_dc3:opmres_x3` = 0, 
          `age_dc4:opmres_x3` = 0, `age_dc5:opmres_x3` = 0, 
          `age_dc2:opmres_x4` = 0, `age_dc3:opmres_x4` = 0, 
          `age_dc4:opmres_x4` = 0, `age_dc5:opmres_x4` = 0, 
          `age_dc2:opmres_x5` = 0, `age_dc3:opmres_x5` = 0, 
          `age_dc4:opmres_x5` = 0, `age_dc5:opmres_x5` = -0.30, 
          `opmres_x2:childx_ca2` = -0.52, `opmres_x3:childx_ca2` = 0, 
          `opmres_x4:childx_ca2` = 0, `opmres_x5:childx_ca2` = 0, 
          `opmres_x2:childx_ca3` = 0, `opmres_x3:childx_ca3` = 0, 
          `opmres_x4:childx_ca3` = 0, `opmres_x5:childx_ca3` = 0, 
          `opmres_x2:childx_ca4` = 0, `opmres_x3:childx_ca4` = 0, 
          `opmres_x4:childx_ca4` = 0, `opmres_x5:childx_ca4` = 0.46)

cat("Inclusion model use logistic model with parameters:",
    "P(I = 1 | X) = logit(X'pbeta)",
    "pbeta: ", sep = "\n")
cat(knitr::kable(PBETA), sep="\n")

cat("Outcome model:",
    "Y_i = 1/p_i + N(0, 1 / p_i^2)", sep = "\n")

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

sloading = model.matrix(SREG, pop)
stopifnot(all(colnames(sloading) == names(PBETA)))
sprob = log_inv(sloading %*% PBETA)

pop[[DEPENDENT]] = 1 / sprob + sapply(sprob, FUN = function(x) {rnorm(1, 0, 1 / x)})
pop_contingency =  pop_contingency %>% arrange(id)
Ntotal = NROW(pop)

cat("Finish preparing population data.", 
    "Simulation start.", sep = "\n")

# ============= Extract marginal information === #

cat("Extract marginals information.\n")

population.margins = lapply(CALIBRATE_MARGIN, FUN = function(x) {xtabs(x, pop)})
Nmargin = margin_vector(population.margins = population.margins)

print_marign(CALIBRATE_MARGIN, population.margins)
J = NROW(pop_contingency) # 100
D = NROW(Nmargin) # 59

cat("Contingency table has: ", J, " non-empty cells", "\n")

# ============= STAN Prepare =================== #

if (PROJECTION_ON) {
  cat("Start complie Stan model: projection.",
      "File: ", PROJECTION_STAN,
      "The warning message is fixed by: `target += log(fabs(0.5 * pow(f[i], -1.5)))`", sep = "\n")
  projection = stan_model(model_name = 'projection_32', 
                          file = PROJECTION_STAN)
}

if (BASIS_ON) {
  cat("Start complie Stan model: basis.",
      "File: ", BASIS_STAN,
      "The warning message is fixed by: `target += log(fabs(0.5 * pow(f[i], -1.5)))`", sep = "\n")
  basis = stan_model(model_name = 'basis_32', 
                     file = BASIS_STAN)
}

if (BAYESRAKING_ON) {
  cat("Start complie Stan model: Bayes raking.",
      "File: ", BAYESRAKING_STAN,
      "The warning message is fixed by: `target += log(fabs(0.5 * pow(f[i], -1.5)))`", sep = "\n")
  bayes_raking = stan_model(model_name = 'bayes_raking', 
                            file = BAYESRAKING_STAN)
}

# ============== Holding quantities ============ #
L_pop = loading_matrix(pop_contingency, SUBGROUP)
popy = template[pop] %>% group_by(id) %>% summarise(ysum = sum(!!sym(DEPENDENT)), Freq = n())
ymean_true = mean(pop[[DEPENDENT]])
ymarginal_true = (L_pop %*% popy$ysum) / (L_pop %*% popy$Freq)

projection_list = basis_list = bayes_raking_list = raking_list = list()
# ============================================== #

for (iter in seq_len(SIMULATION_TIME)) {
  set.seed(SEED + iter)
  cat("Seed for current simulation: ", SEED + iter, "\n")
  Ntotal = NROW(pop)
  sample_data = sampling_from_pop(pop, sprob)
  
  sample_summary = bayes_sample_contingency(sample_data, VARIABLES,
                                            CALIBRATE_MARGIN, SUBGROUP, 
                                            SREG, template)
  
  non_empty_J = sum(sample_summary$ncell != 0)
  tmp_sam = template[sample_data] %>% 
    group_by(id) %>% 
    summarise(ymean = mean(!!sym(DEPENDENT)), 
              yss = var(!!sym(DEPENDENT), na.rm = T) * (n() + 1),
              Freq = n())
  
  y_ave_non_empty = tmp_sam$ymean
  y_id = tmp_sam$id
  y_sum_of_square_non_empty = tmp_sam$yss
  y_sum_of_square_non_empty[is.na(y_sum_of_square_non_empty)] = 0
  y_total = tmp_sam$Freq
  
  ydesign_J = model.matrix(YREG, sample_summary$sample_contingency)[, names(PBETA)[PBETA != 0]]
  ydesign_non_empty = ydesign_J[tmp_sam$id, ]
  py = NCOL(ydesign_J)
  
  sample_summary$pdesign_J = sample_summary$pdesign_J[, names(PBETA)[PBETA != 0]]
  sample_summary$ps = NCOL(sample_summary$pdesign_J)
  
  if (!INDEPENDENT) {
    sclusp <- svydesign(id = ~0, weights = ~Freq, data = sample_summary$sample_contingency)
    rclusp <- rake(sclusp, sample.margins = CALIBRATE_MARGIN, population.margins = population.margins,
                   control = list(maxit = 10, epsilon = 1, verbose = FALSE))
    reg = as.formula(paste0("~", paste0(paste0(VARIABLES, collapse = ' + '))))
    r_table <- svytable(reg, rclusp, round = TRUE)
    r_f_table <- as.data.frame(r_table)
    
    temp <- template[r_f_table] %>% filter(!is.na(id)) %>% arrange(id)
    Ninit <- temp$Freq 
  } else {
    Ninit <- as.vector(Ntotal * exp(t(sample_summary$L) %*% log(Nmargin/Ntotal)))
  }
  
  C <- MASS::Null(t(sample_summary$L_quant))
  d_null = NCOL(C)
  
  data_list = c(sample_summary, 
                list(J = J, Nmargin = as.vector(Nmargin),
                     non_empty_J = non_empty_J, y_ave_non_empty = as.vector(y_ave_non_empty), 
                     y_sum_of_square_non_empty = as.vector(y_sum_of_square_non_empty),
                     y_total = as.vector(y_total), y_id = as.vector(y_id),
                     py = py, ydesign_non_empty = ydesign_non_empty, ydesign_J = ydesign_J,
                     d_null = d_null, C = C, Ninit = as.vector(Ninit)
                     )
                )
  stopifnot(all(rownames(ymarginal_true) == rownames(data_list$L_quant)))

  if (BAYESRAKING_ON) {
    bayes_raking_list[[iter]] = stan_wrapper(bayes_raking, data_list, 
                                             CHAINS, ITER, SEED, 
                                             ymean_true, ymarginal_true)
    }
  
  if (PROJECTION_ON) {
    projection_list[[iter]] = stan_wrapper(projection, data_list, 
                                           CHAINS, ITER, SEED, 
                                           ymean_true, ymarginal_true)
  }
  
  if (BASIS_ON) {
    basis_list[[iter]] = stan_wrapper(basis, data_list, 
                                      CHAINS, ITER, SEED, 
                                      ymean_true, ymarginal_true)
  }
  
  if (RAKING_ON) {

    raking_list[[iter]] = rake_wrapper(sample_data, 
                                       CALIBRATE_MARGIN, population.margins,
                                       DEPENDENT, SUBGROUP, 
                                       ymean_true, ymarginal_true)
  }

  cat("Finish ", iter, "th simulation.\n")
}

sum_plot = c()
cat("True outcome mean: ", ymean_true, '\n')

if (RAKING_ON) {
  raking_sum = summarise_list(raking_list)
  cat("Raking estimate: \n")
  print_summary(raking_sum$ymean)
  raking_sum$ymarginal['Method'] = 'Raking'
  sum_plot = rbind(sum_plot, raking_sum$ymarginal)
}

if (BAYESRAKING_ON) {
  bayes_raking_sum = summarise_list(bayes_raking_list)
  cat("Bayes raking estimate: \n")
  print_summary(bayes_raking_sum$ymean)
  bayes_raking_sum$ymarginal['Method'] = 'bayes_raking'
  sum_plot = rbind(sum_plot, bayes_raking_sum$ymarginal)
}

if (PROJECTION_ON) {
  projection_sum = summarise_list(projection_list)
  cat("Projection estimate: \n")
  print_summary(projection_sum$ymean)
  projection_sum$ymarginal['Method'] = 'projection'
  sum_plot = rbind(sum_plot, projection_sum$ymarginal)
}

if (BASIS_ON) {
  basis_sum = summarise_list(basis_list)
  cat("Basis estimate: \n")
  print_summary(basis_sum$ymean)
  basis_sum$ymarginal['Method'] = 'basis'
  sum_plot = rbind(sum_plot, basis_sum$ymarginal)
}

BiasP = sum_plot %>%
  ggplot(aes(x = Method, y = Margin)) +
  geom_tile(aes(fill = Abs.Bias), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme_pubr()

RMSEP = sum_plot %>%
  ggplot(aes(x = Method, y = Margin)) +
  geom_tile(aes(fill = RMSE), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme_pubr()

CoverageP = sum_plot %>%
  ggplot(aes(x = Method, y = Margin)) +
  geom_tile(aes(fill = Coverage), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_pubr()

SEP = sum_plot %>%
  ggplot(aes(x = Method, y = Margin)) +
  geom_tile(aes(fill = StandardErr), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_pubr()

#print(ggarrange(BiasP, RMSEP, CoverageP, SEP, ncol = 2, nrow = 2))
