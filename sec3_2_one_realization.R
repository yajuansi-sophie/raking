rm(list = ls())
SEED = 1234
set.seed(SEED)

source("utils.R")

library(dplyr)
library(tidyr)
library(data.table)
library(rstan)
library(survey)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

acs_data <- readRDS("Data/acs_simu.rds")

VARIABLES <- c("age_dc", "opmres_x", "childx_ca")
DEPENDENT <- "y" ## Used in simlation outcome name
acs_data = acs_data[, VARIABLES]


CALIBRATE_MARGIN = list( ~age_dc, ~childx_ca, ~opmres_x, ~age_dc+opmres_x, ~childx_ca + opmres_x)
SUBGROUP = c(~age_dc, ~childx_ca, ~opmres_x, ~childx_ca + opmres_x)

SREG = ~ age_dc + opmres_x + childx_ca + age_dc:opmres_x + childx_ca:opmres_x
YREG = SREG


ITERS = 2000
CHAINS = 4
PROJECTION_ON = TRUE
PROJECTION_STAN = "stan/projection_32.stan"
BASIS_ON = TRUE
BASIS_STAN = "stan/basis_32.stan"
BAYESRAKING_ON = TRUE
BAYESRAKING_STAN = "stan/bayes_raking_32.stan"
RAKING_ON = TRUE


cat("\tREAD FIRST:\n")
cat("\tThe following code only implement one realization of simulation in section 3.2.\n")

invisible(readline(prompt="Press [enter] to continue. "))

cat("\tUse variables: \n\t")
cat(VARIABLES, sep = ", ")
cat("\n\tUse additional marginals:")
cat("\n\tage_dc:opmres_x, childx_ca:opmres_x")
cat("\n\n")
cat("\tFor Stan, use", CHAINS, "chains, and", ITERS, "iterations. seed:", SEED, "\n")

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

# ============== Holding quantities ============ #
L_pop = loading_matrix(pop_contingency, SUBGROUP)
popy = template[pop] %>% group_by(id) %>% summarise(ysum = sum(!!sym(DEPENDENT)), Freq = n())
ymean_true = mean(pop[[DEPENDENT]])
ymarginal_true = (L_pop %*% popy$ysum) / (L_pop %*% popy$Freq)

# =========================================== #

# 3. Create sample to test method

## Sample from population
selected = sapply(sprob, FUN = function(x) {rbinom(1, 1, x) == 1})
sam = pop[selected, ]

sam_contingency = table(sam[, VARIABLES]) %>% 
  as.data.frame() %>% 
  template[.] %>%
  filter(!is.na(id))

## Sample information

ncell = sam_contingency$Freq # Observed cell sizes
L = loading_matrix(sam_contingency, CALIBRATE_MARGIN)
stopifnot(all(L == loading_matrix(pop_contingency, CALIBRATE_MARGIN)))
stopifnot(D == NROW(L))

## Inclusion model matrix

pdesign_J = model.matrix(SREG, sam_contingency)[,names(PBETA)[PBETA != 0]]
ps = NCOL(pdesign_J)

### Model for outcome

non_empty_J = sum(ncell != 0)

tmp_sam = template[sam] %>% 
  group_by(id) %>% 
  summarise(ymean = mean(!!sym(DEPENDENT)), 
            yss = var(!!sym(DEPENDENT), na.rm = T) * (n() + 1),
            Freq = n())

y_ave_non_empty = tmp_sam$ymean
y_id = tmp_sam$id
y_sum_of_square_non_empty = tmp_sam$yss
y_sum_of_square_non_empty[is.na(y_sum_of_square_non_empty)] = 0
y_total = tmp_sam$Freq

ydesign_J = model.matrix(YREG, pop_contingency)[,names(PBETA)[PBETA != 0]]
ydesign_non_empty = ydesign_J[tmp_sam$id, ]
py = NCOL(ydesign_J)

# ============= Y information =============== #

L_quant = loading_matrix(pop_contingency, SUBGROUP)
D_quant = NROW(L_quant)

popy = template[pop] %>% group_by(id) %>% summarise(ysum = sum(!!sym(DEPENDENT)), Freq = n())
ymean_true = mean(pop[[DEPENDENT]])
ymarginal_true = (L_quant %*% popy$ysum) / (L_quant %*% popy$Freq)

# ================= Summary Utilites ================= #

sum_margin_plot = c()


# ============= Original Raking ================ #

if (RAKING_ON) {
  
  cat("\tStart use survey::rake\n")
  
  ptm = proc.time()
  design = svydesign(id = ~0, probs = NULL, data = sam)
  
  rclus = rake(design, sample.margins = CALIBRATE_MARGIN, population.margins = population.margins)
  
  raking_time = proc.time() - ptm
  ymean_raking = as.data.frame(svymean(as.formula(paste0('~', DEPENDENT)), rclus))
  ymean_raking['2.5%'] = ymean_raking[1] - 1.96 * ymean_raking[2]
  ymean_raking['97.5%'] = ymean_raking[1] + 1.96 * ymean_raking[2]
  class(ymean_raking) = "numeric"
  names(ymean_raking) = c("mean", "sd", "2.5%", "97.5%")
  
  ymarginal_raking = data.frame()
  for (i in SUBGROUP) {
    vars = all.vars(i)
    test = svyby(as.formula(paste0('~', DEPENDENT)), i, rclus, svymean)
    N = NROW(test)
    tmpname = rep(NA, N)
    for (j in 1:N) {
      varname = rep(N, length(vars))
      for (k in 1:length(vars)) varname[k] = paste0(vars[k], test[j, vars[k]])
      tmpname[j] = paste(varname, collapse = ":")
    }
    rownames(test) = tmpname
    test[vars] = NULL
    ymarginal_raking = rbind(ymarginal_raking, test)
  }
  
  ymarginal_raking['2.5%'] = ymarginal_raking[, DEPENDENT] - 1.96 * ymarginal_raking[, "se"]
  ymarginal_raking['97.5%'] = ymarginal_raking[, DEPENDENT] + 1.96 * ymarginal_raking[, "se"]
  colnames(ymarginal_raking) = c("mean", "sd", "2.5%", "97.5%")
  
  
  ###
  
  ymarginal_raking['2.5%'] = with(ymarginal_raking, mean - 1.96 * sd)
  ymarginal_raking['97.5%'] = with(ymarginal_raking, mean + 1.96 * sd)
  rownames(ymarginal_raking) = rownames(L_quant)
  sum_margin_plot <- rbind(sum_margin_plot, 
                           summary_tmp(ymarginal_raking, ymarginal_true, "raking"))
  cat("\tFinish original rake. Total time: ", raking_time[3])
  cat("\n\n")  
  
}

# ============= Choose the initial ============= #

if (PROJECTION_ON || BASIS_ON) {
  cat("\tSince the scripts will perform PROJECTION or BASIS method.\n",
      "\tFor initial N_0 in section 2.4 margin prior\n",
      "\tPlease choose: 1. independent (May fail to start sampling), 2. raking estimator.\n")
  independent = readline(prompt="1. indenpent, 2. raking estimator: ")
  while (!independent %in% c("1", "2")) {
    cat("\tWrong! Please enter 1 or 2. \n")
    independent = readline(prompt="\t1(indenpent), 2(raking estimator): ")
  }
  if (as.integer(independent) == 1) {
    Ninit <- as.vector(Ntotal * exp(t(L_quant) %*% log(Nmargin/Ntotal)))
  } else {
    sclusp <- svydesign(id = ~0, weights = ~Freq, data = sam_contingency)
    var.formula = population.margins = sample.margins = list()
    
    for (i in seq_along(VARIABLES) ) {
      var.formula[[i]] = sample.margins[[i]] = as.formula(paste0('~', VARIABLES[i])) 
      population.margins[[i]] = xtabs(sample.margins[[i]], pop)
    }
    var.formula[[4]] = sample.margins[[4]] = ~ opmres_x + age_dc
    var.formula[[5]] = sample.margins[[5]] = ~ opmres_x + childx_ca
    population.margins[[4]] = xtabs(sample.margins[[4]], pop)
    population.margins[[5]] = xtabs(sample.margins[[5]], pop)
    
    rclusp <- rake(sclusp, sample.margins = sample.margins, population.margins = population.margins,
                    control = list(maxit = 10, epsilon = 1, verbose = FALSE))
    reg = as.formula(paste0("~", paste0(paste0(VARIABLES, collapse = ' + '))))
    r_table <- svytable(reg, rclusp, round = TRUE)
    r_f_table <- as.data.frame(r_table)

    temp <- template[r_f_table] %>% filter(!is.na(id)) %>% arrange(id)
    Ninit <- temp$Freq  
  }
}


C <- MASS::Null(t(L_quant))
d_null = NCOL(C)


# ============= Stan information =============== #

data_list = list(D = D, J = J, L = L, Nmargin = as.vector(Nmargin),
                 ncell = as.vector(ncell),
                 ps = ps, pdesign_J = pdesign_J, 
                 non_empty_J = non_empty_J, 
                 y_ave_non_empty = as.vector(y_ave_non_empty), 
                 y_sum_of_square_non_empty = as.vector(y_sum_of_square_non_empty),
                 y_total = as.vector(y_total), y_id = as.vector(y_id),
                 py = py, ydesign_non_empty = ydesign_non_empty, ydesign_J = ydesign_J,
                 D_quant = D_quant, L_quant = L_quant,
                 d_null = d_null, C = C, Ninit = as.vector(Ninit))

# ============= Projection  ==================== #

if (PROJECTION_ON) {
  cat("\tStart complie Stan model: projection.\n",
      "\tFile: ", PROJECTION_STAN, "\n",
      "\tThe warning message is fixed by: `target += log(fabs(0.5 * pow(f[i], -1.5)))`\n",
      "\tIt may take a while. \n\n")
  
  projection = stan_model(model_name = 'projection_32', 
                          file = PROJECTION_STAN)
  
  cat("\tStart sampling with iteration: ", ITERS, ", number of chains: ", CHAINS, 
      ", seed: ", SEED, "\n")
  
  ptm <- proc.time()
  projection_fit = sampling(projection, data = data_list, chains = CHAINS,
                         iter = ITERS, seed = SEED, open_progress = FALSE,
                         show_messages = FALSE)
  projection_time <- proc.time() - ptm
  
  ymean_projection = summary_wrap(projection_fit, "ymean")
  ymarginal_projection = summary_wrap(projection_fit, "ymarginals")
  rownames(ymarginal_projection) = rownames(L_quant)
  
  sum_margin_plot <- rbind(sum_margin_plot, 
                           summary_tmp(ymarginal_projection, ymarginal_true, "projection"))
  
  cat("\tFinish projection. Total time (without complie the Stan code):", projection_time[3])
  cat("\n\n")
}


# ============================================== #

# =================== Basis ==================== #

if (BASIS_ON) {
  cat("\tStart complie Stan model: basis.\n",
      "\tFile: ", BASIS_STAN,"\n",
      "\tThe warning message is fixed by: `target += log(fabs(0.5 * pow(f[i], -1.5)))`\n",
      "\tIt may take a while. \n\n")
  
  basis = stan_model(model_name = 'basis_32', 
                     file = BASIS_STAN)
  
  cat("\tStart sampling with iteration: ", ITERS, ", number of chains: ", CHAINS, 
      ", seed: ", SEED, "\n")
  
  ptm <- proc.time()
  basis_fit = sampling(basis, data = data_list, chains = CHAINS,
                       iter = ITERS, seed = SEED, open_progress = FALSE,
                       show_messages = FALSE)
  basis_time <- proc.time() - ptm
  
  ymean_basis = summary_wrap(basis_fit, "ymean")
  ymarginal_basis = summary_wrap(basis_fit, "ymarginals")
  rownames(ymarginal_basis) = rownames(L_quant)
  
  sum_margin_plot <- rbind(sum_margin_plot, 
                           summary_tmp(ymarginal_basis, ymarginal_true, "basis"))
  
  cat("\tFinish basis. Total time (without complie the Stan code):", basis_time[3])
  cat("\n\n")
}

# ============================================== #

# ============= Bayes raking =================== #

if (BAYESRAKING_ON) {
  cat("\tStart complie Stan model: bayes_raking.\n",
      "\tFile:", BAYESRAKING_STAN,"\n",
      "\tThe warning message is fixed by: `target += log(fabs(0.5 * pow(f[i], -1.5)))`\n",
      "\tIt may take a while. \n\n")
  
  bayes_raking = stan_model(model_name = 'bayes_raking_32', 
                            file = BAYESRAKING_STAN)
  
  cat("\tStart sampling with iteration: ", ITERS, ", number of chains: ", CHAINS, 
      ", seed: ", SEED, "\n")
  
  ptm <- proc.time()
  bayes_raking_fit = sampling(bayes_raking, data = data_list, chains = CHAINS,
                              iter = ITERS, seed = SEED, open_progress = FALSE,
                              show_messages = FALSE)
  bayes_raking_time <- proc.time() - ptm
  
  ymean_bayes_raking = summary_wrap(bayes_raking_fit, "ymean")
  ymarginal_bayes_raking = summary_wrap(bayes_raking_fit, "ymarginals")
  rownames(ymarginal_bayes_raking) = rownames(L_quant)
  
  sum_margin_plot <- rbind(sum_margin_plot, 
                           summary_tmp(ymarginal_bayes_raking, ymarginal_true, "bayes-model"))
  
  cat("\tFinish bayes_raking. Total time (without complie the Stan code):", bayes_raking_time[3])
  cat("\n\n")
}

# ============================================== #


cat("\tTrue outcome mean: ", ymean_true, '\n')
if (PROJECTION_ON) {
  cat("\tPROJECTION overall outcome: ", 
      "\n\testimation: ", ymean_projection[1], 
      "\n\tsd: ", ymean_projection[2],
      "\n\t95% CI: (", ymean_projection[3], ",", ymean_projection[4], ")\n")
}

if (BASIS_ON) {
  cat("\tBasis overall outcome: ", 
      "\n\testimation: ", ymean_basis[1], 
      "\n\tsd: ", ymean_basis[2],
      "\n\t95% CI: (", ymean_basis[3], ",", ymean_basis[4], ")\n")
}

if (BAYESRAKING_ON) {
  cat("\tBAYES RAKING overall outcome: ", 
      "\n\testimation: ", ymean_bayes_raking[1], 
      "\n\tsd: ", ymean_bayes_raking[2],
      "\n\t95% CI: (", ymean_bayes_raking[3], ",", ymean_bayes_raking[4], ")\n")
}

if (RAKING_ON) {
  cat("\tRAKING overall outcome: ", 
      "\n\testimation: ", ymean_raking[1], 
      "\n\tsd: ", ymean_raking[2], 
      "\n\t95% CI: (", ymean_raking[3], ",", ymean_raking[4], ")\n")  
}

BiasP = sum_margin_plot %>%
  ggplot(aes(x = Method, y = Margin)) +
  geom_tile(aes(fill = abs(Bias)), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")

MSEP = sum_margin_plot %>%
  ggplot(aes(x = Method, y = Margin)) +
  geom_tile(aes(fill = SquareErr), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")

CoverageP = sum_margin_plot %>%
  ggplot(aes(x = Method, y = Margin)) +
  geom_tile(aes(fill = Coverage), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue") 

SEP = sum_margin_plot %>%
  ggplot(aes(x = Method, y = Margin)) +
  geom_tile(aes(fill = StandardErr), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")

BiasP
MSEP
CoverageP
SEP
