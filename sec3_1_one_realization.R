rm(list = ls())

source("utils.R")

library(dplyr)
library(tidyr)
library(data.table)
library(rstan)
library(survey)
library(ggplot2)
library(knitr) 

rstan_options(auto_write = FALSE)
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

ITER = 1500
CHAINS = 4
STANFILE <- "stan/bayes_raking.stan"

SIMULATION_TIME = 1

# ====================================== # 

cat("\tREAD FIRST:\n")
cat("\tThe following code implement the simulation in section 3.1.\n")

invisible(readline(prompt="\tPress [enter] to continue. "))

cat("\tUse variables: \n")
cat(VARIABLES, sep = ", ", "\n")
cat("\n\n")

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


# ============= STAN Prepare =================== #

cat("Start complie Stan:", STANFILE, 
    "It may take a while", sep='\n')

bayes_raking = stan_model(model_name = 'bayes_raking', 
                          file = STANFILE)

# ============================================== #
# =============== Population prepare =========== #


cat("Extract marginals information.\n")

pop_contingency =  pop_contingency %>% arrange(id)

## pop is our population  data, pop_contingency is the poplation contingency table

# 2. Get marginals distribution
population.margins = lapply(CALIBRATE_MARGIN, FUN = function(x) {xtabs(x, pop)})

Nmargin = margin_vector(population.margins = population.margins)

J = NROW(pop_contingency) # 986
D = NROW(Nmargin) # 21
  
L_pop = loading_matrix(pop_contingency, SUBGROUP)
popy = template[pop] %>% group_by(id) %>% summarise(ysum = sum(!!sym(DEPENDENT)), Freq = n())
ymean_true = mean(pop[[DEPENDENT]])
ymarginal_true = (L_pop %*% popy$ysum) / (L_pop %*% popy$Freq)

# =========================================== #

# =========================================== #

# 3. Create sample to test method

## Sample from population
  
selected = sapply(sprob, FUN = function(x) {rbinom(1, 1, x) == 1})
sam = pop[selected, ]

sam_contingency = table(sam[, VARIABLES]) %>% 
  as.data.frame() %>% 
  template[.] %>%
  filter(!is.na(id))

L = loading_matrix(sam_contingency, CALIBRATE_MARGIN)

## Sample information

ncell = sam_contingency$Freq # Observed cell sizes

## Inclusion model matrix (TRUE model)

cat("We use true model for selection probability\n")
pdesign_J = model.matrix(SREG, pop_contingency)[, names(PBETA)[PBETA != 0]]
ps = NCOL(pdesign_J)

### Model for outcome

non_empty_J = sum(ncell != 0)
tmp_sam = template[sam] %>% 
  group_by(id) %>% 
  summarise(ysum = sum(!!sym(DEPENDENT)), Freq = n())
y_success = tmp_sam$ysum
y_total = tmp_sam$Freq

cat("We use true model for outcome\n")
ydesign_J = model.matrix(YREG, pop_contingency)[, names(YBETA)[YBETA != 0]]
ydesign_non_empty = ydesign_J[tmp_sam$id, ]
py = NCOL(ydesign_J)
  
  ## Loading matrix used in Figure 1
  ## Need to append interaction between age and opmres
  
L_quant = loading_matrix(pop_contingency, SUBGROUP)
D_quant = NROW(L_quant)

# ============= Y information =============== #

popy = template[pop] %>% group_by(id) %>% summarise(ysum = sum(!!sym(DEPENDENT)), Freq = n())
ymean_true = mean(pop[[DEPENDENT]])
ymarginal_true = (L_quant %*% popy$ysum) / (L_quant %*% popy$Freq)

# ============= Stan information =============== #

data_list = list(D = D, J = J, L = L, Nmargin = as.vector(Nmargin),
                  ncell = as.vector(ncell),
                  ps = ps, pdesign_J = pdesign_J, 
                  non_empty_J = non_empty_J, y_success = as.vector(y_success), y_total = as.vector(y_total), 
                  py = py, ydesign_non_empty = ydesign_non_empty, ydesign_J = ydesign_J,
                  D_quant = D_quant, L_quant = L_quant)
  
# ============================================== #
  
cat("   Start sampling with iteration: ", ITER, ", number of chains: ", CHAINS, 
    ", seed: ", SEED, "\n")

ptm <- proc.time()
braking_fit = sampling(bayes_raking, data = data_list, chains = CHAINS,
                        iter = ITER, seed = SEED, open_progress = FALSE,
                        show_messages = FALSE)
bayes_time <- proc.time() - ptm

ymean_bayes_raking = summary_wrap(braking_fit, "ymean")
ymarginal_bayes_raking = summary_wrap(braking_fit, "ymarginals")
rownames(ymarginal_bayes_raking) = rownames(L_quant)

cat("  Finish bayes raking. Total time (without complie the Stan code):", bayes_time[3])
cat("\n\n")
  
  # ============================================== #
  
# ============= Original Raking ================ #

cat("   Start use survey::rake\n")

ptm = proc.time()
design = svydesign(id = ~0, probs = NULL, data = sam)
rclus = rake(design, sample.margins = CALIBRATE_MARGIN, population.margins = population.margins)
raking_time = proc.time() - ptm

ymean_raking = as.data.frame(svymean(as.formula(paste0('~', DEPENDENT)), rclus))
ymean_raking['2.5%'] = ymean_raking[1] - 1.96 * ymean_raking[2]
ymean_raking['97.5%'] = ymean_raking[1] + 1.96 * ymean_raking[2]
class(ymean_raking) = class(ymean_bayes_raking)
names(ymean_raking) = names(ymean_bayes_raking)

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

cat("  Finish original rake. Total time: ", raking_time[3])
cat("\n\n")

cat("   True outcome mean: ", ymean_true, '\n')
cat("  BAYES RAKING overall outcome: ", 
    "\n\testimation: ", ymean_bayes_raking[1], 
    "\n\tsd: ", ymean_bayes_raking[2],
    "\n\t95% CI: (", ymean_bayes_raking[3], ",", ymean_bayes_raking[4], ")\n")
cat("  RAKING overall outcome: ", 
    "\n\testimation: ", ymean_raking[1], 
    "\n\tsd: ", ymean_raking[2], 
    "\n\t95% CI: (", ymean_raking[3], ",", ymean_raking[4], ")\n")

sum1 <- summary_tmp(ymarginal_bayes_raking, ymarginal_true, "bayes_raking")
sum2 <- summary_tmp(ymarginal_raking, ymarginal_true, "raking")

p1 = rbind(sum1, sum2) %>%
  gather(Quantities, Value, -Method, -Margin) %>%
  filter(Quantities != "Coverage") %>%
  ggplot(aes(x = Method, y = Value)) +
  geom_violin() + facet_wrap(~Quantities, scales = "free_y") +
  labs(caption = "!!!! Result only based on 1 simulation")
print(p1)
