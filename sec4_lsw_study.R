rm(list = ls())
SEED = 1234
set.seed(SEED)

source('utils.R')

library(dplyr)
library(data.table)
library(rstan)
library(survey)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

acs_data <- readRDS("Data/acs.rds")
lsw_data <- readRDS("Data/lsw.rds")

VARIABLES = c("age_dc", "educat", "sex", "race_dc", "opmres_x", "childx_ca")
DEPENDENT = 'anyhard'
WEIGHT = "perwt"

MAIN = list( ~age_dc, ~educat, ~sex, ~race_dc, ~opmres_x, ~childx_ca)
INTERACTION = list(~age_dc + opmres_x, ~childx_ca + opmres_x) 
SREG = ~ age_dc + educat + sex + race_dc + opmres_x + childx_ca + opmres_x:age_dc + opmres_x:childx_ca
YREG = ~ age_dc + educat + sex + race_dc + opmres_x + childx_ca

STANFILE = "bayes_raking_horseshoe.stan"
ITER = 500
CHAINS = 4


# =============== Population marginals =========== #

acs_design = svydesign(id = ~1, weights = ~perwt, data = acs_data)

sample.margins = c(MAIN, INTERACTION)
population.margins = lapply(sample.margins, FUN = function(x) {svytable(x, acs_design)})
acs_contingency = svytable(as.formula(paste0("~", paste0(VARIABLES, collapse = "+"))), acs_design) %>%
  as.data.frame()

Nmargin = margin_vector(population.margins = population.margins)
Ntotal = sum(acs_contingency$Freq)

if (Ntotal > 5e5) {
  
  warning(paste("The population size is", Ntotal,", which is larger than 50000. it may affect inclusion model parameters estimation.\nDo you want to re-scale the population total cell size?"))
  rescale = readline(prompt = "It will scale the population size with out affect the marginal distribution (Y/N): ")
  if (rescale == "Y") {
    Nmargin = as.integer(50000 / Ntotal * Nmargin)
    cat("The total population cell size is rescale to: ", 50000)
  } else {
    warning("The current prior setting may affect bayes-raking convergence!")
  }
}

# =============== LSW sample contingency ========= #

lsw_data_use <- lsw_data[, c(VARIABLES, DEPENDENT)] %>% 
  filter(anyhard != "NA") %>%
  mutate(anyhard = as.integer(anyhard)) %>%
  as.data.frame
lsw_data_use[, VARIABLES] = as.data.frame(sapply(lsw_data_use[, VARIABLES], factor))
lsw_data_use = lsw_data_use[, c(VARIABLES, "anyhard")]

lsw_contingency = table(lsw_data_use[, VARIABLES]) %>%
  as.data.frame() %>%
  mutate(id = 1:n())

## Here update L just in case lsw_contingency has different order as acs_contingency
L = loading_matrix(lsw_contingency, sample.margins)
D = NROW(L)
J = NROW(lsw_contingency)
ncell = lsw_contingency$Freq

template <- lsw_contingency %>% as.data.table
setkeyv(template, VARIABLES)

# =============== LSW selection models =========== #
pdesign_J = model.matrix(SREG, lsw_contingency)
ps = NCOL(pdesign_J)

# =============== LSW outcome models ============= #

ch_fo <- paste(as.character(YREG)[c(2,1,3)], collapse = " ")

cat("\t LSW outcome model use logistic regression with formula:\n",
    strwrap(ch_fo), "\n")

tmp_lsw <- template[lsw_data_use] %>% 
  group_by(id) %>%
  summarise(ysum = sum(anyhard), Freq = n())
non_empty_J = sum(ncell != 0)
y_success = tmp_lsw$ysum
y_total = tmp_lsw$Freq

ydesign_J = model.matrix(YREG, lsw_contingency)
ydesign_non_empty = ydesign_J[tmp_lsw$id, ]
py = NCOL(ydesign_J)

L_quant = loading_matrix(lsw_contingency, sample.margins = sample.margins)
D_quant = NROW(L_quant)

# =============== Bayes Raking =================== #
data_list = list(D = D, J = J, L = L, Nmargin = as.vector(Nmargin),
                 Ntotal = Ntotal, ncell = as.vector(ncell),
                 ps = ps, pdesign_J = pdesign_J, 
                 non_empty_J = non_empty_J, y_success = as.vector(y_success), 
                 y_total = as.vector(y_total), py = py, 
                 ydesign_non_empty = ydesign_non_empty, ydesign_J = ydesign_J,
                 D_quant = D_quant, L_quant = L_quant)

cat("   Start complie Stan model: bayes_raking. It may take a while\n")
bayes_raking = stan_model(model_name = 'bayes_raking', 
                          file = STANFILE)

cat("   Start sampling with iteration: ", ITER, ", number of chains: ", CHAINS, 
    ", seed: ", SEED, "\n")

ptm <- proc.time()
braking_fit = sampling(bayes_raking, data = data_list, chains = CHAINS,
                       iter = ITER, seed = SEED, open_progress = FALSE,
                       show_messages = FALSE)
bayes_time <- proc.time() - ptm
print(braking_fit, "ymean")
ymean_bayes_raking = summary_wrap(braking_fit, "ymean")
ymarginal_bayes_raking = summary_wrap(braking_fit, "ymarginals")
rownames(ymarginal_bayes_raking) = rownames(L_quant)
cat("  Finish bayes raking. Total time (without complie the Stan code):", bayes_time[3])
cat("\n\n")


# =========Original Raking ========================= #

ptm = proc.time()
design = svydesign(id = ~0, probs = NULL, data = lsw_data_use)

rclus = rake(design, sample.margins = sample.margins, population.margins = population.margins,
             control = list(maxit = 50))
raking_time = proc.time() - ptm

ymean_raking = as.data.frame(svymean(as.formula(paste0('~', DEPENDENT)), rclus))
ymean_raking['2.5%'] = ymean_raking[1] - 1.96 * ymean_raking[2]
ymean_raking['97.5%'] = ymean_raking[1] + 1.96 * ymean_raking[2]
class(ymean_raking) = "numeric"
names(ymean_raking) = c("mean", "sd", "2.5%", "97.5%")

ymarginal_raking = data.frame()
for (i in sample.margins) {
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

# =============== Plot =========================== #

cat("\tBAYES RAKING overall outcome: ", 
    "\n\testimation: ", ymean_bayes_raking[1], 
    "\n\tsd: ", ymean_bayes_raking[2],
    "\n\t95% CI: (", ymean_bayes_raking[3], ",", ymean_bayes_raking[4], ")\n")

cat("\tRAKING overall outcome: ", 
    "\n\testimation: ", ymean_raking[1], 
    "\n\tsd: ", ymean_raking[2], 
    "\n\t95% CI: (", ymean_raking[3], ",", ymean_raking[4], ")\n") 

ymarginal_bayes_raking = as.data.frame(ymarginal_bayes_raking)
ymarginal_bayes_raking['Method'] = 'Bayes Raking'
ymarginal_bayes_raking['Margin'] = rownames(ymarginal_bayes_raking)
ymarginal_raking['Method'] = 'Raking'
ymarginal_raking['Margin'] = rownames(ymarginal_raking)

p1 = rbind(ymarginal_bayes_raking, ymarginal_raking) %>%
  filter(grepl(":", Margin)) %>%
  ggplot(aes(x = Margin, y = mean, color = Method, shape = Method)) + 
  geom_point() + coord_flip() +
  geom_errorbar(aes(ymin = mean - 1.96*sd, ymax = mean + 1.96 * sd)) +
  theme_pubclean() + labs(x = NULL, y = NULL) + 
  theme(legend.position="right")

print(p1)