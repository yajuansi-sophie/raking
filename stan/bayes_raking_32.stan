/*
bayes_raking.stan with generate quantities
 */
data {
  // dimension parameters following paper section 2.1
  
  int <lower=0> D; // Number of marginals
  int <lower=0> J; // Number of cells include empty cells
  
  matrix[D, J] L; // Loading matrix (binary)
  int <lower=0> Nmargin[D]; // Known marignals
  int <lower=0> Ntotal; // Total population size: Can be calculated based on marginals distribution

  // Observed data: cell size related
  int <lower=0> ncell[J]; // Observed cells sizes
  
  // Models for inclusion mechanism
  int <lower=0> ps; // Number of parameters in inclusion probability model section 2.4
  matrix[J, ps] pdesign_J;
  
  // Observed data: cell outcome related
  int <lower=0> non_empty_J; // Number of non-empty cells
  int <lower=1> y_id[non_empty_J]; // helper id to match the location
  vector[non_empty_J] y_ave_non_empty; // Sample cells' mean
  vector[non_empty_J] y_sum_of_square_non_empty; // Sample cell's sum of square
  vector[non_empty_J] y_total; // Sample cells size 
  
  // Models for outcome
  int py;
  matrix[non_empty_J, py] ydesign_non_empty;
  
  // Models to predict outcome: Since some of the cells are empty, but we can still estimate the potential outcome
  matrix[J, py] ydesign_J;
  
  // Generated Quantities:
  int <lower=0> D_quant; // Dimention of marginals 
  matrix[D_quant, J] L_quant; // Loading matrix
}

parameters {
  vector[ps] pbeta;
  vector<lower=0>[J] Nhat;
}

transformed parameters {
  vector[J] pselect;
  pselect = inv_logit(pdesign_J * pbeta);
}

model {
  // parameters related with model
  vector[non_empty_J] f;
  for (i in 1:non_empty_J) {
    f[i] =  1 / pselect[y_id[i]];
  }
  
  // Model
  Nmargin ~ poisson(L * Nhat); // Margins prior
  for(i in 1:J) {
    if(ncell[i] == 0)
      Nhat[i] ~ cauchy(20, 3);
  }
  ncell ~ poisson(Nhat .* pselect); // Inlcudesion mechanism
  
  for(i in 1:non_empty_J) {
    if (y_total[i] > 1) {
      y_ave_non_empty[i] ~ normal(f[i], f[i] / sqrt(y_total[i])); // Survey outcome
      (y_sum_of_square_non_empty[i] / pow(f[i], 2)) ~ chi_square(y_total[i] - 1);
      target += log(fabs(0.5 * pow(f[i], -1.5)));
    } else {
      y_ave_non_empty[i] ~ normal(f[i], f[i]); // Survey outcome
    }
  }
    
  
}

generated quantities {

  real ymean; // estimator of ymean
  vector[D_quant] ymarginals; // estimator of marginals outcome
  
  vector[J] f_pred;
  
  f_pred = 1 ./ pselect;
  
  ymean = sum(f_pred .* Nhat) / sum(Nhat);
  ymarginals = L_quant * (f_pred .* Nhat) ./ (L_quant * Nhat);
}

