/*
bayes_raking.stan with generate quantities
 */
data {
  // dimension parameters following paper section 2.1
  
  int <lower=0> D; // Number of marginals
  int <lower=0> J; // Number of cells include empty cells
  
  matrix[D, J] L; // Loading matrix (binary)
  int <lower=0> Nmargin[D]; // Known marignals

  // Observed data: cell size related
  int <lower=0> ncell[J]; // Observed cells sizes
  
  // Marigns prior
  int <lower=1> d_null; // deminsion of null space
  vector<lower=0>[J] Ninit; // N_0 in section 2.1, equation 6, 7
  matrix[J, d_null] C; // basis of Null(L)
  
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
  matrix[non_empty_J, ps] ydesign_non_empty;
  
  // Models to predict outcome: Since some of the cells are empty, but we can still estimate the potential outcome
  matrix[J, ps] ydesign_J;
  
  // Generated Quantities:
  int <lower=0> D_quant; // Dimention of marginals 
  matrix[D_quant, J] L_quant; // Loading matrix
}

transformed data {
  matrix[J, J] P; // P in section 2.1
  P = C * inverse(C' * C) * C';
}

parameters {
  vector[ps] pbeta;
  vector[J] v;
}

transformed parameters {
  vector<lower=0>[J] Nhat;
  vector[J] pselect;
  
  Nhat = Ninit + P * v;
  pselect = 1 ./ (1 + exp( - pdesign_J * pbeta));
}

model {
  // parameters related with model
  vector[non_empty_J] f;
  
  for (i in 1:non_empty_J) {
    f[i] =  1 / pselect[y_id[i]];
  }
  // Model
  
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

