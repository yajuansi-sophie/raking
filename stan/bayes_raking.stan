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
  
  // Models for inclusion mechanism
  int <lower=0> ps; // Number of parameters in inclusion probability model section 2.4
  matrix[J, ps] pdesign_J;
  
  // Observed data: cell outcome related
  int <lower=0> non_empty_J; // Number of non-empty cells
  int y_success[non_empty_J]; // Sample cells' mean
  int y_total[non_empty_J]; // Sample cells size 
  
  // Models for outcome
  int <lower=0> py; //Number of parameters
  //matrix[non_empty_J, py] ydesign_non_empty;
  row_vector[py] ydesign_non_empty[non_empty_J];
  
  // Models to predict outcome: Since some of the cells are empty, but we can still estimate the potential outcome
  matrix[J, py] ydesign_J;
  
  // Generated Quantities:
  int <lower=0> D_quant; // Dimention of marginals 
  matrix[D_quant, J] L_quant; // Loading matrix
}

parameters {
  vector[ps] pbeta;
  vector[py] ybeta;
  vector<lower=0>[J] Nhat;
}

model {
  // parameters related with model
  
  vector[J] pselect;
  vector[non_empty_J] f;
  
  pselect = inv_logit(pdesign_J * pbeta);
  // pselect = 1 ./ (1 + exp( - pdesign_J * pbeta));
  //f = 1 ./ (1 + exp( - ydesign_non_empty * ybeta));
  
  // Model
  Nmargin ~ poisson(L * Nhat); // Margins prior
  for(i in 1:J){
    if(ncell[i] == 0)
      Nhat[i] ~ cauchy(10, 3);
  }
  ncell ~ poisson(Nhat .* pselect); // Inlcudesion mechanism
  
  for(i in 1:non_empty_J)
    y_success[i] ~ binomial(y_total[i], inv_logit(ydesign_non_empty[i] * ybeta)); // Survey outcome
  
}

generated quantities {

  real ymean; // estimator of ymean
  vector[D_quant] ymarginals; // estimator of marginals outcome
  
  vector[J] f_pred;

  // f_pred = 1 ./ (1 + exp( - ydesign_J * ybeta));
  f_pred = inv_logit(ydesign_J * ybeta);
  
  ymean = sum(f_pred .* Nhat) / sum(Nhat);
  ymarginals = L_quant * (f_pred .* Nhat) ./ (L_quant * Nhat);
}

