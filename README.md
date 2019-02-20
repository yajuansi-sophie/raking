# Code to reproduce paper

[Bayes-raking: Bayesian Finite Population Inference with Known Margins](https://arxiv.org/abs/1901.02117)



# layouts
  
1. Data:

  - [`acs_simu.rds`](Data/acs_simu.rds): Data used for simulation, which mimic ACS
  - [`lsw.rds`](Data/lsw.rds): Data used for real data analysis. **EXCLUDED IN THIS REPO**

2. Scripts:

    - [utils.R](utils.R) some functions used in simulation.
      - `loading_matrix` used to generate matrix $L$, describe in section 2.1 eq(5)
      - `margin_vector` used to generate vector $\overrightarrow{N}_{..}$
      - `bayes_sample_contingency` a small wrapper to genenar a sample contingency table information for bayes raking. Used similar grammar as R package [`survey::rake`](https://www.rdocumentation.org/packages/survey/versions/3.35/topics/rake)
      - NO NEED TO READ

    a. Section 3.1
      - [`sec3_1_simu.R`](sec3_1_simu.R) TIME CONSUMING
      - [`sec3_1_one_realization.R`](sec3_1_one_realization.R) With detail step how to generate loading matrix L.
      - [`stan/bayes_raking.stan`](stan/bayse_raking.stan) Bayes raking + Binary outcome

    b. Section 3.2

      - [`sec3_2_simu.R`](sec3_2_simu.R) TIME CONSUMING
      - [`sec3_2_one_realization.R`](sec3_2_one_realization.R) With detail step how to generate loading matrix L.
      - [`bayes_raking_32.stan`](stan/bayes_raking_32.stan) Bayes raking + continue outcome
      - [`basis_32.stan`](stan/basis_32.stan) Basis + continue outcome
      - [`projection_32.stan`](stan/projection_32.stan) Projection + continue outcome

    c. Section 4 (Real data, require Data/lsw.rds)

      - [`sec4_lsw_study.R`](sec4_lsw_study.R) !!**REQUIRE** lsw.rds required
      - [`bayes_raking_horseshoe.stan`](stan/bayes_raking_horseshoe.stan)
