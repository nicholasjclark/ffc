// Stan model code for GP dynamic factor forecasting models
functions {
  vector rep_each(vector x, int K) {
    int N = rows(x);
    vector[N * K] y;

    for (n in 1 : N) {
     // Compute start and end positions for the block
      int start = (n - 1) * K + 1;
      int end = n * K;

      // Assign the same value to a block of y
      y[start : end] = rep_vector(x[n], K);
    }

   return y;
  }

  /* Spectral density GP eigenvalues*/
  /* see Riutort-Mayol et al 2023 for details (https://doi.org/10.1007/s11222-022-10167-2)*/
  real lambda_gp(real L, int m) {
    return square((m * pi()) / (2 * L));
  }

  /* Spectral density GP eigenfunctions*/
  /* see Riutort-Mayol et al 2023 for details (https://doi.org/10.1007/s11222-022-10167-2)*/
  vector phi_SE(real L, int m, vector x) {
    return inv_sqrt(L) * sin((m * pi()) / (2 * L) * (x + L));

  }

  /* Spectral density squared exponential Gaussian Process*/
  /* see Riutort-Mayol et al 2023 for details (https://doi.org/10.1007/s11222-022-10167-2)*/
  vector spd_SE(real alpha, real rho, vector w) {
    return square(alpha) * sqrt(2 * pi()) * rho * exp(-0.5 * square(rho * w));
  }
}

data {
  int<lower=0> n; // number of timepoints per series
  int<lower=0> K; // number of dynamic factors
  int<lower=0> n_series; // number of series
  vector[n_series] sample_sd; // known sampling SDs
  int<lower=0> M; // number of nonzero lower-triangular factor loadings
  int<lower=0> n_nonmissing; // number of nonmissing observations
  vector[n_nonmissing] flat_ys; // flattened nonmissing observations
  array[n_nonmissing] int<lower=0> obs_ind; // indices of nonmissing observations
  int<lower=1> family; // 1 = normal, 2 = student-t
  row_vector[n_series] alpha; // series intercepts
  vector[1] beta; // beta coefficient (1) to use in id_glm functions
}

transformed data {
 // Create centered time vector and GP basis functions
 vector[n] times = linspaced_vector(n, 1, n);
 vector[n] times_cent = times - mean(times);
 real boundary = 1.25 * (max(times_cent) - min(times_cent));
 int num_gp_basis = min(25, n);
 matrix[n, num_gp_basis] gp_phi;
 vector[num_gp_basis] lambda_gp_vals;

 // Compute GP basis functions and lambdas
 for (m in 1 : num_gp_basis) {
   gp_phi[, m] = phi_SE(boundary, m, times_cent);
   lambda_gp_vals[m] = sqrt(lambda_gp(boundary, m));
 }
}

parameters {
  // gp length scale parameters
  vector<lower=0>[K] rho_gp;

  // gp coefficient weights
  matrix[num_gp_basis, K] b_gp;

  // factor loading lower triangle
  vector[M] L;

  // observation df for student-t
  real<lower=0> nu;

  // observation variances
  vector<lower=0>[n_series] sigma_obs;
}

transformed parameters {
  // dynamic factor loading matrix
  matrix[n_series, K] Lambda = rep_matrix(0, n_series, K);

  // factor loadings, with constraints
  {
    int index = 0;
    for (j in 1 : K) {
      for (s in j : n_series) {
        index += 1;
        Lambda[s, j] = L[index];
      }
    }
  }

  // derived latent factors
  matrix[num_gp_basis, K] SPD_beta;
  for (j in 1 : K) {
    SPD_beta[, j] = sqrt(spd_SE(1.0, rho_gp[j], lambda_gp_vals)) .* b_gp[, j];
  }
  matrix[n, K] LV = gp_phi * SPD_beta;

  // derived series-level trends
  matrix[n, n_series] trend = LV * Lambda' + rep_matrix(alpha', n)';
}

model {
  // prior for observation df parameter
  nu ~ gamma(2, 0.1);

  // priors for gp parameters
  to_vector(b_gp) ~ std_normal();
  rho_gp ~ inv_gamma(1.499007, 5.670433);

  // priors for factor loading coefficients
  L ~ double_exponential(0, 0.5);

  // priors for observation error
  sigma_obs ~ exponential(4);

  {
    // likelihood functions
    // Sampling SDs for vectorized likelihood calculations
    vector[n_nonmissing] flat_sigma_obs = rep_each(sigma_obs, n)[obs_ind];
    if(family == 1) {
      matrix[n_nonmissing, 1] flat_trends;
      flat_trends[1 : n_nonmissing, 1] = to_vector(trend)[obs_ind];
      flat_ys ~ normal_id_glm(flat_trends, 0.0, beta, flat_sigma_obs);
    }

    if(family == 2) {
      vector[n_nonmissing] flat_trends = to_vector(trend)[obs_ind];
      flat_ys ~ student_t(nu, flat_trends, flat_sigma_obs);
    }
  }
}

generated quantities {
  array[n, n_series] real ypred;
  vector[n] nu_vec;
  vector[n_series] sigma_combined;
  matrix[n, n_series] sigma_obs_vec;

  // Precompute combined standard deviations
  for (s in 1 : n_series)
    sigma_combined[s] = sqrt(square(sigma_obs[s]) + square(sample_sd[s]));

  // Broadcast to full matrix
  for (s in 1 : n_series)
    sigma_obs_vec[, s] = rep_vector(sigma_combined[s], n);

  // Posterior predictions
  if (family == 1) {
    for (s in 1 : n_series)
      ypred[, s] = normal_rng(trend[, s], sigma_obs_vec[, s]);
  } else if (family == 2) {
    nu_vec = rep_vector(nu, n);
    for (s in 1 : n_series)
      ypred[, s] = student_t_rng(nu_vec, trend[, s], sigma_obs_vec[, s]);
  }
}
