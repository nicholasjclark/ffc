// Stan model code for autoregressive dynamic factor forecasting models
functions {
  vector rep_each(vector x, int K) {
    int N = rows(x);
    vector[N * K] y;
    int pos = 1;
    for (n in 1 : N) {
      for (k in 1 : K) {
        y[pos] = x[n];
        pos += 1;
      }
    }
    return y;
  }
}
data {
  int<lower=0> n; // number of timepoints per series
  int<lower=0> n_lv; // number of dynamic factors
  int<lower=0> n_series; // number of series
  vector[n_series] sample_var; // known sampling variance
  int<lower=0> M; // number of nonzero lower-triangular factor loadings
  int<lower=0> n_nonmissing; // number of nonmissing observations
  vector[n_nonmissing] flat_ys; // flattened nonmissing observations
  array[n_nonmissing] int<lower=0> obs_ind; // indices of nonmissing observations
  int<lower=1> family; // 1 = normal, 2 = student-t
  vector[2] prior_alpha; // prior intercept control parameters
  vector[3] prior_sigma; // prior obs error control parameters
  vector[3] prior_ar; // prior ar coefficient control parameters
  int<lower=1> P; // AR order
  vector[1] beta; // beta coefficient (1) to use in id_glm functions
}
transformed data {

}
parameters {
  // dynamic factor AR terms
  array[P, n_lv] real<lower=-1, upper=1> ar;

  // dynamic factor process errors
  matrix[n, n_lv] LV_z;

  // factor loading lower triangle
  vector[M] L;

  // series-level intercepts
  vector[n_series] alpha;

  // series-level observation errors
  vector<lower=0>[n_series] sigma_obs;

  // observation df for student-t
  real<lower=0> nu;

}
transformed parameters {
  // trends and dynamic factor loading matrix
  matrix[n, n_series] trend;
  matrix[n_series, n_lv] Lambda = rep_matrix(0, n_series, n_lv);

  // factor loadings, with constraints
  {
    int index;
    index = 0;
    for (j in 1 : n_lv) {
      for (s in j : n_series) {
        index = index + 1;
        Lambda[s, j] = L[index];
      }
    }
  }

  // derived latent factors
  matrix[n, n_lv] LV  = rep_matrix(0, n, n_lv);
  LV[1, 1 : n_lv] = LV_z[1, 1 : n_lv];
  for (j in 1 : n_lv) {
    for (p in 1 : P) {
      for (t in (P + 1) : n){
        LV[t, j] += ar[p, j] * LV[t - p, j] + LV_z[t, j];
      }
    }
  }

  // derived series-level trends
  for (i in 1 : n) {
    for (s in 1 : n_series) {
      trend[i, s] = alpha[s] + dot_product(Lambda[s,  : ], LV[i,  : ]);
    }
  }
}
model {
  // priors for series-level intercepts
  target += normal_lpdf(alpha | prior_alpha[1], prior_alpha[2]);

  // prior for observation df parameter
  nu ~ gamma(2, 0.1);

  // priors for observation error parameters
  if (prior_sigma[3] == 1) target += normal_lpdf(sigma_obs | prior_sigma[1], prior_sigma[2]);
  else if (prior_sigma[3] == 2) target += beta_lpdf(sigma_obs | prior_sigma[1], prior_sigma[2]);
  else if (prior_sigma[3] == 3) target += inv_gamma_lpdf(sigma_obs | prior_sigma[1], prior_sigma[2]);
  else if (prior_sigma[3] == 4) target += exponential_lpdf(sigma_obs | prior_sigma[1]);

  // priors for AR parameters
  for (j in 1 : n_lv) {
    for (p in 1 : P) {
      if (prior_ar[3] == 1) target += std_normal_lpdf(ar[p, j]);
      else if (prior_ar[3] == 2) target += normal_lpdf(ar[p, j] | prior_ar[1], prior_ar[2]);
    }
  }

  // priors for factor loading coefficients
  L ~ std_normal();

  // factor process errors
  to_vector(LV_z) ~ std_normal();

  {
    // likelihood functions
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
  // adjusted observation error
  vector<lower=0>[n_series] sigma_obs_adj;
  for (s in 1 : n_series){
    sigma_obs_adj[s] = sqrt(square(sigma_obs[s]) + sample_var[s]);
  }

  // Posterior predictions
  array[n, n_series] real ypred;
  matrix[n, n_series] sigma_obs_vec;
  if(family == 1) {
    for (s in 1 : n_series) {
      sigma_obs_vec[1 : n, s] = rep_vector(sigma_obs_adj[s], n);
      ypred[1 : n, s] = normal_rng(trend[1 : n, s], sigma_obs_vec[1 : n, s]);
    }
  }

  if(family == 2) {
    vector[n] nu_vec = rep_vector(nu, n);
    for (s in 1 : n_series) {
      sigma_obs_vec[1 : n, s] = rep_vector(sigma_obs_adj[s], n);
      ypred[1 : n, s] = student_t_rng(nu_vec, trend[1 : n, s],
                                      sigma_obs_vec[1 : n, s]);
    }
  }
}
