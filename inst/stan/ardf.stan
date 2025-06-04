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
  int<lower=0> K; // number of dynamic factors
  int<lower=0> n_series; // number of series
  vector[n_series] sample_sd; // known sampling SDs
  int<lower=0> M; // number of nonzero lower-triangular factor loadings
  int<lower=0> n_nonmissing; // number of nonmissing observations
  vector[n_nonmissing] flat_ys; // flattened nonmissing observations
  array[n_nonmissing] int<lower=0> obs_ind; // indices of nonmissing observations
  int<lower=1> family; // 1 = normal, 2 = student-t
  vector[n_series] alpha; // series intercepts
  vector[3] prior_ar; // prior ar coefficient control parameters
  int<lower=1> P; // AR order
  vector[1] beta; // beta coefficient (1) to use in id_glm functions
}
transformed data {

}
parameters {
  // dynamic factor AR terms
  array[P, K] real<lower=-1, upper=1> ar;

  // dynamic factor process errors
  matrix[n, K] LV_z;

  // factor loading lower triangle
  vector[M] L;

  // observation df for student-t
  real<lower=0> nu;

  // observation variances
  vector<lower=0>[n_series] sigma_obs;

}
transformed parameters {
  // trends and dynamic factor loading matrix
  matrix[n, n_series] trend;
  matrix[n_series, K] Lambda = rep_matrix(0, n_series, K);

  // factor loadings, with constraints
  {
    int index;
    index = 0;
    for (j in 1 : K) {
      for (s in j : n_series) {
        index = index + 1;
        Lambda[s, j] = L[index];
      }
    }
  }

  // derived latent factors
  matrix[n, K] LV  = rep_matrix(0, n, K);
  LV[1, 1 : K] = LV_z[1, 1 : K];
  for (j in 1 : K) {
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

  // prior for observation df parameter
  nu ~ gamma(2, 0.1);

  // priors for AR parameters
  for (j in 1 : K) {
    for (p in 1 : P) {
      if (prior_ar[3] == 1) target += std_normal_lpdf(ar[p, j]);
      else if (prior_ar[3] == 2) target += normal_lpdf(ar[p, j] | prior_ar[1], prior_ar[2]);
    }
  }

  // priors for factor loading coefficients
  L ~ double_exponential(0, 0.5);

  // factor process errors
  to_vector(LV_z) ~ std_normal();

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
  // Posterior predictions
  array[n, n_series] real ypred;
  matrix[n, n_series] sigma_obs_vec;
  for (s in 1 : n_series) {
    sigma_obs_vec[1 : n, s] = rep_vector(sqrt(sigma_obs[s] ^ 2 + sample_sd[s] ^ 2), n);
  }

  if(family == 1) {
    for (s in 1 : n_series) {
      ypred[1 : n, s] = normal_rng(trend[1 : n, s],
                                   sigma_obs_vec[1 : n, s]);
    }
  }

  if(family == 2) {
    vector[n] nu_vec = rep_vector(nu, n);
    for (s in 1 : n_series) {
      ypred[1 : n, s] = student_t_rng(nu_vec, trend[1 : n, s],
                                      sigma_obs_vec[1 : n, s]);
    }
  }
}
