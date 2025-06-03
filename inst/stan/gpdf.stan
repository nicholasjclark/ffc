// Stan model code for GP dynamic factor forecasting models
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
  /* Spectral density GP eigenvalues*/
  /* see Riutort-Mayol et al 2023 for details (https://doi.org/10.1007/s11222-022-10167-2)*/
  real lambda_gp(real L, int m) {
    real lam;
    lam = ((m * pi()) / (2 * L)) ^ 2;
    return lam;
  }
  /* Spectral density GP eigenfunctions*/
  /* see Riutort-Mayol et al 2023 for details (https://doi.org/10.1007/s11222-022-10167-2)*/
  vector phi_SE(real L, int m, vector x) {
    vector[rows(x)] fi;
    fi = 1 / sqrt(L) * sin(m * pi() / (2 * L) * (x + L));
    return fi;
  }
  /* Spectral density squared exponential Gaussian Process*/
  /* see Riutort-Mayol et al 2023 for details (https://doi.org/10.1007/s11222-022-10167-2)*/
  real spd_SE(real alpha, real rho, real w) {
    real S;
    S = (alpha ^ 2) * sqrt(2 * pi()) * rho * exp(-0.5 * (rho ^ 2) * (w ^ 2));
    return S;
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
  vector[2] prior_alpha; // prior intercept control parameters
  vector[1] beta; // beta coefficient (1) to use in id_glm functions
}
transformed data {
  // gp phi
  vector<lower=1>[n] times;
  vector[n] times_cent;
  real mean_times;
  real<lower=0> boundary;
  int<lower=1> num_gp_basis;
  num_gp_basis = min(25, n);
  matrix[n, num_gp_basis] gp_phi;
  for (t in 1 : n) {
    times[t] = t;
  }
  mean_times = mean(times);
  times_cent = times - mean_times;
  boundary = (5.0 / 4) * (max(times_cent) - min(times_cent));
  for (m in 1 : num_gp_basis) {
    gp_phi[ : , m] = phi_SE(boundary, m, times_cent);
  }

  # Sampling SDs for vectorized likelihood calculations
  vector[n_nonmissing] flat_sigma_obs = rep_each(sample_sd, n)[obs_ind];

  matrix[n, n_series] sigma_obs_vec;
  for (s in 1 : n_series) {
    sigma_obs_vec[1 : n, s] = rep_vector(sample_sd[s], n);
  }
}
parameters {
  // gp length scale parameters
  vector<lower=0>[K] rho_gp;

  // gp coefficient weights
  matrix[num_gp_basis, K] b_gp;

  // factor loading lower triangle
  vector[M] L;

  // series-level intercepts
  vector[n_series] alpha;

  // observation df for student-t
  real<lower=0> nu;

}
transformed parameters {
  // trends and dynamic factor loading matrix
  matrix[n, n_series] trend;
  matrix[n_series, K] Lambda = rep_matrix(0, n_series, K);

  // gp spectral densities
  matrix[n, K] LV;
  matrix[num_gp_basis, K] diag_SPD;
  matrix[num_gp_basis, K] SPD_beta;

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
  for (m in 1 : num_gp_basis) {
    for (j in 1 : K) {
      diag_SPD[m, j] = sqrt(spd_SE(1.0, rho_gp[j],
                            sqrt(lambda_gp(boundary, m))));
    }
  }
  SPD_beta = diag_SPD .* b_gp;
  LV = gp_phi * SPD_beta;

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

  // priors for gp parameters
  to_vector(b_gp) ~ std_normal();
  rho_gp ~ inv_gamma(1.499007, 5.670433);

  // priors for factor loading coefficients
  L ~ student_t(3, 0, 1.5);

  {
    // likelihood functions
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

  if(family == 1) {
    for (s in 1 : n_series) {
      ypred[1 : n, s] = normal_rng(trend[1 : n, s], sigma_obs_vec[1 : n, s]);
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
