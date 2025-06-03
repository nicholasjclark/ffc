// Stan model code for Vector autoregressive dynamic factor forecasting models
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

  /* Function to compute the matrix square root */
  /* see Heaps 2022 for details (https://doi.org/10.1080/10618600.2022.2079648)*/
  matrix sqrtm(matrix A) {
    int m = rows(A);
    vector[m] root_root_evals = sqrt(sqrt(eigenvalues_sym(A)));
    matrix[m, m] evecs = eigenvectors_sym(A);
    matrix[m, m] eprod = diag_post_multiply(evecs, root_root_evals);
    return tcrossprod(eprod);
  }

  /* Function to transform P_real to P */
  /* see Heaps 2022 for details (https://doi.org/10.1080/10618600.2022.2079648)*/
  matrix P_realtoP(matrix P_real) {
    int m = rows(P_real);
    matrix[m, m] B = tcrossprod(P_real);
    for (i in 1 : m) {
      B[i, i] += 1.0;
    }
    return mdivide_left_spd(sqrtm(B), P_real);
  }

  /* Function to perform the reverse mapping*/
  /* see Heaps 2022 for details (https://doi.org/10.1080/10618600.2022.2079648)*/
  array[,] matrix rev_mapping(array[] matrix P, matrix Sigma) {
    int p = size(P);
    int m = rows(Sigma);
    array[p, p] matrix[m, m] phi_for;
    array[p, p] matrix[m, m] phi_rev;
    array[p + 1] matrix[m, m] Sigma_for;
    array[p + 1] matrix[m, m] Sigma_rev;
    matrix[m, m] S_for;
    matrix[m, m] S_rev;
    array[p + 1] matrix[m, m] S_for_list;
    array[p + 1] matrix[m, m] Gamma_trans;
    array[2, p] matrix[m, m] phiGamma;

    // Step 1:
    Sigma_for[p + 1] = Sigma;
    S_for_list[p + 1] = sqrtm(Sigma);
    for (s in 1 : p) {
      // In this block of code S_rev is B^{-1} and S_for is a working matrix
      S_for = -tcrossprod(P[p - s + 1]);
      for (i in 1 : m) {
        S_for[i, i] += 1.0;
      }
      S_rev = sqrtm(S_for);
      S_for_list[p - s + 1] = mdivide_right_spd(mdivide_left_spd(S_rev,
                                                                 sqrtm(
                                                                 quad_form_sym(
                                                                 Sigma_for[
                                                                 p - s + 2],
                                                                 S_rev))),
                                                S_rev);
      Sigma_for[p - s + 1] = tcrossprod(S_for_list[p - s + 1]);
    }

    // Step 2:
    Sigma_rev[1] = Sigma_for[1];
    Gamma_trans[1] = Sigma_for[1];
    for (s in 0 : (p - 1)) {
      S_for = S_for_list[s + 1];
      S_rev = sqrtm(Sigma_rev[s + 1]);
      phi_for[s + 1, s + 1] = mdivide_right_spd(S_for * P[s + 1], S_rev);
      phi_rev[s + 1, s + 1] = mdivide_right_spd(S_rev * P[s + 1]', S_for);
      Gamma_trans[s + 2] = phi_for[s + 1, s + 1] * Sigma_rev[s + 1];
      if (s >= 1) {
        for (k in 1 : s) {
          phi_for[s + 1, k] = phi_for[s, k]
                              - phi_for[s + 1, s + 1] * phi_rev[s, s - k + 1];
          phi_rev[s + 1, k] = phi_rev[s, k]
                              - phi_rev[s + 1, s + 1] * phi_for[s, s - k + 1];
        }
        for (k in 1 : s) {
          Gamma_trans[s + 2] = Gamma_trans[s + 2]
                               + phi_for[s, k] * Gamma_trans[s + 2 - k];
        }
      }
      Sigma_rev[s + 2] = Sigma_rev[s + 1]
                         - quad_form_sym(Sigma_for[s + 1],
                                         phi_rev[s + 1, s + 1]');
    }
    for (i in 1 : p) {
      phiGamma[1, i] = phi_for[p, i];
    }
    for (i in 1 : p) {
      phiGamma[2, i] = Gamma_trans[i]';
    }
    return phiGamma;
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
  int<lower=1> P; // AR order
  vector[1] beta; // beta coefficient (1) to use in id_glm functions
}
transformed data {
  # Sampling SDs for vectorized likelihood calculations
  vector[n_nonmissing] flat_sigma_obs = rep_each(sample_sd, n)[obs_ind];

  matrix[n, n_series] sigma_obs_vec;
  for (s in 1 : n_series) {
    sigma_obs_vec[1 : n, s] = rep_vector(sample_sd[s], n);
  }

  # Error variances (fixed to 1)
  vector[K] sigma = rep_vector(1.0, K);

  # Zero vector
  vector[K] zero_vec = rep_vector(0.0, K);
}
parameters {
  // unconstrained VAR partial autocorrelations
  array[P] matrix[K, K] P_real;

  // partial autocorrelation hyperparameters
  vector[2] Pmu;
  vector<lower=0>[2] Pomega;

  // factor error correlations
  cholesky_factor_corr[K] L_Omega;

  // raw factors
  array[n] vector[K] LV_raw;

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

  // constrained VAR matrices
  array[P] matrix[K, K] A;

  // LKJ form of factor covariance matrix
  matrix[K, K] L_Sigma = diag_pre_multiply(sigma, L_Omega);

  // computed error covariance matrix
  cov_matrix[K] Sigma = multiply_lower_tri_self_transpose(L_Sigma);

  // initial factor covariance
  cov_matrix[P * K] Gamma;

  // factor estimates in matrix-form
  matrix[n, K] LV;
  for (i in 1 : n) {
    LV[i, 1 : K] = to_row_vector(LV_raw[i]);
  }

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

  // derived stationary VAR latent factors
  {
    array[P] matrix[K, K] Ptilde;
    array[2, P] matrix[K, K] phiGamma;
    for (i in 1 : P) Ptilde[i] = P_realtoP(P_real[i]);
    phiGamma = rev_mapping(Ptilde, Sigma);
    A = phiGamma[1];
    for(i in 1 : P) {
      for(j in 1 : P) {
        if(i <= j) Gamma[((i-1)*K+1):(i*K), ((j-1)*K+1):(j*K)] = phiGamma[2, j-i+1];
        else Gamma[((i-1)*K+1):(i*K), ((j-1)*K+1):(j*K)] = phiGamma[2, i-j+1]';
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
  // factor mean parameters
  array[n] vector[K] mu;

  // priors for series-level intercepts
  target += normal_lpdf(alpha | prior_alpha[1], prior_alpha[2]);

  // prior for observation df parameter
  nu ~ gamma(2, 0.1);

  // LKJ error correlation prior
  L_Omega ~ lkj_corr_cholesky(2);

  // unconstrained partial autocorrelations
  Pmu ~ normal(0.0, 0.675);
  Pomega ~ gamma(1.365, 0.071);

  for (s in 1 : P) {
    diagonal(P_real[s]) ~ normal(Pmu[1], 1 / sqrt(Pomega[1]));
     for(i in 1 : K) {
       for(j in 1 : K) {
        if(i != j) P_real[s, i, j] ~ normal(Pmu[2], 1 / sqrt(Pomega[2]));
      }
    }
  }

  // factor means
  for (t in 1 : n) mu[t] = zero_vec;
  for (p in 1 : P) {
    for (t in (P + 1) : n){
      mu[t] += A[p] * LV_raw[t - p];
    }
  }

  // stochastic factors (contemporaneously correlated)
  for (t in 1 : P) {
    LV_raw[t] ~ multi_normal(zero_vec, Gamma[((t-1)*K+1):(t*K), ((t-1)*K+1):(t*K)]);
  }
  for (t in (P + 1) : n) {
    LV_raw[t] ~ multi_normal_cholesky(mu[t], L_Sigma);
  }

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
