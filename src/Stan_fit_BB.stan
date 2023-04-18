data {
  // Dimensions
  int<lower=2> num_days;   //
  int<lower=1> num_groups; //

  // Data
  real t[num_days];     // 
  int pos[num_days];   // 
  int tot[num_days];   // 
  int day_to_group[num_days];

  // Dispersion prior hyperparameters
  real<lower=0> m_rho;
  real<lower=0> sig_rho;

  // Day-of-the-week effect prior hyper-parameters
  real<lower=0> a_w;
  real<lower=0> b_w;

  // GP prior hyper-parameters
  vector[2] log_theta_0; // [log_range_0, log_sigma_0]
  cov_matrix[2] B;
}

transformed data {
  matrix[num_groups, num_groups] A = diag_matrix(rep_vector(1, num_groups));
  matrix[num_groups, num_groups-1] A_qr;
  for (i in 1: num_groups-1) A[num_groups,i] = -1;
  A[num_groups, num_groups] = 0;
  A_qr = qr_Q(A)[ , 1:(num_groups-1)];
}

parameters {
  // Dispersion
  real<lower=0, upper=1> rho;

  // GP
  vector[num_days] x_t;
  vector[2] log_theta_x; // [log range, log sigma]

  // Day-of-the-week effect
  vector[num_groups - 1] wraw_d; // vector[num_groups] w_d;
  real<lower=0> tau_w;
}

transformed parameters {
  real<lower=0> sig_w = pow(sqrt(tau_w), -1);
  vector[num_groups] w_d = sig_w * (A_qr * wraw_d);

  real<lower=0> range_x = exp(log_theta_x[1]);
  real<lower=0> sig2_x = square(exp(log_theta_x[2]));
  real<lower=0> tau_x = pow(sig2_x, -1);
}

model{
  // Dispersion
  logit(rho) ~ normal(m_rho, sig_rho);

  // Day-of-the-week effect (Stan p20)
  tau_w ~ gamma(a_w, b_w);
  wraw_d ~ normal(0, inv(sqrt(1 - inv(num_groups)))); // w_d ~ normal(0, sig_w);

  // GP
  matrix[num_days, num_days] L_K;
  matrix[num_days, num_days] K;
  for (i in 1:(num_days - 1)) {
    K[i, i] = sig2_x + 0.0001;
    for (j in (i + 1):num_days) {
      K[i, j] = sig2_x*(1 + sqrt(12)*abs(t[i] - t[j])/range_x)*exp(-sqrt(12)*abs(t[i] - t[j])/range_x);
      K[j, i] = K[i, j];
    }
  }
  K[num_days, num_days] = sig2_x + 0.0001;
  L_K = cholesky_decompose(K);

  vector[num_days] mu_x = rep_vector(0, num_days);
  log_theta_x ~ multi_normal(log_theta_0, B);
  x_t ~ multi_normal_cholesky(mu_x, L_K);

  // Observations
  // See p.272
  vector[num_days] w_d_long;
  vector[num_days] alpha;
  vector[num_days] beta;
  for (n in 1:num_days) {
    w_d_long[n] = w_d[day_to_group[n]];
    alpha[n] = inv_logit(x_t[n] + w_d_long[n])*(1/rho - 1);
    beta[n] = (1 - inv_logit(x_t[n] + w_d_long[n]))*(1/rho - 1);
  }

  pos ~ beta_binomial(tot, alpha, beta);
}
