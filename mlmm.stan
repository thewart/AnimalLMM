data {
  int N;    //number of animals
  int D;    //number of data dimensions
  int P;    //number of (between-animal) predictors
  int V;    //number of variance components
  int K[V]; //number of levels per variance component
  
  vector[D] Y[N];  //dependent variables
  vector[P] X[N];  //predictors
  int Z[V,N];   //levels
}

transformed data {
  int start[V+1];
  
  start[1] = 1;
  start[V+1] = sum(K)+1;
  for (v in 2:V) start[v] = start[v-1] + K[v-1];
}

parameters {
  cholesky_factor_corr[D] L_u[V];   //correlation matrix for random effects
  cholesky_factor_corr[D] L_eps;    //correlation matrix for observation noise
  vector<lower=0>[D] sigma_u[V];     //random effect std 
  vector<lower=0>[D] sigma_eps;          //observation std
  matrix[D,P] beta;                //regression coefficients
  matrix[D,sum(K)] u_raw;          //raw random effects
  vector[D] alpha;                  //intercept
}

transformed parameters {
  vector[D] eta[N];
  
  for (i in 1:N) eta[i] = alpha + beta*X[i];
  for (v in 1:V) {
    matrix[D,K[v]] u;
    u = diag_pre_multiply(sigma_u[v],L_u[v]) * u_raw[,start[v]:(start[v+1]-1)];
    for (i in 1:N) eta[i] = eta[i] + u[:,Z[v,i]];
  }
}

model {
  real u_tmp[D,sum(K)];
  real beta_tmp[D,P];
  
  Y ~ multi_normal_cholesky(eta,diag_pre_multiply(sigma_eps,L_eps));
  #Y ~ multi_normal_cholesky(eta,diag_matrix(sigma_eps.*sigma_eps));
  for (v in 1:V) {
    L_u[v] ~ lkj_corr_cholesky(2.0);
    sigma_u[v] ~ cauchy(0.0,1.0);
  }
  sigma_eps ~ cauchy(0.0,1.0);
  
  to_vector(u_raw) ~ normal(0.0,1.0);
  to_vector(beta) ~ student_t(4,0.0,1.0);
  alpha ~ normal(0,1);
}

generated quantities {
  corr_matrix[D] rho_eps;
  corr_matrix[D] rho_u[V];
  
  rho_eps = L_eps*L_eps';
  for (v in 1:V) rho_u[v] = L_u[v] * L_u[v]';
}
