// functions {
//   matrix kronecker_prod(matrix A, matrix B) {
//     int m = rows(A);
//     int n = cols(A);
//     int p = rows(B);
//     int q = cols(B);
//     matrix[m*p,n*q] C;
//     for (i in 1:m)
//       for (j in 1:n)
//         for (k in 1:p)
//           for (l in 1:q)
//             C[p*(i-1)+k,q*(j-1)+l] = A[i,j]*B[k,l];
//     return C;
//   }
// }

data {
  int N;    //number of animals
  int D;    //number of data dimensions
  int P;    //number of (between-animal) predictors
  int V;    //number of variance components
  int K[V]; //number of levels per variance component
  
  vector[D] Y[N];  //dependent variables
  vector[P] X[N];  //predictors
  int Z[V,N];   //levels
  matrix[N,N] L_A; //cholesky-like kinship matrix
}

transformed data {
  cov_matrix[N] A;
  int start[V+1];

  start[1] = 1;
  start[V+1] = sum(K)+1;
  for (v in 2:V) start[v] = start[v-1] + K[v-1];
  A = L_A * L_A';
}

parameters {
  cholesky_factor_corr[D] L_u[V];   //correlation matrix for random effects
  cholesky_factor_corr[D] L_eps;    //correlation matrix for observation noise
  cholesky_factor_corr[D] L_g;      //genetic correlation matrix 
  vector<lower=0>[D] sigma_u[V];    //random effect std 
  vector<lower=0>[D] sigma_eps;     //observation std
  vector<lower=0>[D] sigma_g;       //genetic effect std
  matrix[D,P] beta;                 //regression coefficients
  matrix[D,sum(K)] u_raw;           //raw random effects
  matrix[D,N] g_raw;                //raw genetic effects
  vector[D] alpha;                  //intercept
}

transformed parameters {
  vector[D] eta_Bu[N];
  vector[D] eta[N];
  
  for (i in 1:N) eta_Bu[i] = alpha + beta*X[i];
  for (v in 1:V) {
    matrix[D,K[v]] u;
    u = diag_pre_multiply(sigma_u[v],L_u[v]) * u_raw[:,start[v]:(start[v+1]-1)];
    for (i in 1:N) eta_Bu[i] = eta_Bu[i] + u[:,Z[v,i]];
  }
  
  {  matrix[D,N] g;
    g = diag_pre_multiply(sigma_g,L_g) * g_raw * L_A';
    for (i in 1:N) eta[i] = eta_Bu[i] + g[:,i];
  }
}

model {

  Y ~ multi_normal_cholesky(eta,diag_pre_multiply(sigma_eps,L_eps));
  for (v in 1:V) {
    L_u[v] ~ lkj_corr_cholesky(2.0);
    sigma_u[v] ~ cauchy(0.0,1.0);
  }
  L_eps ~ lkj_corr_cholesky(1.0);
  L_g ~ lkj_corr_cholesky(1.0);
  sigma_eps ~ cauchy(0.0,1.0);
  sigma_g ~ cauchy(0.0,1.0);
  
  to_vector(u_raw) ~ normal(0.0,1.0);
  to_vector(g_raw) ~ normal(0.0,1.0);
  to_vector(beta) ~ student_t(4,0.0,1.0);
  alpha ~ normal(0,1);
}

generated quantities {
  corr_matrix[D] rho_eps;
  corr_matrix[D] rho_u[V];
  corr_matrix[D] rho_g;
  real h2[D];
  // cov_matrix[D*N] Sigma;
  // real log_lik;
  
  rho_eps = L_eps*L_eps';
  for (v in 1:V) rho_u[v] = L_u[v]*L_u[v]';
  rho_g = L_g*L_g';
  
  for (d in 1:D) h2[d] = pow(sigma_g[d],2)/(pow(sigma_g[d],2) + pow(sigma_eps[d],2));
  
  // Sigma = kronecker_prod(A,quad_form_diag(rho_g,sigma_g)) + 
  //   kronecker_prod(diag(rep_vector(1,N)),quad_form_diag(rho_eps,sigma_eps));
  // 
  // for (i in 1:N) {
  //   cov_matrix[D*(N-1)] Sigma_y_y;
  //   matrix[D*(N-1),D] Sigma_y_yrep;
  //   cov_matrix[D,D] Sigma_yrep_yrep;
  //   vector[D*(N-1)] yresid;
  // 
  //   log_lik = 0;
  //   {
  //     int y_ind[D*(N-1)];
  //     int yrep_ind[D];
  //     
  //     for (j in 1:N) {
  //       for (d in 1:D) {
  //         if (i==j) {
  //           yrep_ind[d] = (j-1)*D+d;
  //         } else if (j<i) {
  //           y_ind[(j-1)*D+d] = (j-1)*D+d;
  //           yresid[(j-1)*D+d] = Y[j][d]-eta_Bu[j][d];
  //         } else if (j>i) {
  //           y_ind[(j-2)*D+d] = (j-1)*D+d;
  //           yresid[(j-2)*D+d] = Y[j][d]-eta_Bu[j][d];
  //         }
  //       }
  //     }
  //     
  //     Sigma_yrep_yrep = Sigma[yrep_ind[1]:yrep_ind[D],yrep_ind[1]:yrep_ind[D]];
  //     Sigma_y_yrep = Sigma[y_ind,yrep_ind];
  //     Sigma_y_y = Sigma[y_ind,y_ind];
  //     Sigma_y_y = inverse_spd(Sigma_y_y);
  //     
  //     vector[D] mu_yrep_c_y = pre_eta[i] + Sigma_y_yrep' * Sigma_y_y * yresid;
  //     cov_matrix[D] Sigma_yrep_c_y = Sigma_yrep_yrep - quad_form_sym(Sigma_y_y,Sigma_y_yrep);
  //     
  //     log_lik = log_lik + multi_normal_lpdf(Y[i],mu_yrep_c_y,Sigma_yrep_c_y);
  //   }
  // }

}
