data {
  int N;    //number of animals
  int P;    //number of (between-animal) predictors
  #int V;    //number of variance components
  #int K[V]; //number of levels per variance component
  int R;    //number of biological replicates
  
  vector[N] Y;  //dependent variables
  matrix[N,P] X;  //predictors
  #int Z[V,N];   //levels
  matrix[N,R] L_A; //cholesky-like kinship matrix
}

transformed data {
  // int start[V+1];
  // 
  // start[1] = 1;
  // start[V+1] = sum(K)+1;
  // for (v in 2:V) start[v] = start[v-1] + K[v-1];
}

parameters {
  #real<lower=0> sigma_u[V];    //random effect std 
  real<lower=0> sigma_eps;     //observation std
  real<lower=0> sigma_g;       //genetic effect std
  vector[P] beta;                 //regression coefficients
  #vector[sum(K)] u_raw;           //raw random effects
  vector[R] g_raw;                //raw genetic effects
  real alpha;                  //intercept
}

transformed parameters {
  vector[N] eta;
  
  eta = alpha + X*beta + sigma_g * L_A * g_raw;

  // for (v in 1:V) {
  //   vector[K[v]] u;
  //   u = sigma_u[v] * u_raw[start[v]:(start[v+1]-1)];
  //   for (i in 1:N) eta[i] = eta[i] + u[Z[v,i]];
  // }
}

model {

  Y ~ normal(eta,sigma_eps);
  // for (v in 1:V) {
  //   sigma_u[v] ~ cauchy(0.0,1.0);
  // }
  sigma_eps ~ cauchy(0.0,1.0);
  sigma_g ~ cauchy(0.0,1.0);
  
  #u_raw ~ normal(0.0,1.0);
  g_raw ~ normal(0.0,1.0);
  beta ~ student_t(4,0.0,1.0);
  alpha ~ normal(0,1);
}

generated quantities {
  real h2;
  h2 = pow(sigma_g,2)/(pow(sigma_g,2) + pow(sigma_eps,2));
}
