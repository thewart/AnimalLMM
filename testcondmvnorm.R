N <- 3
D <- 2

A <- MCMCpack::riwish(N+1,diag(N))
A <- diag(1/sqrt(diag(A))) %*% A %*% diag(1/sqrt(diag(A)))

rho <- MCMCpack::riwish(D+1,diag(D))
rho <- diag(1/sqrt(diag(rho))) %*% rho %*% diag(1/sqrt(diag(rho)))

Sigma <- kronecker(A,rho)

y_ind <- vector("numeric",length=D*(N-1));
yrep_ind <- vector("numeric",length=D);
i <- 1
for (j in 1:N) {
  for (d in 1:D) {
    if (i==j) {
      yrep_ind[d] = (j-1)*D+d;
    } else if (j<i) {
      y_ind[(j-1)*D+d] = (j-1)*D+d;
    } else if (j>i) {
      y_ind[(j-2)*D+d] = (j-1)*D+d;
    }
  }
}

Sigma_yrep_yrep = Sigma[yrep_ind[1]:yrep_ind[D],yrep_ind[1]:yrep_ind[D]];
Sigma_y_yrep = Sigma[y_ind,yrep_ind];
Sigma_y_y = Sigma[y_ind,y_ind];

library(mvtnorm)
y <- rmvnorm(1,rep(0,nrow(Sigma)),Sigma)
dmvnorm(y,mean = rep(0,nrow(Sigma)),sigma = Sigma,log = T)

foo <- t(Sigma_y_yrep) %*% solve(Sigma_y_y) %*% y[3:6]
moo <- Sigma_yrep_yrep - t(Sigma_y_yrep) %*% solve(Sigma_y_y) %*% Sigma_y_yrep
dmvnorm(y[1:2],foo,moo,log=T) + dmvnorm(y[3:6],rep(0,4),sigma = Sigma[3:6,3:6],log = T)
