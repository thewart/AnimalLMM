library(rstanarm)
library(rstan)
options(mc.cores = parallel::detectCores())
P <- 4
N <- 500
D <- 2
K <- c(5,10,20)
X <- rnorm(P*N) %>% matrix(nrow=P)
sigma_z <- matrix(c(0.5,0.1,0.2,1.0,0.33,0.33),2)
sigma <- c(0.5,1.0)

L <- matrix(c(1.0,0.5,0.5,1.0),2) %>% chol() %>% t()
L_Z <- list(matrix(c(1.0,0.0,0.0,1.0),2) %>% chol() %>% t(),
            matrix(c(1.0,-0.7,-0.7,1.0),2) %>% chol() %>% t(),
            matrix(c(1.0,0.9,0.9,1.0),2) %>% chol() %>% t())
u <- lapply(1:length(K),function(x) diag(sigma_z[,x]) %*% L_Z[[x]] %*% matrix(rnorm(2*K[x]),2))

Z <- lapply(K,function(x) sample(1:x,N,replace=T)) 

beta <- rnorm(D*P) %>% matrix(nrow=D)
Y <- beta %*% X + diag(sigma) %*% L %*% matrix(rnorm(D*N),nrow=D)
for (k in 1:length(K)) Y <- Y + u[[k]][,Z[[k]]]

dat <- cbind(data.table(y=Y[2,]),t(X),as.data.table(Z))
setnames(dat,c("y","x1","x2","x3","x4","l1","l2","l3"))
dat[,y:=(y-mean(y))/sd(y)]
fit <- stan_lmer(y~x1+x2+x3+x4+(1|l1) + (1|l2) + (1|l3),dat,iter=1000)

standat <- list(N=N,D=D,P=P,K=K,V=length(K),
                Y=t(Y),X=t(X),Z=matrix(unlist(Z),nrow=length(K),byrow = T))

mvlmm <- stan_model("~/code/multivarlmm/mlmm.stan")
stanfit <- sampling(mvlmm,data=standat,chains=4,iter=1000,
                    pars=c("alpha","beta","sigma_eps","rho_eps","sigma_u","rho_u"))

#heritability
sigma_g <- c(1.0,0.5)
L_g <- matrix(c(1.0,0.75,0.75,1.0),2) %>% chol() %>% t()
A <- MCMCpack::riwish(N+1,diag(N))
A <- diag(1/sqrt(diag(A))) %*% A %*% diag(1/sqrt(diag(A)))
L_A <- chol(A) %>% t()
A <- L_A%*%t(L_A)
Y <- Y + diag(sigma_g) %*% L_g %*% matrix(rnorm(D*N),nrow=D) %*% t(L_A)

dat <- cbind(data.table(y=Y[2,]),t(X),as.data.table(Z))
setnames(dat,c("y","x1","x2","x3","x4","l1","l2","l3"))
dat$ID <- 1:N
rownames(A) <- colnames(A) <- 1:N
source("~/Dropbox/tools/cheetahmm/cheetahmm.R")
cheetah.mm(y ~ x1+x2+x3+x4 + (1|l1) + (1|l2) + (1|l3) + (1|ID),data=dat,pedigree = list(ID=A))

standat$Y <- t(Y)
standat$L_A <- t(chol(A))

mvlmm <- stan_model("~/code/multivarlmm/mlmm_animal.stan")
stanfit <- sampling(mvlmm,data=standat,chains=1,iter=200,
                    pars=c("alpha","beta","sigma_eps","sigma_u","sigma_g","rho_eps","rho_u","rho_g"))
