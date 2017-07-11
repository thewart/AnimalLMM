library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

geno <- "none"
source("~/code/multivarlmm/morphdatprep.R")
mvlmm <- stan_model("~/code/multivarlmm/mlmm.stan")
stanfit <- sampling(mvlmm,data=standat,chains=4,iter=1000,
                    pars=c("alpha","beta","sigma_eps","sigma_u","rho_eps","rho_u"))


geno <- "ped"
source("~/code/multivarlmm/morphdatprep.R")
mvlmm <- stan_model("~/code/multivarlmm/mlmm_animal.stan")
load("~/Dropbox/Pedigree and Life-History Data/pedigreeKW2016.RData")
A <- 2*kinship2::kinship(bigped$id,dadid=bigped$sire,momid=bigped$dam)
A <- A[match(morphdat$ID,rownames(A)),match(morphdat$ID,rownames(A))]
standat$L_A <- t(chol(A))
stanfit <- sampling(mvlmm,data=standat,chains=,iter=500,warmup=250,
                    pars=c("alpha","beta","sigma_eps","sigma_u","sigma_g","rho_eps","rho_u","rho_g","h2"))

geno <- "snp"
#source('~/code/snparray/processraw.R')
#G <- as.data.table(cdat$X)
G <- fread("~/analysis/SNPannotation/SNPmaster_qc.csv")
library(pcaMethods)
source("~/code/snparray/SNPpreproc.R")
source("~/code/multivarlmm/morphdatprep.R")
G <- G[ID %in% morphdat$ID]
setkey(G,"ID")
G <- redund(G)
Gpca <- pca(t(as.matrix(G[,-1])),method = "ppca",nPcs = 135,scale = "uv",center = T)# %>% completeObs()
standat$L_A <- loadings(Gpca) %*% diag(sDev(Gpca))
standat$L_A <- cbind(L_A,matrix(0,nrow(L_A),nrow(L_A)-ncol(L_A)))
mvlmm <- stan_model("~/code/multivarlmm/mlmm_animal.stan")
stanfit_G <- sampling(mvlmm,data=standat,chains=4,iter=2000,warmup=500,
                    pars=c("alpha","beta","sigma_eps","sigma_u","sigma_g","rho_eps","rho_u","rho_g"))

A <- A[match(morphdat$ID,rownames(A)),match(morphdat$ID,rownames(A))]
standat$L_A <- t(chol(A))
stanfit_A <- sampling(mvlmm,data=standat,chains=4,iter=2000,warmup=500,
                      pars=c("alpha","beta","sigma_eps","sigma_u","sigma_g","rho_eps","rho_u","rho_g"))


guh <- stan_hist(stanfit,pars="h2")$data
names(guh)[1] <- "heritability"
ggplot(guh,aes(x=heritability)) + geom_density() + facet_wrap(~parameter,scale="free") + geom_vline(xintercept=0,size=0) 
