
rm(list=ls())
# install.packages("lars")
library(MASS) # for mvrnorm function
library(glmnet)
# read source files
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # setting the workspace to source file location


####################################################
# 2. Data generation given at Peterson et al. (2016)
####################################################

# (a) first setting
n = 100
p = 150




####################################################
# generation of covariate matrix X: first and second settings 
set.seed(12)
Omega0 = matrix(0, p,p)


for(j in 1:(p-1)){
  Omega0[j,j+1] = Omega0[j+1,j] = 0.3
}


diag(Omega0) = 1
library(MASS)
X = mvrnorm(n,rep(0, p),Sigma = solve(Omega0))
G0 = 1*(Omega0 != 0)
diag(G0) = 0




##################################################
# generation of coefficient vector \beta
beta0 = matrix(0, nrow=p, ncol=1)

##############################################################

t = 10

beta0[,1] = c(rep(3,t), rep(0,p-t))


gam0.loc = 1:t
gam0.len = t
gamma0 = c(rep(1,gam0.len), rep(0, p-gam0.len))


####################################################
# Generate a data vector Y : n X 1
####################################################

set.seed(12)
Y = X %*% beta0 + rnorm(n = n, mean = rep(0,n), sd = 1)
Y = rnorm(n = n, mean = X %*% beta0, sd = 1)
Y = ifelse(Y<=0, 0, 1)



####################################################
# hyperparameters
####################################################
a.hp = 2.75 #2.75 tau 2 beta all 3 gamma 6 iter 200 200
b.hp = 0.5
qn = 0.005

tau2 = 1


r = 1e-4
s = 1e-8



# niter = 10000
# nburn = .2*niter
niter = 2000
nburn = 2000
nadap = 0

init.gam = rep(0, p)

sd = 12

set.seed(sd)
time = Sys.time()
source("joint_selection_new.R")
res = BSSC(Y, X, qn, init.A=NULL, init.gam=NULL, Rj=NULL, a.hp, b.hp, tau2, r, s, niter, nburn)
elapsed = Sys.time() - time
elapsed

# variable selection result
post.ind = as.numeric(colMeans(res$gamma.mat[2000:4000,])>0.5)
summary.joint(gamma0, post.ind)

# inverse covariance estimation result
post.Omega = apply(res$Omega.mat[,,-(1:nburn)], c(1, 2), mean, na.rm = TRUE)
post.Omega[which(abs(post.Omega) < 10^(-4))] = 0
EvaluationNorm(Omega0, post.Omega)

# inverse graph selection result
post.G = apply(res$G.mat[,,-(1:nburn)], c(1, 2), mean, na.rm = TRUE)
post.G[which(post.G > 0.5)] = 1
post.G[which(post.G < 0.5)] = 0
post.G[upper.tri(post.G)] = 0
Evaluation.G(G0, post.G)


