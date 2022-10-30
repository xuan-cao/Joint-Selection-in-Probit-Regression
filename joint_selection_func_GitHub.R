library(lars)
library(truncnorm)
# Y : n X 1 data 
# X : n X p covariate matrix
################################# need to sample beta!!!
BSSC <- function(Y, X, qn, init.A=NULL, init.gam=NULL, Rj=NULL, a.hp, b.hp, tau2, r, s, niter, nburn){
   
   n = nrow(X)
   p = ncol(X)
   S = 1/n*t(X)%*%X
   
   
   # set Rj value
   if(is.null(Rj)) Rj = floor( n/(log(n, base=10)) )
   dj = NULL
   
   res = list()
   
   G.mat = array(0, dim = c(p, p, niter + nburn))
   Omega.mat = array(0, dim = c(p, p, niter + nburn))

   
   # set initial value for Omega and G
   library(glasso)
   Vs = var(X)
   a = glasso(Vs, rho=.055)
   aa = glasso(Vs,rho=.1, w.init=a$w, wi.init=a$wi)
   glassoaa = aa$wi
   Omega.mat[,,1] = glassoaa
   diag(glassoaa) = 0
   G.mat[,,1] = 1*(abs(glassoaa) > 0)

   # set initial value for gamma
   gamma.mat = matrix(0, nrow=niter + nburn, ncol=p)
   gamma.mat[1,] = rep(0, p)
   beta = rep(0,p)
   
   # set initial value for Z
   Z = rep(1,n)*ifelse(Y == 1, 1, -1)
   
   # MCMC sampling
   for(i in 2:(niter+nburn)){
      

      for(j in 1:n) {
      # updating Z
      temp = t(X[j,]) %*% beta
      if(Y[j] > 0) Z[j] = rtruncnorm(1, a=0, mean = temp, sd = 1)
      else Z[j] = rtruncnorm(1, b=0, mean = temp, sd = 1)
      }
   
      gamma = gamma.mat[i-1,]
      
      for (g in 1:p) {
      # update gamma
      betag = beta[g]
      Xg = X[, g]
      beta_g = beta[-g]
      X_g = X[, -g]
      invSigmag = t(Xg)%*%Xg + 1/tau2
      Sigmag = 1/invSigmag

      mug = Sigmag*t(Xg)%*%(Z - X_g%*%beta_g)
      temp = 1 + sqrt(Sigmag/tau2)*exp(0.5*mug^2*invSigmag - a.hp + 2*b.hp*t(gamma)%*%G.mat[,,i-1][,g])
      temp = 1 - 1/temp
      gamma[g] = rbinom(n = 1, size = 1, prob = temp)
      
      
      # update beta
      if(gamma[g] == 0) beta[g] = 0
      else
         beta[g] = rnorm(n = 1, mean = mug, sd = sqrt(Sigmag))
      }
      
      gamma.mat[i,] = gamma
   
      
      # Gibbs for G and Omega
      Omega = Omega.mat[,,i-1]
      G = matrix(0,p,p)
      for(j in 1:(p-1)){
         for(k in (j+1):p){
            lambda = rgamma(1,shape = r + 0.5, rate = 0.5*Omega[j,k]*Omega[j,k] + s)
            #print(lambda)
            a = S[j,j] + S[k,k] + lambda/n
            b = t(Omega[j,])%*%S[k,] + t(Omega[,k])%*%S[j,] - Omega[j,k]*(S[j,j] + S[k,k])
            logc = log(qn) - log(1 - qn) + 0.5*log(lambda/(n*a)) + n*b^2/(2*a) + 2*b.hp*gamma[j]*gamma[k]
            
      

            # sample G_jk
            if(logc < 10) G[j,k] = G[k,j] = rbinom(n = 1, size = 1, prob = 1 - 1/(1+exp(logc)))
            else G[j,k] = G[k,j] = rbinom(n = 1, size = 1, prob = 1 - 1/(1+exp(10)))
            
            # sample Omega_jk
            if(G[j,k] == 0) Omega[j,k] = Omega[k,j] = 0
            else Omega[j,k] = Omega[k,j] = rnorm(n = 1, mean = -b/a, sd = sqrt(1/(n*a)))
         }
         bj = t(Omega[j,])%*%S[j,] - Omega[j,j]*S[j,j]
         Omega[j,j] = (sqrt(bj^2 + 4*S[j,j]) - bj)/(2*S[j,j])
      }
      bp = t(Omega[p,])%*%S[p,] - Omega[p,p]*S[p,p]
      Omega[p,p] = (sqrt(bp^2 + 4*S[p,p]) - bp)/(2*S[p,p])
      
      G.mat[,,i] = G
      Omega.mat[,,i] = Omega
      
      
      if(i %% 100 == 0)
      cat(i,"th iteration is completed. . . . . .\n")
      }# end of (i in 2:(niter+nburn)) for loop
   
   
   res = list(G.mat = G.mat, Omega.mat = Omega.mat, gamma.mat = gamma.mat)
   
   return(res)
   
   
   }
   
  
      


###################################################################
# Auxiliary functions
###################################################################



summary.joint <- function(gamma0, gam.est){
  
  TP = length(intersect(which(gamma0 != 0), which(gam.est >= 0.5)))
  TN = p - sum(gamma0!=0) - length(intersect(which(gamma0 == 0), which(gam.est >= 0.5)))
  FP = length(intersect(which(gamma0 == 0), which(gam.est >= 0.5)))
  FN = sum(gamma0!=0) - length(intersect(which(gamma0 != 0), which(gam.est >= 0.5)))
  
  Sen = TP/(TP+FN)
  Spe = TN/(TN+FP)
  MCC = (TP*TN - FP*FN)/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) )
  # Sen; Spe; MCC
  char = paste(c("Number of errors:", FP+FN, ",   Sensitivity:", round(Sen, digits = 4), ",   Specificity:", round(Spe, digits = 4), ",   MCC:",round(MCC, digits = 4)), collapse = "")
  plot(gam.est>0.5, main=char, ylab="gamma_i")
  lines(gamma0, col=2)
  
  num.pos = sum(gam.est>0.5)
  loc.pos = which(gam.est>0.5)
  
  cat("The number of positive: " , num.pos, "\n")
  cat("The locations of positive: " , loc.pos, "\n")
  
  res = list(TP=TP, TN=TN, FP=FP, FN=FN, Sen=Sen, Spe=Spe, MCC=MCC)
  return(res)
}

Joint_auc <- function(gamma0, inc.prob){
   
   thres = inc.prob[order(inc.prob)] # thresholds for ROC curve
   thres[-1] = thres[-1] + abs(min(diff(thres)))/2
   thres[1] = thres[1] - 0.001
   
   TPR = rep(0, length(thres))
   FPR = rep(0, length(thres))
   for(i in 1:length(thres)){
      gamma.est = (inc.prob > thres[i])
      TPR[i] = sum(gamma0*gamma.est)/sum(gamma0)
      FPR[i] = sum(gamma.est-gamma0 == 1)/sum(1-gamma0)
   }
   
   TPR = TPR[order(TPR)]
   FPR = FPR[order(FPR)]
   
   # inputs already sorted, best scores first 
   dFPR <- c(diff(FPR), 0)
   dTPR <- c(diff(TPR), 0)
   sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}



EvaluationNorm <- function(Omega1, Omega2){
   diff = Omega1 - Omega2
   Onenorm = norm(diff, type = "O")/norm(Omega1, type = "O")
   Fnorm = norm(diff, type = "F")/norm(Omega1, type = "F")
   Specnorm = norm(diff, type = "2")/norm(Omega1, type = "2")
   infnorm = max(abs(diff))/max(abs((Omega1)))
   return(list(Onenorm=Onenorm,Specnorm=Specnorm,Fnorm=Fnorm,infnorm=infnorm))
}



Evaluation.G <- function(Adj1, Adj2){
   true.index <- which(Adj1==1)
   false.index <- setdiff(which(lower.tri(Adj1)),true.index)
   positive.index <- which(Adj2==1)
   negative.index <- setdiff(which(lower.tri(Adj2)),positive.index)
   
   TP <- length(intersect(true.index,positive.index))
   FP <- length(intersect(false.index,positive.index))
   FN <- length(intersect(true.index,negative.index))
   TN <- length(intersect(false.index,negative.index))
   
   MCC.denom <- sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)
   if(MCC.denom==0) MCC.denom <- 1
   MCC <- (TP*TN-FP*FN)/MCC.denom
   if((TN+FP)==0) MCC <- 1
   
   Precision <- TP/(TP+FP)
   if((TP+FP)==0) Precision <- 1
   Recall <- TP/(TP+FN)
   if((TP+FN)==0) Recall <- 1
   Sensitivity <- Recall
   Specific <- TN/(TN+FP)
   if((TN+FP)==0) Specific <- 1
   
   
   return(list(Precision=Precision,Recall=Recall,Sensitivity=Sensitivity,Specific=Specific,MCC=MCC,TP=TP,FP=FP,TN=TN,FN=FN))
}





