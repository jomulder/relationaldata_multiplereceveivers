
library(mvtnorm)
library(truncnorm)
library(MCMCpack)
library(bayesplot)
library(dplyr)
library(parallel)
library(combinat)
library(gtools)
library(coda)


# helper function
list.rem <- function(list,element){
  list[[element]] <- NULL
  return(list)
}
# generate data with multiple receivers with p nodes and q latent variables
generate_MR <- function(p=8,q=2,sigma2b=1,sigma2L=1,beta=.5*(-1:1),nvec=20,Win=T,rhoX=.3,mu_c=NA,tau_c=NA,sd_c=NA){
  # number of nodes is p
  # dim(u_r) = q (dimension of latent variables)
  # sigma2b is variance of random effects popularity
  # sigma2L is latent variable variance
  # beta are fixed effects
  # nvec is number of messages send by each node
  # Win=T implies the inclusion of message specific latent variables
  # rhoX is the correlation between the predictor variables
  
  k <- length(beta)
  #generate list of covariate matrices X_s for all nodes, s=1,...,p
  
  nvec <- rpois(p,nvec)
  #nvec[nvec<=1] <- 2
  
  X <- lapply(1:p,function(s){
    if(nvec[s]>0){
      covm <- diag(k)*(1-rhoX) + rhoX
      do.call(rbind,
              lapply(1:nvec[s],function(i){
                Xmat <- rmvnorm(p,sigma=covm)
                Xmat[s,] <- NA
                Xmat
              })
      )
    }
  })
  
  # latent receveiver effects
  #  V <- matrix(0,nrow=p,ncol=q)
  V <- matrix(rnorm(q*p,sd=sqrt(sigma2L)),nrow=p,ncol=q)
  # V <- matrix(round(runif(p*q))*2-1,ncol=q)
  # while(min(apply(V,2,sd))==0){
  #   V <- matrix(round(runif(p*q))*2-1,ncol=q)
  # }
  
  # latent sender effects
  #  V <- V - rep(1,p)%*%t(apply(V,2,mean))
  #  U <- matrix(0,nrow=p,ncol=q) 
  U <- matrix(rnorm(q*p,sd=sqrt(sigma2L)),nrow=p,ncol=q)
  # U <- matrix(round(runif(p*q))*2-1,ncol=q)
  # while(min(apply(U,2,sd))==0){
  #   U <- matrix(round(runif(p*q))*2-1,ncol=q)
  # }
  
  # latent message effect
  W <- lapply(1:p,function(s){
    if(nvec[s]>0){
      Ws <- matrix(rnorm(nvec[s]*q,sd=sqrt(sigma2L)),nrow=nvec[s],ncol=q)*Win
      Ws[1,] <- 0
      Ws
    }
  })
  
  #  U <- U - rep(1,p)%*%t(apply(U,2,mean))
  bvec <- rnorm(p,mean=0,sd=sqrt(sigma2b)) # random receiver effects
  #  bvec <- round(runif(p))-.5 # random receiver effects
  #bvec <- bvec - mean(bvec)
  
  Z <- lapply(1:p,function(s){
    if(nvec[s]>0){
      matrix(rnorm(nvec[s]*p),nrow=nvec[s]) + t(matrix(c(X[[s]]%*%beta),ncol=nvec[s])) +
        rep(1,nvec[s])%*%t(bvec + U%*%V[s,]) + W[[s]]%*%t(U)
    }
  })
  
  # Draw thresholds c
  minZ <- min(unlist(Z),na.rm=T)
  maxZ <- max(unlist(Z),na.rm=T)
  if(is.na(mu_c)){ #use default based on the scale of Z
    mu_c <- maxZ
  }
  if(is.na(sd_c)){ #use default based on the scale of Z
    sd_c <- (maxZ-minZ)/15
  }
  if(is.na(tau_c)){ #use default based on the scale of Z
    tau_c <- 0
  }
  c <- lapply(1:p,function(s){
    mu_c_s <- rnorm(1,mu_c,sd=tau_c)
    if(nvec[s]>0){
      maxZ_s <- apply(Z[[s]],1,max,na.rm=T)
      rtruncnorm(nvec[s],a=-Inf,b=maxZ_s,mean=mu_c_s,sd=sd_c)
    }
  })
  Y <- lapply(1:p,function(s){
    if(nvec[s]>0){
      Ys <- matrix(as.integer(Z[[s]] > c[[s]]%*%t(rep(1,p))),ncol=p)
      Ys[,s] <- NA
      return(Ys)
    }
  })
  m <- lapply(1:p,function(s){
    if(nvec[s]>0){
      apply(Y[[s]],1,sum,na.rm=T)
    }
  })
  
  return(list(Y=Y,X=X,beta=beta,bvec=bvec,V=V,U=U,W=W,Z=Z,c=c,mu_c=mu_c,sd_c=sd_c,m=m))
}

# Gibbs sampler where the argument q is the dimension of the latent variables
Gibbs_sampler_MR <- function(Y,X,q=1,samsize=3e3,burnin=3e3,store=2,cores=1,whichrep=1:4,shape0c=1,scale0c=1,sd2=.5,shape0b=2,scale0b=2,startingvalues=NA){
  
  #prior hyperparemeters for distribution of c ~ N(mu_c,sd_c)
  #prior shape and scale for sd_c are shape0c and scale0c respectivevly
  #sd for random walk of mu_c and var_c are sd1 and sd2 respectively
  #prior hyperparameters for random effects variance for receiver popularity b, shape0b and scale0b
  
  p <- ncol(Y[[1]])
  k <- ncol(X[[1]])
  nvec <- unlist(lapply(1:length(Y),function(s){
    if(!is.null(Y[[s]])){
      nrow(Y[[s]])
    }else{ 0 }
  })) # number of messages per sender
  send0 <- (1:p)[nvec==0] #which nodes send no messages
  rec0 <- apply(do.call(rbind,Y),2,mean)==0 #which nodes receive no messages
  N <- sum(nvec) #total number of receivers
  mvec <- lapply(1:length(Y),function(s){
    if(!is.null(Y[[s]])){
      if(nrow(Y[[s]])>1){
        apply(Y[[s]][,-s],1,sum)
      }else{
        sum(Y[[s]][,-s])
      }
    }
  })
  
  #pre-process design matrices
  select1 <- unlist(lapply(1:p,function(s){
    vec1 <- rep(TRUE,p)
    vec1[s] <- FALSE
    rep(vec1,nvec[s])
  }))
  select1 <- which(select1==FALSE)
  Xstack <- do.call(rbind,X)[-select1,]
  tXX <- t(Xstack)%*%Xstack
  tXXi <- solve(tXX)
  #use bvec for initial value of beta as bvec are intercepts
  mean0c <- 0
  bvec <- mean0c - qnorm(1-apply(do.call(rbind,Y),2,mean,na.rm=TRUE))
  if(sum(bvec== -Inf)>0){
    bvec[which(bvec==-Inf)] <- min(bvec[-which(bvec==-Inf)]) - .2
  }
  tXy <- Reduce('+',
                lapply(which(nvec>0),function(s){
                  if(nvec[s] > 0){
                    t(X[[s]][-(p*(0:(nvec[s]-1))+s),]) %*% c(t(Y[[s]] - rep(1,nvec[s])%*%t(bvec))[-s,] )
                  }
                }) )
  estBeta <- c(solve(tXX) %*% tXy)
  remove(tXy)
  
  # initial values
  if(!is.na(startingvalues[1])){
    welk1 <- which(names(startingvalues)=="beta")
    if(length(welk1)==1){
      beta <- estBeta
    }else{
      beta <- rep(0,ncol(X[[1]]))
    }
    welk1 <- which(names(startingvalues)=="V")
    if(length(welk1)==1){
      V <- startingvalues[[welk1]]
    }else{
      V <- matrix(0,ncol=q,nrow=p)
    }
    welk1 <- which(names(startingvalues)=="U")
    if(length(welk1)==1){
      U <- startingvalues[[welk1]]
    }else{
      U <- matrix(0,ncol=q,nrow=p)
    }
    welk1 <- which(names(startingvalues)=="sigma2b")
    if(length(welk1)==1){
      sigma2b <- startingvalues[[welk1]]
    }else{
      sigma2b <- 1
    }
    welk1 <- which(names(startingvalues)=="Z")
    if(length(welk1)==1){
      Z <- startingvalues[[welk1]]
    }else{
      Z <- Y
    }
    welk1 <- which(names(startingvalues)=="W")
    if(length(welk1)==1){
      W <- startingvalues[[welk1]]
    }else{
      W <- lapply(1:p,function(s){ # latent message effects
        matrix(0,nrow=nvec[s],ncol=q)
      })
    }
    welk1 <- which(names(startingvalues)=="c")
    if(length(welk1)==1){
      c <- startingvalues[[welk1]]
    }else{
      c <- lapply(1:length(Y),function(s){
        if(!is.null(Y[[s]])){
          rep(0,nrow(Y[[s]]))
        }
      })
    }
    welk1 <- which(names(startingvalues)=="mu_c")
    if(length(welk1)==1){
      mu_c <- startingvalues[[welk1]]
    }else{
      mu_c <- mean0c
    }
    welk1 <- which(names(startingvalues)=="var_c")
    if(length(welk1)==1){
      var_c <- startingvalues[[welk1]]
    }else{
      var_c <- scale0c/shape0c
    }
    sd_c <- sqrt(var_c)
    welk1 <- which(names(startingvalues)=="bvec")
    if(length(welk1)==1){
      bvec <- startingvalues[[welk1]]
    }else{
      bvec <- mean0c - qnorm(1-apply(do.call(rbind,Y),2,mean,na.rm=TRUE))
      if(sum(bvec== -Inf)>0){
        bvec[which(bvec==-Inf)] <- min(bvec[-which(bvec==-Inf)]) - .2
      }
    }
  }else{
    beta <- estBeta #fixed effects
    V <- matrix(0,ncol=q,nrow=p) # latent sender effects
    U <- matrix(0,ncol=q,nrow=p) # latent receiver effects
    W <- lapply(1:p,function(s){ # latent message effects
      matrix(0,nrow=nvec[s],ncol=q)
    })
    bvec <- mean0c - qnorm(1-apply(do.call(rbind,Y),2,mean,na.rm=TRUE))
    if(sum(bvec== -Inf)>0){
      bvec[which(bvec==-Inf)] <- min(bvec[-which(bvec==-Inf)]) - .2
    }
    sigma2b <- 1
    Z <- Y
    c <- lapply(1:length(Y),function(s){
      if(!is.null(Y[[s]])){
        rep(0,nrow(Y[[s]]))
      }
    })
    mu_c <- mean0c
    var_c <- scale0c/shape0c
    sd_c <- sqrt(var_c)
  }
  
  store_beta <- matrix(0,nrow=samsize/store,ncol=length(beta))
  store_V <- array(0,dim=c(samsize/store,nrow(V),ncol(V)))
  store_U <- array(0,dim=c(samsize/store,nrow(V),ncol(V)))
  welkw <- which(nvec==max(nvec))[1]
  store_W1 <- array(0,dim=c(samsize/store,2,ncol(W[[welkw]])))
  store_bvec <- matrix(0,nrow=samsize/store,ncol=length(bvec))
  store_Z1 <- array(0,dim=c(samsize/store,2,2))
  store_Yrep <- lapply(1:length(whichrep),function(wi){
    array(0,dim=c(samsize/store,p))
  })
  store_c <- store_m <- array(0,dim=c(samsize/store,sum(nvec)))
  store_muvarc <- matrix(0,nrow=samsize/store,ncol=2)
  store_sigma2b <- matrix(1,nrow=samsize/store,ncol=1)
  
  Xbeta <- lapply(1:length(X),function(s){
    if(nvec[s]>0){
      t(matrix(X[[s]]%*%beta,ncol=nvec[s]))
    }
  })
  
  #g in Zellner's g prior for beta
  g <- N
  
  transitivity <- 1
  
  print("burn-in")
  
  pb = txtProgressBar(min = 0, max = burnin, initial = 0)
  
  #burn-in period
  for(ss in 1:burnin){
    
    #draw latent Z's
    Z <- mclapply(1:length(Y),function(s){
      if(nvec[s]>0){
        mean_s <- rep(1,nvec[s])%*%t(c(bvec + U%*%V[s,])) + Xbeta[[s]] +
          W[[s]]%*%t(U)
        lb <- (c[[s]]+1e3)*Y[[s]]-1e3
        ub <- c[[s]]*(1-Y[[s]]) + 1e3*Y[[s]]
        Zs_can <- matrix(rtruncnorm(nvec[s]*p,a=lb,b=ub,mean=mean_s),nrow=nvec[s]) #when a uniform improper distribution would be used for c.
        Zs_can_max <- apply(Zs_can,1,max,na.rm=T)
        Zs_max <- apply(Z[[s]],1,max,na.rm=T)
        welk <- runif(nvec[s])<pnorm(Zs_max,mean=mu_c,sd=sd_c)/pnorm(Zs_can_max,mean=mu_c,sd=sd_c)
        welk[which(pnorm(Zs_can_max,mean=mu_c,sd=sd_c)==0)] <- T
        Zs_can[!welk,] <- Z[[s]][!welk,]
        return(Zs_can)
      }
    },mc.cores=cores)
    
    #draw threshold values c
    lb <- mclapply(1:length(Y),function(s){ #these are also needed when sampling mu_c and cd_c
      if(nvec[s]>0){
        apply(Z[[s]]*(Y[[s]]==0) - 1e3*(Y[[s]]==1),1,max,na.rm=T)
      }
    },mc.cores=cores)
    ub <- mclapply(1:length(Y),function(s){ #these are also needed when sampling mu_c and cd_c
      if(nvec[s]>0){
        apply(Z[[s]]*(Y[[s]]==1) + 1e3*(Y[[s]]==0),1,min,na.rm=T)
      }
    },mc.cores=cores)
    c <- mclapply(1:length(Y),function(s){
      if(nvec[s]>0){
        rtruncnorm(nvec[s],a=lb[[s]],b=ub[[s]],mean=mu_c,sd=sd_c)
      }
    },mc.cores=cores)
    
    #draw sd sd_c of population N(mu_c,sd_c) for c using a random walk
    var_c_can <- rtruncnorm(1,a=0,b=Inf,mean=var_c,sd=sd2)
    sd_c_can <- sqrt(var_c_can)
    diffs <- pnorm(unlist(ub),mean=mu_c,sd=sd_c) - pnorm(unlist(lb),mean=mu_c,sd=sd_c)
    diffs_can <- pnorm(unlist(ub),mean=mu_c,sd=sd_c_can) - pnorm(unlist(lb),mean=mu_c,sd=sd_c_can)
    #only these provide information for mu_c & sd_c
    welk <- which((diffs>0)*(diffs_can>0)==1)
    R_MH <- exp(  sum(dnorm(unlist(c)[welk],mean=mu_c,sd=sd_c_can,log=T)) - 
                    sum(dnorm(unlist(c)[welk],mean=mu_c,sd=sd_c,log=T)) + 
                    sum(log(diffs[welk])) - sum(log(diffs_can[welk])) + 
                    log(dinvgamma(var_c_can,shape=shape0c,scale=scale0c)) -
                    log(dinvgamma(var_c,shape=shape0c,scale=scale0c))) * #IG(shape0c,scale0c) prior
      pnorm(var_c/sd2)/pnorm(var_c_can/sd2) #correction for sampling from a nonsymmetric truncated normal proposal distribution
    var_c <- ifelse(runif(1)<R_MH,var_c_can,var_c)
    sd_c <- sqrt(var_c)
    
    #draw fixed effects beta
    Lstack <- unlist(lapply(1:p,function(s){
      if(nvec[s]>0){
        rep(U%*%V[s,],nvec[s]) + c(U%*%t(W[[s]]))
      }
    }))
    ZAstack <- (c(t(do.call(rbind,Z))) - rep(bvec,sum(nvec)) - Lstack)[-select1]
    #Use Zellner's g prior
    covBeta <- tXXi*g/(g+1)
    meanBeta <- c(tXXi%*%t(Xstack)%*%ZAstack)*g/(g+1)
    beta <- c(rmvnorm(1,mean=meanBeta,sigma=covBeta))
    
    Xbeta <- lapply(1:length(X),function(s){
      if(nvec[s]>0){
        t(matrix(X[[s]]%*%beta,ncol=nvec[s]))
      }
    })
    
    # draw popularity effects b and latent receiver effects U
    bU <- t(matrix(unlist(mclapply(1:p,function(r){
      
      ZA_r <- unlist(lapply((1:p)[-r],function(s){
        Z[[s]][,r] - Xbeta[[s]][,r] 
      }))
      
      XA_r <- cbind(rep(1,sum(nvec)-nvec[r]),V[rep((1:p)[-r],times=nvec[-r]),] + 
                      do.call(rbind,list.rem(W,r)))
      Sigmau <- solve(t(XA_r)%*%XA_r+diag(c(1/sigma2b,rep(1,q))))
      meanu <- Sigmau%*%t(XA_r)%*%ZA_r
      c(rmvnorm(1,mean=meanu,sigma=Sigmau))
      
    },mc.cores=cores)),nrow=q+1))
    bvec <- bU[,1]
    U <- as.matrix(bU[,-1])
    
    # draw random effects variance
    sigma2b = rinvgamma(1,shape=(p-1)/2+shape0b,scale=sum(bvec[-1]**2)/2+scale0b)
    
    # draw latent sender effects V
    V <- t(matrix(unlist(mclapply(1:p,function(s){
      if(nvec[s]>0){
        ZA_s <- c(t((Z[[s]] - Xbeta[[s]])[,-s])) -
          rep(bvec[-s],nvec[s]) -
          c((U%*%t(W[[s]]))[-s,])
        XA_s <- do.call(rbind,replicate(nvec[s],as.matrix(U[-s,]),simplify=F))
        Sigmav <- solve(t(XA_s)%*%XA_s+diag(q))
        meanv <- Sigmav%*%t(XA_s)%*%ZA_s
        c(rmvnorm(1,mean=meanv,sigma=Sigmav))
      }else{rnorm(q)}
    },mc.cores=cores)),nrow=q))
    
    # draw latent message effects W
    W <- mclapply(1:p,function(s){
      if(nvec[s]>0){
        Sigma_Ws <- solve(t(U[-s,])%*%U[-s,] + diag(q))
        Ws <- rmvnorm(nvec[s],sigma=Sigma_Ws) +
          ((Z[[s]] - Xbeta[[s]] - rep(1,nvec[s])%*%t(bvec + U%*%V[s,]))[,-s]) %*%
          U[-s,] %*% Sigma_Ws
        Ws[1,] <- 0
        Ws
      }
    },mc.cores=cores)
    
    setTxtProgressBar(pb,ss)
  }
  cat("\n")
  
  postburnin <- list(Z=Z,c=c,lb=lb,ub=ub,mu_c=mu_c,sd_c=sd_c,beta=beta,bvec=bvec,sigma2b=sigma2b,U=U,V=V,W=W)

  cat("\n")
  print(Sys.time())
  cat("\n")
  print("posterior sampling")
  
  pb = txtProgressBar(min = 0, max = samsize, initial = 0)
  
  teller <- 0
  
  for(ss in 1:samsize){
    
    #draw latent Z's
    Z <- mclapply(1:length(Y),function(s){
      if(nvec[s]>0){
        mean_s <- rep(1,nvec[s])%*%t(c(bvec + U%*%V[s,])) + Xbeta[[s]] +
          W[[s]]%*%t(U)
        lb <- (c[[s]]+1e3)*Y[[s]]-1e3
        ub <- c[[s]]*(1-Y[[s]]) + 1e3*Y[[s]]
        Zs_can <- matrix(rtruncnorm(nvec[s]*p,a=lb,b=ub,mean=mean_s),nrow=nvec[s]) #when a uniform improper distribution would be used for c.
        Zs_can_max <- apply(Zs_can,1,max,na.rm=T)
        Zs_max <- apply(Z[[s]],1,max,na.rm=T)
        welk <- runif(nvec[s])<pnorm(Zs_max,mean=mu_c,sd=sd_c)/pnorm(Zs_can_max,mean=mu_c,sd=sd_c)
        welk[which(pnorm(Zs_can_max,mean=mu_c,sd=sd_c)==0)] <- T
        Zs_can[!welk,] <- Z[[s]][!welk,]
        return(Zs_can)
      }
    },mc.cores=cores)
    
    #draw threshold values c
    lb <- mclapply(1:length(Y),function(s){ #these are also needed when sampling mu_c and cd_c
      if(nvec[s]>0){
        apply(Z[[s]]*(Y[[s]]==0) - 1e3*(Y[[s]]==1),1,max,na.rm=T)
      }
    },mc.cores=cores)
    ub <- mclapply(1:length(Y),function(s){ #these are also needed when sampling mu_c and cd_c
      if(nvec[s]>0){
        apply(Z[[s]]*(Y[[s]]==1) + 1e3*(Y[[s]]==0),1,min,na.rm=T)
      }
    },mc.cores=cores)
    c <- mclapply(1:length(Y),function(s){
      if(nvec[s]>0){
        rtruncnorm(nvec[s],a=lb[[s]],b=ub[[s]],mean=mu_c,sd=sd_c)
      }
    },mc.cores=cores)
    
    #draw sd sd_c of population N(mu_c,sd_c) for c using a random walk
    var_c_can <- rtruncnorm(1,a=0,b=Inf,mean=var_c,sd=sd2)
    sd_c_can <- sqrt(var_c_can)
    diffs <- pnorm(unlist(ub),mean=mu_c,sd=sd_c) - pnorm(unlist(lb),mean=mu_c,sd=sd_c)
    diffs_can <- pnorm(unlist(ub),mean=mu_c,sd=sd_c_can) - pnorm(unlist(lb),mean=mu_c,sd=sd_c_can)
    #only these provide information for mu_c & sd_c
    welk <- which((diffs>0)*(diffs_can>0)==1)
    R_MH <- exp(  sum(dnorm(unlist(c)[welk],mean=mu_c,sd=sd_c_can,log=T)) - 
                    sum(dnorm(unlist(c)[welk],mean=mu_c,sd=sd_c,log=T)) + 
                    sum(log(diffs[welk])) - sum(log(diffs_can[welk])) + 
                    log(dinvgamma(var_c_can,shape=shape0c,scale=scale0c)) -
                    log(dinvgamma(var_c,shape=shape0c,scale=scale0c))) * #IG(shape0c,scale0c) prior
      pnorm(var_c/sd2)/pnorm(var_c_can/sd2) #correction for sampling from a nonsymmetric truncated normal proposal distribution
    var_c <- ifelse(runif(1)<R_MH,var_c_can,var_c)
    sd_c <- sqrt(var_c)
    
    #draw fixed effects beta
    Lstack <- unlist(lapply(1:p,function(s){
      if(nvec[s]>0){
        rep(U%*%V[s,],nvec[s]) + c(U%*%t(W[[s]]))
      }
    }))
    ZAstack <- (c(t(do.call(rbind,Z))) - rep(bvec,sum(nvec)) - Lstack)[-select1]
    #Use Zellner's g prior
    covBeta <- tXXi*g/(g+1)
    meanBeta <- c(tXXi%*%t(Xstack)%*%ZAstack)*g/(g+1)
    beta <- c(rmvnorm(1,mean=meanBeta,sigma=covBeta))
    
    Xbeta <- lapply(1:length(X),function(s){
      if(nvec[s]>0){
        t(matrix(X[[s]]%*%beta,ncol=nvec[s]))
      }
    })
    
    # draw popularity effects b and latent receiver effects U
    bU <- t(matrix(unlist(mclapply(1:p,function(r){
      
      ZA_r <- unlist(lapply((1:p)[-r],function(s){
        Z[[s]][,r] - Xbeta[[s]][,r] #- bvec[r]
      }))
      
      XA_r <- cbind(rep(1,sum(nvec)-nvec[r]),V[rep((1:p)[-r],times=nvec[-r]),] + 
                      do.call(rbind,list.rem(W,r)))
      Sigmau <- solve(t(XA_r)%*%XA_r+diag(c(1/sigma2b,rep(1,q))))
      meanu <- Sigmau%*%t(XA_r)%*%ZA_r
      c(rmvnorm(1,mean=meanu,sigma=Sigmau))
      
    },mc.cores=cores)),nrow=q+1))
    bvec <- bU[,1]
    U <- as.matrix(bU[,-1])

    # draw random effects variance
    sigma2b = rinvgamma(1,shape=(p-1)/2+shape0b,scale=sum(bvec[-1]**2)/2+scale0b)
    
    # draw latent sender effects V
    V <- t(matrix(unlist(mclapply(1:p,function(s){
      if(nvec[s]>0){
        ZA_s <- c(t((Z[[s]] - Xbeta[[s]])[,-s])) -
          rep(bvec[-s],nvec[s]) -
          c((U%*%t(W[[s]]))[-s,])
        XA_s <- do.call(rbind,replicate(nvec[s],as.matrix(U[-s,]),simplify=F))
        Sigmav <- solve(t(XA_s)%*%XA_s+diag(q))
        meanv <- Sigmav%*%t(XA_s)%*%ZA_s
        c(rmvnorm(1,mean=meanv,sigma=Sigmav))
      }else{rnorm(q)}
    },mc.cores=cores)),nrow=q))
    
    # draw latent message effects W
    W <- mclapply(1:p,function(s){
      if(nvec[s]>0){
        Sigma_Ws <- solve(t(U[-s,])%*%U[-s,] + diag(q))
        Ws <- rmvnorm(nvec[s],sigma=Sigma_Ws) +
          ((Z[[s]] - Xbeta[[s]] - rep(1,nvec[s])%*%t(bvec + U%*%V[s,]))[,-s]) %*%
          U[-s,] %*% Sigma_Ws
        Ws[1,] <- 0
        Ws
      }
    },mc.cores=cores)
    
    if(ceiling(ss/store)==ss/store){
      #generate replicated data for ppc's
      Yrep <- lapply(1:length(Y),function(s){
        if(nvec[s]>0){
          Zs <- Xbeta[[s]] + rep(1,nvec[s])%*%(bvec + t(V[s,])%*%t(U)) +
            W[[s]]%*%t(U) + matrix(rnorm(nvec[s]*p),nrow=nvec[s])
          Zs[,s] <- NA
          Zs_max <- apply(Zs,1,max,na.rm=T)
          c_s <- rtruncnorm(nvec[s],a=-Inf,b=Zs_max,mean=mu_c,sd=sd_c)
          matrix(as.integer(Zs > c_s%*%t(rep(1,p))),nrow=nvec[s])
        }
      })
      
      recpref <- lapply(1:length(whichrep),function(wi){
        apply(Yrep[[whichrep[wi]]],2,sum) / nrow(Yrep[[whichrep[wi]]])
      })
      
      #number of receivers
      mvec <- unlist(lapply(1:length(Y),function(s){
        if(nvec[s]>0){
          apply(Yrep[[s]],1,sum,na.rm=T)
        }
      }))

      teller <- teller + 1
      store_beta[teller,] <- beta
      store_V[teller,,] <- V
      store_U[teller,,] <- U
      store_W1[teller,,] <- W[[welkw]][1:2,]
      store_bvec[teller,] <- bvec
      store_Z1[teller,,] <- Z[[1]][1:2,c(5,8)]
      for(wi in 1:length(whichrep)){
        store_Yrep[[wi]][teller,] <- recpref[[wi]]
      }
      store_c[teller,] <- unlist(c)
      store_muvarc[teller,] <- c(mu_c,sd_c)
      store_m[teller,] <- mvec
      store_sigma2b[teller,1] <- sigma2b
      
      if(ceiling(teller/10)==teller/10){
        par(mfrow=c(4,2))
        plot(1:teller,store_beta[1:teller,1],"l",xlab="iterations",ylab="beta1")
        plot(1:teller,store_bvec[1:teller,1],"l",xlab="iterations",ylab="b1")
        plot(1:teller,store_V[1:teller,1,1],"l",xlab="iterations",ylab="V11")
        plot(1:teller,store_U[1:teller,1,1],"l",xlab="iterations",ylab="U11")
        welk <- 1
        recpref <- apply(Y[[welk]],2,sum) / nrow(Y[[welk]])
        welkpop <- order(recpref,decreasing=TRUE)[1:10]
        plot(1:length(welkpop),recpref[welkpop],"l",ylab="recpop1",ylim=c(0,1))
        lb1 <- max(teller-100,1) # only plot for last 100 replicated datasets
        for(rep in lb1:teller){
          lines(store_Yrep[[welk]][rep,welkpop],col="grey")
        }
        lines(recpref[welkpop],col=2)
        welk <- 2
        recpref <- apply(Y[[welk]],2,sum) / nrow(Y[[welk]])
        welkpop <- order(recpref,decreasing=TRUE)[1:10]
        plot(1:length(welkpop),recpref[welkpop],"l",ylab="recpop2",ylim=c(0,1))
        for(rep in lb1:teller){
          lines(store_Yrep[[welk]][rep,welkpop],col="grey")
        }
        lines(recpref[welkpop],col=2)
        welk <- 3
        recpref <- apply(Y[[welk]],2,sum) / nrow(Y[[welk]])
        welkpop <- order(recpref,decreasing=TRUE)[1:10]
        plot(1:length(welkpop),recpref[welkpop],"l",ylab="recpop3",ylim=c(0,1))
        for(rep in lb1:teller){
          lines(store_Yrep[[welk]][rep,welkpop],col="grey")
        }
        lines(recpref[welkpop],col=2)
        plot(1:teller,store_muvarc[1:teller,2],"l",xlab="iterations",ylab="sd_c")
      }
      
    }
    
    setTxtProgressBar(pb,ss)
  }
  
  postmcmc <- list(beta=beta,Z=Z,c=c,lb=lb,ub=ub,mu_c=mu_c,var_c=var_c,bvec=bvec,sigma2b=sigma2b,U=U,V=V,W=W)
  
  return(list(beta=store_beta,bvec=store_bvec,V=store_V,
              U=store_U,W1=store_W1,
              Z1=store_Z1,
              c=store_c,
              muvarc=store_muvarc,
              m=store_m,
              Yrep=store_Yrep,
              sigma2b=store_sigma2b,
              postburnin=postburnin,
              postmcmc=postmcmc
  ))
}
# Gibbs sampler in case of no latent variables
Gibbs_sampler_MR0 <- function(Y,X,samsize=3e3,burnin=3e3,store=2,cores=1,whichrep=1:4,shape0c=1,scale0c=1,sd2=.5,shape0b=2,scale0b=2,startingvalues=NA){
  
  #prior hyperparemeters for distribution of c ~ N(mu_c,sd_c)
  #prior mean and sd for mu_c are mean0c and sd0c respectively
  #prior shape and scale for sd_c are shape0c and scale0c respectivevly
  #sd for random walk of mu_c and var_c are sd1 and sd2 respectively
  #prior hyperparameters for random effects variance for receiver popularity b, shape0b and scale0b
  
  p <- ncol(Y[[1]])
  k <- ncol(X[[1]])
  nvec <- unlist(lapply(1:length(Y),function(s){
    if(!is.null(Y[[s]])){
      nrow(Y[[s]])
    }else{ 0 }
  })) # number of messages per sender
  send0 <- (1:p)[nvec==0] #which nodes send no messages
  rec0 <- apply(do.call(rbind,Y),2,mean)==0 #which nodes receive no messages
  N <- sum(nvec) #total number of receivers
  mvec <- lapply(1:length(Y),function(s){
    if(!is.null(Y[[s]])){
      if(nrow(Y[[s]])>1){
        apply(Y[[s]][,-s],1,sum)
      }else{
        sum(Y[[s]][,-s])
      }
    }
  })
  select1 <- unlist(lapply(1:p,function(s){
    vec1 <- rep(T,p)
    vec1[s] <- F
    rep(vec1,nvec[s])
  }))
  select1 <- which(select1==FALSE)
  Xstack <- do.call(rbind,X)[-select1,]
  tXX <- t(Xstack)%*%Xstack
  tXXi <- solve(tXX)
  #use bvec for initial value of beta as bvec are intercepts
  mean0c <- 0
  bvec <- mean0c - qnorm(1-apply(do.call(rbind,Y),2,mean,na.rm=TRUE))
  if(sum(bvec== -Inf)>0){
    bvec[which(bvec==-Inf)] <- min(bvec[-which(bvec==-Inf)]) - .2
  }
  
  tXy <- Reduce('+',
                lapply(1:p,function(s){
                  if(nvec[s] > 0){
                    t(X[[s]][-((0:(nvec[s]-1))*p+s),]) %*% c(t(Y[[s]] - rep(1,nvec[s])%*%t(bvec))[-s,])
                  }
                }) )
  estBeta <- c(solve(tXX) %*% tXy)
  remove(tXy)
  
  # initial values
  mean0c <- 0
  if(!is.na(startingvalues[1])){
    welk1 <- which(names(startingvalues)=="beta")
    if(length(welk1)==1){
      beta <- startingvalues[[welk1]]
    }else{
      beta <- estBeta
    }
    welk1 <- which(names(startingvalues)=="bvec")
    if(length(welk1)==1){
      bvec <- startingvalues[[welk1]]
    }else{
      bvec <- mean0c - qnorm(1-apply(do.call(rbind,Y),2,mean,na.rm=TRUE))
      if(sum(bvec== -Inf)>0){
        bvec[which(bvec==-Inf)] <- min(bvec[-which(bvec==-Inf)]) - .2
      }
    }
    welk1 <- which(names(startingvalues)=="sigma2b")
    if(length(welk1)==1){
      sigma2b <- startingvalues[[welk1]]
    }else{
      sigma2b <- 1
    }
    welk1 <- which(names(startingvalues)=="Z")
    if(length(welk1)==1){
      Z <- startingvalues[[welk1]]
    }else{
      Z <- Y
    }
    welk1 <- which(names(startingvalues)=="mu_c")
    if(length(welk1)==1){
      mu_c <- startingvalues[[welk1]]
    }else{
      mu_c <- mean0c
    }
    welk1 <- which(names(startingvalues)=="c")
    if(length(welk1)==1){
      c <- startingvalues[[welk1]]
    }else{
      c <- lapply(1:length(Y),function(s){
        if(!is.null(Y[[s]])){
          rep(0,nrow(Y[[s]]))
        }
      })
    }
    welk1 <- which(names(startingvalues)=="var_c")
    if(length(welk1)==1){
      var_c <- startingvalues[[welk1]]
    }else{
      var_c <- scale0c/shape0c
    }
    sd_c <- sqrt(var_c)
  }else{
    beta <- estBeta #fixed effects
    bvec <- mean0c - qnorm(1-apply(do.call(rbind,Y),2,mean,na.rm=TRUE))
    if(sum(bvec== -Inf)>0){
      bvec[which(bvec==-Inf)] <- min(bvec[-which(bvec==-Inf)]) - .2
    }
    sigma2b <- 1
    Z <- Y
    c <- lapply(1:length(Y),function(s){
      if(!is.null(Y[[s]])){
        rep(0,nrow(Y[[s]]))
      }
    })
    mu_c <- mean0c
    var_c <- scale0c/shape0c
    sd_c <- sqrt(var_c)
  }
  
  store_beta <- matrix(0,nrow=samsize/store,ncol=length(beta))
  store_bvec <- store_prep <- matrix(0,nrow=samsize/store,ncol=length(bvec))
  store_Z1 <- array(0,dim=c(samsize/store,2,2))
  store_Yrep <- lapply(1:length(whichrep),function(wi){
    array(0,dim=c(samsize/store,p))
  })
  store_c <- store_m <- array(0,dim=c(samsize/store,sum(nvec)))
  store_muvarc <- matrix(0,nrow=samsize/store,ncol=2)
  store_sigma2b <- store_gLS <- matrix(1,nrow=samsize/store,ncol=1)
  
  Xbeta <- lapply(1:length(X),function(s){
    if(nvec[s]>0){
      t(matrix(X[[s]]%*%beta,ncol=nvec[s]))
    }
  })
  
  #g in Zellner's g prior for beta
  g <- N

  transitivity <- 1
  
  print("burn-in")
  
  pb = txtProgressBar(min = 0, max = burnin, initial = 0)
  
  #burn-in period
  for(ss in 1:burnin){
    
    #draw latent Z's
    Z <- mclapply(1:length(Y),function(s){
      if(nvec[s]>0){
        mean_s <- rep(1,nvec[s])%*%t(bvec) + Xbeta[[s]]
        lb <- (c[[s]]+1e3)*Y[[s]]-1e3
        ub <- c[[s]]*(1-Y[[s]]) + 1e3*Y[[s]]
        Zs_can <- matrix(rtruncnorm(nvec[s]*p,a=lb,b=ub,mean=mean_s),nrow=nvec[s]) #when a uniform improper distribution would be used for c.
        Zs_can_max <- apply(Zs_can,1,max,na.rm=T)
        Zs_max <- apply(Z[[s]],1,max,na.rm=T)
        welk <- runif(nvec[s])<pnorm(Zs_max,mean=mu_c,sd=sd_c)/pnorm(Zs_can_max,mean=mu_c,sd=sd_c)
        welk[which(pnorm(Zs_can_max,mean=mu_c,sd=sd_c)==0)] <- T
        Zs_can[!welk,] <- Z[[s]][!welk,]
        return(Zs_can)
      }
    },mc.cores=cores)
    
    lb <- mclapply(1:length(Y),function(s){ #these are also needed when sampling mu_c and cd_c
      if(nvec[s]>0){
        apply(Z[[s]]*(Y[[s]]==0) - 1e3*(Y[[s]]==1),1,max,na.rm=T)
      }
    },mc.cores=cores)
    ub <- mclapply(1:length(Y),function(s){ #these are also needed when sampling mu_c and cd_c
      if(nvec[s]>0){
        apply(Z[[s]]*(Y[[s]]==1) + 1e3*(Y[[s]]==0),1,min,na.rm=T)
      }
    },mc.cores=cores)
    c <- mclapply(1:length(Y),function(s){
      if(nvec[s]>0){
        rtruncnorm(nvec[s],a=lb[[s]],b=ub[[s]],mean=mu_c,sd=sd_c)
      }
    },mc.cores=cores)
    
    #draw sd sd_c of population N(mu_c,sd_c) for c using a random walk
    var_c_can <- rtruncnorm(1,a=0,b=Inf,mean=var_c,sd=sd2)
    sd_c_can <- sqrt(var_c_can)
    diffs <- pnorm(unlist(ub),mean=mu_c,sd=sd_c) - pnorm(unlist(lb),mean=mu_c,sd=sd_c)
    diffs_can <- pnorm(unlist(ub),mean=mu_c,sd=sd_c_can) - pnorm(unlist(lb),mean=mu_c,sd=sd_c_can)
    #only these provide information for mu_c & sd_c
    welk <- which((diffs>0)*(diffs_can>0)==1)
    R_MH <- exp( sum(dnorm(unlist(c)[welk],mean=mu_c,sd=sd_c_can,log=T)) - 
                   sum(dnorm(unlist(c)[welk],mean=mu_c,sd=sd_c,log=T)) + 
                   sum(log(diffs[welk])) - sum(log(diffs_can[welk])) + 
                   log(dinvgamma(var_c_can,shape=shape0c,scale=scale0c)) -
                   log(dinvgamma(var_c,shape=shape0c,scale=scale0c))) * #IG(shape0c,scale0c) prior
      pnorm(var_c/sd2)/pnorm(var_c_can/sd2) #correction for sampling from a nonsymmetric truncated normal proposal distribution
    var_c <- ifelse(runif(1) < R_MH, var_c_can, var_c)
    sd_c <- sqrt(var_c)
    
    #draw fixed effects beta
    ZAstack <- (c(t(do.call(rbind,Z))) - rep(bvec,sum(nvec)))[-select1]
    #Use Zellner's g prior
    covBeta <- tXXi*g/(g+1)
    meanBeta <- c(tXXi%*%t(Xstack)%*%ZAstack)*g/(g+1)
    beta <- c(rmvnorm(1,mean=meanBeta,sigma=covBeta))
    Xbeta <- lapply(1:length(X),function(s){
      if(nvec[s]>0){
        t(matrix(X[[s]]%*%beta,ncol=nvec[s]))
      }
    })
    
    # draw (random) effects bvec
    ZAstack <- do.call(rbind,Z) - do.call(rbind,Xbeta) # t(matrix(XstackNA%*%beta,nrow=p))
    meanpostb <- (apply(ZAstack,2,mean,na.rm=TRUE) * (rep(N,p) - nvec)) /
      (rep(N,p) - nvec + 1/sigma2b) #posterior mean
    sig2postb <- 1/(rep(N,p) - nvec + 1/sigma2b)
    bvec <- rnorm(p,mean=meanpostb,sd=sqrt(sig2postb))

    # draw random effects variance
    sigma2b = rinvgamma(1,shape=(p-1)/2+shape0b,scale=sum(bvec[-1]**2)/2+scale0b)
    
    setTxtProgressBar(pb,ss)
  }
  cat("\n")
  
  postburnin <- list(Z=Z,c=c,lb=lb,ub=ub,mu_c=mu_c,sd_c=sd_c,beta=beta,bvec=bvec,sigma2b=sigma2b)

  cat("\n")
  print(Sys.time())
  cat("\n")
  print("posterior sampling")
  
  pb = txtProgressBar(min = 0, max = samsize, initial = 0)
  
  teller <- 0
  
  for(ss in 1:samsize){
    
    #draw latent Z's
    Z <- mclapply(1:length(Y),function(s){
      if(nvec[s]>0){
        mean_s <- rep(1,nvec[s])%*%t(bvec) + Xbeta[[s]]
        lb <- (c[[s]]+1e3)*Y[[s]]-1e3
        ub <- c[[s]]*(1-Y[[s]]) + 1e3*Y[[s]]
        Zs_can <- matrix(rtruncnorm(nvec[s]*p,a=lb,b=ub,mean=mean_s),nrow=nvec[s])
        Zs_can_max <- apply(Zs_can,1,max,na.rm=T)
        Zs_max <- apply(Z[[s]],1,max,na.rm=T)
        welk <- runif(nvec[s])<pnorm(Zs_max,mean=mu_c,sd=sd_c)/pnorm(Zs_can_max,mean=mu_c,sd=sd_c)
        welk[which(pnorm(Zs_can_max,mean=mu_c,sd=sd_c)==0)] <- T
        Zs_can[!welk,] <- Z[[s]][!welk,]
        return(Zs_can)
      }
    },mc.cores=cores)
    
    #draw threshold values c
    lb <- mclapply(1:length(Y),function(s){ #these are also needed when sampling mu_c and cd_c
      if(nvec[s]>0){
        apply(Z[[s]]*(Y[[s]]==0) - 1e3*(Y[[s]]==1),1,max,na.rm=T)
      }
    },mc.cores=cores)
    ub <- mclapply(1:length(Y),function(s){ #these are also needed when sampling mu_c and cd_c
      if(nvec[s]>0){
        apply(Z[[s]]*(Y[[s]]==1) + 1e3*(Y[[s]]==0),1,min,na.rm=T)
      }
    },mc.cores=cores)
    c <- mclapply(1:length(Y),function(s){
      if(nvec[s]>0){
        rtruncnorm(nvec[s],a=lb[[s]],b=ub[[s]],mean=mu_c,sd=sd_c)
      }
    },mc.cores=cores)
    
    #draw sd sd_c of population N(mu_c,sd_c) for c using a random walk
    var_c_can <- rtruncnorm(1,a=0,b=Inf,mean=var_c,sd=sd2)
    sd_c_can <- sqrt(var_c_can)
    diffs <- pnorm(unlist(ub),mean=mu_c,sd=sd_c) - pnorm(unlist(lb),mean=mu_c,sd=sd_c)
    diffs_can <- pnorm(unlist(ub),mean=mu_c,sd=sd_c_can) - pnorm(unlist(lb),mean=mu_c,sd=sd_c_can)
    #only these provide information for mu_c & sd_c
    welk <- which((diffs>0)*(diffs_can>0)==1)
    R_MH <- exp(  sum(dnorm(unlist(c)[welk],mean=mu_c,sd=sd_c_can,log=T)) - 
                    sum(dnorm(unlist(c)[welk],mean=mu_c,sd=sd_c,log=T)) + 
                    sum(log(diffs[welk])) - sum(log(diffs_can[welk])) + 
                    log(dinvgamma(var_c_can,shape=shape0c,scale=scale0c)) -
                    log(dinvgamma(var_c,shape=shape0c,scale=scale0c))) * #IG(shape0c,scale0c) prior
      pnorm(var_c/sd2)/pnorm(var_c_can/sd2) #correction for sampling from a nonsymmetric truncated normal proposal distribution
    var_c <- ifelse(runif(1)<R_MH,var_c_can,var_c)
    sd_c <- sqrt(var_c)
    
    #draw fixed effects beta
    ZAstack <- (c(t(do.call(rbind,Z))) - rep(bvec,sum(nvec)))[-select1]
    
    #Use Zellner's g prior
    covBeta <- tXXi*g/(g+1)
    meanBeta <- c(tXXi%*%t(Xstack)%*%ZAstack)*g/(g+1)
    beta <- c(rmvnorm(1,mean=meanBeta,sigma=covBeta))
    #beta <- rep(0,length(beta))
    Xbeta <- lapply(1:length(X),function(s){
      if(nvec[s]>0){
        t(matrix(X[[s]]%*%beta,ncol=nvec[s]))
      }
    })
    
    # draw (random) effects bvec
    ZAstack <- do.call(rbind,Z) - do.call(rbind,Xbeta) # t(matrix(XstackNA%*%beta,nrow=p))
    meanpostb <- (apply(ZAstack,2,mean,na.rm=TRUE) * (rep(N,p) - nvec)) /
      (rep(N,p) - nvec + 1/sigma2b) #posterior mean
    sig2postb <- 1/(rep(N,p) - nvec + 1/sigma2b)
    bvec <- rnorm(p,mean=meanpostb,sd=sqrt(sig2postb))
    
    # draw random effects variance
    sigma2b = rinvgamma(1,shape=(p-1)/2+shape0b,scale=sum(bvec[-1]**2)/2+scale0b)
    
    if(ceiling(ss/store)==ss/store){
      #generate replicated data for ppc's
      Yrep <- lapply(1:length(Y),function(s){
        if(nvec[s]>0){
          Zs <- rep(1,nvec[s])%*%t(bvec) + Xbeta[[s]] +
            matrix(rnorm(nvec[s]*p),nrow=nvec[s])
          Zs[,s] <- NA
          Zs_max <- apply(Zs,1,max,na.rm=T)
          c_s <- rtruncnorm(nvec[s],a=-Inf,b=Zs_max,mean=mu_c,sd=sd_c)
          matrix(as.integer(Zs > c_s%*%t(rep(1,p))),nrow=nvec[s])
        }
      })
      
      recpref <- lapply(1:length(whichrep),function(wi){
        apply(Yrep[[whichrep[wi]]],2,sum) / nrow(Yrep[[whichrep[wi]]])
      })
      
      #number of receivers
      mvec <- unlist(lapply(1:length(Y),function(s){
        if(nvec[s]>0){
          apply(Yrep[[s]],1,sum,na.rm=T)
        }
      }))
      
      prep <- apply(do.call(rbind,Yrep),2,mean,na.rm=TRUE)
      
      teller <- teller + 1
      store_beta[teller,] <- beta
      store_bvec[teller,] <- bvec
      store_prep[teller,] <- prep
      store_Z1[teller,1:2,1:2] <- Z[[1]][1:2,2:3]
      for(wi in 1:length(whichrep)){
        store_Yrep[[wi]][teller,] <- recpref[[wi]]
      }
      store_c[teller,] <- unlist(c)
      store_muvarc[teller,] <- c(mu_c,sd_c)
      store_m[teller,] <- mvec
      store_sigma2b[teller,1] <- sigma2b

      if(ceiling(teller/10)==teller/10){
        par(mfrow=c(3,2))
        plot(1:teller,store_beta[1:teller,1],"l",xlab="iterations",ylab="beta1")
        plot(1:teller,store_bvec[1:teller,1],"l",xlab="iterations",ylab="b1")
        welk <- 1
        recpref <- apply(Y[[welk]],2,sum) / nrow(Y[[welk]])
        welkpop <- order(recpref,decreasing=TRUE)[1:10]
        plot(1:length(welkpop),recpref[welkpop],"l",ylab="recpop1",ylim=c(0,1))
        lb1 <- max(teller-100,1)
        for(rep in lb1:teller){
          lines(store_Yrep[[welk]][rep,welkpop],col="grey")
        }
        lines(recpref[welkpop],col=2)
        welk <- 2
        recpref <- apply(Y[[welk]],2,sum) / nrow(Y[[welk]])
        welkpop <- order(recpref,decreasing=TRUE)[1:10]
        plot(1:length(welkpop),recpref[welkpop],"l",ylab="recpop2",ylim=c(0,1))
        for(rep in lb1:teller){
          lines(store_Yrep[[welk]][rep,welkpop],col="grey")
        }
        lines(recpref[welkpop],col=2)
        welk <- 3
        recpref <- apply(Y[[welk]],2,sum) / nrow(Y[[welk]])
        welkpop <- order(recpref,decreasing=TRUE)[1:10]
        plot(1:length(welkpop),recpref[welkpop],"l",ylab="recpop3",ylim=c(0,1))
        for(rep in lb1:teller){
          lines(store_Yrep[[welk]][rep,welkpop],col="grey")
        }
        lines(recpref[welkpop],col=2)
        plot(1:teller,store_muvarc[1:teller,2],"l",xlab="iterations",ylab="sd_c")
      }
      
    }
    
    setTxtProgressBar(pb,ss)
  }
  
  postmcmc <- list(beta=beta,Z=Z,c=c,lb=lb,ub=ub,mu_c=mu_c,var_c=var_c,bvec=bvec,sigma2b=sigma2b)
  
  return(list(beta=store_beta,bvec=store_bvec,
              c=store_c,
              muvarc=store_muvarc,
              m=store_m,
              Yrep=store_Yrep,
              sigma2b=store_sigma2b,
              postburnin=postburnin,
              postmcmc=postmcmc
  ))
}

# Example run
# Generate relational events with multiple receivers
set.seed(1234)
#genData <- generate_MR(p=20,q=1,sigma2b=1**2,sigma2L=1**2,beta=c(.5,0,-1),nvec=40,Win=T,rhoX=.3)
genData <- generate_MR(p=10,q=1,sigma2b=1**2,sigma2L=8**2,beta=(-2:2)/2,nvec=10,Win=T,rhoX=.3,mu_c=10,sd_c=2)

Y <- genData$Y
X <- genData$X

#fit with 0 latent variables
out0 <- Gibbs_sampler_MR0(Y,X,samsize=3e3,burnin=1e3,store=2,cores=1,whichrep=1:4,shape0c=10,scale0c=1,sd2=.5,shape0b=2,scale0b=2,startingvalues=NA)
#fit with 1 latent variable
out1 <- Gibbs_sampler_MR(Y,X,q=1,samsize=3e3,burnin=1e3,store=2,cores=1,whichrep=1:5,shape0c=10,scale0c=1,sd2=.5,shape0b=2,scale0b=2,startingvalues=NA)

mcmc_trace(data.frame(out1$beta))
mcmc_trace(data.frame(out1$bvec[,1:10]))
mcmc_trace(data.frame(out1$U[,1:10,1]))
mcmc_trace(data.frame(out1$V[,1:10,1]))

mcmc_areas_ridges(data.frame(out1$beta))
mcmc_areas_ridges(data.frame(out1$bvec[,1:10]))
mcmc_areas_ridges(data.frame(out1$U[,1:10,1]))
mcmc_areas_ridges(data.frame(out1$V[,1:10,1]))

#posterior predictive check of receiver popularity
welk <- 4
recpref <- apply(Y[[welk]],2,sum) / nrow(Y[[welk]])
plot(1:ncol(Y[[1]]),recpref,"l",col=2,ylim=c(0,1))
for(rep in 1:nrow(out1$Yrep[[welk]])){
  lines(out1$Yrep[[welk]][rep,],col="grey")
}
lines(recpref,col=2)

