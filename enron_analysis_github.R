
library(dplyr)
library(parallel)
library(bayesplot)

load("enron.rda")

# helper function
checkreceiver <- function(listrec,reci){
  check <- function(vec){ sum(vec==reci) }
  unlist(lapply(listrec,check))
}

# Create list Y of binary matrices of receivers per sender
p <- nrow(enron.actors) #total number of nodes
n <- nrow(enron.messages) #total number of messages
ns <- rep(0,p) # number of messages of sender s
numrec <- rep(0,n) # number of receivers per message
Y <- Xendo <- lapply(1:p,function(i){NULL})
pb = txtProgressBar(min = 0, max = n, initial = 0)
for(i in 1:n){
  sender <- enron.messages[i,]$sender
  receiver <- unlist(enron.messages[i,]$receiver)
  numrec[i] <- length(receiver)
  receivers <- rep(0,p)
  receivers[sender] <- NA
  receivers[receiver] <- 1
  Y[[sender]] <- rbind(Y[[sender]],receivers)
  ns[sender] <- ns[sender] + 1
  row.names(Y[[sender]]) <- NULL
  
  time <- enron.messages[i,]$time
  
  Xendo[[sender]] <- rbind(Xendo[[sender]],t(matrix(unlist(lapply(1:p,function(r){
    timediff <- (time - enron.messages$time[1:(i-1)]) / 60
    timemin <- min(which(timediff < 21.3*24))

    if(r != sender){
      
	timediff1 <- timediff[timemin:(i-1)] < .5
	timediff2 <- (timediff[timemin:(i-1)] > .5) * (timediff[timemin:(i-1)] < 1.1) == 1
	timediff3 <- (timediff[timemin:(i-1)] > 1.1) * (timediff[timemin:(i-1)] < 2) == 1
      timediff4 <- (timediff[timemin:(i-1)] > 2) * (timediff[timemin:(i-1)] < 8) == 1
      timediff5 <- (timediff[timemin:(i-1)] > 8) * (timediff[timemin:(i-1)] < 1.3*24) == 1
      timediff6 <- (timediff[timemin:(i-1)] > 1.3*24) * (timediff[timemin:(i-1)] < 5.3*24) == 1
      timediff7 <- (timediff[timemin:(i-1)] > 5.3*24) * (timediff[timemin:(i-1)] < 21.3*24) == 1

      samesender <- enron.messages$sender[timemin:(i-1)]==sender
      samereceiver <- checkreceiver(listrec=enron.messages$receiver[timemin:(i-1)],r)
      inertia <- c(sum(samesender * samereceiver * timediff1),
                   sum(samesender * samereceiver * timediff2),
                   sum(samesender * samereceiver * timediff3),
                   sum(samesender * samereceiver * timediff4),
                   sum(samesender * samereceiver * timediff5),
                   sum(samesender * samereceiver * timediff6),
                   sum(samesender * samereceiver * timediff7))

      receiverissender <- enron.messages$sender[timemin:(i-1)] == r
      senderisreceiver <- checkreceiver(listrec=enron.messages$receiver[timemin:(i-1)],sender)
      reciprocity <- c(sum(receiverissender * senderisreceiver * timediff1),
                       sum(receiverissender * senderisreceiver * timediff2),
                       sum(receiverissender * senderisreceiver * timediff3),
                       sum(receiverissender * senderisreceiver * timediff4),
                       sum(receiverissender * senderisreceiver * timediff5),
                       sum(receiverissender * senderisreceiver * timediff6),
                       sum(receiverissender * senderisreceiver * timediff7))
      return(c(inertia,reciprocity))
    }else return(rep(NA,14))
  })),nrow=14)))

  setTxtProgressBar(pb,i)
}

enron.actors$long.department[which(is.na(enron.actors$long.department))] <- c("X","XX","XXX") #three missings result in FALSE dyadic cov for equal long.dept
enron.actors$trading <- 
  grepl("Trading",as.character(enron.actors$title)) + grepl("Trader",as.character(enron.actors$title))
# Let X[[s]] be a matrix of covariates of sender s and all receivers
Xexo <- lapply(1:p,function(s){
  
  sender <- c(enron.actors$department[s]=="Legal",
    enron.actors$department[s]=="Trading",
    enron.actors$seniority[s]=="Junior",
    enron.actors$gender[s]=="Female")
  
  Xs <- t(matrix(unlist(lapply(1:p,function(r){
    receiver <- c(enron.actors$department[r]=="Legal",
      enron.actors$department[r]=="Trading",
      enron.actors$seniority[r]=="Junior",
      enron.actors$gender[r]=="Female")
    c(sender %*% t(receiver))
  })),nrow=16))
  Xs <- cbind(Xs,as.integer(enron.actors$title[s] == enron.actors$title))
  Xs <- cbind(Xs,as.integer(enron.actors$long.department[s] == enron.actors$long.department))

  Xs[s,] <- NA
  return(Xs)
})
#combine Xendo and Xexo, and transform counts in Xendo using log(count+1)

X <- lapply(1:length(Xexo),function(s){
  if(ns[s]>0){
    cbind(log(log(Xendo[[s]]+1)+1),
      do.call(rbind,replicate(ns[s],Xexo[[s]],simplify=FALSE))
    )
  }
})

#standardize X
Xmean <- apply(do.call(rbind,X),2,mean,na.rm=TRUE)
Xsd <- apply(do.call(rbind,X),2,sd,na.rm=TRUE)
X <- lapply(1:p,function(s){
  if(ns[s]>0){
    Xs <- X[[s]] - rep(1,length=nrow(X[[s]]))%*%t(Xmean)
    Xs <- Xs %*% diag(1/Xsd)
    return(Xs)
  }
})


out0s <- Gibbs_sampler_MR30(Y,X,samsize=2e4,burnin=1e5,store=100,cores=1,whichrep=1:30,
	mean0c=0,shape0c=200,scale0c=0.5,sd2=.0004,shape0b=4,scale0b=4,startingvalues=starting)

out1s <- Gibbs_sampler_MR3(Y,X,q=1,samsize=2e4,burnin=1e5,store=100,cores=1,whichrep=1:30,
	mean0c=0,shape0c=200,scale0c=15,sd2=.007,shape0b=4,scale0b=4,startingvalues=starting)

out2s <- Gibbs_sampler_MR3(Y,X,q=2,samsize=2e4,burnin=1e5,store=100,cores=1,whichrep=1:30,
	mean0c=0,shape0c=200,scale0c=15,sd2=.007,shape0b=4,scale0b=4,startingvalues=starting)


