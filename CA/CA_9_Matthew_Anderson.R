##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# CA #9
# April 1, 2019
##--------------------------------------------------------

# clear environment
rm(list = ls())


# y from lam=4 distribution (group 1) 1000 iterations
g1<-replicate(1000,rpois(50,4))

# y from lam=3.2 distribution (group 2) 1000 iterations
g2<-replicate(1000,rpois(50,3.2))

# combined into response vector y
y<-as.matrix((rbind(g1,g2)))

# Design matrix X 
X<-cbind(    rep(1,nrow(y)), c(rep(1,nrow(g1)), rep(0,nrow(g2)))    ) 


# b_1 is inital guess, y is response vector, x is design matrix, tol is tolerence of convergence 
risk<-function(b_1, y, x, tol, stop){
  
  # initializing counter at 0
  counter <- 0
  
  # initialize beta 
  beta <- b_1
  
  while(counter < stop){
    
    # log poisson model mean
    pi<-exp(X%*%beta)
    
    # log poisson model variance 
    pi2<-exp(X%*%beta)
    
    # score vector U
    U<- t(X)%*%( y-pi)
    
    
    # Covariance matrix 
    j11<-sum(X[,1]*X[,1] * pi )
    j12<-sum(X[,1]*X[,2] * pi )
 
    j21<-sum(X[,2]*X[,1] * pi )
    j22<-sum(X[,2]*X[,2] * pi )
   

    
    
    J<-matrix(c(j11,j12,j21,j22), nrow=2, ncol = 2, byrow=T)
    
    # score equation
    beta_new<-beta + solve(J)%*%U
    
    # convergence of updates 
    b_norm<-sqrt(sum((beta_new - beta)^2))
    
    # when to stop if tol reached 
    if ( b_norm <= tol ){
      break
    }
    
    #updating each estimate for beta 
    beta = beta_new
    counter = counter + 1
    
  }
  
  if (counter == stop){
    print('Finish')
  }
  
  # output of the function, estimates, and covariance matrix 
  output<-list('beta' = beta, 'Covariance' = solve(J))
  
  return(output)
  
}

# function with initial guess, y, x, tol and stop size 
estimate<-risk(b_1 = c(2,2), y=y[,1], x=X, tol = 10^-6, stop=1000) 

# beta estimates
estimate$beta

# covaraince matrix
estimate$Covariance

# applying this estimate to every random vector of y simulated in the 10000 simulations 
results <- apply(y,2,risk, b_1=c(1,1), x=X,tol=10^-6,stop=1000)

# extracting the estimate of b1
b1<-sapply(lapply(results,"[[",1),"[",2)

# covaraince matrix entry [2,2] extracted and stored into a vector 
cov1 <- sapply(lapply(results, "[[",2), "[",2,2)

# confidence level 
alpha<-.05  

# confidence interval upper and lower bound 
ci_upper<- exp(b1 +qnorm((1-alpha/2))*sqrt(cov1) ) 
ci_lower<- exp(b1- qnorm((1-alpha/2))*sqrt((cov1) ) )

# Confidence Intervals (first 6 printed)
head(cbind(ci_lower,ci_upper))  

# true relative risk pi_1/pi_2
true_RR<-rep(4/3.2, length(b1))

# estimate of RR from GLM 
RR_hat<-exp(b1)


### Summaries

# Average Bias
sum(RR_hat-true_RR)/length(b1)

# MSE
sum((RR_hat-true_RR)^2)/length(b1)

# Mean Absolute Error 
sum(abs(RR_hat-true_RR))/ length(b1)
  
# average interval width for RR
mean(ci_upper-ci_lower)

# coverage probablity pretty good 
length(which(true_RR[1]>ci_lower & true_RR[1]< ci_upper))/ncol(y)

# empirical power not the best 
length(which(ci_lower>1 |ci_upper<1))/ncol(y)
  
  
  
  
  
  
  
