##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# CA #8
# March 24, 2019
##--------------------------------------------------------

#-----------------Question 3------------------------------

# y from p=.4 distribution (group 1) 1000 iterations
g1<-replicate(1000,rbinom(50,1,.4))

# y from p=.25 distribution (group 2) 1000 iterations
g2<-replicate(1000,rbinom(50,1,.25))

# combined into response vector y
y<-rbind(g1,g2)

# Design matrix. First collumn all 1's second collumn group factor 1,0
X<-cbind(    rep(1,nrow(y)), c(rep(1,nrow(g1)), rep(0,nrow(g2)))    )  


#####################  logit link function

# b_1 is inital guess, y is response vector, x is design matrix, tol is tolerence of convergence 
logit<-function(b_1, y, x, tol, stop){
  
  # initializing counter at 0
  counter <- 0
  
  # initialize beta 
  beta <- b_1
  
  while(counter < stop){
    
    # logit model mean
    pi<-exp(X%*%beta)/(1+exp(X%*%beta))
    
    # logit model variance 
    pi2<-exp(X%*%beta)/  ((1+exp(X%*%beta))^2)
    
    # score vector U
    U<- t(X)%*%( y-pi)

    
    # Covariance matrix 
    J11<- (t(X[,1])^2%*%(exp(X%*%beta)/(1 + exp(X%*%beta))^2))
    J12<-sum(X[,1]*X[,2] * (exp(X%*%beta)/(1 + exp(X%*%beta))^2))
    J21<-sum(X[,1]*X[,2] * (exp(X%*%beta)/(1 + exp(X%*%beta))^2))
    J22<-t(X[,2])^2 %*%(exp(X%*%beta)/(1 + exp(X%*%beta))^2)
    J<-matrix(c(J11, J12,J21, J22), nrow=2, ncol = 2, byrow=T)
    
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
logit_estimate<-logit(b_1 = c(0.5,0.5), y=y[,10], x=X, tol = 10^-6, stop=1000) 

# beta estimates
logit_estimate$beta

# covaraince matrix
logit_estimate$Covariance

# applying this estimate to every random vector of y simulated in the 10000 simulations 
results <- apply(y,2,logit, b_1=c(0.5,0.5), x=X,tol=10^-6,stop=1000)

# extracting the estimate of b1
b1<-sapply(lapply(results,"[[",1),"[",2)

# covaraince matrix entry [2,2] extracted and stored into a vector 
cov1 <- sapply(lapply(results, "[[",2), "[",2,2)

# confidence level 
alpha<-.05

# True Odds Ratio parameter (.4/.6 )/ (.25/ .75) = 2
OR<- 2 

# confidence interval upper and lower bound 
ci_upper<- exp(b1 +qnorm((1-alpha/2))*sqrt(cov1) ) 
ci_lower<- exp(b1- qnorm((1-alpha/2))*sqrt((cov1) ) )

# Confidence Intervals (first 6 printed)
head(cbind(ci_lower,ci_upper))

# Width Confidence Interval for Odds Ratio
width_OR<-mean(ci_upper-ci_lower)


# counting coverage probablity 
coverage_probability<-c()

for (i in 1:length(ci_lower)){
  if (OR>ci_lower[i] & OR< ci_upper[i]){
    coverage_probability[i]<-1
  }
  else { coverage_probability[i]<-0}
}



# counting power
power_count<-c()

for (i in 1:length(ci_lower)){
  if (ci_lower[i] >1 |  ci_upper[i]<1 )  {
    power_count[i]<-1
    }
 else { power_count[i]<-0}
}

# Summaries

# average interval width
width_OR

# coverage probablity 
sum(coverage_probability)/length(ci_lower)

# empirical power
sum(power_count)/length(ci_lower)










