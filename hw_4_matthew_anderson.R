##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# HW #4
# March 28, 2019
##--------------------------------------------------------

#-----------------Question 2------------------------------
# a) 

# response vector y
y<-c( rep(1,334), rep(0,227), rep(1,309), rep(0,279)  )

# collumns of design matrix
x_1<-rep(1, sum(334,227,309,279) )
x_2<-c(rep(1,sum(334,227)), rep(0,sum(309,279)))

# Desing matrix, X
X<-cbind( x_1, x_2  )


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
    
    
    # Covariance matrix components 
    J11<- (t(X[,1])^2%*%(exp(X%*%beta)/(1 + exp(X%*%beta))^2))
    J12<-sum(X[,1]*X[,2] * (exp(X%*%beta)/(1 + exp(X%*%beta))^2))
    J21<-sum(X[,1]*X[,2] * (exp(X%*%beta)/(1 + exp(X%*%beta))^2))
    J22<-t(X[,2])^2 %*%(exp(X%*%beta)/(1 + exp(X%*%beta))^2)
    
    # Covariance matrix, J
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
  
  # convergence greater than stop 
  if (counter == stop){
    print('Finish')
  }
  
  # output of the function, estimates, and covariance matrix 
  output<-list('beta' = beta, 'Covariance' = solve(J))
  

  return(output)
  
}




# function with initial guess, y, x, tol and stop size 
logit_estimate<-logit(b_1 = c(1,1), y=y, x=X, tol = 10^-6, stop=1000) 



# beta estimates
logit_estimate$beta

# covaraince matrix
logit_estimate$Covariance

# Standard error 
diag(sqrt(logit_estimate$Covariance))

# CI
lb_0<-exp(logit_estimate$beta[1,]-qnorm(1-.05/2)*sqrt(logit_estimate$Covariance[1,1]))
ub_0<-exp(logit_estimate$beta[1,]+qnorm(1-.05/2)*sqrt(logit_estimate$Covariance[1,1]))

lb_1<-exp(logit_estimate$beta[2,]-qnorm(1-.05/2)*sqrt(logit_estimate$Covariance[2,2]))
ub_1<-exp(logit_estimate$beta[2,]+qnorm(1-.05/2)*sqrt(logit_estimate$Covariance[2,2]))

# CI intercept
c(lb_0,ub_0)

# CI beta_1
c(lb_1, ub_1)


# b)

# vector of mean estimate pi_hat
pi_hat<- exp(X%*%as.vector(logit_estimate$beta))/(1+exp(X%*%as.vector(logit_estimate$beta)))

# binding y and pi_hat so I can use the values of Y as an index 
mat<-cbind(y,pi_hat)

# Vector of pi_hat were y==0
c1<-mat[which(y==0),2]

# Vector of pi_hat were y==1
c2<-mat[which(y==1),2]

# Deviance test statistic
sum(2*sum(log(1/(1-c1))), 2*sum(log(1/c2)))

# Chi-squared test goodness of fit p-value
pchisq(sum(2*sum(log(1/(1-c1))), 2*sum(log(1/c2))), length(y)-1,lower.tail = FALSE)


#-----------------Question 3------------------------------
# clear global environment 
rm(list=ls())


# data

# combined  response for Rural and Urban in  vector y
y<-c(65,2,65,5,52,4,310,36,98,7,159,10,175,22,877,102,41,5,117,7,137,16,477,63,11,0,35,6,39,8,167,33)

# number of trials for each response
m<-c(317,20,476,33,486,40,3259,316,486,31,1004,81,1355,122,7660,724,223,18,539,39,697,68,3442,344,40,3,148,16,214,25,1019,114)

# factors type, age, urban rural (ur) treating them as numeric for this problem
type<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4)
age<-c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4)


# Design matrix.First collumn is 1's for the intercept, then factors:type, age, ur
X<-cbind(rep(1,length(y)), type, age  )  


####################################### 
# a) logistic binomial regression

logit_binomial<-function(b_1, y, x, tol, stop){
  
  # initializing counter at 0
  counter <- 0
  
  # initialize beta 
  beta <- b_1
  
  while(counter < stop){
    
    # logit model mean
    pi<-m* exp(X%*%beta)/(1+exp(X%*%beta))
    
    # logit model variance 
    pi2<- m* exp(X%*%beta)/  ((1+exp(X%*%beta))^2)
    
    # score vector U
    U<- t(X)%*%( y-pi)
    
    
    # Covariance matrix 
    j11<-sum(X[,1]*X[,1] * pi2  )
    j12<-sum(X[,1]*X[,2] * pi2  )
    j13<-sum(X[,1]*X[,3] * pi2  )
    j21<-sum(X[,2]*X[,1] * pi2  )
    j22<-sum(X[,2]*X[,2] * pi2  )
    j23<-sum(X[,2]*X[,3] * pi2  )
    j31<-sum(X[,3]*X[,1] * pi2  )
    j32<-sum(X[,3]*X[,2] * pi2  )
    j33<-sum(X[,3]*X[,3] * pi2  )
    
    J<-matrix(c(j11,j12,j13,j21,j22,j23,j31,j32,j33), nrow=3, ncol = 3, byrow=T)
    
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
logit_binomial<-logit_binomial(b_1 = c(.1,.1,.1), y=y, x=X, tol = 10^-5, stop=1000) 

# beta estimates
logit_binomial$beta

# covaraince matrix
logit_binomial$Covariance


################################# 
# b) no offset term log link function for poisson linear model 

# b_1 is inital guess, y is response vector, x is design matrix, tol is tolerence of convergence 
log_poisson<-function(b_1, y, x, tol, stop){
  
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
    j11<-sum(X[,1]*X[,1] * (exp(X%*%beta) ) )
    j12<-sum(X[,1]*X[,2] * (exp(X%*%beta) ) )
    j13<-sum(X[,1]*X[,3] * (exp(X%*%beta) ) )
    j21<-sum(X[,2]*X[,1] * (exp(X%*%beta) ) )
    j22<-sum(X[,2]*X[,2] * (exp(X%*%beta) ) )
    j23<-sum(X[,2]*X[,3] * (exp(X%*%beta) ) )
    j31<-sum(X[,3]*X[,1] * (exp(X%*%beta) ) )
    j32<-sum(X[,3]*X[,2] * (exp(X%*%beta) ) )
    j33<-sum(X[,3]*X[,3] * (exp(X%*%beta) ) )
    
    J<-matrix(c(j11,j12,j13,j21,j22,j23,j31,j32,j33), nrow=3, ncol = 3, byrow=T)
    
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
log_poisson<-log_poisson(b_1 = c(2,2,2), y=y, x=X, tol = 10^-6, stop=1000) 

# beta estimates
log_poisson$beta

# covaraince matrix
log_poisson$Covariance

###############################################
# c) offset term included poisson 


# b_1 is inital guess, y is response vector, x is design matrix, tol is tolerence of convergence 
log_poisson_off<-function(b_1, y, x, tol, stop){
  
  # initializing counter at 0
  counter <- 0
  
  # initialize beta 
  beta <- b_1
  
  while(counter < stop){
    
    # log poisson model mean
    pi<- m*exp(X%*%beta)
    
    # log poisson model variance 
    pi2<-m*exp(X%*%beta)
    
    # score vector U
    U<- t(X)%*%( y-pi)
    
    
    # Covariance matrix 
    j11<-sum(X[,1]*X[,1] * m*(exp(X%*%beta) ) )
    j12<-sum(X[,1]*X[,2] * m*(exp(X%*%beta) ) )
    j13<-sum(X[,1]*X[,3] * m*(exp(X%*%beta) ) )
    j21<-sum(X[,2]*X[,1] * m*(exp(X%*%beta) ) )
    j22<-sum(X[,2]*X[,2] * m*(exp(X%*%beta) ) )
    j23<-sum(X[,2]*X[,3] * m*(exp(X%*%beta) ) )
    j31<-sum(X[,3]*X[,1] * m*(exp(X%*%beta) ) )
    j32<-sum(X[,3]*X[,2] * m*(exp(X%*%beta) ) )
    j33<-sum(X[,3]*X[,3] * m*(exp(X%*%beta) ) )
    
    J<-matrix(c(j11,j12,j13,j21,j22,j23,j31,j32,j33), nrow=3, ncol = 3, byrow=T)
    
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
log_poisson_off<-log_poisson_off(b_1 = c(2,2,2), y=y, x=X, tol = 10^-6, stop=1000) 

# beta estimates
log_poisson_off$beta

# covaraince matrix
log_poisson_off$Covariance








