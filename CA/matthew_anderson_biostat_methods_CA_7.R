##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# Hw #7 
# March 17, 2019
##--------------------------------------------------------
#-----------------Question 1------------------------------

# First 15 entries of  x and y 
i1<-1:15
y1<-rep(1,length(i1))
x1<-rep(1,length(i1))

#  entries 16 to 20 of  x and y 
i2<-16:20
y2<-rep(0,length(i2))
x2<-rep(1,length(i2))

#  entries 21 to 25 of  x and y 
i3<-21:25
y3<-rep(1,length(i3))
x3<-rep(0,length(i3))

# entries 26 to 40
i4<-26:40
y4<-rep(0,length(i4))
x4<-rep(0,length(i4))

# Response vector y
y<-c(y1,y2,y3,y4)

# data x
x<-c(x1,x2,x3,x4)

# Desing Matrix X
X<-cbind(rep(1,length(x)),x)



# a) A logit link function.

#b1 is inital guess, y is response vector, x is design matrix, tol is tolerence of convergence 
logit<-function(b1, y, x, tol, stop){
  
  counter <- 0
  
  # initialize beta 
  beta <- b1
  
  while(counter < stop){
                          # logit model 
    W<-diag(as.numeric(exp(X%*%beta)/(1+exp(X%*%beta))))
    
    # Covariance matrix 
    J<-t(X)%*%W%*%X
    
    # Z vector 
    z <- X%*%beta + solve(W)%*%(y-exp(X%*%beta)/(1+exp(X%*%beta)))
    
    # score equation
    beta_new<-solve(J)%*%t(X)%*%W%*%z
    
    # convergence of updates 
    b_norm<-sqrt(sum((beta_new - beta)^2))
    
    # when to stop if tol reached 
  if ( b_norm < tol ){
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
logit_estimate<-logit(b1 = c(1,1), y=y, x=X, tol = 10^-6, stop=1000)

#beta estimates
logit_estimate$beta

# covaraince matrix
logit_estimate$Covariance

# b) A probit link function

#b2 is inital guess, y is response vector, x is design matrix, tol is tolerence of convergence 
probit<-function(b2, y, x, tol, stop){
  
  counter <- 0
  
  # initialize beta 
  beta <- b2
  
  while(counter < stop){
                                    # Probit Model
    W<-diag( as.numeric(( 1/ ( (pnorm(X%*%beta)*(1-pnorm(X%*%beta)) )    ))*(dnorm(X%*%beta)^2)  )) 
    
    # Covariance matrix 
    J<-t(X)%*%W%*%X
    
    # Z vector
    z <- X%*%beta + solve(W)%*%(y-pnorm(X%*%beta) )
    
    # score equation
    beta_new<-solve(J)%*%t(X)%*%W%*%z
    
    # convergence of updates
    b_norm<-sqrt(sum((beta_new - beta)^2))
    
    # when to stop if tol reached 
    if ( b_norm < tol ){
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
probit_estimate<-probit(b2 = c(1,1), y=y, x=X, tol = 10^-6, stop=1000)

# beta estimate
probit_estimate$beta

#Covariance 
probit_estimate$Covariance


# c) Compare parameter estimates the results from these two models

logit_estimate$beta
probit_estimate$beta
# same direction in slope and intecept, the slope of the logit is larger slightly 

logit_estimate$Covariance
probit_estimate$Covariance
# varaince matrix is very similar between the two 