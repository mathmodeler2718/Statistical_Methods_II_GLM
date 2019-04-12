
##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# Computing Assignment 6
# February 25, 2019
##--------------------------------------------------------

#-----------------Question 1------------------------------

#-----------------
# (a)
# the purpose is to estimate $\theta with grid search for our sample 
# Y=c(1,3,4,6,2,5,9) where m=15 for all subjects and y_i~ Binomial(m,$\theta)

# vector of Y
Y<- c(1,3,4,6,2,5,9)

# trials m
m<-15

#sequence of probability ranging from 0,1 by .001 increments  
delta<-.001
min_theta<-0+delta
max_theta<-1-delta

# vector of sequential theta values by .001 increments
theta_vec<-seq(min_theta,max_theta,by=delta)

# log likelihood function that should be maximized 
log_L_theta_x<-log( sum( choose(m,Y) ) ) + sum(Y)*log(theta_vec) + sum(m-Y)*log(1-theta_vec) 

# visualization of the likelihood
#plot(theta_vec, exp(log_L_theta_x) ) 

# visualization of the log likelihood 
#plot(theta_vec,log_L_theta_x)

# maximum y_value on the graph
y_val<-max(log_L_theta_x)

# x value or the mle
mle<-theta_vec[which.max(log_L_theta_x)]
mle

#------------
# (b)
# the purpose is to approximate $\theta using Newton-Raphson

#initial guess
theta<- .5

#empty vector
theta_vec<-c()

# initial guess is stored as the first element of the vector
theta_vec[1]<-theta

#number of iterations
k<-10

#tolerance
tol<-.001

# Newton-Raphson method to calculate root to optimize the log likelihood function
for (i in 2:k ){
  
  theta_vec[i]<- theta_vec[i-1]- sum(Y/theta_vec[i-1] - (m-Y)/(1-theta_vec[i-1]) ) / sum(-Y/theta_vec[i-1]^2 - (m-Y)/(1-theta_vec[i-1])^2) 
  
  if ( abs(theta_vec[i-1]-theta_vec[i]) < tol) {
    
    root.approx <- tail(theta_vec[i], n=1)
    res <- list('root approximation' = root.approx, 'iterations' = k)
  }
  
}

# root approximation answer
res


#-----------------Question 2------------------------------

#Fisher's Exact Test 

# power function
power_fish.test<-function(p_1,p_2,n1,n2,alpha,iterations){
  
  # Fisher's exact test function
  fish.test<-function( y1,n1,y2,n2,alternative=c("p.1.less.than.p.2", "p.1.greater.than.p.2", "two.sided") ){
    
    alternative <- match.arg(alternative)
    
    #min between max of 0 and n1+(y1+y2)-(n1+m2)
    min_n11<-max(0,n1 + (y1+y2) - (n1+n2))
    
    #max between min of n1, total successes
    max_n11<- min(n1,(y1+y2))
    
    #condition for pi_1 < pi_2
    if(alternative=="p.1.less.than.p.2"){
      
      #range for x from min to y1 (less than)  and cumulative sum of probability 
      x<-min_n11:y1
      p_value<-sum(dhyper(x,n1,n2,y1+y2))
      
    } 
    
    #condition for pi_1 > pi_2
    else if(alternative=="p.1.greater.than.p.2"){
      
      #range for x from y1 to max (greater than) and cumulative sum of probability 
      x<-y1:max_n11
      p_value<- sum(dhyper(x,n1,n2,y1+y2))
      
    }
    
    #condition for two sided case
    else if (alternative=="two.sided"){
      
      #range for x from min to the max
      x1<-min_n11:max_n11
      
      
      
      #indicator function when the probability of x1 is less than y1
      I<-dhyper(x1,n1,n2,y1+y2)<=dhyper(y1,n1,n2,y1+y2)
      
      #cumulative sum of probability
      p_value<-sum(dhyper(x1,n1,n2,y1+y2)*I)
      
    }
    
    # returns the p-value of the function
    return(p_value)
  }
  
  # generating sample from rbinom for pi 1 and pi 2
  vec_1<-rbinom(iterations,n1,p_1)  
  vec_2<-rbinom(iterations,n2,p_2)
  
  # empty vector 
  v<-c()
  
  # loop to calculate p values to be stored in vector v
  for (i in 1:iterations){
    
    v[i]<-fish.test(vec_1[i], n1, vec_2[i], n2, alternative = "two.sided")
    
  }
  # emperical power 
  return(length(v[v<alpha])/iterations)
  
}

# power calculations for 3 examples

power_fish.test(.3, .2, 30, 30, .05, 1000)
power_fish.test(.3, .2, 30, 40, .05, 1000)
power_fish.test(.3, .2, 40, 30, .05, 1000)



















