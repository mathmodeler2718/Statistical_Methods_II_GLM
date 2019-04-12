##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# Hw 1
# January 29, 2018
##--------------------------------------------------------

#-----------------Question 3------------------------------

#Function for binomial CDF
# inputs arguments are a sample size (m), success probability (p), lower bound to accumulate probability (l), and an upper bound (u)

binom_cdf<-function(m,p,l,u){
  
  #empty vector to store individual probabilities 
  prob<-c()
  
  #range of values of y (if less than problem input 0 for l )
  y<-l:u
  
  #for loop to calculate and store probabilities according to the binomial pmf
  for (i in 1:length(y)){
    prob[i]<-choose(m, y[i])*p^y[i] * (1-p)^(m-y[i])
  }
  
  #returns cumulative sum 
  return(sum(prob))
  
}


#a since less than 0 is the lower bound 
binom_cdf(20,.5,0,4)

#b
binom_cdf(15,.7,9,15)

#c
binom_cdf(19,.17,12,14)


#-----------------Question 4------------------------------

# case a simulation 

# initializing confidence level (a), total iterations (1000), sample size per sample (30), pi null(pi_o)
a<-.05
m<-1000
n<-30
pi_o<-.1

#simulating from binomial distribution and pi_hat calculated
x<-rbinom(m,n,pi_o)
pi_hat<- x/n


# Wald test

#test statistic
Z_wald<- (pi_hat-pi_o) / (sqrt(pi_hat*(1-pi_hat)/n))

# power calculation (only keeping test statiscs  greater than alpha) equilivent to 1-beta
power_wald<-length(Z_wald[Z_wald>=qnorm(a)])/m



# Score test

#test statistic
Z_score<- (pi_hat-pi_o) / (sqrt(pi_o*(1-pi_o)/n))

# power calculation (only keeping test statiscs greater than alpha
power_score<-length(Z_score[ Z_score>=qnorm(a)])/m

#clear environment 
rm(list = ls())




####################################################################################################################

# case b simulation 

# initializing confidence level (a), total iterations (1000), sample size per sample (30), pi null(pi_o)
a<-.05
m<-1000
n<-30
pi_o<-.1

#simulating from binomial distribution and pi_hat calculated
x<-rbinom(m,n,pi_o)
pi_hat<- x/n

# Wald test

#test statistic
Z_wald<- (pi_hat-pi_o) / (sqrt(pi_hat*(1-pi_hat)/n))

#type 1 error calculation, counting how many times the test statistic is less than -1.95 
type_1_wald<-length(Z_wald[Z_wald<=qnorm(a)])/m



# Score test

#test statistic
Z_score<- (pi_hat-pi_o) / (sqrt(pi_o*(1-pi_o)/n))

#type 1 error calculation, counting how many times the test statistic is less than -1.95 
type_1_score<-length(Z_score[Z_score<=qnorm(a)])/m



#-----------------Question 5------------------------------

#Function for Poisson CDF
# inputs arguments are rate (lam), lower bound to accumulate probability (l), and an upper bound (u)

poisson_cdf<-function(lam,l,u){
  
  #empty vector to store individual probabilities 
  prob<-c()
  
  #range of values of y 
  y<-l:u
  
  #for loop to calculate and store probabilities according to the Poisson pmf
  for (i in 1:length(y)){
    prob[i]<-exp(-lam)*lam^y[i]/factorial(y[i])
  }
  
  #returns cumulative sum 
  return(sum(prob))
  
}


#a since less than 0 is the lower bound 
poisson_cdf(1,0,4)

#b
poisson_cdf(3,3,3)

#c
1-poisson_cdf(2,0,3)

#same values as accumulating dpois

