##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# Computing Assignment 2
# January 28, 2018
##--------------------------------------------------------

#-----------------Question 1------------------------------


#Binomial test function for all three cases

btest<-function( y,n,Ho,alternative=c("less.than", "greater.than", "two.sided") ){
  
  alternative <- match.arg(alternative)
  
  #condition for less than case
  if(alternative=="less.than"){
    
    #range for x from 0 to y (left than) and cumulative sum of probability 
    x<-0:y
    p_value<-sum(dbinom(x,n,Ho))
    
  } 
  
  #condition for greater than case
  else if(alternative=="greater.than"){
    
    #range for x from y to n (greater than) and cumulative sum of probability 
    x<-y:n
    p_value<-sum(dbinom(x,n,Ho))
    
  }
  
  #condition for two sided case
  else if (alternative=="two.sided"){
    
    #range for x from 0 to n (greater than) 
    x<-0:n
    
    #indicator function when the probability of x is less than y
    I<-dbinom(x,n,Ho)<=dbinom(y,n,Ho)
    
    #and cumulative sum of probability 
    p_value<-sum(dbinom(x,n,Ho)*I)
    
  }
  
  
  return(p_value)
}

#-----------------Question 2------------------------------

#functions for when n=12. y=2, pi=.3 in all three cases
btest(2,12,.3,"two.sided")
btest(2,12,.3,"greater.than")
btest(2,12,.3,"less.than")



#-----------------Question 3------------------------------
#------a 

#successes
y<-7

#sample size
n<-20

#estimate of p
p_hat<-y/n

#parameter values of the uninformative prior
shape_1<-1
shape_2<-1

#values of alpha and beta for the posterior distribution calculated in class
alpha<-y+shape_1
beta_b<-n-y+shape_2

# lower and upper quantiles
q1<-.025
q2<-.975

#inverse beta distribution function 
qbeta(q1,alpha,beta_b)
qbeta(q2,alpha,beta_b)

#expected value of the distribution
mean_beta<-(y+shape_1)/(n+shape_1+shape_2)


#------b
#successes
y<-7

#sample size
n<-20

#estimate of p
p_hat<-y/n

##parameter values of the informative prior
shape_1<-7
shape_2<-5

#values of alpha and beta for the posterior distribution calculated in class
alpha<-y+shape_1
beta_b<-n-y+shape_2

#inverse beta distribution function 
qbeta(q1,alpha,beta_b)
qbeta(q2,alpha,beta_b)

#expected value of the distribution
mean_beta<-(y+shape_1)/(n+shape_1+shape_2)


#------------------------Comments-----------------------------------

#uninformative prior

mean_beta
# .3636

qbeta(q1,alpha,beta_b)
# .1810716

qbeta(q2,alpha,beta_b)
#.5696755

#length of interval 
qbeta(q2,alpha,beta_b)-qbeta(q1,alpha,beta_b)
# 0.3886039


#informative prior

mean_beta
#0.4375

qbeta(q1,alpha,beta_b)
# 0.273165

qbeta(q2,alpha,beta_b)
# 0.6092408

#length of interval 
qbeta(q2,alpha,beta_b)-qbeta(q1,alpha,beta_b)
#0.3360758

# the uninformative prior has a better point estimate but a  slightly  longer interval while the informative prior has a 
# smaller interval but further point estimate 


#-----------------Question 4------------------------------

#sequence of probability ranging from 0,1 by .00001 increments  
min_p_o<-0
max_p_o<-1
delta<-.00001
p_vec<-seq(min_p_o,max_p_o,by=delta)

#successes
y<-7

#sample size
n<-20

#estimate of p
p_hat<-y/n

#lower quantile and upper quantile
q1<-.025
q2<-.975

#empty vector
confidence_vec<-c()

# Z score statistic W from Wald's test
z_w<-(p_hat-p_vec)/(sqrt(p_hat*(1-p_hat)/n))

#vector of true false logic where w is greater than the lower quantile and w is less than the upper quantile 
tf<-z_w>qnorm(q1,0,1) & z_w<qnorm(q2)

#The confidence limits will be the min and max of the product of these two vectors p_vec and tf
confidence_limits<-p_vec*tf

#lower bound of the confidence interval 
min(confidence_limits[confidence_limits>0])

#upper bound of the confidence interval
max(confidence_limits[confidence_limits>0])


#from the standard equation, we verify that they are the same 
alpha_level<-.05
z_star<-qnorm(1-alpha_level/2)

#lower then upper bound of the interval 
p_hat- z_star* sqrt(p_hat*(1-p_hat)/n)
p_hat+ z_star* sqrt(p_hat*(1-p_hat)/n)




