##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# Exam 1
# March 6, 2019
##--------------------------------------------------------

#-----------------Question 1------------------------------

#-------------- score interval ------------------------------

# data
data<-c(10,12,8,15,21,11,14,6,7,22)

#sample_size
n<-length(data)

#alpha level 
alpha_level<-.05

# expanding grid to seach through all values of lamda starting at 2. I used 100 as an upper bound, thought this was reasonable for our data 
all_lamda<-2:100

#mle
lamda_hat<-mean(data)


#lower quantile and upper quantile 
q1<-alpha_level/2
q2<- 1-(alpha_level/2)

#-------------- score interval ------------------------------

#vector of lamda that satisify abs(z)<= z_(1-alpha/2)
bounds<-all_lamda[abs( (log(lamda_hat)) - (log(all_lamda)) ) / (sqrt(1/(n*all_lamda)))  <= qnorm(q2)]

# The 95% confidence interval for the log of lamda_hat for the score interval 
c(log(min(bounds)),log(max(bounds)))

#-------------- walds interval ------------------------------

#empty vectors before inversion 
bounds<-c()

#vector of lamda that satisify abs(z)<= z_(1-alpha-2), notice the denominator has lamda_hat instead of all the values of lamda 
bounds<-all_lamda[abs( (log(lamda_hat)) - (log(all_lamda)) ) / (sqrt(1/(n*lamda_hat)))  <= qnorm(q2)]

# The 95% confidence interval for the log of lamda_hat 
c(log(min(bounds)),log(max(bounds)))


#-----------------Question 2------------------------------

#Fisher's Exact Test 

# set seed for random number generator
set.seed(1234789)
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

# sample sizes in increments of 20 (for perfect sizes of n)
n<-seq(20,2000,20)

#empty vector to store power calcuations through all values of n
power<-c()

# takes a minute 
for (i in n){
power[i]<-power_fish.test(.2, .25, i, i, .05, 100)
}

# returns the sample size with the desired power
which(power==.8)


#-----------------Question 3------------------------------

# Response y
y<-c(25,32,26,19,22,29,31,25,23,18,22,15,26,31,27,27)

# sample size 
n<-length(y)

# data
x<-c(13,17,15,11,10,16,19,14,12,11,10,9,17,21,19,18)

# intecept of 1's
ones<-rep(1,length(x))

# desing matrix X
X<-cbind(ones,x)

# Slope Estimate $\beta^hat using conventional least squares estimator
beta_hat<-  sum( ( X[,2]-mean(X[,2]) ) * ( y-mean(y) )) / (sum( ( X[,2]-mean(X[,2]) )^2 ))



# Initial Guess
b_vec<-t(t(c(7,1)))

# Estimate of the nuisance parameter $/sigma^{2} the variance 
sig_hat<- sum( (y-t(X[,2])*beta_hat)/ (n-2) )

# elements of covaraince matrix 
t1<-sum(X[,1]^2) 
t2<-sum(X[,1]*X[,2])
t4<-sum(X[,2]^2) 

T_matrix<- (1/sig_hat) *matrix(c(t1,t2,t2,t4), nrow = 2)


#initial guess
BB<-b_vec
k<-20
# Scoring for Newton Raphson  U=t(cbind(sum( (y-X%*%b_vec) *X[,1] ) / sig_hat^2, sum( (y-X%*%b_vec) *X[,2] ) / sig_hat^2))
for (i in 1:k){
  
  BB<-BB+solve(T_matrix)%*%(t(cbind(sum( (y-X%*%BB) *X[,1] ) / sig_hat^2,sum( (y-X%*%BB) *X[,2] ) / sig_hat^2)))
}

# estimates for Beta's
BB

#covariance matrix 
T_matrix


#-----------------Question 4 extra credit------------------------------

stirling<-function(n,k){
  1/(sqrt(2*pi*k)) * (n*exp(1)/k)^k
}

# 61 digets of 9 choose 5

stirling(9999999999999999999999999999999999999999999999999999999999999,5)

