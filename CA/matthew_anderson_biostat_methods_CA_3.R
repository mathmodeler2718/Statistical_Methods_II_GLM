##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# CA #3
# Febuary 3, 2019
##--------------------------------------------------------

#-----------------Question 1------------------------------

#score interval function 
score_interval<-function(n,y,alpha_level){
  
  #sequence of probability ranging from 0,1 by .00001 increments  
  delta<-.001
  min_p_o<-0+delta
  max_p_o<-1-delta
  
  p_vec<-seq(min_p_o,max_p_o,by=delta)
  
  #estimate of p
  p_hat<-y/n
  
  #lower quantile and upper quantile
  q1<-alpha_level/2
  q2<- 1-(alpha_level/2)
  
  #empty vector
  confidence_vec<-c()
  
  # Z score statistic from score test 
  z_score<-(p_hat-p_vec)/(sqrt(p_vec*(1-p_vec)/n))
  
  #vector of true false logic where is greater than the lower quantile and w is less than the upper quantile 
  tf<-(z_score>qnorm(q1,0,1) & z_score<qnorm(q2))
  
  #The confidence limits will be the min and max of the product of these two vectors p_vec and tf
  confidence_limits<-p_vec*tf
  
  #lower bound of the confidence interval 
  min(confidence_limits[confidence_limits>0])
  
  #upper bound of the confidence interval
  max(confidence_limits[confidence_limits>0])
  
  
  return(c(min(confidence_limits[confidence_limits>0]),max(confidence_limits[confidence_limits>0])))
  
}

#score interval output 
score_interval(30,18,.05)





# wald interval function
wald_interval<-function(n,y,alpha_level){
  
  #sequence of probability ranging from 0,1 by .00001 increments  
  delta<-.00001
  min_p_o<-0+delta
  max_p_o<-1-delta
  delta<-.00001
  p_vec<-seq(min_p_o,max_p_o,by=delta)
  
  #estimate of p
  p_hat<-y/n
  
  #lower quantile and upper quantile
  q1<-alpha_level/2
  q2<- 1-(alpha_level/2)
  
  #empty vector
  confidence_vec<-c()
  
  # Z wald statistic from wald test 
  z_wald<-(p_hat-p_vec)/(sqrt(p_hat*(1-p_hat)/n))
  
  #vector of true false logic where w is greater than the lower quantile and w is less than the upper quantile 
  tf<-z_wald>qnorm(q1,0,1) & z_wald<qnorm(q2)
  
  #The confidence limits will be the min and max of the product of these two vectors p_vec and tf
  confidence_limits<-p_vec*tf
  
  #lower bound of the confidence interval 
  min(confidence_limits[confidence_limits>0])
  
  #upper bound of the confidence interval
  max(confidence_limits[confidence_limits>0])
  
  
  return(c(min(confidence_limits[confidence_limits>0]),max(confidence_limits[confidence_limits>0])))
  
}

# wald interval output 
wald_interval(30,18,.05)


#-----------------Question 2------------------------------

#initialization  iterations(1000), sample size(30), parameter PI (.3), alpha level (.05)
iterations<-1000
n<-30
PI<-.3
alpha_level<-.05

# this is generating random success from a binomial distribution, binomial(30, .3)
# y will be a vector of 1000 that samples from this distribution
y<-rbinom(iterations,n,PI)



# score interval

#empty vectors
s<-c()
length_s<-c()
tf_score<-c()

# this loop is going to make a list and store a CI in each list, then I will assign a true or false to each interval that contains .3
for (i in 1:iterations){
  s[[i]]<-score_interval(n,y[i],alpha_level)
  length_s[i]<-(s[[i]][2]-s[[i]][1])
  tf_score[i]<- (s[[i]][1] <=PI & s[[i]][2] >= PI )
  }

#adding up all of the true intervals that contain .3 and dividing by 1000 will give us the coverage probability
coverage_probability_score<-sum(tf_score)/iterations

#here I am averaging the length of each interval 
average_length_score<-mean(length_s)


# wald interval 

w<-c()
length_w<-c()
tf_wald<-c()

# this loop is going to make a list and store a CI in each list, then I will assign a true or false to each interval that contains .3
for (i in 1:iterations){
  w[[i]]<-wald_interval(n,y[i],alpha_level)
  length_w[i]<-(w[[i]][2]-w[[i]][1])
  tf_wald[i]<- (w[[i]][1] <=PI & w[[i]][2] >= PI )
}

#adding up all of the true intervals that contain .3 and dividing by 1000 will give us the coverage probabilty
coverage_probability_wald<-sum(tf_wald)/iterations

#here I am averaging the length of each interval 
average_length_wald<-mean(length_w)

# results of experiment 
coverage_probability_score
coverage_probability_wald
average_length_score
average_length_wald






