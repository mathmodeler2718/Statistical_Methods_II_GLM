##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# CA #10
# April 15, 2019
##--------------------------------------------------------

###------Question 1

# vector of relapse where 1=relapse, 0=no relapse
d<-c(1,0,1,1,1,0,1,0,1,1,1,0,0,1,1,1)

# vector of number of subjects at risk after t_{j-1}
r<-c(16:1)

#cumulative product 
s<-cumprod(1-d/r)

# survival plot 
plot(r,rev(s),"s")

# print table of r, d, and s at each time, t_i
cbind(r,d,s)

###------------- Question 2

# a)
# time vector
time<-c(10,11,11,15,19,21,24,28,29,29,30,31,32,32,35,36)

# our estimate of lamda, lamda_hat
lamda_hat<-sum(d)/sum(time)
lamda_hat

# the standard error of lamda_hat
se<-lamda_hat/sqrt(sum(d))
se

# b) 

# estimate of survival function at each t assuming exponential distribution
exp(-lamda_hat*time)

# answer from the product moment to compare 
s

# the estimates are quite different 
# constant hazard assumption is violated so we should not assume an exponential distribution  