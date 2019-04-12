##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# CA #4
# Febuary 11, 2019
##--------------------------------------------------------

#-----------------Question 1------------------------------
# a

#sample sizes and success counts
n1<-40
n2<-40
y1<-20
y2<-16

#mle estimate, p_hat
p1_hat<-y1/n1
p2_hat<-y2/n2

#sequence of probability ranging from 0,1 by .001 increments  
delta<-.001
min_p_o<-0+delta
max_p_o<-1-delta
p1_o<-seq(min_p_o,max_p_o,by=delta)
p2_o<-seq(min_p_o,max_p_o,by=delta)

#numerator of the LRT, the null density
lamda_num<- p1_o^y1 * (1-p1_o)^(n1-y1) * p1_o^y2 * (1-p1_o)^(n2-y2)

#denominator of the LRT, the MLE density
lamda_den<- p1_hat^y1 * (1-p1_hat)^(n1-y1) * p2_hat^y2 * (1-p2_hat)^(n2-y2)

#  lambda, ratio lamda_numerator and lamda_denominator 
lamda<-lamda_num/lamda_den

#test statistic
t<- -2*log(max(lamda))
t

#test statisc using built in dbinom function
-2*log(max(dbinom(y1,n1,p1_o)*dbinom(y2,n2,p2_o)) / max(dbinom(y1,n1,p1_hat)*dbinom(y2,n2,p2_hat))  )

#p value from chi squared distribution with 1 degree of freedom
pchisq(t,1,lower.tail = F)

#########################b 

# creating our table 
table<-matrix(c(y1,n1-y1,y2,n2-y2),2,2)
table

# table of empty values to stor expected counts (warning is ok on line 53)
expected_counts<-matrix(, nrow = 2, ncol = 2)

# loop for storing expected counts under the null hypothesis into the table, row total * column total / total total 
for ( i in 1:length(table[1,])   ){
  
  for ( j in 1:length(table[,1]) ){
    
    expected_counts[i,j]<-sum(table[i,])*sum(table[,j])/ sum(table)
  }
}

# table of our observed counts
observed_counts<-c(y1,n1-y1,y2, n2-y2)

# chi squred test statistic 
t2<- sum((observed_counts-expected_counts)^2 / expected_counts)
t2

#p value from chi squared distribution with 1 degree of freedom
pchisq(t2,1,lower.tail = F)

#-----------------Question 2------------------------------

#alpha level 
alpha_level<-.05

# expanding grid to seach through all values of pi_01 and pi_02
all_pairs<-expand.grid(p1_o,p2_o)

# values of p1_o and p2_o stored as individual vectors 
p1_o<-all_pairs[,1]
p2_o<-all_pairs[,2]

# vector of differences
difference<-p1_o-p2_o


#lower quantile and upper quantile 
q1<-alpha_level/2
q2<- 1-(alpha_level/2)

#empty vectors before loop
bounds<-c()

# Z score statistic from score test  takes bout 2 min ( still learning how to write more efficent code)
for (i in 1:length(all_pairs[,1])){
  if ( abs( (p1_hat-p2_hat) - (all_pairs[i,1] -all_pairs[i,2])) / (sqrt(all_pairs[i,1]*(1-all_pairs[i,1])/n1+all_pairs[i,2]*(1-all_pairs[i,2])/n2))  < qnorm(q2) )
  {
    bounds<-c(bounds,difference[i])
  }
}

# The 95% confidence interval 
c(min(bounds),max(bounds))












