##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# Compuing Assignmnet 1
# January 21, 2018
##--------------------------------------------------------


#-----------------Question 1------------------------------

#here i am initializing our n, 7, and our starting position for the factorial calculation fac=1
n<-7
fac<-1

#for loop to calculate the multiplication of each i from 1:7
for (i in 1:n){
  fac<-fac*i
}
fac


#-----------------Question 2------------------------------

#I am initializing the starting position at 1 and multiplying each additional i value to the total calculation using a for loop for a,b, c where c=a-b

combination<-function(a,b){
  
  fac_a<-1
  
  for (i in 1:a){
   fac_a<-fac_a*i
  }
  
  fac_b<-1
  
  for (i in 1:b){
    fac_b<-fac_b*i
  }

  
  fac_c<-1
  
  for (i in 1:(a-b)){
    fac_c<-fac_c*i
    
  }
  
  print(fac_a/(fac_b*fac_c))
  
}

#a
combination(10,3)

#b
combination(15,7)


#-----------------Question 3------------------------------

#mean calculation

#choose delta, min y, max y
delta<-.01
min_y<- -2
max_y<- 8
v<-(max_y-min_y)/delta

#empty vector a and and initializing  j at 0
a<-c()
j<-0


#this loop will make a vector a and each entry of a will have the area of the rectangle 
  for (i in seq(min_y,max_y,by=delta)){
  
 
    a[j]<- (i)*dnorm(i,3,sqrt(2))*delta
    j=j+1
    }

#summing each rectangle will give us the expected value of the random variable which is 3
sum(a)


# the variance is the same except x=(x-3)^2 since we want the variance and 3 is the expected value of x 

delta<-.01
min_y<- -2
max_y<- 8
v<-(max_y-min_y)/delta
a<-c()
j<-0


for (i in seq(min_y,max_y,by=delta)){
  
  
  a[j]<- (i-3)^2*dnorm(i,3,sqrt(2))*delta
  j=j+1
}

#summing will give us the variance which is 2
sum(a)


#-----------------Question 4------------------------------
#sample size is 30, 1 trial, and p is .4
n<-30
trials<-1
p<-.4

#random binomial function 
success<-rbinom(n,trials,p)

#total successes/n
sum(success)/length(success)





#-----------------Question 5------------------------------

#s, n, trials, p initialized
s<-10
n<-30
trials<-1
p<-.4

#p_hat is an empty vector to fill all 10 estimates of p
p_hat<-c()

#mat is an empty matrix to fill each row as an observation and each column as a sample 
mat <- matrix(, nrow = n, ncol = s)

#filing the matrix with random number generator from binomial distribution. Each column is the sample 
for(j in 1:s){
  mat[, j] <- rbinom(n,trials,p)
}

#calculating p_hat for each column by summing the column and dividing by 30
for (i in 1:s){
  p_hat[i]<-sum(mat[,i]/length(mat[,i]))
}

#vector of 10 p_hat
p_hat

#mean p_hat
mean(p_hat)

#standard error of p_hat
sd(p_hat)


