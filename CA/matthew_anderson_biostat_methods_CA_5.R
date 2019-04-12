##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# Computing Assignment 5
# Febuary 18, 2019
##--------------------------------------------------------

#-----------------Question 1------------------------------


#Fisher's Exact Test 

# function fish.test takes input arguments, y1, n1, y2, n2 and specifies an alternative hypotheisis
fish.test<-function( y1,n1,y2,n2,alternative=c("p.1.less.than.p.2", "p.1.greater.than.p.2", "two.sided") ){
  
  alternative <- match.arg(alternative)
  
  #min between max of 0 and n1+(y1+y2)-(n1+m2)
  min_n11<-max(0,n1 + (y1+y2) - (n1+n2))
  
  #max between min of n1, total sucesses
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


# first calculation
fish.test(2,10,7,13,alternative = "p.1.less.than.p.2")

# second calculation
fish.test(6,11,2,8,alternative = "p.1.greater.than.p.2")

# third calculation
fish.test(1,9,5,10,alternative = "p.1.less.than.p.2")












