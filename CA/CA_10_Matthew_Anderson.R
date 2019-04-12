##--------------------------------------------------------
# Matthew Anderson
# Biostat Methods II
# CA #10
# April 12, 2019
##--------------------------------------------------------

# vector of Relapse where 1=relpase, 0=no relapse
d<-c(1,0,1,1,1,0,1,0,1,1,1,0,0,1,1,1)

# vector of number of subjects at risk after t(_{j-1})
r<-c(16:1)

#cumulative product 
s<-cumprod(1-d/r)

# survival plot 
plot(r,rev(s))

# r, d, and s at each time, t_i
cbind(r,d,s)



























