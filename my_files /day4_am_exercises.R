# Q2 2N=100, f0=10/100=0.1

# a) f1=f0=0.1
dbinom(x=10,size=100,prob=0.1)

# b)f1=0
dbinom(x=0,size=100,prob=0.1)

# c)f1=1
dbinom(x=100,size=100,prob=0.1)

# d)f1<f0
pbinom(9,size=100,prob=0.1)

# e)f1>f0
1-pbinom(10,size=100,prob=0.1)

# Q3 plot prob distribution of Plot the probability distribution on the allele 
#frequency for a current allele frequency f=0.1 and 2N=10, 2N=100, 2N=1000 
#and 2N=10000. What is the effect of population size?
par(mfrow=c(2,2))
for(twoN in c(10, 100, 1000, 10000)){
  p <- 0:twoN
  d <- dbinom(p, size = twoN, prob = 0.1)
  plot(p, d, type = 'b')
}

# Simulating neutral allele frequency trajectories

# Q4 write a function simulateWF() that simulates a Wright-Fisher allele 
# frequency trajectory. The function should take the following arguments: 
  # initial allele frequency (f), 
  # population size (2N), 
  # number of generations (G). 
# The function should return the allele frequency (between 0 and 1) for each generation as a vector.
simulateWF <- function (twoN, f, G) {
  p <- numeric(G+1)
  p[1] <- f
  for(i in 1:G){
    p[i+1] <- rbinom(1, size = twoN, prob = p[i])/twoN
  }
  return(p)
}

# Q5 use the function simulateWF() to simulate 1000 trajectories with:
    # 2N=100 and f=0.1 for G= 1000 
  # and plot the results
trajectories <- replicate(1000, simulateWF(twoN = 100, f = 0.1, G = 1000))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))

print(paste("Allele was lost in", sum(trajectories[1000,] == 0), "/", ncol(trajectories), "cases."))


# Q6 Use your function simulateWF() to study fixation probability of a new 
# mutation under different population sizes (G) (ensure G is large enough).
twoN=100
trajectories <- replicate(1000, simulateWF(twoN = 1000, f = 1/twoN, G = 1000))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))

print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))

# Simulating selection and genetic drift

# Q7 Write a function simulateWFWithSelection() that simulates allele trajectories under both genetic drift and 
# viability selection. Similar to your function simulateWF() you wrote above, it should take as input:
    #i) the population size 2N, 
    #ii) the initial allele frequency f, 
    #iii) the number G of generations to simulate and 
    #iv) also a vector v of viabilities for the genotypes AA, Aa and aa.
# It should then return the allele frequency (between 0 and 1) for each generation as a vector. 
# In each generation, your function should apply selection to alter the allele frequency,
# and then use binomial sampling to simulate genetic drift in that modified allele frequency.
simulateWFWithSelection <- function(twoN, f, G, v){
  p <- numeric(G + 1)
  p[1] <- f
  for(i in 1:G){
    # selection
    fA <- p[i]
    fa <- 1-fA
    fPrime <- (v[1]*fA*fA + v[2]*fA*fa)/(v[1]*fA*fA + v[2]*2*fA*fa + v[3]*fa*fa);
    # drift
    p[i+1] <- rbinom(1, size = twoN, prob = fPrime) / twoN
  }
  return(p)
}

# Q8 Use your function simulateWFWithSelection() to simulate 1000 trajectories with:
    # 2N=100,  
    #f=0.1  
    # genic selection with viabilities v=(1,1−s,(1−s)2) for G=1000 and plot them in one plot. 
# Start with s=0.01. 
s <- 0.01
trajectories <- replicate(100, simulateWFWithSelection(twoN = 1000, f = 0.1, G = 1000, v=c(1,1-s,(1-s)^2)))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))
print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))

# Q9 Set 2N=10^6 and dominant v=(1,1,1−s) with s=0.05.
s <- 0.05
trajectories <- replicate(100, simulateWFWithSelection(twoN = 10^6, f = 0.1, G = 1000, v=c(1,1,1-s)))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))
print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))
