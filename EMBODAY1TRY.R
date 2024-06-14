x<-5
y<-2

mean(x,y)

#trial 123


#trying for Day 2 at EMBO POPGEN

#Genetic Drift 

#new exercise =1

# N= population size
#fA Allel frequency 
#rbinom (n, size, frequecnyc)/size
# rbinom (1, 2*N, 0.5)/2*N
#n= generation 


N <- 100
fA <- 0.5

rbinom(1,2*N,0.5)/(2*N)

#frequency of next generation = 0.525


rbinom(1,2*N,0.525)/(2*N)

#frequency of next generation = 0.48

rbinom(1,2*N,0.48)/(2*N)

#frequency of next generation = 0.485

rbinom(1,2*N,0.485)/(2*N)

#frequency of next generation = 0.46


# I did for multiple generation with different frequencies got at the end of each answer
# now i will take the frequency constant but change the number of generation 

N <- 100
fA <- 0.5

rbinom(1,2*N,0.5)/(2*N)

#frequency of next generation = 0.525

rbinom(2,2*N,0.5)/(2*N)

#frequency of 2 generation = 0.535, 0.510

rbinom(3,2*N,0.5)/(2*N)

#frequency of 3 generation = 0.505, 0.495, 0.460

rbinom(4,2*N,0.5)/(2*N)

#frequency of 4 generation = 0.490 0.520 0.555 0.480

rbinom(100,2*N,0.5)/(2*N)
#frequency of 4 generation = 0.485 0.495 0.570 0.505 0.510 0.475 0.505 0.460 0.485 0.490 0.530 0.500 0.545 0.505 0.435 0.530 0.520 0.475 0.570 0.565 0.445 0.450 0.475 0.540 0.490 0.500 0.490 0.430 0.510 0.515 0.555 0.500 0.480 0.565 0.545 0.480 0.450 0.490 0.535 0.520 0.500 0.520 0.485 0.530 0.525 0.495 0.450 0.520 0.570 0.485 0.520 0.480 0.500 0.495 0.505 0.510 0.485 0.445 0.475 0.435 0.515 0.510 0.525 0.460 0.490 0.545 0.550 0.520 0.505 0.485 0.530 0.500 0.540 0.480 0.520 0.485 0.495 0.460 0.510 0.415 0.515 0.450 0.540 0.530 0.475

frequency<- rbinom(100,2*N,0.5)/(2*N)

print(frequency)
plot(frequency)

# new exercise=2

fA <-rep(NA,100)

# at t=0

fA [1] <- 0.50
for (t in 1:99) fA [t+1] <- rbinom(1,2*N, fA [t])/(2*N)

plot (x=1:100, y=fA, type = "l", ylim = c(0,1), lwd=2)


# did this for multiple generation with frequencies from previous generation

# new exercise =3

#two different populations
#populations small=red
#population large= blue 

NR <- 60

fA <-rep(NA,100)

# at t=0

fA [1] <- 0.50
for (t in 1:99) fA [t+1] <- rbinom(1,2*NR, fA [t])/(2*NR)

plot (x=1:100, y=fA, type = "l", ylim = c(0,1), lwd=2)

NB <- 100000
fA <-rep(NA,100)

# at t=0

fA [1] <- 0.50
for (t in 1:99) fA [t+1] <- rbinom(1,2*NB, fA [t])/(2*NB)


plot (x=1:100, y=fA, type = "l", ylim = c(0,1), lwd=2)


#using velentias code

# large population
N = 500

# initial allele frequency for allele A
fA = 0.5

# number of generations to simulate
n_gen = 100

# table to store the allele frequency over time
large = data.frame(Generation = 1:n_gen, Allele_frequency = fA)

for(i in 1:n_gen){
  
  fA = rbinom(1, 2*N, fA) / (2*N)
  large[i, 2] = fA 
  
}

# small population
N = 10

# table to store the allele frequency over time
small = data.frame(Generation = 1:n_gen, Allele_frequency = fA)

for(i in 1:n_gen){
  
  fA = rbinom(1, 2*N, fA) / (2*N)
  small[i, 2] = fA 
  
}
small$Type = "small"
large$Type = "large"
data = rbind(large, small)

ggplot(data = data, aes(x = Generation, y = Allele_frequency, col = Type)) + geom_line() + geom_point() + theme_bw() + ylab("Allele frequency") + xlab("Generation") + ggtitle("Initial allele frequency of 0.5")

#calculations 
1/(1+1/2+1/3)

#example file by Matteo

import msprime

# Simulate 10 diploid samples under the coalescent with recombination on a 10kb region.
ts = msprime.sim_ancestry(
  samples=10,
  recombination_rate=1e-8, # as in humans
  sequence_length=10_000,
  population_size=10_000,
  random_seed=1234)

# we can add mutations
mutated_ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=1234)
print(mutated_ts.tables.sites)

for variant in mutated_ts.variants():
  print(variant)

#no didnt work had to do in python, so ran in embo serve, worked fine. 

#Day 3

install.packages("reticulate")

msprime <- reticulate::import("msprime")

#trying task 2 day 3 afternon session 

sims <- read.csv("mosquito-task2.csv", head=T)

# check prior distributions
x11()
par(mfrow=c(2,2))
hist(sims$N1)
hist(sims$N2)
hist(sims$T_split)
hist(sims$MigRate)

# remove simulations with NaN for some summary stats!

# find useful summary stats which correlate with T_split
cor(sims$Fst, sims$T_split)
cor(sims$dxy, sims$T_split)
cor(sims$segsites1, sims$T_split)
cor(sims$segsites2, sims$T_split)
cor(sims$pi1, sims$T_split)
cor(sims$pi2, sims$T_split)
cor(sims$tajima1, sims$T_split)
cor(sims$tajima2, sims$T_split)

# load observed summary stats
obs <- read.csv("mosquito-observed.csv", head=T)

# check if simulated retained summary stats contain the observed one
quantile(sims$Fst); cat(obs$Fst)
quantile(sims$segsites1); cat(obs$segsites1)
quantile(sims$segsites2); cat(obs$segsites2)

# merge obs with retained sims to scale them
sumstats <- scale(rbind(obs[c(1,3,4)],sims[,c(5,7,8)]))

library(abc)

est <- abc(target=sumstats[1,], param=sims$T_split, sumstat=sumstats[-1,], tol=0.05, method="rejection")

# check distances in the acceptance region
hist(est$dist)
abline(v=max(est$dist[which(est$region)]), lty=2)

# posterior
x11()
par(mfrow=c(2,1))
hist(est$unadj.values, freq=FALSE, xlim=range(sims$T_split), col=rgb(0,0,1,1/4), main="Posterior probability", xlab="Split time")
# MAP
map <- mean(est$unadj.values)
abline(v=map, lty=2)

# confidence intervals
hpd <- quantile(x=est$unadj.values, probs=c(0.025,0.975))
abline(v=hpd, lty=3)

# prior
hist(sims$T_split, freq=FALSE, xlim=range(sims$T_split), col=rgb(1,0,0,1/4))



#day 4

library(relater)
library(ggplot2)

install.packages(relater)

test1

