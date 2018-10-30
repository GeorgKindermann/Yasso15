library(Rcpp)
sourceCpp("useInR.cc")

theta <- c(0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-7.90E-05,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26)
init <- rep(0,5)
infall <- c(0.5,0.1,0.1,0.2,0) #0..Ethanol sol.(waxes), 1..Water sol.(sugar), 2..Acid soluble (Celluloses), 3..Insoluble (lignin), 4..Humus
time <- 1   #Time to run
avgT <- 10  #Temp annual average [C]
sumP <- 600 #Precip annual summ [mm]
ampT <- 12  #Amplitude (max. difference of month averages) [C]
diam <- 2   #size [cm]
leach <- 0  #Leaching

for(year in 0:9) {
  res <- yasso(theta, time, avgT, sumP, ampT, init, infall, diam, leach)
  print(c(year,res));
  init <- res
}
res <- yassoSpinn(theta, avgT, sumP, ampT, infall, diam, leach)
print(res)






library(Matrix)
library(expm)
library(expoRkit)

A <- t(matrix(c(-0.639507,2.81383,0.078307,0.0504443,0,0.633112,-6.39507,0.0263111,0.000603138,0,0.000390099,0.00306963,-0.313228,0.00361883,0,0.00049242,0.639507,0.203598,-0.0548308,0,0.00294173,0.0294173,0.00144085,0.000252222,-0.00183784),nrow=5,ncol=5))
z1 <- matrix(c(0.5,0.1,0.1,0.2,0), ncol=1)

expm(A) %*% z1
expAtv(A,z1)
expv(A, z1)
