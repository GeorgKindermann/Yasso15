source("https://raw.githubusercontent.com/GeorgKindermann/Yasso15/master/y15_subroutine.r")
#source("y15_subroutine.r")  #In case you have it local on your computer


theta <- c(0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-7.90E-05,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26)
init <- rep(0,5)
infall <- c(0.5,0.1,0.1,0.2,0) #0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
time <- 1   #Time to run
avgT <- 10  #Temp annual average [C]
sumP <- 600 #Precip annual summ [mm]
ampT <- 12  #Amplitude (max. difference of month averages / 2) [C]
diam <- 2   #size [cm]
leach <- 0  #Leaching

for(year in 0:9) {
  res <- yasso.getNextTimestep(theta, avgT, sumP, ampT, diam, leach, init, infall, time)
  print(c(year,res));
  init <- res
}
res <- yasso.getSpin(theta, avgT, sumP, ampT, diam, leach, infall)
print(res)
