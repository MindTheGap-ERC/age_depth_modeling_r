######################################################################################################
## This R script is supplied in context of the GLOMAR short course "Age Models and Geochronology: An Introduction to Different Age-depth Modelling Approaches"
## It provides numerical recipes for some widely-adopted techniques for basic age-depth modeling 
## by dr. David De Vleeschouwer and dr. Christian Zeeden
######################################################################################################

rm(list=ls()) # Empty global environment
install.packages(c("astrochron", "Hmisc", "Bchron", "devtools"))
library(astrochron) # Load some useful R libraries that must be installed
library(Hmisc)
library(Bchron)
library(devtools)

install_github("robintrayler/modifiedBChron")
library(modifiedBChron)

######################################################################################################
# Let's take a look at the dataset.
######################################################################################################
dates = matrix(c(0,9,34,46,66,72,0,44.8,116,129,523,797.3,1,10,20,40,120,140), nrow = 6, ncol = 3)
colnames(dates) <- c("Depth (m)", "Age (ka)", "Uncertainty (2sigma, kyr)")
dates = as.data.frame(dates)
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim=c(0,1000),xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}


######################################################################################################
# [1] Linear interpolation
######################################################################################################
depths_i=seq(0,100,1)
ages_li=approx(dates$`Depth (m)`, dates$`Age (ka)`, xout = depths_i)
ages_li=as.data.frame(ages_li)
colnames(ages_li) <- c("Depth (m)", "Age (ka)")
# ages_li2 <- linterp(dates[,1:2], dt=1) ### The same result can be obtained through the 'linterp' {astrochron} function
dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim=c(0,1000), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_li$`Age (ka)`, ages_li$`Depth (m)`, col = "grey")

######################################################################################################
# [2] Linear interpolation and extrapolation
######################################################################################################
ages_lie=approxExtrap(dates$`Depth (m)`, dates$`Age (ka)`, xout = depths_i)
ages_lie=as.data.frame(ages_lie)
colnames(ages_lie) <- c("Depth (m)", "Age (ka)")
# ages_lie2 <- tune(cb(depths_i, depths_i), dates[,1:2], extrapolate=T) ### The same result can be obtained through the 'linterp' {astrochron} function
dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, max(ages_lie$`Age (ka)`, na.rm = T)), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_lie$`Age (ka)`, ages_lie$`Depth (m)`, col = "grey")

######################################################################################################
# [3] Linear interpolation and extrapolation with outlier
######################################################################################################
ages_lie_outlier=approxExtrap(dates[c(1:3,5,6),1], dates[c(1:3,5,6),2], xout = depths_i)
ages_lie_outlier=as.data.frame(ages_lie_outlier)
colnames(ages_lie_outlier) <- c("Depth (m)", "Age (ka)")
dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, max(ages_lie_outlier$`Age (ka)`, na.rm = T)), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_lie_outlier$`Age (ka)`, ages_lie_outlier$`Depth (m)`, col = "grey")

######################################################################################################
# [4] Linear interpolation and extrapolation with hiatus
######################################################################################################
depths_i_top = seq(0,61,1)
depths_i_bottom = seq(63,100,1)
ages_lie_hiatus_top=approxExtrap(dates[c(1:4),1], dates[c(1:4),2], xout = depths_i_top)
ages_lie_hiatus_bottom=approxExtrap(dates[c(5:6),1], dates[c(5:6),2], xout = depths_i_bottom)
ages_lie_hiatus_top=as.data.frame(ages_lie_hiatus_top)
ages_lie_hiatus_bottom=as.data.frame(ages_lie_hiatus_bottom)
ages_lie_hiatus=rbind(ages_lie_hiatus_top,ages_lie_hiatus_bottom)
colnames(ages_lie_hiatus) <- c("Depth (m)", "Age (ka)")
rm(ages_lie_hiatus_top,ages_lie_hiatus_bottom)
dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, max(ages_lie_hiatus$`Age (ka)`, na.rm = T)), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_lie_hiatus$`Age (ka)`, ages_lie_hiatus$`Depth (m)`, col = "grey")
abline(h = 62)
text(1500,64, "Hiatus / Hardground")

######################################################################################################
# [5] Linear regression with hiatus
######################################################################################################
lreg_top = lm(`Age (ka)`~`Depth (m)`, data = as.data.frame(dates[c(1:4),c(1:2)]))
summary(lreg_top)
depths_i_top = as.data.frame(depths_i_top)
colnames(depths_i_top) <- c("Depth (m)")
ages_lreg_top = predict(lreg_top, depths_i_top, interval = "confidence")

lreg_bottom= lm(`Age (ka)`~`Depth (m)`, data = as.data.frame(dates[c(5:6),c(1:2)]))
summary(lreg_bottom)
ages_lreg_bottom = lreg_bottom$coefficients[1]+depths_i_bottom*lreg_bottom$coefficients[2]

dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, max(ages_lreg_bottom)), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_lreg_top[,1], depths_i_top[,1], col = "grey")
lines(ages_lreg_top[,2], depths_i_top[,1])
lines(ages_lreg_top[,3], depths_i_top[,1])
points(ages_lreg_bottom, depths_i_bottom, col = "grey")
abline(h = 62)
text(1500,64, "Hiatus / Hardground")
text(1000,30, "You can calculate confidence levels on a regression. But should you?")
text(1000,33, "Confidence bands do not consider the uncertainties on the input dates")
text(1000,36, "These uncertainties have little geologic meaning!", col = "red")

######################################################################################################
# [6] Polynomial regression
######################################################################################################
order = 4
preg= lm(dates$`Age (ka)` ~ poly(dates$`Depth (m)`,order, raw = T))
summary(preg)

ages_preg = matrix(0,nrow = 1, ncol = length(depths_i))
for (o in 1:order+1){
  ages_preg = ages_preg + preg$coefficients[o]*depths_i^(o-1)
  }

dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, max(ages_lreg_bottom)), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_preg, depths_i, col = "grey")
text(1200,30,"Careful!!! Depending on the dataset", col = "red")
text(1200,33,"and the order of the polynomial regression", col = "red")
text(1200,36,"one can get age inversions and/or", col = "red")
text(1200,39,"extremely high or low sedimentation rates", col = "red")

######################################################################################################
# [7] Linear interpolation and extrapolation with analytic uncertainty calculation! 
# CAREFUL CAREFUL CAREFUL: This "naive" error propagation is an over-optimistic approach to quantifying uncertainty!!! 
# This is because it assumes constant sedimentation rates in-between dated levels: This is of course a false assumption for most geologic archives!
######################################################################################################
dyi=c() # dyi is the uncertainty on the modeled age at depth level i
for (i in 1:72){
  temp1=sort(c(depths_i[i], dates$`Depth (m)`))
  temp2=which(temp1==depths_i[i])
  x1=temp1[max(temp2)-1] # x1 is the dated level above depth level i
  x2=temp1[max(temp2)+1] #x2 is the dated level above depth level i
  dy1=dates$`Uncertainty (2sigma, kyr)`[max(temp2)-1] # dy1 is the uncertainty on the age y1 at dated level x1
  dy2=dates$`Uncertainty (2sigma, kyr)`[max(temp2)] # dy2 is the uncertainty on the age y2 at dated level x2
  dyi[i]=sqrt(((1-(depths_i[i]-x1)/(x2-x1))*dy1)^2+(((depths_i[i]-x1)/(x2-x1))*dy2)^2) 
}
for (i in 73:101){
  temp1=sort(c(depths_i[i], dates$`Depth (m)`))
  temp2=which(temp1==depths_i[i])
  x1=temp1[max(temp2)-2]
  x2=temp1[max(temp2)-1]
  dy1=dates$`Uncertainty (2sigma, kyr)`[max(temp2)-2]
  dy2=dates$`Uncertainty (2sigma, kyr)`[max(temp2)-1]
  dyi[i]=sqrt(((1+(depths_i[i]-x2)/(x2-x1))*dy2)^2+(((-depths_i[i]+x2)/(x2-x1))*dy1)^2)
}

dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, 2500), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_lie$`Age (ka)`, ages_lie$`Depth (m)`, col = "grey")
lines(ages_lie$`Age (ka)`+dyi, ages_lie$`Depth (m)`) # Upper confidence level
lines(ages_lie$`Age (ka)`-dyi, ages_lie$`Depth (m)`) # Lower confidence level
text(1200,30,"Careful!!! This approach is over-optimistic!", col = "red")
text(1200,33,"Uncertainty decreases when moving stratigraphically away from dated levels", col = "red")
text(1200,36,"This approach wrongly assumes constant sedimentation rates between dated levels!", col = "red")

######################################################################################################
# [8] Linear interpolation and extrapolation with Monte-Carlo based uncertainty calculation! 
# CAREFUL CAREFUL CAREFUL: The Monte-Carlo approach is equally "naive" as the analytic uncertainty calculation 
# It also assumes constant sedimentation rates in-between dated levels: This is of course a false assumption for most geologic archives!
######################################################################################################
Nsimulations = 1000
ages_lie_MonteCarlo=matrix(data = NA, nrow = 101, ncol = Nsimulations )
for (i in 1:Nsimulations){
  age_depth_MonteCarlo=c()
  for (j in 1:6){
  age_depth_MonteCarlo[j]=rnorm(1, mean = dates$`Age (ka)`[j], sd = dates$`Uncertainty (2sigma, kyr)`[j]/2)}
    temp=approxExtrap(dates$`Depth (m)`, age_depth_MonteCarlo, xout = depths_i)$y
  ages_lie_MonteCarlo[,i]=temp
}

ages_lie_MC=matrix(data = NA, nrow = 101, ncol = 3)
colnames(ages_lie_MC) <- c("2.5% CL", "Median Age", "97.5% CL")
for (k in 1:101){
ages_lie_MC[k,]=quantile(ages_lie_MonteCarlo[k,], probs = c(0.025,0.5,0.975))}

dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, 2500), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_lie_MC[,2], depths_i, col = "grey")
lines(ages_lie_MC[,1], depths_i) # Upper confidence level
lines(ages_lie_MC[,3], depths_i) # Lower confidence level
text(1200,30,"Careful!!! This approach is over-optimistic!", col = "red")
text(1200,33,"Uncertainty decreases when moving stratigraphically away from dated levels", col = "red")
text(1200,36,"This approach wrongly assumes constant sedimentation rates between dated levels!", col = "red")

######################################################################################################
# [9] Linear interpolation and extrapolation with Monte-Carlo based uncertainty calculation. 
# But before interpolation, we verify whether the sampled ages (age_depth_MonteCarlo) are monotonically increasing. Samples with age inversions are ignored.
# CAREFUL: This Monte-Carlo approach is less over-optimistic than [8], but it still assumes constant sedimentation rates in-between dated levels!
######################################################################################################
ages_lie_MonteCarlo=matrix(data = NA, nrow = 101, ncol = Nsimulations )
for (i in 1:Nsimulations){
  age_depth_MonteCarlo=c()
  for (j in 1:6){
    age_depth_MonteCarlo[j]=rnorm(1, mean = dates$`Age (ka)`[j], sd = dates$`Uncertainty (2sigma, kyr)`[j]/2)}
  if (length(which(diff(age_depth_MonteCarlo)<0))>0){ages_lie_MonteCarlo[,i] = rep(NA, 101)}
  else{temp=approxExtrap(dates$`Depth (m)`, age_depth_MonteCarlo, xout = depths_i)$y
  ages_lie_MonteCarlo[,i]=temp}
}

ages_lie_MC=matrix(data = NA, nrow = 101, ncol = 3)
colnames(ages_lie_MC) <- c("2.5% CL", "Median Age", "97.5% CL")
for (k in 1:101){
  ages_lie_MC[k,]=quantile(ages_lie_MonteCarlo[k,], probs = c(0.025,0.5,0.975), na.rm = T)}

dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, 2500), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_lie_MC[,2], depths_i, col = "grey")
lines(ages_lie_MC[,1], depths_i) # Upper confidence level
lines(ages_lie_MC[,3], depths_i) # Lower confidence level
text(1200,30,"Careful!!! This approach is still somewhat over-optimistic!", col = "red")
text(1200,33,"Uncertainty still decrease when moving stratigraphically away from dated levels", col = "red")
text(1200,36,"This approach still assumes constant sedimentation rates between dated levels!", col = "red")

######################################################################################################
# [10] Linear interpolation and extrapolation with interpolation in-between confidence levels
# This is a simple solution to alleviate the "over-optimistic" character demonstrated in the two previous blocs [7] & [8]
######################################################################################################
dyi=c() # dyi is the uncertainty on the modeled age at depth level i
dyi=approx(dates$`Depth (m)`, dates$`Uncertainty (2sigma, kyr)`, xout = depths_i)$y # Simple linear interpolation of the uncertainty in-between dated levels
for (i in 73:101){ # In the extrapolation part, we keep using the analytic approach because it does not suffer from being over-optimistic.
  temp1=sort(c(depths_i[i], dates$`Depth (m)`))
  temp2=which(temp1==depths_i[i])
  x1=temp1[max(temp2)-2]
  x2=temp1[max(temp2)-1]
  dy1=dates$`Uncertainty (2sigma, kyr)`[max(temp2)-2]
  dy2=dates$`Uncertainty (2sigma, kyr)`[max(temp2)-1]
  dyi[i]=sqrt(((1+(depths_i[i]-x2)/(x2-x1))*dy2)^2+(((-depths_i[i]+x2)/(x2-x1))*dy1)^2)
}

dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, 2500), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_lie$`Age (ka)`, ages_lie$`Depth (m)`, col = "grey")
lines(ages_lie$`Age (ka)`+dyi, ages_lie$`Depth (m)`) # Upper confidence level
lines(ages_lie$`Age (ka)`-dyi, ages_lie$`Depth (m)`) # Lower confidence level


######################################################################################################
# [11] Using a Poisson distribution to simulate break-points in-between dated levels where sedimentation rate changes.
######################################################################################################
Nsimulations = 1000
ages_lie_MonteCarlo_Poisson=matrix(data = NA, nrow = 101, ncol = Nsimulations )

for (i in 1:Nsimulations){
  age_depth_MonteCarlo=c()
  depth_Poisson=c()
  age_Poisson=c()
  
  for (j in 1:6){
  age_depth_MonteCarlo[j]=rnorm(1, mean = dates$`Age (ka)`[j], sd = dates$`Uncertainty (2sigma, kyr)`[j]/2)
  }
  if (length(which(diff(age_depth_MonteCarlo)<0))>0){ages_lie_MonteCarlo_Poisson[,i] = rep(NA, 101)}
  
  else{
  for (k in 1:5){
    thickness = dates$`Depth (m)`[k+1]-dates$`Depth (m)`[k]
    lambda = thickness / 20
    N_breaks = rpois(1, lambda = lambda)
    temp1=runif(N_breaks, min = dates$`Depth (m)`[k], max = dates$`Depth (m)`[k+1]) # temp1 represents a depth-level at which a break-point is introduced. temp1 is sampled from a uniform distribution between the two adjacent dated levels
    temp2=runif(N_breaks, min = age_depth_MonteCarlo[k], max = age_depth_MonteCarlo[k+1]) # temp2 represents the age of the introduced break-point. temp2 is sampled from a uniform distribution between the two adjacent dates
    depth_Poisson = c(depth_Poisson, temp1) 
    age_Poisson = c(age_Poisson, temp2)
  }
  depths_sim = c(dates$`Depth (m)`, depth_Poisson) # We now concatenate the depths of the dated levels and the introduced break-points for this specific Monte-Carlo simulation
  ages_sim = c(age_depth_MonteCarlo, age_Poisson) # We now concatenate the ages of the dated levels and the introduced break-points for this specific Monte-Carlo simulation
  temp=approxExtrap(depths_sim, ages_sim , xout = depths_i)$y # Temp is the age-depth model (through linear interpolation) of this specific Monte-Carlo simulation
  ages_lie_MonteCarlo_Poisson[,i]=temp # ages_lie_MonteCarlo_Poisson is a large matrix, storing each of the 1000 age-depth models that are generated during Monte-Carlo. 
  }
}

ages_lie_MC_Poisson=matrix(data = NA, nrow = 101, ncol = 3)
colnames(ages_lie_MC_Poisson) <- c("2.5% CL", "Median Age", "97.5% CL")
for (k in 1:101){
  ages_lie_MC_Poisson[k,]=quantile(ages_lie_MonteCarlo_Poisson[k,], probs = c(0.025,0.5,0.975), na.rm = T)}

dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, 2500), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_lie_MC_Poisson[,2], depths_i, col = "grey")
lines(ages_lie_MC_Poisson[,1], depths_i) # Upper confidence level
lines(ages_lie_MC_Poisson[,3], depths_i) # Lower confidence level
text(1200,30,"Uncertainties strongly depend on the adopted Poisson model", col = "red")
text(1200,33,"Here, lambda is set at 20.", col = "red")
text(1200,36,"Here, this means that sedimentation rates are expected to change every 20 meters.", col = "red")

######################################################################################################
# [12] Bchron
######################################################################################################

output = Bchronology(
  ages = dates$`Age (ka)`,
  ageSds = dates$`Uncertainty (2sigma, kyr)`/2,
  positions = dates$`Depth (m)`,
  calCurves = rep("normal", 6),
  predictPositions = depths_i)

ages_Bchron=matrix(data = NA, nrow = 101, ncol = 3)
colnames(ages_Bchron) <- c("2.5% CL", "Median Age", "97.5% CL")
for (k in 1:101){
  ages_Bchron[k,]=quantile(output$thetaPredict[,k], probs = c(0.025,0.5,0.975))}

dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, 2500), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_Bchron[,2], depths_i, col = "grey")
lines(ages_Bchron[,1], depths_i) # Upper confidence level
lines(ages_Bchron[,3], depths_i) # Lower confidence level

######################################################################################################
# [13] Modified Bchron
######################################################################################################


output2 = ageModel(
  ages = dates$`Age (ka)`,
  ageSds = dates$`Uncertainty (2sigma, kyr)`/2,
  positions = -dates$`Depth (m)`,
  distTypes = rep("G", 6),
  predictPositions = -depths_i,
  MC = 5000,
  burn = 1000,
  ids = c("1","2","3","4","5","6"))

modelPlot(output2, scale = 0.1, xlim = c(0,2500), ylim = c(-100,0), xaxs = "i", yaxs = "i")
points(ages_Bchron[,2], -depths_i, col = "gray40")
lines(ages_Bchron[,1], -depths_i, col = "gray40") # Upper confidence level
lines(ages_Bchron[,3], -depths_i, col = "gray40") # Lower confidence level

ages_Bchron2=matrix(data = NA, nrow = 101, ncol = 3)
colnames(ages_Bchron2) <- c("2.5% CL", "Median Age", "97.5% CL")
for (k in 1:101){
  ages_Bchron2[k,]=quantile(output2$model[k,], probs = c(0.025,0.5,0.975))}

dev.off()
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, 2500), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_Bchron2[,2], depths_i, col = "grey")
lines(ages_Bchron2[,1], depths_i) # Upper confidence level
lines(ages_Bchron2[,3], depths_i) # Lower confidence level
points(ages_Bchron[,2], depths_i, col = "sienna1")
lines(ages_Bchron[,1], depths_i, col = "sienna3") # Upper confidence level
lines(ages_Bchron[,3], depths_i, col = "sienna3") # Lower confidence level







######################################################################################################
# [13] Error-free time-axis: Transfer uncertainty in the proxy domain. 
# This code is using the the Bchron age-depth model and uncertainties that have been determined in [11]. Please execute block [11] before running block [13]
######################################################################################################
proxy_series=ar1(401, dt = 0.25, rho = 0.99)
dev.off()
plot(proxy_series[,2], proxy_series[,1], ylim = c(100,0), xlim = c(-3,3), xaxs = "i", yaxs = "i", type = "l", col = "black", pch = 19, xlab = "Proxy (unit)", ylab = "Depth (m)")

proxy_timeseries=tune(proxy_series, cbind(depths_i,ages_Bchron[,2]), extrapolate = T)
dev.off()
plot(proxy_timeseries[,2], proxy_timeseries[,1], ylim = c(1000,0), xlim = c(-3,3), xaxs = "i", yaxs = "i", type = "l", col = "black", pch = 19, xlab = "Proxy (unit)", ylab = "Time (ka)")
for (j in 1:6){
  points(-2.8, dates[j,2], pch = 19, col = "blue")
  segments(-2.8,dates[j,2]-dates[j,3],-2.8,dates[j,2]+dates[j,3], col = "blue", lwd = 2)}


time_i = seq(0, max(proxy_timeseries[,1]), 5)
proxy_lowerCL=c()
proxy_upperCL=c()
for (i in 1:length(time_i)){
  idx = which(ages_Bchron[,1] < time_i[i] & ages_Bchron[,3] > time_i[i])
  depth_range = depths_i[idx]
  idx2 = which(proxy_series[,1] >= min(depth_range) & proxy_series[,1]<= max(depth_range))
  proxy_range = proxy_series[idx2,2]
  proxy_lowerCL[i] = min(proxy_range)
  proxy_upperCL[i] = max(proxy_range)
}

dev.off()
plot(proxy_timeseries[,2], proxy_timeseries[,1], ylim = c(1000,0), xlim = c(-3,3), xaxs = "i", yaxs = "i", type = "l", col = "black", pch = 19, xlab = "Proxy (unit)", ylab = "Time (ka)")
polygon(c(proxy_lowerCL, rev(proxy_upperCL)), c(time_i, rev(time_i)),col = "grey70", border = NA)
lines(proxy_timeseries[,2], proxy_timeseries[,1], lwd = 2)

######################################################################################################
# [14] Error-free time-axis: Transfer uncertainty in the proxy domain. 
# This code is using the age-depth model (ages_lie) and uncertainty (dyi) that have been determined in [10]. Please execute block [10] before running block [13]!
######################################################################################################
proxy_series=ar1(401, dt = 0.25, rho = 0.99)
proxy_timeseries=tune(proxy_series, ages_lie, extrapolate = T)
dev.off()
plot(proxy_timeseries[,2], proxy_timeseries[,1], ylim = c(2000,0), xlim = c(-3,3), xaxs = "i", yaxs = "i", type = "l", col = "black", pch = 19, xlab = "Proxy (unit)", ylab = "Time (ka)")
for (j in 1:6){
  points(-2.8, dates[j,2], pch = 19, col = "blue")
  segments(-2.8,dates[j,2]-dates[j,3],-2.8,dates[j,2]+dates[j,3], col = "blue", lwd = 2)}


time_i = seq(0, max(proxy_timeseries[,1]), 2)
proxy_lowerCL=c()
proxy_upperCL=c()
for (i in 1:length(time_i)){
  idx = which(ages_lie$`Age (ka)`-dyi < time_i[i] & ages_lie$`Age (ka)`+dyi > time_i[i])
  depth_range = depths_i[idx]
  idx2 = which(proxy_series[,1] >= min(depth_range) & proxy_series[,1]<= max(depth_range))
  proxy_range = proxy_series[idx2,2]
  proxy_lowerCL[i] = min(proxy_range)
  proxy_upperCL[i] = max(proxy_range)
}

dev.off()
plot(proxy_timeseries[,2], proxy_timeseries[,1], ylim = c(2000,0), xlim = c(-3,3), xaxs = "i", yaxs = "i", type = "l", col = "black", pch = 19, xlab = "Proxy (unit)", ylab = "Time (ka)")
polygon(c(proxy_lowerCL, rev(proxy_upperCL)), c(time_i, rev(time_i)),col = "grey70", border = NA)
lines(proxy_timeseries[,2], proxy_timeseries[,1], lwd = 2)
