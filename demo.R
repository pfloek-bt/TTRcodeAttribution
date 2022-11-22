# load the necessary functions
source("include.R")
# load an example data set
MyData<-readRDS(file="MyData_SA-BFA-BNP.rds")

# plot the vegetation index from the example data set
plot(MyData$date,MyData$y,type="l")

# simulate the model using parameters previously estimated
sim <- ModelABU(MyData$pars,MyData)

# add a line to show the simulated data
lines(MyData$date,sim,col=2)
