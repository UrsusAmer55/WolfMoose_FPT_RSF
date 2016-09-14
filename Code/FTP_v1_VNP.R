library(adehabitatLT)


#read in data with interp points
M3spALLI3<-readRDS("E:/Moose/StateSpaceOutput/Combo_AllMooseCovar_072816.Rda")
str(M3spALLI3)
#needed to make traj in ade
#EndTimeL,IDcollar, X, Y
unique(M3spALLI3$IDcollar)
testM<-M3spALLI3[M3spALLI3$IDcollar=="73_99211"|M3spALLI3$IDcollar=="4_99216"|M3spALLI3$IDcollar=="13_101216",]

testMspr<-testM[testM$season=="Spring",]

testM<-testMspr

length(unique(testM$EndTimeL))
head(testM)
names(M3spALLI3)
#remove the duplicates
X<-M3spALLI3[!duplicated(M3spALLI3[,16]),]
X<-testM[!duplicated(testM[,16]),]
VNPm <- as.ltraj(xy = X[,c("X","Y")], date = X$EndTimeL, id = X$IDcollar)
VNPm 
plot(VNPm)
###covert to a data frame
VNPmDF<-ld(VNPm)

#function dl() converts back to traj
is.regular(VNPm)
head(rec(VNPm)[[2]])

plotltr(VNPm, "dt/3600")
head(VNPm)
###if greater than 22 mins
22*60
foo <- function(dt) {
  return(dt> (3600))
}

VNPm2 <- cutltraj(VNPm, "foo(dt)", nextr = TRUE)
VNPm2
plotltr(VNPm2, "dt/3600")
is.sd(VNPm2)


fpt2M<-fpt(VNPm,  seq(300,5000, length=30), units = c( "hours"))

fpt2M
varout<-varlogfpt(fpt2M, graph = TRUE)
meanout<-meanfpt(fpt2M, graph = TRUE)
## S3 method for class 'fipati'
attr(fpt2M, "radii")

plot(fpt2M, scale=1200, warn = FALSE)

###covert to a data frame
VNPmDF<-ld(VNPm)
str(VNPmDF)
??gBuffer

install.packages("rgeos")
library(rgeos)
library(sp)

# points from scratch
VNPmDF$ptid<-1:nrow(VNPmDF)
VNPmDF<-VNPmDF[1:100,]
coords = cbind(VNPmDF$x, VNPmDF$y)
sp = SpatialPoints(coords)
# make spatial data frame
VNPmDFsp = SpatialPointsDataFrame(coords, VNPmDF)

?gBuffer
str(VNPmDFsp)
pc1100km <- gBuffer(VNPmDFsp, width=1100, id=VNPmDFsp$ptid )
# Add data, and write to shapefile
pc100km <- SpatialPolygonsDataFrame( pc100km, data=pc100km@data )
writeOGR( pc100km, "pc100km", "pc100km", driver="ESRI Shapefile" )



