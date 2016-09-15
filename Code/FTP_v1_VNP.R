library(adehabitatLT)
library(bcpa)
library(intervals)
require(waddle)
library(chron)
library(spatstat)
library(raster)
library(ggplot2)
library(gtools)
library(sp)  # classes for spatial data
library(raster)  # grids, rasters
library(rasterVis)  # raster visualisation
library(maptools)
library(rgeos)
library(geosphere)
library(rgdal)
library(chron)

M1<-read.csv("E:/Moose/Code_Projects/VNP/Output/AllMooseLocs.csv",header=TRUE)

###remove: Collar_ID: 47695, 31188, not from VNP study
M2<-M1[c(M1$Collar_ID!=31188),]
unique(M1$Uniq_ID)

M2$Time <- chron(times=as.character(M2$Time))
M2$Date<- as.POSIXct(M2$Date, format= "%m/%d/%Y ")
M2$dtm <- as.factor(paste(M2$Date, M2$Time, sep = " "))
M2$dtPm<-strptime(M2$dtm, "%Y-%m-%d %H:%M:%S")
# ####convert to local time
M2$dtPm <- as.POSIXct(M2$dtPm, tz="GMT")
M2$dtPmL<-format(M2$dtPm, tz="America/Chicago",usetz=TRUE)
M2$dtPmL<-M2$EndTimeL

###remove locs with large HDOP
M3<-M2[c(M2$HDOP<15),]


####convert to local time
M3$dtPm <- as.POSIXct(M3$dtPm, tz="GMT")
M3$dtPmL<-format(M3$dtPm, tz="America/Chicago",usetz=TRUE)
M3$dtPmL<-as.POSIXct(M3$dtPmL)

M3$dtPmDiff<-c(0,diff(M3$dtPm))


M3$Uniq_IDN<-as.numeric(M3$Uniq_ID)
M3<-M3[order(M3$Uniq_IDN,M3$dtPmL),]
M3$obsID<-1:nrow(M3)

###rename things for the waddle interpolate package
M3$Time<-as.POSIXct(M3$dtPmL)
M3$X<-M3$easting
M3$Y<-M3$northing
M3$id<-M3$IDcollar
M3$IDcollar<-paste(M3$Uniq_ID,M3$Collar_ID,sep="_")


M3<-M3[c(-26,-27,-28,-29,-30,-31,-32,-33)]

###fix issue of 38 and 39 being mixed up
M3$IDflip<-M3$Uniq_ID
M3$IDflip[M3$Uniq_ID==38] <- 39
M3$IDflip[M3$Uniq_ID==39] <- 38
M3$Uniq_ID<-M3$IDflip

M3_8<-M3[M3$Uniq_ID==8,]
M3_8<-M3_8[M3_8$Collar_ID==99210,]
unique(M3_8$Collar_ID)
M3_8$dtPmLDiff<-c(0,diff(M3_8$dtPmL))
M3_8_50<-M3_8[M3_8$dtPmLDiff>62,]


###remove 1st fix from 45 (remove first few days) & 73 (weird first fix time)
M373<-M3[M3$Uniq_ID==73,]
M373<-M373[M373$dtPmL>="2012-03-15 10:30:19",]
M3<-M3[M3$Uniq_ID!=73,]
M3<-rbind(M3,M373)

M345<-M3[M3$Uniq_ID==45,]
head(M345,50)
M345<-M345[M345$Uniq_ID==45&M345$dtPmL>="2011-01-28 08:36:15",]
M3<-M3[M3$Uniq_ID!=45,]
M3<-rbind(M3,M345)


###what's up with the really fast move rates of 41
M341<-M3[M3$Uniq_ID==41,]
par(mfrow=c(1,1))
hist(M341$steplength,freq=FALSE)
hist(M345$steplength,freq=FALSE)
hist(M373$steplength,freq=FALSE)

M341<-M341[order(M341$dtPmL),]

XY<- project(cbind( M341$Longitude,M341$Latitude), "+proj=utm +zone=15 ellps=WGS84")
plot(XY)
M341$Y_test<-XY[,2]
M341$X_test<-XY[,1]
tail(sort(M341$X_test))

###look at outliers
M341<-M341[M341$X_test<=553024,]
plot(M341$X_test,M341$Y_test)

# Distances
nobs<-nrow(M341)  
M341$dist_test<-c(0,sqrt((M341$X_test[-1]-M341$X_test[-nobs])^2+(M341$Y_test[-1]-M341$Y_test[-nobs])^2))

par(mfrow=c(1,1))
plot(M341$dtPmL,M341$steplength)
hist(M341$dist_test)

###looks fixed to me - replace original data with this
M341$steplength<-c(NA,sqrt((M341$X_test[-1]-M341$X_test[-nobs])^2+(M341$Y_test[-1]-M341$Y_test[-nobs])^2))
head(M341)
M341$northing<-M341$Y_test
M341$easting<-M341$X_test

M341<-M341[c(1:32)]
head(M341)
M3<-M3[M3$Uniq_ID!=41,]

names(M341)
M341<-M341[,1:31]
M3<-rbind(M3,M341)



####2010 had 15 min intervals....
M3$Time<-as.POSIXct(M3$dtPm)
M3$X<-M3$easting
M3$Y<-M3$northing
M3$id<-M3$IDcollar
M3$IDcollar<-paste(M3$Uniq_ID,M3$Collar_ID,sep="_")
M3<-M3[order(M3$Uniq_IDN,M3$dtPmL),]

M3_15<-M3[M3$dtPm<"2011-01-24 00:00:00",]
M3_15<-M3_15[M3_15$dtPmDiff>=0,]
M3_15$year<-format(M3_15$dtPmL,"%Y")


M3_20<-M3[M3$dtPm>="2011-01-24 00:00:00",]
M3_20<-M3_20[M3_20$dtPmDiff>=0,]
M3_20$year<-format(M3_20$dtPmL,"%Y")


M3spALLI3<-readRDS("E:/Moose/StateSpaceOutput/Combo_AllMooseCovar_072816.Rda")
str(M3spALLI3)

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


fpt2M<-fpt(VNPm,  seq(10,500, length=100), units = c( "hours"))


varout<-varlogfpt(fpt2M, graph = TRUE)
meanout<-meanfpt(fpt2M, graph = TRUE)
## S3 method for class 'fipati'
attr(fpt2M, "radii")

plot(fpt2M, scale=320, warn = FALSE)

###covert to a data frame
VNPmDF<-ld(VNPm)
str(VNPmDF)
??gBuffer


# points from scratch
VNPmDF$ptid<-1:nrow(VNPmDF)
VNPmDF<-VNPmDF[1:1000,]
coords = cbind(VNPmDF$x, VNPmDF$y)
sp = SpatialPoints(coords)
# make spatial data frame
VNPmDFsp = SpatialPointsDataFrame(coords, VNPmDF)

?gBuffer
str(VNPmDFsp)
pc1100km <- gBuffer(VNPmDFsp, width=320,byid=TRUE )

plot(pc1100km)
plot(VNPm[1],add=TRUE)
plot(VNPmDFsp,add=TRUE)
plot(VNPm[1])
plot(pc1100km,add=TRUE)
# Add data, and write to shapefile
pc100km2 <- SpatialPolygonsDataFrame( pc1100km, data=pc1100km@data )
head(pc100km2)
#writeOGR( pc100km, "pc100km", "pc100km", driver="ESRI Shapefile" )
vegR<-raster("E:/Moose/Data/OutsideSpatial/NLCD11VNPclip1.tif") #raster: NLCD habitat


###can use points and buffer with the "buffer" value in extract function
###look into functions "zonal"
#we will want both the distance from point to 
#and the length inside of the following polylines:
#1)roads and/ or snowmob trails
#lakes
#ponds
#lakes and ponds
#habitat edges

#we will also want temp included if important to wolves?


system.time(testEx2<-extract(vegR,VNPmDFsp[1:1000,],buffer=320))
system.time(testEx<-extract(vegR,pc100km2))

head(testEx2)
data.frame(sapply(testEx2, mean))
testout<-unlist(testEx2)
head(testout)

library(reshape2)

head(melt(testEx2),10)
testEx2df<-melt(testEx2)


str(testEx2df)

library(dplyr)

sub<-testEx2df %>%
  group_by('L1','value') %>%
  mutate(count = n())

sub<-aadt %>%
  select(FID_b30021,FID_AADT_T,CURR_VOL) %>%
  group_by(FID_b30021) %>%
  mutate(Rd.count=n_distinct(FID_AADT_T)) %>%
  mutate(maxV=max(CURR_VOL)) %>%
  mutate(meanV=mean(CURR_VOL))

L1<-testEx2df$L1
value<-testEx2df$value
ddply(testEx2df,L1, transform, count = length(value))
y <- xtabs(~ L1 + value, testEx2df)
x <- count(testEx2df, c('L1', 'value'))
y<-as.data.frame(y)
w <- reshape(y, 
             timevar = "value",
             idvar = c("L1"),
             direction = "wide")
head(w)
str(w)
w<-as.data.frame(w)
w2 <- data.frame(matrix(unlist(w), nrow=1000, byrow=T))
head(w2)



