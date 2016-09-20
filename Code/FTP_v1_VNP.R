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

#need to run separately on M115 and M20
X<-M3_20

X<-X[order(X$Uniq_ID,X$dtPmL),]
head(X)

#create datetime/idcollar combo
X$collardate<-paste(X$IDcollar,X$dtPmL,sep="_")

##get rid of duplicates
names(X)

dups<-X[duplicated(X[,34]),]
head(dups,40)
table((dups$IDcollar))
X2<-X[!duplicated(X[,34]),]
table((X2$IDcollar))

VNPm <- as.ltraj(xy = X2[,c("X","Y")], date = X2$dtPmL, id = X2$IDcollar)

VNPm
is.regular(VNPm)
plotltr(VNPm[5], "dt/3600")


##function to identify if time lag greater than 1 day
foo <- function(dt) {
   return(dt> (.5*(3600*24)))
   }

VNPm2 <- cutltraj(VNPm, "foo(dt)", nextr = TRUE)
VNPm2

VNPm2df<-ld(VNPm2)
head(VNPm2df)

###rename columns for Interpolate funcs
VNPm2df$X<-VNPm2df$x
names(VNPm2df)[names(VNPm2df)=="x"] <- "X"
names(VNPm2df)[names(VNPm2df)=="y"] <- "Y"
names(VNPm2df)[names(VNPm2df)=="date"] <- "Time"

UID<-unique(VNPm2df$burst)
samplehold<-NULL
for(i in 1:length(unique(UID))){
  sub<-VNPm2df[VNPm2df$burst==UID[i],]
  Moose.Data <- InterpolatePoints(sub, n = 20,units="min")$Data
  Moose.VT <- GetVT(Moose.Data, units = "min")
  Moose.VT$burstID<-UID[i]
  Moose.VT$ID<-rep(unique(sub$burst),nrow(Moose.VT))
  Moose.VT$Collar<-rep(unique(sub$id),nrow(Moose.VT))
  #Moose.VT<-cbind(Moose.VT,sub$X[c(3:nrow(Moose.VT))],sub$Y[c(3:nrow(Moose.VT))],sub$dtPmL[c(3:nrow(Moose.VT))])
  samplehold<-rbind(samplehold,Moose.VT)
}



###446751 obvs going in - not w/ interp 548959
1-(446751/548959)

VNPmN<-samplehold
str(VNPmN)
###convert T.Posix to Local
VNPmN$T.POSIX <- as.POSIXct(VNPmN$T.POSIX, tz="GMT")
# Moose.VT$T.POSIXL<-format(Moose.VT$T.POSIX, tz="America/Chicago",usetz=TRUE)
# Moose.VT$T.POSIXL<-as.POSIXct(Moose.VT$T.POSIXL)
###for 20 min group
VNPmN$EndTimeL<-VNPmN$T.POSIX+(10*60)
head(VNPmN)

Re(VNPmN$Z.end[1])
VNPmN$Z.endCH<-as.character(VNPmN$Z.end)

step2<-strsplit(VNPmN$Z.endCH, "i", fixed = TRUE)
step2<-as.character(step2)
step3 <-as.numeric(unlist(strsplit(step2, "+", fixed = TRUE)))

step3even <- step3[seq(1, length(step3), 2)]
step3odd <- step3[seq(2, length(step3), 2)]


step4even<-as.data.frame((step3even))
step4odd<-as.data.frame((step3odd))
colnames(step4even)<-"X.ENDloc"
colnames(step4odd)<-"Y.ENDloc"

VNPmN<-cbind(VNPmN,step4even,step4odd)

VNPmN$X<-VNPmN$X.ENDloc #change names of coords
VNPmN$Y<-VNPmN$Y.ENDloc #change names of coords

saveRDS(VNPmN,"E:/Moose/FPT/DataProc/VNPmN_091916.R")
VNPmN<-readRDS("E:/Moose/FPT/DataProc/VNPmN_091916.R")




VNPmN$month<-format(VNPmN$EndTimeL,"%m")
VNPmN$month<-as.numeric(VNPmN$month)
VNPmN$season<-NA
VNPmN$season[VNPmN$month>=4 & VNPmN$month<7]<- "Spring"
VNPmN$season[VNPmN$month>=7 & VNPmN$month<11]<- "Summer"
VNPmN$season[VNPmN$month==11 |VNPmN$month==12 | VNPmN$month==01 | VNPmN$month==02 | VNPmN$month==03 ]<- "Winter"

VNPmNspr<-VNPmN[VNPmN$season=="Spring",]
VNPmNsum<-VNPmN[VNPmN$season=="Summer",]
VNPmNwin<-VNPmN[VNPmN$season=="Winter",]

VNPmNspr$burstID<-droplevels(VNPmNspr$burstID)
VNPmNsum$burstID<-droplevels(VNPmNsum$burstID)
VNPmNwin$burstID<-droplevels(VNPmNwin$burstID)



head(VNPmN)
VNPmNsprtrajINT <- as.ltraj(xy = VNPmNspr[,c("X","Y")], date = VNPmNspr$EndTimeL, id = VNPmNspr$burstID)
VNPmNsumtrajINT <- as.ltraj(xy = VNPmNsum[,c("X","Y")], date = VNPmNsum$EndTimeL, id = VNPmNsum$burstID)
VNPmNwintrajINT <- as.ltraj(xy = VNPmNwin[,c("X","Y")], date = VNPmNwin$EndTimeL, id = VNPmNwin$burstID)

##function to identify if time lag greater than 1 day
foo <- function(dt) {
  return(dt> (.5*(3600*24)))
}

VNPmNsprtrajINT2 <- cutltraj(VNPmNsprtrajINT, "foo(dt)", nextr = TRUE)
VNPmNsprtrajINT2

VNPmNsumtrajINT2 <- cutltraj(VNPmNsumtrajINT, "foo(dt)", nextr = TRUE)
VNPmNsumtrajINT2

VNPmNwintrajINT2 <- cutltraj(VNPmNwintrajINT, "foo(dt)", nextr = TRUE)
VNPmNwintrajINT2

VNPmNsprtrajINTdf<-ld(VNPmNsprtrajINT)
str(VNPmNsprtrajINTdf)

VNPmNsprtrajINT2 
is.regular(VNPmNwintrajINT2)
plot(VNPmNsumtrajINT[48])
str(VNPmNltrajINT)
VNPmNwintrajINT2df<-ld(VNPmNwintrajINT2)
str(VNPmNwintrajINT2df)


totlocs<-with(VNPmNwintrajINT2df, aggregate(dist, list(burst), FUN = function(x) length(x)))
winbursts<-unique(VNPmNwintrajINT2df$burst)




WINfpt<-fpt(VNPmNwintrajINT2,  seq(30,1000, length=100), units = c( "hours"))

varoutWIN<-varlogfpt(WINfpt, graph = FALSE)
head(varoutWIN)
plot(varoutWIN$r30)
names(varoutWIN)
varoutWINdf<-data.frame(varoutWIN)
head(varoutWINdf)
str(varoutWINdf)

varoutWINdf2<-cbind(winbursts,totlocs,varoutWINdf)
head(varoutWINdf2)

#remove bursts with less than 360 aka 5 days
((5*24)*60)/20
varoutWINdf3<-varoutWINdf2[varoutWINdf2$x>=500,]
names(varoutWINdf3)

str(varoutWINdf3)
Winmeans<-aggregate(varoutWINdf3[,4:103],by=list(varoutWINdf3$winbursts), mean)
nobs<-varoutWINdf3[,3]
str(nobs)
head(Winmeans)
str(Winmeans)
unique(Winmeans$Group.1)
WinmeansN<-cbind(Winmeans,nobs)
head(WinmeansN)


###get the weighted mean by individual ID-collar

WinmeansN$IDcollar<-substr(WinmeansN$Group.1,1,9)
WinmeansN$ID<-substr(WinmeansN$Group.1,1,2)
WinmeansN$ID<-as.factor(WinmeansN$ID)
unique(WinmeansN$ID)
WinmeansN$ID<-as.factor(WinmeansN$IDcollar)
unique(WinmeansN$IDcollar)

library(data.table)

dt <- as.data.table(WinmeansN)
head(dt)
names(dt)
colsToKeep = c(names(WinmeansN[,2:102]))


dt2 <- dt[,lapply(.SD,weighted.mean,w=nobs), 
          by = list(IDcollar), .SDcols = colsToKeep]

dfWinID<-data.frame(dt2)
head(dfWinID)
unique(dfWinID$IDcollar)
names(dfWinID)

WinmeansbyCollID<-colMeans()
  ?colMeans



library(matrixStats)
colSds(dfWinID[,2:101])

head(dfWinID)

dfWinID$one<-"1"
?aggregate
aggregate(dfWinID[,2:101], FUN=mean,by=list(dfWinID$one),na.action = na.omit)
head(dfWinID)


with(dfWinID[,2:101], aggregate(dfWinID[,2:101], FUN =  function(x) c( SD = sd(x), MN= mean(x) ) ) )

head(dfWinID)
plot(colMeans(dfWinID[,2:101],na.rm = TRUE))

mean<-aggregate(dfWinID[,2:101], FUN=mean,by=list(dfWinID$one),na.rm=TRUE)
stdev<-aggregate(dfWinID[,2:101], FUN=sd,by=list(dfWinID$one),na.rm=TRUE)
num<-with(dfWinID, aggregate(dfWinID[,2:101], list(one), FUN = function(x) length(x)))

names(mean)
str(mean)

meanSDn<-rbind(mean[,2:101],stdev[,2:101],num[,2:101])
meanSDnt<-t(meanSDn)
head(meanSDnt)
colnames(meanSDnt)<-c("mean","sd","n")
meanSDnt<-data.frame(meanSDnt)
str(meanSDnt)
radii<-as.data.frame(attr(varoutWIN,"radii"))
colnames(radii)<-c("radii")
str(radii)
meanSDntR<-cbind(meanSDnt,radii)
head(meanSDntR)
meanSDntR$UPCI<-meanSDntR$mean+(1.96*(meanSDntR$sd/sqrt(meanSDntR$n)))
meanSDntR$LOCI<-meanSDntR$mean-(1.96*(meanSDntR$sd/sqrt(meanSDntR$n)))  
  

plot(meanSDntR$radii,meanSDntR$mean,ylim=c(0,1.5))
lines(meanSDntR$radii,meanSDntR$mean)
lines(meanSDntR$radii,meanSDntR$UPCI,col="red")
lines(meanSDntR$radii,meanSDntR$LOCI,col="red")

str(varoutWIN)


attributes(varoutWINdf) <- NULL
str(varoutWINdf)
varoutWINdf<-data.frame(unlist(str(varoutWINdf)))

str(VNPmNspr)

# radiit<-t(radii)
# str(radiit)

str(varoutWIN)

varoutWINt<-t(varoutWIN)
head(varoutWINt)

varoutWINr<-rbind(varoutWIN,radii)

###need to find the max and return the corresponding radii



str(radii)
structure(fpt2M320)


varoutWIN[, "max"] <- apply(df[, 2:26], 1, max)

str(varout)
plot(varout$r1)
meanout<-meanfpt(WINfpt, graph = TRUE)




#This was used to read in interp data used in HMM analysis
#M3spALLI3<-readRDS("E:/Moose/StateSpaceOutput/Combo_AllMooseCovar_072816.Rda")
str(M3spALLI3)

str(M3spALLI3)
#needed to make traj in ade
#EndTimeL,IDcollar, X, Y
unique(M3spALLI3$IDcollar)
testM<-M3spALLI3[M3spALLI3$IDcollar=="73_99211"|M3spALLI3$IDcollar=="4_99216"|M3spALLI3$IDcollar=="13_101216",]
testM<-M3spALLI3[M3spALLI3$IDcollar=="73_99211",]


testMspr<-testM[testM$season=="Spring",]

testM<-testMspr

length(unique(testM$EndTimeL))
head(testM)
names(M3spALLI3)
#remove the duplicates
X<-M3spALLI3[!duplicated(M3spALLI3[,16]),]
X<-testM[!duplicated(testM[,16]),]
head(X)
X<-X[order(X$ID,X$T.POSIXL),]

VNPm <- as.ltraj(xy = X[,c("X","Y")], date = X$EndTimeL, id = X$IDcollar)



plot(VNPm)
###covert to a data frame
VNPmDF<-ld(VNPm)

#function dl() converts back to traj
is.regular(VNPm)
`head(rec(VNPm)[[2]])

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

plot(fpt2M, scale=40, warn = FALSE)

###covert to a data frame
VNPmDF<-ld(VNPm)
str(VNPmDF)
??gBuffer


# points from scratch
VNPmDF$ptid<-1:nrow(VNPmDF)
VNPmDF<-VNPmDF[1:400,]
coords = cbind(VNPmDF$x, VNPmDF$y)
sp = SpatialPoints(coords)
# make spatial data frame
VNPmDFsp = SpatialPointsDataFrame(coords, VNPmDF)

?gBuffer
str(VNPmDFsp)
pc1100km <- gBuffer(VNPmDFsp, width=40,byid=TRUE )

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
nlcdtab<-read.csv("E:/Moose/Data/OutsideSpatial/NLCD11_Table.csv",header=TRUE) #table so you can see what raster value = what hab class
nlcdtab[1:2]

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


system.time(testEx2<-extract(vegR,VNPmDFsp[1:400,],buffer=40))
system.time(testEx<-extract(vegR,pc100km2))

head(testEx2)
# data.frame(sapply(testEx2, mean))
# testout<-unlist(testEx2)
# head(testout)

library(reshape2)
testEx2df<-melt(testEx2)


# str(testEx2df)
# 
# library(dplyr)
# testEx2df$value<-as.factor(testEx2df$value)
# testEx2df$L1<-as.factor(testEx2df$L1)
# 
# sub<-testEx2df %>%
#   group_by(L1,value) %>%
#   mutate(count = n())
# 
# head(sub)
# 
# library(plyr)
# 
# str(testEx2df)
# count(testEx2df, vars=c("L1","value"))
# ?aggregate
# 
# wolfave<-aggregate(testEx2df$value,list(testEx2df$L1),length)
# head(wolfave)
# 
# head(w3)
# 
# sub<-aadt %>%
#   select(FID_b30021,FID_AADT_T,CURR_VOL) %>%
#   group_by(FID_b30021) %>%
#   mutate(Rd.count=n_distinct(FID_AADT_T)) %>%
#   mutate(maxV=max(CURR_VOL)) %>%
#   mutate(meanV=mean(CURR_VOL))
# 
# L1<-testEx2df$L1
# value<-testEx2df$value
# ddply(testEx2df,L1, transform, count = length(value))
y <- xtabs(~ L1 + value, testEx2df)
# x <- count(testEx2df, c('L1', 'value'))
y<-as.data.frame(y)
w <- reshape(y, 
             timevar = "value",
             idvar = c("L1"),
             direction = "wide")
head(w)
nlcdtab[1:2]
#colnames(w)<-c("ptid","OpenWater","DeciduousForest","EvergreenForest","MixedForest","ShrubScrub","WoodyWet","HerbWet")
colnames(w)<-c("ptid","DeciduousForest","EvergreenForest","MixedForest","WoodyWet","HerbWet")
head(w)
w$sum<-rowSums(w[2:ncol(w)])
summary(w$sum)

w$PerWater<-w$OpenWater/w$sum
w$PerWDecFor<-w$DeciduousForest/w$sum
w$PerEverFor<-w$EvergreenForest/w$sum
w$PerMixFor<-w$MixedForest/w$sum
w$PerShrub<-w$ShrubScrub/w$sum
w$PerWoodWet<-w$WoodyWet/w$sum
w$PerHerbWet<-w$HerbWet/w$sum

hist(w$PerWDecFor)

##subset to just the percentage and ID
names(w)
w2 <- w[c(1,8:12)]
names(w2)


str(VNPmDFsp@data)

mergetest<-merge(VNPmDFsp@data,w2,by.x="ptid",by.y="ptid")
head(mergetest)

str(X)
VNPm400<-X[1:700,]
VNPm400 <- as.ltraj(xy = VNPm400[,c("X","Y")], date = VNPm400$EndTimeL, id = VNPm400$IDcollar)



fpt2M320<-fpt(VNPm400,  c(40), units = c( "hours"))
plot(fpt2M320, scale=40, warn = TRUE)

?fpt


str(fpt2M320)
attr(fpt2M320$data.frame$r1,"date")
structure(fpt2M320)

dftest<-as.data.frame(fpt2M320[[1]])
head(dftest)
dftest2 <- dftest[!is.na(dftest$r1)]

dates<-as.data.frame(attr(dftest,"date"))
head(dates)
colnames(dates)<-c("fipatiDATE")


fptval<-dftest$r1
head(fptval)


fptvaldate<-cbind(dates,fptval)
head(fptvaldate)

head(VNPm400)
VNPm400DF<-ld(VNPm400)
head(VNPm400DF)


str(mergetest)
mergetest2<-merge(mergetest,fptvaldate,by.x="date",by.y="fipatiDATE")
head(mergetest2)

plot(mergetest2$PerWDecFor,mergetest2$fptval)
plot(mergetest2$PerEverFor,mergetest2$fptval)
plot(mergetest2$PerMixFor,mergetest2$fptval)
plot(mergetest2$PerWoodWet,mergetest2$fptval)
plot(mergetest2$PerHerbWet,mergetest2$fptval)

#attr(fpt2M320,"radii")

fpt2M320["date"]

dfT<-as.data.frame(fpt2M320[1])
plot(dfT$r1)
head(mergetest)

newdf1 <- [['fpt2M320']]

## S3 method for class 'fipati'
test<-as.data.frame(fpt2M320[1])
head(test)
test2 <- test[!is.na(test$r1)]
test2<-test[complete.cases(test),]

test2<-lapply(list, function(x) test[!is.na(test)])

test<-fpt2M320[1]
str(test)
summary(test[1]$r1)
test2<-unlist(test[1])
str(test2)
summary(test2)
head(test2,300)


test<-unlist(fpt2M320[1])
head(test)
