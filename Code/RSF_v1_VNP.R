```{r,echo=FALSE}
library(chron)
library(spatstat)
library(rgdal)
library(maptools)
library(raster)
library(adehabitatHR)
library(adehabitat)
detach("package:adehabitat", unload=TRUE)
library(KernSmooth)
library(ResourceSelection) 

```

###Read in data
```{r}
###orig wolf data
wolf<-read.csv("E:/Moose/Data/Voyageurs_National_Park/AllVNPWolf_MD102815.csv",header=TRUE)
table(as.factor(wolf$WOLFID),wolf$WOLFID)

###2015 Lotek Data
lowolf<-read.csv("E:/Moose/Data/Voyageurs_National_Park/VNPwolfdataMark112515/VNPwolfdataMark112515/Lotek111715.csv",header=TRUE)
###remove misses
lowolfC<-lowolf[lowolf$Latitude!=0,]
table(as.factor(lowolfC$WOLFID),lowolfC$WOLFID)

###2015 Lotek Data
vecwolf<-read.csv("E:/Moose/Data/Voyageurs_National_Park/VNPwolfdataMark112515/VNPwolfdataMark112515/Vect111715.csv",header=TRUE)
###remove misses
vecwolfC<-vecwolf[vecwolf$Latitude!=0,]
table(as.factor(vecwolfC$WOLFID),vecwolfC$WOLFID)

###subset orig data to just the ones not included in updated data
### set up all the attributes so they can all be combined with the same names
wolfO<-wolf
wolfOrigC<-wolfO[wolfO$WOLFID=="V003"|wolfO$WOLFID=="V007"|wolfO$WOLFID=="V008"|wolfO$WOLFID=="V010"|wolfO$WOLFID=="V011"|wolfO$WOLFID=="V012",]
unique(wolfOrigC$WOLFID)
head(wolfOrigC)
wolfvar<-c("OBJECTID","WOLFID","Latitude","Longitude","Altitude","Date","Time","Temp_C")
wolfOrigC<-wolfOrigC[(wolfvar)]
head(wolfOrigC)
wolfOrigC$Time <- chron(times=as.character(wolfOrigC$Time))
wolfOrigC$Date<- as.POSIXlt(wolfOrigC$Date, format= "%m/%d/%Y")
wolfOrigC$dtm <- as.factor(paste(wolfOrigC$Date, wolfOrigC$Time, sep = " "))
wolfOrigC$dtPm<-strptime(wolfOrigC$dtm, "%Y-%m-%d %H:%M:%S")
str(wolfOrigC$dtPm)
wolfOrigC$dtPm <- as.POSIXct(wolfOrigC$dtPm)



###keep all LOtek - change blanks to "V047, and set up attribute names
lowolfC$WOLFID<-as.character(lowolfC$WOLFID)
lowolfC$WOLFID[lowolfC$WOLFID==""] <- "V047"
lowolfC2<-lowolfC
head(lowolfC2)
###remove locations with a ton of error
lowolfC3<-lowolfC2[c(lowolfC2$DOP<15),]
lovars<-c("OBJECTID","WOLFID","Latitude","Longitude","Altitude","DateTimeGMT","Temp")
lowolfC3<-lowolfC3[(lovars)]
head(lowolfC3)
lowolfC3$DateTimeGMT<- as.POSIXlt(lowolfC3$DateTimeGMT, format= "%m/%d/%Y %H:%M")
lowolfC3$dtPm <- as.POSIXct(lowolfC3$DateTimeGMT, tz="GMT")
lowolfC3$dtPm<-format(lowolfC3$dtPm, tz="America/Chicago",usetz=TRUE)

### Vectronic data - do the same as above - different attribute names
###remove locations with a ton of error
vecwolfC2<-vecwolfC[c(vecwolfC$DOP<15),]
vecvars<-c("OBJECTID","WOLFID","Latitude","Longitude","Height","UTC_Date","UTC_Time","TempC")
vecwolfC3<-vecwolfC2[(vecvars)]
head(vecwolfC3)
#vecwolfC3$UTC_Time <- chron(times=as.character(vecwolfC3$UTC_Time))
vecwolfC3$UTC_Date<- as.POSIXlt(vecwolfC3$UTC_Date, format= "%m/%d/%Y ")
vecwolfC3$dtm <- as.factor(paste(vecwolfC3$UTC_Date, vecwolfC3$UTC_Time, sep = " "))
vecwolfC3$dtPm<-strptime(vecwolfC3$dtm, "%Y-%m-%d %H:%M:%S")
vecwolfC3$dtPm <- as.POSIXct(vecwolfC3$dtPm, tz="GMT")
vecwolfC3$dtPm<-format(vecwolfC3$dtPm, tz="America/Chicago",usetz=TRUE)
head(vecwolfC3)
###change names of attributes for combining into single data set
colnames(wolfOrigC) <- c("OBJECTID","WOLFID","Latitude","Longitude","Altitude","Date","Time","Temp","dtm","dtPm")
colnames(lowolfC3) <- c("OBJECTID","WOLFID","Latitude","Longitude","Altitude","DateTime","Temp","dtPm")
colnames(vecwolfC3) <- c("OBJECTID","WOLFID","Latitude","Longitude","Altitude","Date","Time","Temp","dtm","dtPm")

wolfvarF<-c("OBJECTID","WOLFID","Latitude","Longitude","Altitude","dtPm","Temp")
wolfOrigCF<-wolfOrigC[(wolfvarF)]
lowolfC3F<-lowolfC3[(wolfvarF)]
vecwolfC3F<-vecwolfC3[(wolfvarF)]

###combine to 1 dataset
Wolf1215<-rbind(wolfOrigCF,lowolfC3F,vecwolfC3F)
###order the data
Wolf1215 <- Wolf1215[order(Wolf1215$WOLFID,Wolf1215$dtPm),]
###number of locs by WOLFID
table(Wolf1215$WOLFID)

###create time attributes
Wolf1215$jul<-format(Wolf1215$dtPm,"%j")
Wolf1215$month<-format(Wolf1215$dtPm,format="%m")
Wolf1215$month<-as.numeric(Wolf1215$month)

###duplicate check
Wolf1215$DupCheck<-paste(Wolf1215$WOLFID,Wolf1215$jul,Wolf1215$Latitude,Wolf1215$Longitude,sep="")
which(duplicated(Wolf1215$DupCheck))
dups<-Wolf1215[duplicated(Wolf1215$DupCheck), ]

###add seasonality attribute
Wolf1215$season<-NA
Wolf1215$season[Wolf1215$month>=4 & Wolf1215$month<7]<- "Spring"
Wolf1215$season[Wolf1215$month>=7 & Wolf1215$month<=11]<- "Summer"
Wolf1215$season[Wolf1215$month==11 | Wolf1215$month==12 | Wolf1215$month==01 | Wolf1215$month==02 | Wolf1215$month==03 ]<- "Winter"


###calculate the diff (in seconds) between fixes
Wolf1215$time.diff<-c(0,as.numeric(diff(Wolf1215$dtPm)))
### remove the negative times (result from end of 1 data set to the next - we wanted to remove 1 day or so of data anyways after capture)
Wolf1215c<-subset(Wolf1215,Wolf1215$time.diff>0)

###here I am using notes from previous person and visuals to determine who belonged in packs together - this may need to be updated further
### I put the date as 8/5/2013 for V014 for all locs since didnt know it - this will be summer
### may want to exclude 023 and 032 from formal analysis (constant disperse) but still calc HR's
Wolf1215c_4hr<-Wolf1215c
#combine ID to new group/pack ID
Wolf1215c_4hr$PACKID<-NA
### 007 and 030 need to be combined
Wolf1215c_4hr$PACKID[Wolf1215c_4hr$WOLFID=="V030" | Wolf1215c_4hr$WOLFID=="V007"]<- "V0730"
###remove Wolf 28
Wolf1215c_4hr<-Wolf1215c_4hr[c(Wolf1215c_4hr$WOLFID!="V028"),] 
#wolf$WOLFID<-as.factor(wolf$WOLFID)
### combine 24 and 14
Wolf1215c_4hr$PACKID[Wolf1215c_4hr$WOLFID=="V024" | Wolf1215c_4hr$WOLFID=="V014"]<- "V1424"
### V011 7/7/13 to 11/12/13 and (separate territory) 1/30/14 to end of dataset 
Wolf1215c_4hr$PACKID[Wolf1215c_4hr$WOLFID=="V011"&Wolf1215c_4hr$dtPm<="2013-11-12" ] <- "V1101"
Wolf1215c_4hr$PACKID[Wolf1215c_4hr$WOLFID=="V011"&Wolf1215c_4hr$dtPm>"2013-11-12" ] <- "V1102"
### V029 7/22/14 to 11/14/14 and (separate territory) 12/16/14 to end of dataset 
Wolf1215c_4hr$PACKID[Wolf1215c_4hr$WOLFID=="V029"&Wolf1215c_4hr$dtPm<="2014-11-14" ] <- "V2901"
Wolf1215c_4hr$PACKID[Wolf1215c_4hr$WOLFID=="V029"&Wolf1215c_4hr$dtPm>"2014-11-14" ] <- "V2902"
### group V003 prior to 10/9/2012 with V029 main (first) territory
Wolf1215c_4hr$PACKID[Wolf1215c_4hr$WOLFID=="V003"&Wolf1215c_4hr$dtPm<="2012-10-09" ] <- "V2901"
Wolf1215c_4hr$PACKID[Wolf1215c_4hr$WOLFID=="V003"&Wolf1215c_4hr$dtPm>"2012-10-09" ] <- "V003"
### exclude V003 after 10/9/2012 - dispersal
#Wolf1215c_4hr$PACKID<-Wolf1215c_4hr[c(Wolf1215c_4hr$PACKID!="V003"),] 

Wolf1215c_4hr$PACKID<-as.character(Wolf1215c_4hr$PACKID)
Wolf1215c_4hr$WOLFID<-as.character(Wolf1215c_4hr$WOLFID)

Wolf1215c_4hr$PACKID <- ifelse(is.na(Wolf1215c_4hr$PACKID), Wolf1215c_4hr$WOLFID, Wolf1215c_4hr$PACKID)

###remaining locs per pack, wolf, 
table(Wolf1215c_4hr$PACKID)
table(Wolf1215c_4hr$WOLFID)


###combine season and packGRoup
Wolf1215c_4hr$packIDseas<-paste(Wolf1215c_4hr$PACKID,Wolf1215c_4hr$season,sep="_")
unique(Wolf1215c_4hr$packIDseas)
table(Wolf1215c_4hr$packIDseas)

#remove additional packs for various reasons (dispersal, few fixes etc)
### summer
#V023_Summer,V038_Summer,V035_Summer,V1101_Summer
Wolf1215c_4hr<-subset(Wolf1215c_4hr,Wolf1215c_4hr$packIDseas!="V023_Summer"&Wolf1215c_4hr$packIDseas!="V038_Summer"&Wolf1215c_4hr$packIDseas!="V035_Summer"&Wolf1215c_4hr$packIDseas!="V1101_Summer")
table(Wolf1215c_4hr$packIDseas)

### winter
#V023_Winter,V1424_Winter,V1102_Winter,V2902_Winter
Wolf1215c_4hr<-subset(Wolf1215c_4hr,Wolf1215c_4hr$packIDseas!="V023_Winter"&Wolf1215c_4hr$packIDseas!="V1424_Winter"&Wolf1215c_4hr$packIDseas!="V1102_Winter"&Wolf1215c_4hr$packIDseas!="V2902_Winter")
table(Wolf1215c_4hr$packIDseas)

###spring 
Wolf1215c_4hr<-subset(Wolf1215c_4hr,Wolf1215c_4hr$packIDseas!="V035_Spring"&Wolf1215c_4hr$packIDseas!="V1424_Spring"&Wolf1215c_4hr$packIDseas!="V1102_Spring"&Wolf1215c_4hr$packIDseas!="V2902_Spring")
table(Wolf1215c_4hr$packIDseas)

summary(Wolf1215c_4hr$time.diff)


#remove locations further than 12 hours and 90 secs from previous
#remove 39199-37903 (removes 1,296 fixes)
Wolf1215c_4hr2<-Wolf1215c_4hr[Wolf1215c_4hr$time.diff<=43290,]


###save and read in R file
saveRDS(Wolf1215c_4hr2,"E:/Moose/FPT/DataProc/Wolf_RSF_1003316.Rda")
Wolf1215c_4hr<-readRDS("E:/Moose/FPT/DataProc/Wolf_RSF_1003316.Rda")

table(Wolf1215c_4hr$packIDseas)

###project wolf locs into UTM from lat/long
Wolf1215c_4hr$XY<- project(cbind(Wolf1215c_4hr$Longitude, Wolf1215c_4hr$Latitude), "+proj=utm +zone=15 ellps=WGS84")
Wolf1215c_4hr$Y<-Wolf1215c_4hr$XY[,2]
Wolf1215c_4hr$X<-Wolf1215c_4hr$XY[,1]

Wolf1215c_4hr_SP<-SpatialPointsDataFrame(coords=Wolf1215c_4hr[c("X","Y")],proj4string= CRS("+proj=utm +zone=15 ellps=WGS84"),data=Wolf1215c_4hr)


```{r cache=TRUE}
use<-Wolf1215c_4hr[Wolf1215c_4hr$packIDseas=="V008_Summer" ,] # sp class
str(use)

use<-SpatialPoints(use[c("X","Y")])
plot(use)
str(use)
nuse<-nrow(use@coords)
nuse
```
# Lets generate some random points from within the `area used by Apollo.
# One option is to use an overly smoothed KDE of the UD (e.g., using KDE with h = href)
# to define availability

```{r cache=TRUE}
# UD estimate
UD.apollo<-kernelUD(use, h="href")

#UD.apollo <- kernelUD(SpatialPointsDataFrame(use[c("X","Y")], data=data.frame(id=as.factor(use$packIDseas))),h = "href", kern = c("bivnorm"))
image(UD.apollo)

# Outer 95% contour
ver <- getverticeshr(UD.apollo, 95)
plot(ver)
points(use) # Not big enough to encompass all points

```
# 
# Generate random points using the csr function in splancs library, but it would be 'best' to generate points on a uniform grid.
# 
# 
# Note: you can change nav to vary the ratio of available:used points to see what effect it has 
# on the regression parameters.

```{r }
nav<-1 # number of availble points per used point
navail<-nuse*nav
#  avail<-csr(ver@polygons[[1]]@Polygons[[1]]@coords,navail)
avail<-spsample(ver, n=navail, 'regular')
plot(ver)
points(avail, col="red")
points(use)
```

# Drop used points that fall outside off the 95% KDE. 

```{r }
kdecontour<-ver@polygons[[1]]@Polygons[[1]]@coords
inkde<-point.in.polygon(use@coords[,1],use@coords[,2], kdecontour[,1], kdecontour[,2])
use<-use[inkde==1,]
```
# Lets consider NDVI, Veg, distance to road, distance to river, with data from Apollo.  
# We need to extract covariate values for all of the available points

```{r }
avail<-avail@coords
avails<-data.frame(LONG_X=avail[,1], LAT_Y=avail[,2]) 

vegR<-raster("E:/Moose/Data/OutsideSpatial/NLCD11VNPclip1.tif") #raster: NLCD habitat
nlcdtab<-read.csv("E:/Moose/Data/OutsideSpatial/NLCD11_Table.csv",header=TRUE) #table so you can see what raster value = what hab class
nlcdtab[1:2]

lake30<-raster("E:/Moose/Data/Voyageurs_National_Park/lakedist30")
lakepnd30<-raster("E:/Moose/Data/Voyageurs_National_Park/lakepnddist30")

####distance to water body, at least .5HA for larger region - use this one!
waterboddist<-raster("E:/Moose/Data/OutsideSpatial/All_MN_Water/distwat_5ha2")


##snowmobile trails

snowmob<-readOGR("E:/Moose/Data/OutsideSpatial/NewSteveSnowmobileVNP",layer="All_DNR_VNPtrails") #shapfefiles: snowmobile trails
snowmob.psp = as.psp(snowmob)


##roads
roads<-readOGR("E:/Moose/Data/OutsideSpatial/Roads/Koochiching",layer="Koochiching_Roads")
roads.psp = as.psp(roads)



avails$veg<-extract(vegR, avails[,1:2])
avails$waterboddist<-extract(waterboddist, avails[,1:2])
#avails$lakepnddist<-extract(lakepnd30, avails[,1:2])

avails$distsnowmob<-nncross(avail, snowmob.psp)$dist 
avails$distroad<-nncross(avail, roads.psp)$dist 



used<-data.frame(use)
str(used)
used$veg<-extract(vegR, used[,1:2])
used$waterboddist<-extract(waterboddist, used[,1:2])
#used$lakepnddist<-extract(lakepnd30, used[,1:2])
str(use)

used$distsnowmob<-nncross(used, snowmob.psp)$dist 
used$distroad<-nncross(used, roads.psp)$dist 

use<-used




#avails$distroad<-nncross(avail, roads.pp)$dist
```
# 
# Now, create the response (y = 0 for available, 1 for used points).  For used points,
# let's create a data.frame called 'used' from the hyena dataset (which contains the used point coordinates, time/date,
# and also covariate values).  

```{r }
avails$y<-0 # response
used<-use
used<-used[inkde==1,]
used$y<-1
colnames(used)[1]<-"LONG_X"    
colnames(used)[2]<-"LAT_Y"    


used<-used[,names(avails)]  # only keep necessary columns...
```

```{r cache=TRUE}


# Merge the use and availability data sets for fitting logistic regression models  
# Create weights (weight availability points using a large number)

```{r }
alldat<-rbind(used,avails)
alldat$wt<-1
alldat$wt[alldat$y==0]<-10000
unique(alldat$y)
```

# Create an offset = log(area/$n_a$) so that the intercept in the logistic regression models is comparable to that
# of the point process model (See Warton and Shepherd (2010))

```{r cache=TRUE}
kde.owin<-owin(poly=list(x=rev(kdecontour[duplicated(kdecontour)!=T,1]), y=rev(kdecontour[duplicated(kdecontour)!=T,2])))
hrarea<-area.owin(kde.owin)
alldat$off<-log(hrarea/navail)
```
# Before fitting models, lets look at some plots to look at distribution of 
# covariates at used and available locations.

```{r fig.height=10, fig.width=10}
```{r cache=TRUE}
par(mfrow=c(2,2))
uriv<-bkde(used$lakepnddist, range.x=c(0,20000), gridsize=1000)
ariv<-bkde(avails$lakepnddist, range.x=c(0,20000), gridsize=1000)
plot(uriv, type="l", xlab="Distrince to nearest lake/pond", ylab="Density")
lines(ariv, lty=2)
legend("topright", lty=c(1,2), bty="n", c("Used","Available"))
uroad<-bkde(used$distsnowmob, range.x=c(0,5000), gridsize=500)
aroad<-bkde(avails$distsnowmob, range.x=c(0,5000), gridsize=500)
plot(uroad, type="l", xlab="Distrince to snowmob trail", ylab="Density")
lines(aroad, lty=2)
plot(uriv$x, uriv$y/ariv$y, type="l", xlab="Distance to lake/pond", ylab="Use/Avail")
plot(uroad$x, uroad$y/aroad$y, type="l", xlab="Distance to nearest snowmob trail", ylab="Use/Avail")
```
Can also look at a table of categorical values

```{r }
par(mfrow=c(1,1))
vegtab<-table(alldat$y, alldat$veg)
prop.table(vegtab,1)
plot(vegtab)
nlcdtab[1:2]
```
# Fit logistic regression models with and without weights  

```{r }
names(alldat)
fit1<-glm(y~lakepnddist+distsnowmob+as.factor(veg)+offset(off), data=alldat, family=binomial())
summary(fit1)

fit1w<-glm(y~lakepnddist+distsnowmob+as.factor(veg)+offset(off), weights=wt,data=alldat, family=binomial())
summary(fit1w)
```
# Fit ppp model using the same use and availability points

```{r }
used.ppp<-ppp(used[,1], used[,2],window=as.owin(kde.owin))
Q <- quadscheme(data=used.ppp, dummy=list(x=avails[,1], y=avails[,2]))
df <- data.frame(lakepnddist=c(used$lakepnddist, avails$lakepnddist), veg=c(used$veg, avails$veg),
                 distsnowmob=c(used$distsnowmob, avails$distsnowmob))

ppfit<- ppm(Q, ~ distsnowmob+lakepnddist+as.factor(veg), Poisson(), covariates=df)
summary(ppfit)

# Lastly, rsf fit from Lele & Keim 

```{r }

rsffit<-rsf(y~distsnowmob+lakepnddist+as.factor(veg)+offset(off), data=alldat, m=0)
str(rsffit)
summary(rsffit)

coef(rsffit)
rsffit$std.error[1]
cfH<-data.frame(summary(rsffit)$Std.Error)

```

# Now, compare the coefficients from the different fitted models

```{r }
coefs<-cbind(coef(fit1), coef(fit1w), coef(ppfit), c(NA,coef(rsffit)))
colnames(coefs)<-c("LR no wts", "LR wt", "PPM", "rsf") 
rownames(coefs)<-names(coef(fit1))             
outcoefs<-print(format(coefs,digits=3, scientific=T), quote=F)



```

# Possible diagnostics using output of ppm model  
Some diagnostic plots

```{r fig.width=10, fig.height=10}
diagnose.ppm(ppfit)
te<-quadrat.test(ppfit)
te
plot(te)
```





###project wolf locs into UTM from lat/long
Wolf1215c_4hr$XY<- project(cbind(Wolf1215c_4hr$Longitude, Wolf1215c_4hr$Latitude), "+proj=utm +zone=15 ellps=WGS84")
Wolf1215c_4hr$Y<-Wolf1215c_4hr$XY[,2]
Wolf1215c_4hr$X<-Wolf1215c_4hr$XY[,1]

str(Wolf1215c_4hr)


###remove pack-seasons that don't have at least 30 fixes
Wolf1215c_4hrC<-Wolf1215c_4hr[Wolf1215c_4hr$packIDseas!="V025_Summer"&Wolf1215c_4hr$packIDseas!="V033_Spring"&Wolf1215c_4hr$packIDseas!="V038_Winter"&Wolf1215c_4hr$packIDseas!="V047_Summer"&Wolf1215c_4hr$packIDseas!="V1102_Summer",]


str(Wolf1215c_4hrC_win$packIDseas)


unique(Wolf1215c_4hrC_win$packIDseas)

###create kernel UD home ranges by pack-season
Wolf1215c_4hrC_spr<-Wolf1215c_4hrC[Wolf1215c_4hrC$season=="Spring",]
Wolf1215c_4hrC_sum<-Wolf1215c_4hrC[Wolf1215c_4hrC$season=="Summer",]
Wolf1215c_4hrC_win<-Wolf1215c_4hrC[Wolf1215c_4hrC$season=="Winter",]

str(Wolf1215c_4hrC_spr)
unique(Wolf1215c_4hrC_spr$PACKID)

str(Wolf1215c_4hrC_spr_SP)

Wolf1215c_4hrC_spr<-subset(Wolf1215c_4hrC_spr, select=-c(X_UTM))
Wolf1215c_4hrC_sum<-subset(Wolf1215c_4hrC_sum, select=-c(X_UTM))
Wolf1215c_4hrC_win<-subset(Wolf1215c_4hrC_win, select=-c(X_UTM))

