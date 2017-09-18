# Brice Randolph
# UCLA Statistics MS Thesis Code
# last edit: June 29 2017

# Data acquisition: 
# 1. Visit https://gis.water.ca.gov/app/gicima/
# 2. Under 'Select Layer Group', choose Fall 2016 Depth
# 3. Select 'Download'

# Load necessary packages:
#library(formatR)
library(readr)
library(FRK)
library(sp)
library(ggplot2)
library(ggmap)
library(geoR)
library(gstat)
library(INLA)
library(foreign)
library(rgdal)
library(gridExtra)
######
h20<-read.dbf("~/Downloads/Fall_2016_Depth_Points/F2016_DBGS_Points_20170327_102311.dbf")
#WSEL = groundwater surface elevation
#DGBS = groundwater depth below ground surface
#########
# Preprocessing 
h20<-h20[,c("Site_Code","WSEL","Latitude","Longitude","DGBS")]
h20=h20[unique(h20$Site_Code),] #remove duplicate reading from same site
h20<-h20[,c("Latitude","Longitude","DGBS")]

coordinates(h20) = ~Longitude+Latitude
proj4string(h20)=CRS("+proj=longlat") # Fixes distance calculations
h20<-remove.duplicates(h20)
h20Edit<-h20[h20$DGBS>1,] # Working with log transformed data with values > 1


# checking normality assumptions 
hist(log(h20Edit$DGBS),main = "Log Depth to Groundwater",xlab="")#Appears normal
lmodel<-lm(log(h20Edit$DGBS)~(Latitude+1),data=h20Edit)
h20DF<-data.frame(h20Edit)

summary(lmodel)
par(mfrow=c(1,2))
hist(h20Edit$DGBS,main="Histogram of depth to groundwater",xlab = "depth to groundwater (ft)")
hist(log(h20Edit$DGBS),main="Histogram of log(depth to groundwater)",xlab = "log(depth to groundwater)")
### looking at normality in detrended values
hist(lmodel$residuals,main="Histogram of detrended data",xlab="residuals")
qqnorm(lmodel$residuals);qqline(lmodel$residuals, col = 2)
###############
# Procedure for creating single train/test set to evaluate prediction
# Code at end contains the procedure for replicating this 10 different times
set.seed(1)
trainIndices<-sample(1:length(h20Edit),length(h20Edit)/4,replace = FALSE)
test<-h20Edit[trainIndices,]
train<-h20Edit[-trainIndices,]
#####
# Creating the Basic Aerial Units(BAUs) 
# Essentially discretizes the domain of interest

# @ arguments:  2D plane, BAU cellsize, grid (not hex),
# data around which to create BAUs, border buffer factor, 
#
GridBAUs1 <- auto_BAUs(manifold = plane(), 
                       cellsize = c(0.2,0.2), 
                       type = "grid", 
                       data = h20Edit, 
                       convex=-0.05, 
                       nonconvex_hull=FALSE)
h20Pts<-h20Edit
h20Pts$DGBS<-NULL
GridBAUs2<-BAUs_from_points(h20Pts)  
#this (GridBAUs2) has terrible performance 
#compared to the other BAU methods and is left out of the report
GridBAUs3<- auto_BAUs(manifold = plane(), 
                      cellsize = c(0.1,0.1),
                      type = "grid",
                      data = h20Edit,
                      convex=-0.05, 
                      nonconvex_hull=FALSE)
GridBAUs4<- auto_BAUs(manifold = plane(), 
                      cellsize = c(0.05,0.05), 
                      type = "grid", 
                      data = h20Edit,
                      convex=-0.05, 
                      nonconvex_hull=FALSE)

plot(GridBAUs1)
plot(GridBAUs2)
plot(GridBAUs3)
plot(GridBAUs4)

GridBAUs1$fs <- 1 # fine-scale variation at BAU level
GridBAUs2$fs <- 1
GridBAUs3$fs <- 1
GridBAUs4$fs <- 1

#Types of basis functions used: "bisquare", "Gaussian","exp", "Matern32"
G <- auto_basis(manifold = plane(), 
                data=h20Pts, 
                nres = 2, 
                type = "Matern32", 
                regular = 0) 
show_basis(G) + 
  coord_fixed() + 
  xlab("Longitude") + 
  ylab("Latitude") 
## Note: show basis assumes spherical distance functions when plotting
###
G2 <- auto_basis(manifold = plane(), 
                 data=h20Pts,
                 nres = 3, 
                 type = "bisquare", 
                 regular = 0) 
show_basis(G2) + 
  coord_fixed() + 
  xlab("(b)") + 
  ylab("Latitude") 
## Note: show basis assumes spherical distance functions when plotting
###

f1 <- log(DGBS) ~ 1+Latitude # formula for SRE model
f2 <- log(DGBS) ~ 1 # formula for SRE model
# @ parameters: formula, list of datasets, BAUs, basis functions, 
# estimation measurement error
S1 <- SRE(f = f1, 
          data = list(train), 
          BAUs = GridBAUs1, 
          basis = G, 
          est_error = TRUE,
          average_in_BAU = FALSE) 
# @ parameters: model, max. num EM iterations, tolerance at which EM is
# assumed to have converged, bool print log-likelihood at each iteration
S1 <- SRE.fit(SRE_model = S1, 
              n_EM = 100,
              tol = 0.1, 
              print_lik=TRUE) 
GridBAUs1 <- SRE.predict(SRE_model = S1, 
                         obs_fs = FALSE)


MSE_FRK1<-mean((log(test$DGBS)-over(test,GridBAUs1)$mu)^2) #returns the GridBAUs1 entry with closest pixel
MSE_FRK1
head(over(test,GridBAUs1)$Latitude)
head(test$Latitude)# sanity check

################ second FRK BAU setting
S2 <- SRE(f = f1, 
          data = list(train), 
          BAUs = GridBAUs2, 
          basis = G, 
          est_error = TRUE,
          average_in_BAU = FALSE) 

S2 <- SRE.fit(SRE_model = S2, 
              n_EM = 100, 
              tol = 0.01,
              print_lik=TRUE) 
GridBAUs2 <- SRE.predict(SRE_model = S2, 
                         obs_fs = FALSE)

MSE_FRK2<-mean((log(test$DGBS)-over(test,GridBAUs2)$mu)^2) #returns the GridBAUs1 entry with closest pixel
MSE_FRK2

################ third FRK BAU setting
S3 <- SRE(f = f1, 
          data = list(train), 
          BAUs = GridBAUs3, 
          basis = G, 
          est_error = TRUE, 
          average_in_BAU = FALSE) 

S3 <- SRE.fit(SRE_model = S3, 
              n_EM = 100, 
              tol = 0.01, 
              print_lik=TRUE) 
GridBAUs3 <- SRE.predict(SRE_model = S3, 
                         obs_fs = FALSE)


MSE_FRK3<-mean((log(test$DGBS)-over(test,GridBAUs3)$mu)^2) #returns the GridBAUs1 entry with closest pixel
MSE_FRK3
################################# Other kriging methods

### variogram for ordinary kriging
train_g <- gstat(id="log_dist", formula = log(DGBS)~1, data = train)
vg<-variogram(train_g)
vgRobust<-variogram(train_g,cressie=TRUE) # robust variogram calculation
plot(vg)
plot(vgRobust)
lm(vgRobust$gamma~vgRobust$dist) # sigma2e = 0.664
v.fit <- fit.variogram(vgRobust, vgm(1.5,"Sph",300,0.5))
v.fit
plot(vg,v.fit)
##### variogram for universal kriging
train_g_U <- gstat(id="log_dist", formula = log(DGBS)~1+Latitude, data = train)
vg_U<-variogram(train_g_U)
dir.vgm<-variogram(train_g_U, alpha=c(0,45,90,135))
plot(dir.vgm)
v.fit_U <- fit.variogram(vg_U, vgm(1.5,"Sph",400,0.5))
plot(vg_U,v.fit_U)


### ordinary and universal kriging
OK <- krige(id="logDist", formula=log(DGBS)~1, train, newdata=test, model = v.fit)
MSE_OK<-mean((log(test$DGBS)-OK$logDist.pred)^2) #returns the GridBAUs1 entry with closest pixel

UK <- krige(id="logDist", formula=log(DGBS)~1+Latitude, train, newdata=test, model = v.fit_U)
MSE_UK<-mean((log(test$DGBS)-UK$logDist.pred)^2) 

#UK2 is used to produce the analogous kriging map of california using OK and UK as opposed to FRK  
UK2 <- krige(id="logDist", formula=log(DGBS)~1+Latitude, train, newdata=GridBAUs3, model = v.fit_U)


#compare MSE of single trial:  Actual comparisons 
# are done after replicating the procedure 10 times (shown at end of code)
MSE_FRK1
MSE_FRK2
MSE_FRK3
MSE_OK
MSE_UK
MSE_OKL
###
###    Plotting the predictions and standard errors


BAUs_df <- as(GridBAUs3,"data.frame") # need to coerce the BAUs to a data frame for ggplot2 

g_grid2FRK <- ggplot() + 
  geom_tile(data=BAUs_df , 
            aes(Longitude,Latitude,fill=mu), 
            colour="light grey") + 
  scale_fill_distiller(palette="Spectral", 
                       name="pred.") + 
  coord_fixed() + 
  xlab("Longitude") + ylab("Latitude") + 
  theme_bw()
# the following can be added to plot data over the map
#geom_point(data=data.frame(h20Edit), # Plot data
#           aes(Longitude,Latitude,fill=log(DGBS)), # Colour <-> log(zinc)
#           colour="black", # point outer colour
#           pch=21, size=1)  # size of point
g_grid2FRK # to view the map

# Similar to above but with standard errors
g2 <- ggplot() + 
  geom_tile(data=BAUs_df,
            aes(Longitude,Latitude,fill=sqrt(var)),
            colour="light grey") +
  scale_fill_distiller(palette="BrBG",name = "s.e.",
                       guide = guide_legend(title="se")) +
  coord_fixed() +
  xlab("Longitude") + ylab("Latitude") + theme_bw()
g2
grid.arrange(g1, g2, ncol=2)
#############################################
# plotting the data over california
myLocation <- c(lon = -119, lat = 36)
###myMap = watercolor for just showing locations
###myMap2 = toner for contrasting color grad
myMap <- get_map(location=myLocation,
                 source="stamen", maptype="watercolor", crop=FALSE,zoom = 6)
myMap2 <- get_map(location=myLocation,
                  source="stamen", maptype="toner", crop=FALSE,zoom = 6)
ggmap(myMap2)+
  geom_point(aes(x = Longitude, y = Latitude,colour=log(DGBS)), size=1,data = data.frame(h20Edit),
             alpha = .5)+scale_colour_gradientn(colours = rainbow(4))
#### 

#  Plotting universal kriging predictions and standard errors 

UK_df <- as(UK2,"data.frame") #logDist.pred logDist.var
gUKp <- ggplot() + 
  geom_tile(data=UK_df , 
            aes(Longitude,Latitude,fill=logDist.pred), 
            colour="light grey") + 
  scale_fill_distiller(palette="Spectral", 
                       name="pred.") + 
  coord_fixed() + 
  xlab("Longitude") + ylab("Latitude") + 
  theme_bw()
gUKp

gUKse <- ggplot() + 
  geom_tile(data=UK_df,
            aes(Longitude,Latitude,fill=sqrt(logDist.var)),
            colour="light grey") +
  scale_fill_distiller(palette="BrBG",name = "s.e.",
                       guide = guide_legend(title="se")) +
  coord_fixed() +
  xlab("Longitude") + ylab("Latitude") + theme_bw()
gUKse
grid.arrange(g_grid2FRK, gUKp, ncol=2)

######################## 
# Replicating the procedure on 10 different test sets
# Seed changes 10 times 

mseListOK<-list()
mseListUK<-list()

mseListFRK1<-list()
time<-list()

for(i in 1:10){
  set.seed(i)
  trainIndices<-sample(1:length(h20Edit),length(h20Edit)/4,replace = FALSE)
  testi<-h20Edit[trainIndices,]
  traini<-h20Edit[-trainIndices,]
  
  start.timei<-Sys.time()
  S1i <- SRE(f = f1, 
             data = list(traini), 
             BAUs = GridBAUs1, 
             basis = G, 
             est_error = TRUE, 
             average_in_BAU = FALSE) 
  
  S1i <- SRE.fit(SRE_model = S1i, 
                 n_EM = 100, 
                 tol = 0.1, 
                 print_lik=TRUE) 
  GridBAUs1i <- SRE.predict(SRE_model = S1i, 
                            obs_fs = FALSE)
  
  end.timei<-Sys.time()
  time.takeni<-end.timei-start.timei
  MSE_FRK1i<-mean((log(testi$DGBS)-over(testi,GridBAUs1i)$mu)^2) #returns the GridBAUs1 entry with closest pixel
  
  
  time<-c(time,time.takeni)
  mseListFRK1<-c(mseListFRK1,MSE_FRK1i)
}
# results of first FRK procedure (larger BAUs)
mean(unlist(mseListFRK1))
sd(unlist(mseListFRK1))
mean(unlist(time))
sd(unlist(time))


####
mseListFRK3<-list()
timeFRK3<-list()

for(i in 1:10){
  set.seed(i)
  trainIndices<-sample(1:length(h20Edit),length(h20Edit)/4,replace = FALSE)
  testi<-h20Edit[trainIndices,]
  traini<-h20Edit[-trainIndices,]
  
  
  start.timei<-Sys.time()
  S3i <- SRE(f = f1, 
             data = list(traini), 
             BAUs = GridBAUs3, 
             basis = G2, 
             est_error = TRUE, 
             average_in_BAU = FALSE) 
  
  S3i <- SRE.fit(SRE_model = S3i, 
                 n_EM = 100, 
                 tol = 1, 
                 print_lik=TRUE) 
  GridBAUs3i <- SRE.predict(SRE_model = S3i, 
                            obs_fs = FALSE)
  end.timei<-Sys.time()
  time.takeni<-end.timei-start.timei
  MSE_FRK3i<-mean((log(testi$DGBS)-over(testi,GridBAUs3i)$mu)^2) #returns the GridBAUs1 entry with closest pixel
  
  
  timeFRK3<-c(timeFRK3,time.takeni)
  mseListFRK3<-c(mseListFRK3,MSE_FRK3i)
}
# results of finer grid model FRK
mean(unlist(timeFRK3))
sd(unlist(timeFRK3))
mean(unlist(mseListFRK3))
sd(unlist(mseListFRK3))


####
mseListFRK1<-list()
time<-list()

# Universal kriging and Ordinary kriging MSE are computed as before, using
# MSE_UK<-mean((log(test$DGBS)-UK$logDist.pred)^2) 
# for the different test sets

for(i in 1:10){
  set.seed(i)
  trainIndices<-sample(1:length(h20Edit),length(h20Edit)/4,replace = FALSE)
  testi<-h20Edit[trainIndices,]
  traini<-h20Edit[-trainIndices,]
  
  start.timei<-Sys.time()
  ##### variogram for universal kriging
  train_g_U <- gstat(id="log_dist", formula = log(DGBS)~1+Latitude, data = traini)
  vg_U<-variogram(train_g_U)
  v.fit_U <- fit.variogram(vg_U, vgm(1.5,"Sph",400,0.5))
  
  UK2 <- krige(id="logDist", formula=log(DGBS)~1+Latitude, traini, newdata=GridBAUs3, model = v.fit_U)
  
  end.timei<-Sys.time()
  time.takeni<-end.timei-start.timei
  time<-c(time,time.takeni)
}

mean(unlist(time)) # avg time = 32 sec. sd 19.6
sd(unlist(time))

# formatting for latex
#tidy_source(width.cutoff = 50)
