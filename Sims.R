################################################################################
####################################################Carnivore Monitoring script#
#################################################### M.E. - June 2018          #
################################################################################
setwd("C:/Users/s91r448/Documents/GitHub/")

library(raster)
library(rgdal)
library(rgeos)
library(Rcpp)


#sourceCpp('./src/SPACE.cpp')
source('./rSPACE/rSPACE/R/createGrid.R')
source('./rSPACE/rSPACE/R/buildUseLayers.R')
source('./scr/sample_ind.R')
source('./scr/use_surface.R')


## Parameters ##################################################################  
  Fisher<-list(
    N               = 400,                # Initial population size
    lmda            = 0.933,              # Population growth rate
    n_yrs           = 10,                 # Maximum number of years in simulation
    n_visits        = 5,                  # Maximum number of visits per year
    grid_size       = 25,                 # Cell size in grid
    MFratio         = c(0.6, 0.4),        # Ratio of types of individuals
    buffer          = c(3.75, 6.30),      # Distance between individual center locations
    moveDist        = c(2.5, 4.2),        # Movement radius
    moveDistQ       = c(0.25, 0.25),        # Proportion of time in radius
    maxDistQ        = c(0.25, 0.25),       # Truncate movements above 1 SD
    habitat.cutoff  = 0.5,                # Minimum habitat value required for individual center locations
    turnover        = 0.35                # Turnover rate 
    )  


  Marten<-list(
    N               = 400,                # Initial population size
    lmda            = 0.933,              # Population growth rate
    n_yrs           = 10,                 # Maximum number of years in simulation
    n_visits        = 5,                  # Maximum number of visits per year
    grid_size       = 6.25,               # Cell size in grid
    MFratio         = c(0.64, 0.36),      # Ratio of types of individuals
    buffer          = c(1.62, 1.16),      # Distance between individual center locations
    moveDist        = c(0.81, 0.58),      # Movement radius
    moveDistQ       = c(0.25, 0.25),        # Proportion of time in radius
    maxDistQ        = c(0.25, 0.25),         # Truncate movements above 1 SD
    habitat.cutoff  = 0.5,                # Minimum habitat value required for individual center locations
    turnover        = 0.25,               # Turnover rate
    sample.cutoff=0 
    )      
    
# Landscape ####################################################################
WolverineHabitat<-raster('C:/Users/s91r448/Documents/GitHub/rSPACE/rSPACE/inst/external/WolvHabitat_Bitterroot.tif')
  WolverineHabitat<-projectRaster(WolverineHabitat, res=100, crs=CRS('+proj=utm +ellps=WGS84 +zone=12 +units=m'))
  ext<-cbind(bbox(WolverineHabitat)[,1], bbox(WolverineHabitat)[,1]+c(50*10^3, 200*10^3))
map2<-raster(ext=extent(ext), res=100, crs=CRS('+proj=utm +ellps=WGS84 +zone=12 +units=m'), val=1) 
  rm(ext, WolverineHabitat)                            

# Marten grid
grd_poly<-makeGrid(map2, 6.25*10^6, type='square')
  grd<-rasterize(grd_poly, map2, field=1:length(grd_poly))
  grd<-getValues(grd)

  plot(map2); plot(grd_poly, add=T)
# Add Fisher grid labels
  #To do -----------------------------------------------------------------------

# Place individuals ############################################################ 
#Test
sm<-sample(ncell(map2))
system.time({                                       
(USE<-new_sample(coordinates(map2), 3000, 6.32, rep(0,ncell(map2)), sm-1    ))
}) 
   
sum(USE)      # Number of individuals placed
 dev.set(2)     # Plot of individuals
plot(map2); plot(SpatialPoints(coordinates(map2)[USE==1,]), add=T,pch='.')

xyw<-coordinates(map2)[USE,] 
MX<-matrix(0,nrow(xyw),nrow(xyw))
for(i in 1:nrow(xyw)){
 for(j in 1:nrow(xyw)){
    MX[i,j]= sqrt((xyw[i,1]-xyw[j,1])^2+(xyw[i,2]-xyw[j,2])^2)}}
 dev.set(3)
 hist(apply(MX,2,function(x) min(x[x>0])))
 
 (which(apply(MX,2,function(x) min(x[x>0]))<6000)->wch)
    MX[wch,wch]
 

# Create encounters ############################################################     
xyzg<-cbind(coordinates(map2), rep(1,ncell(map2)), grd)
sd_xy<-solveSD(Marten$moveDistQ[1], Marten$moveDist[1], map2)


system.time({     
MX <-use_surfaceC(coordinates(map2)[USE,],                              #XY coordinates of wolverines
                 xyzg,                                                 #Map matrix with coords, habitat, and grid
                 n_grid=length(unique(grd))+1,                           #Number of unique levels for grid 
                                                                       #(0 = excluded areas, +1 to count if no excluded areas)
                 sd_x=sd_xy[1], sd_y=sd_xy[2], trunc_cutoff=qnorm(1.7/2))
})


# Loop #########################################################################
nRuns<-5
system.time({
spp<-list('Marten',1) # Fisher or Marten - And the index for individualtype to use
P<-lapply(get(spp[[1]]), function(x) ifelse(length(x)==1, x[1],x[spp[[2]]]))
 P$N<-P$N*P$MFratio                  

for(rn in 1:nRuns){     
for(yr in 0:P$n_yrs){
  if(yr==0){
     popN<-P$N[1]
     loc<-new_sample(xyzg[,1:2], popN, P$buffer, rep(0,nrow(xyzg)), sample(nrow(xyzg)-1))
     P$Nout[yr+1]<-popN<-sum(loc)
     message(sum(loc), ' fit at start')
  } else { 
     wch<-which(loc)                     #Mortality
     wch<-sample(wch, round((1-0.933)*popN))
     loc[wch]<-F
     
     turn<-round(P$turnover*length(wch))
     wch<-which(loc)                     #Turnover
     wch<-sample(wch, turn)
     loc[wch]<-F
     
     loc<-new_sample(xyzg[,1:2], turn, P$buffer, as.numeric(loc), sample(nrow(xyzg)-1))
     P$Nout[yr+1]<-popN<-sum(loc)
  }
 
  sd_xy<-solveSD(P$moveDistQ[1], P$moveDist[1], map2) 
  
  yr_mx<-use_surfaceC(xyzg[loc,], 
                      xyzg, 
                      n_grid=length(unique(c(grd,0))), 
                      sd_x=sd_xy[1], sd_y=sd_xy[2], qnorm((P$maxDistQ+1)/2))
    yr_ch<-sapply(yr_mx[,2], function(x) paste0(rbinom(n=5, size=1, prob=x),collapse=''))
  out<-data.frame(yr_mx[,1], yr_mx[,2], yr_ch)
       names(out)<-paste(c('N','enc','ch'),yr,sep='_')
  if(yr>0){
     OUT<-data.frame(OUT,out)     
 } else {  OUT<-out }   
}
fname<-paste0('rSPACE_', paste0(spp, collapse=''),'_x',rn,'.txt')
write.table(OUT, file=fname)
cat(gsub('\\.txt|rSPACE_','',fname), P$Nout,'\n', file='rSPACE_nTotal.txt', append=T)  
}
  })
                   
   
     
             