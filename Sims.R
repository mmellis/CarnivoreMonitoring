library(raster)
library(rgdal)
library(rgeos)
library(Rcpp)

sourceCpp('./src/SPACE.cpp')

landscape<-raster(ext=extent(c(0, 50, 0, 200)), res=0.1, vals=1)
  grd<-make_grid(as(extent(landscape), 'SpatialPolygons'), type='square', cell_area=(2.5*2.5), drop=T)
  #gr<-getValues(rasterize(grd, landscape))

  plot(landscape)
  plot(grd, add=T)
  plot(grd25, add=T, border='blue')
  plot(grd)
  plot(x, add=T, border='red')  
  plot(gIntersection(grd,x,byid=T) , border='blue')
  
  Fisher<-list(
    N               = 250,                # Initial population size
    lmda            = 0.933,              # Population growth rate
    n_yrs           = 10,                 # Maximum number of years in simulation
    n_visits        = 5,                  # Maximum number of visits per year
    grid_size       = 25,                 # Cell size in grid
    MFratio         = c(0.6, 0.4),        # Ratio of types of individuals
    buffer          = c(3.75, 6.30),      # Distance between individual center locations
    moveDist        = c(2.5, 4.2),        # Movement radius
    moveDistQ       = c(0.9, 0.9),        # Proportion of time in radius
    maxDistQ        = c(0.68,0.68),       # Truncate movements above 1 SD
    habitat.cutoff  = 0.5,                # Minimum habitat value required for individual center locations
    turnover        = 0.35                # Turnover rate 
    )  


  Marten<-list(
    N               = 250,                # Initial population size
    lmda            = 0.933,              # Population growth rate
    n_yrs           = 10,                 # Maximum number of years in simulation
    n_visits        = 5,                  # Maximum number of visits per year
    grid_size       = 6.25,               # Cell size in grid
    MFratio         = c(0.36, 0.64),      # Ratio of types of individuals
    buffer          = c(1.16, 1.62),      # Distance between individual center locations
    moveDist        = c(0.58, 0.81),      # Movement radius
    moveDistQ       = c(0.5, 0.5),        # Proportion of time in radius
    maxDistQ        = c(0.98,0.98),         # Truncate movements above 1 SD
    habitat.cutoff  = 0.5,                # Minimum habitat value required for individual center locations
    turnover        = 0.25,               # Turnover rate
    sample.cutoff=0 
    )      
    
    # Place individuals...saving locations over 10 years
WolverineHabitat<-raster('C:/Users/martha.ellis/Documents/GitHub/rSPACE/rSPACE/inst/external/WolvHabitat_Bitterroot.tif')
    ext<-cbind(bbox(WolverineHabitat)[,1], bbox(WolverineHabitat)[,1]+c(50*10^3, 200*10^3))
map2<-raster(ext=extent(ext), res=100, crs=CRS('+proj=utm +ellps=WGS84 +zone=12 +units=m'), val=1) 

#grd<-makeGrid(map2, 6.25*10^6, type='square')
#  grd<-rasterize(grd, map2, field=1:length(grd))
#  grd<-getValues(grd)
#Example1<-encounter.history(map=map2, Parameters=Marten, showSteps=T, grid_layer=grd, n_cells=1599)
  
  system.time({  
(USE<-new_sample(coordinates(map2), 100, 3.62, rep(0,ncell(map2)), sample(ncell(map2))-1   ))
   })

sum(USE) 
dev.new()   
plot(map2); plot(SpatialPoints(coordinates(map2)[USE==1,]), add=T)
     # Loop over individuals, build useLayer for each individual, calc n & prob pres for each cell  
     
     
     
             