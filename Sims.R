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

# Fisher grid ##################################################################
gF<-makeGrid(map2, 25*10^6, type='square')  #Fisher
  names(gF)<-paste0('V',1:length(gF))

# Marten grid ################################################################## 
gM<-local({  
     ptsF<-lapply(1:length(gF), function(x){
            rbind(gF[x]@polygons[[1]]@Polygons[[1]]@coords[1:4,],
                   gF[x]@polygons[[1]]@Polygons[[1]]@labpt)}) 
      
      ptsM<-lapply(1:length(ptsF), function(x){
        ord<-c(1,2,5,4,1)
        pts<-ptsF[[x]]
        pts<-expand.grid(x=sort(unique(pts[,1])), y=rev(sort(unique(pts[,2]))))
        newpts<-list(pts[ord,], pts[ord+1,], pts[ord+3,], pts[ord+4,])
        names(newpts)<-paste0('V',x,'.',1:4)   
        return(newpts)})
      ptsM<-unlist(ptsM, recursive=F,use.names=T)
      
     Srl<-lapply(1:length(ptsM), function(x) {
         Polygons(srl=list(Polygon(coords=ptsM[[x]])), 
         ID=paste(names(ptsM[x]),x,sep='.'))
       })
  
  return(SpatialPolygons(Srl, proj4string=CRS(proj4string(gF)))) 
})

# Reduced areas ################################################################
reduceArea_polygon<-function(gPoly, fc=0.5){
  editCoords<-lapply(1:length(gPoly), function(x){
                org=gPoly[x]@polygons[[1]]@Polygons[[1]]@coords
                 midX<-(max(org[,1])-min(org[,1]))*fc + min(org[,1])
                 midY<-(max(org[,2])-min(org[,2]))*fc + min(org[,2])
                org[org==max(org[,1])]<-midX
                org[org==max(org[,2])]<-midY
                return(org)  })
  Srl<-lapply(1:length(editCoords), function(x){
        Polygons(srl=list(Polygon(coords=unname(editCoords[[x]]))),
        ID=names(gPoly)[x]) })
  return(SpatialPolygons(Srl, proj4string=CRS(proj4string(gPoly)))) }

gM_B<-gM  
gM_C<-reduceArea_polygon(gM,fc=0.5)  
gM_D<-reduceArea_polygon(gM,fc=0.25)

gF_A<-gF
gF_B<-reduceArea_polygon(gF,fc=0.5)
gF_C<-reduceArea_polygon(gF,fc=0.25) 
  
  rF<-rasterize(gF, map2, field=1:length(gF), background=0)
    rF<-addLayer(rF, rasterize(gF_B, map2, field=1:length(gF), background=0),
                     rasterize(gF_C, map2, field=1:length(gF), background=0))
  rM<-rasterize(gM, map2, field=1:length(gM), background=0)
    rM<-addLayer(rM, rasterize(gM_C, map2, field=1:length(gM), background=0),
                     rasterize(gM_D, map2, field=1:length(gM), background=0))
  
write.table(data.frame(coordinates(map2), hab=rep(1, ncells(map2)), 
   Fisher=getValues(rF), Marten=getValues(rM)), file='grid.txt', row.names=F)


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
system.time({     
MX <-use_surfaceC(coordinates(map2)[USE,],                              #XY coordinates of wolverines
                 xyzg,                                                 #Map matrix with coords, habitat, and grid
                 n_grid=length(unique(c(0,grd)))+1,                           #Number of unique levels for grid 
                                                                       #(0 = excluded areas, +1 to count if no excluded areas)
                 sd_x=sd_xy[1], sd_y=sd_xy[2], trunc_cutoff=qnorm(1.7/2))
})


# Loop #########################################################################
nRuns<-5

xyzg<-as.matrix(read.table('grid.txt',header=T)[,c('x','y','hab','Marten.layer.1')])

spp<-list('Marten',1) # Fisher or Marten - And the index for individualtype to use
P<-lapply(get(spp[[1]]), function(x) ifelse(length(x)==1, x[1],x[spp[[2]]]))
 P$N<-P$N*P$MFratio  
 P$MoveP<-local({  
            sd_xy<-solveSD(P$moveDistQ[1], P$moveDist[1], map2)
            return(c(sd_xy,                               # sd_x and sd_y
                     qnorm((1+P$maxDistQ)/2),             # Truncation 1
                     P$maxD2*1000 / mean(sd_xy),          # Second distance converted to SDs
                     0.95)                                # Relative prob loss in tail
                     ) })   

if(spp[[1]] == 'Marten'){
   grd_names<-paste(rep(1:400,each=4), rep(1:4, times=400), 1:1600, sep='.')
} else { grd_names<-1:400 }

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
   
  yr_mx<-use_surfaceC(xyzg[loc,], 
                      xyzg, 
                      n_grid=length(unique(c(xyzg[,4],0))), 
                      MoveP=P$MoveP)
    yr_ch<-sapply(yr_mx[,2], function(x) paste0(rbinom(n=5, size=1, prob=x),collapse=''))
  out<-data.frame(yr_mx[,1], yr_mx[,2], yr_ch)
       names(out)<-paste(c('N','enc','ch'),yr,sep='_')
  if(yr>0){
     OUT<-data.frame(OUT,out)     
 } else {  OUT<-out }   
}

if( 0 %in% unique(c(xyzg[,4])) )
   OUT<-OUT[-1,]
      
fname<-paste0('rSPACE_Tails_', paste0(spp, collapse=''),'_x',rn,'.txt')
write.table(data.frame(grd=grd_names,OUT), file=fname, row.names=F)
cat(gsub('\\.txt|rSPACE_','',fname), P$Nout,'\n', file='rSPACE_Eff50_nTotal.txt', append=T)  
}
  
                   
   
     
             