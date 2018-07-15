### To run: ####################################################################
##   1.) Set working directory
 #setwd("C:/Users/s91r448/Documents/WeaselSims/")
 setwd("C:/Users/martha.ellis/Documents/WeaselSims") 
##   2.) Make sure required packages are installed
  pck<-c('rgeos','sp','raster','Rcpp')
  pck<-pck[!(pck %in% installed.packages()[,"Package"])]
  if(length(pck)) install.packages(pck)
  rm(pck)
##   3.) Make sure R scripts are in working directory
  stopifnot(file.exists(c('./ScriptDir/grid.txt', 
                          './ScriptDir/UniformLandscape.tif',
                          './ScriptDir/scr/createGrid.R',
                          './ScriptDir/scr/buildUseLayers.R',
                          './ScriptDir/scr/sample_ind.R',
                          './ScriptDir/scr/use_surface.R')))
##   4.) Set scenario number to run
  sc=7
  nRuns=50    
##   5.) Execute entire file to R console  
################################################################################


#### Scenarios to simulatate for Weasel Team ###################################
Scenarios<-expand.grid(N=c(150,250,400), lmda=c(0.933,0.978), ESA=c(25,6.25,1.56,0.39))
  output_dir<-paste0('./Scenario',sc, '_1kTails')
  dir.create(output_dir)
  dir.create(paste(output_dir, 'output', sep='/')) 

## Parameters ##################################################################  
  Fisher<-list(
    N               = Scenarios$N[sc],    # Initial population size
    lmda            = Scenarios$lmda[sc], # Population growth rate
    n_yrs           = 10,                 # Maximum number of years in simulation
    n_visits        = 10,                 # Maximum number of visits per year
    grid_size       = 25,                 # Cell size in grid
    MFratio         = c(0.6, 0.4),        # Ratio of types of individuals
    buffer          = c(3.75, 6.30),      # Distance between individual center locations
    moveDist        = c(2.5, 4.2),        # Movement radius
    moveDistQ       = c(0.25, 0.25),      # Proportion of time in radius
    maxDistQ        = c(0.25, 0.25),      # Truncate movements above 1 SD
    habitat.cutoff  = 0.5,                # Minimum habitat value required for individual center locations
    turnover        = 0.35                # Turnover rate 
    )  
        
  Marten<-list(
    N               = Scenarios$N[sc],    # Initial population size
    lmda            = Scenarios$lmda[sc], # Population growth rate
    n_yrs           = 10,                 # Maximum number of years in simulation
    n_visits        = 10,                 # Maximum number of visits per year
    grid_size       = 6.25,               # Cell size in grid
    MFratio         = c(0.64, 0.36),      # Ratio of types of individuals
    buffer          = c(1.62, 1.16),      # Distance between individual center locations
    moveDist        = c(0.81, 0.58),      # Movement radius
    moveDistQ       = c(0.25, 0.25),      # Proportion of time in radius
    maxDistQ        = c(0.25, 0.25),      # Truncate movements above 1 SD
    habitat.cutoff  = 0.5,                # Minimum habitat value required for individual center locations
    turnover        = 0.25,               # Turnover rate
    sample.cutoff=0 
    )  

# Libraries ####################################################################
library(raster)
library(rgdal)
library(rgeos)
library(Rcpp)

source('./ScriptDir/scr/createGrid.R')
source('./ScriptDir/scr/buildUseLayers.R')
source('./ScriptDir/scr/sample_ind.R')
source('./ScriptDir/scr/use_surface.R')

# Landscape ####################################################################
MAP<-raster('./ScriptDir/UniformLandscape.tif')                                                   

# Grid #########################################################################
grd<-read.table('./ScriptDir/grid.txt', header=T)
  names(grd)[3:8]<-c(paste('Fisher', c(25,6.25,1.56), sep='_'),
                     paste('Marten', c(6.25,1.56, 0.39), sep='_'))                


# Species setup ################################################################
SPP<-list(list('Marten',1), list('Marten',2), list('Fisher',1), list('Fisher',2))
  if(Scenarios$ESA[sc]==25){
   SPP<-SPP[c(3,4)]
  } else if(Scenarios$ESA[sc]==0.39)
   SPP<-SPP[c(1,2)]
   
################################################################################
for(sp in 1:length(SPP)){   ################################ LOOP OVER SPECIES #
  spp<-SPP[[sp]] # Fisher or Marten - index for individualtype to use

  # Grid
  xyzg<-subset(grd, select=c('x','y','hab',paste(spp[[1]], Scenarios$ESA[sc], sep='_')))
  xyzg<-as.matrix(xyzg)
  if(spp[[1]] == 'Marten'){
     grd_names<-paste(rep(1:400,each=4), rep(1:4, times=400), 1:1600, sep='.')
  } else  grd_names = 1:400 

  # Parameters
  P<-lapply(get(spp[[1]]), function(x) ifelse(length(x)==1, x[1],x[spp[[2]]]))
   P$N<-P$N*P$MFratio  
   #P$maxD2<-P$moveDist[1] #Proportional = repeat moveDist  
   P$maxD2<- 1            #1#km
   P$MoveP<-local({  
              sd_xy<-solveSD(P$moveDistQ[1], P$moveDist[1], MAP)
              return(c(sd_xy,                               # sd_x and sd_y
                       qnorm((1+P$maxDistQ)/2),             # Truncation 1
                       P$maxD2*1000 / mean(sd_xy),          # Second distance converted to SDs
                       0.95)                                # Relative prob loss in tail
                       ) })   

  
  for(rn in 1:nRuns){    ################################ LOOP OVER REPLICATES #  
  for(yr in 0:P$n_yrs){  ################################ LOOP OVER YEARS ######
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
      yr_ch<-sapply(yr_mx[,2], function(x) paste0(rbinom(n=P$n_visits, size=1, prob=x),collapse=''))
    out<-data.frame(yr_mx[,1], yr_mx[,2], yr_ch)
         names(out)<-paste(c('N','enc','ch'),yr,sep='_')
    if(yr>0){
       OUT<-data.frame(OUT,out)     
   } else {  OUT<-out }   
  }                                                                # END YEARS #
        
  fname<-paste0('rSPACE_sc',sc,'_', paste0(spp, collapse=''),'_x',rn,'.txt')
  write.table(data.frame(grd=grd_names,OUT[-1,]), file=paste0(output_dir,'/',fname), row.names=F)
  cat(gsub('\\.txt|rSPACE_','',fname), P$Nout,'\n', file=paste0(output_dir, '/output/rSPACE_sc', sc, '_nTotal.txt'), append=T)  
}}                                                   # END REPLICATES, SPECIES #

