### To run: ####################################################################
##   0.) Download/install Program MARK

##   1.) Set working directory
 #setwd("C:/Users/s91r448/Documents/WeaselSims/")
 setwd("C:/Users/martha.ellis/Documents/WeaselSims") 

##   2.) Make sure required packages are installed
  pck<-c('RMark')
  pck<-pck[!(pck %in% installed.packages()[,"Package"])]
  if(length(pck)) install.packages(pck)
  rm(pck)

##   3.) Make sure R scripts are in working directory
  stopifnot(file.exists(c('./ScriptDir/scr/Test_Samples.R', 
                          './ScriptDir/scr/AnalysisFun.R',
                          './ScriptDir/scr/WeaselFun.R',
                          './ScriptDir/scr/WeaselerFun.R')))

##   4.) Set scenario number to run
   #(a) Normal runs
  nRuns=50 
  scenarios_to_test=paste0('./Scenario', 1:6) 
  species_to_test = lapply(1:length(scenarios_to_test), function(x) c('Fisher','Marten')) 
      PList<-list(n_yrs=11,
                n_visits=10,
                n_visit_test=c(3,5,10),                    
                detP_test = c(0.2, 0.7),                              
                grid_sample=c(0.05,0.25,0.5,0.75,0.95), 
                alt_model=c(0))
  FUN<-list("WeaselFun", T)
  addtext='' 
 
   #(b) RevGrid
  nRuns=50 
  scenarios_to_test=paste0('./Scenario', c(2,8,14,20), '_RevGrid')[1] 
  species_to_test = lapply(1:length(scenarios_to_test), function(x) c('Fisher','Marten')) 
      PList<-list(n_yrs=11,
                n_visits=10,
                n_visit_test=c(3,5,10),                    
                detP_test = c(0.2, 0.7),                              
                grid_sample=c(0.05,0.25,0.5,0.75,0.95), 
                alt_model=c(0))
  FUN<-list("WeaselFun", T)
  addtext=''  
    
   #(c) AltM runs
  nRuns=50 
  scenarios_to_test=paste0('./Scenario', c(8,11,14,17)) 
  species_to_test = list('Fisher','Fisher','Marten','Marten')
     PList<-list(n_yrs=11,
                n_visits=10,
                n_visit_test=c(5),                    
                detP_test = c(0.2, 0.7),                              
                grid_sample=c(0.05,0.25,0.5,0.75,0.95), 
                alt_model=c(0:3))
  FUN<-list("WeaselFun2", T)
  addtext='_altM' 
  
   #(d) AltM=5 runs
  nRuns=50 
  scenarios_to_test=paste0('./Scenario', c(8,14)) 
  species_to_test = list('Fisher','Marten')
     PList<-list(n_yrs=11,
                n_visits=10,
                n_visit_test=c(5),                    
                detP_test = c(0.7),                              
                grid_sample=c(1.0), 
                alt_model=c(5)) 
  FUN<-list("WeaselFun3", F) 
  addtext='_altM5'   
##   5.) Execute entire file to R console  

################################################################################

# Libraries ####################################################################
library(RMark)

source('./ScriptDir/scr/Test_Samples.R')
source('./ScriptDir/scr/AnalysisFun.R')
source('./ScriptDir/scr/WeaselFun.R')
source('./ScriptDir/scr/WeaselerFun.R')

MX<-matrix(rbinom(p=0.2, n=1600*11, size=1),1600,11)
################################################################################ 
for(i in 1:length(scenarios_to_test)){   
  sc=scenarios_to_test[i]
   scN<-gsub('^.*Scenario|_.*$','', sc)
  SPP<-species_to_test[[i]]
  SPP<-SPP[sapply(SPP, function(x) any(grepl(x, dir(path=sc))))]
  stopifnot(length(SPP)>0)

            
  for(j in 1:length(SPP)){
    sp<-SPP[j]
    
    iFile<-paste0('rSPACE_sc',scN,'_',sp,'_x')
    oFile<-paste0(sp,'_Scenario', scN,'_results', addtext, '.txt')
    
    set.seed(1)                                                     
                 
    # 1.) Collate encounter histories across individual types
    for(rn in 1:nRuns){
      outputfiles<-dir(path=sc, pattern=paste0(sp,'._x',rn), full.names=T) 
      dta1<-read.table(outputfiles[1], header=T,colClasses="character")[,seq(1,34,by=3)]
        ch1<-strsplit(apply(dta1[,-1],1,paste, collapse=''),'')                              
      dta2<-read.table(outputfiles[2], header=T,colClasses="character")[,seq(1,34,by=3)]
        ch2<-strsplit(apply(dta2[,-1],1,paste, collapse=''),'')
      
      dta<-data.frame(grd=dta1$grd, V1=1:nrow(dta1), V2=NA, ch=sapply(1:nrow(dta1), 
             function(x) paste(as.numeric(grepl('1',ch1[[x]]) | grepl('1',ch2[[x]])),collapse='')))
      write.table(dta, file=gsub(paste0(sp,'.'), sp, outputfiles[1]), row.names=F, col.names=F)
      rm(dta1,dta2,ch1,ch2,dta,outputfiles)
    }
    
    # 2.) Analyze collated data files
    testReplicates(sc, PList, 
                     base.name=iFile, 
                     function_name=FUN[[1]], jttr=FUN[[2]],
                     results.file=oFile,
                     skipConfirm=T, overwrite=T, 
                     sample_matrix=MX, n_runs=nRuns) 

  }}
  
## 3.) Check out results
# library(dplyr)
# library(ggplot2)
#
### Copy from above 
#sc='./Scenario7'
#sp<-c('Fisher','Marten')[1]
#nRuns=50 
#iFile<-paste0('rSPACE_sc',sub('^..Scenario','',sc),'_',sp,'_x')
#oFile<-paste0(sp,gsub('\\./','_',sc),'_results.txt')
#
#
# Scenarios<-expand.grid(N=c(150,250,400), lmda=c(0.933,0.978), ESA=c(25,6.25,1.56,0.39))
# Sc<-Scenarios[as.numeric(sub('^..Scenario','',sc)),]
# CI = qnorm(0.9)
#
# dta<-read.table(paste0(sc,'/output/',oFile), header=T)
#  dta<-dta %>% mutate(ind=as.numeric(-trend > CI*trendSE))
#
#dev.new(width=10, height=7)
#ggplot(dta, aes(x=n_grid, y=ind, 
#   colour=factor(n_visits), linetype=factor(detP)))+
#   stat_smooth(se=F, size=0.5)+
#   stat_smooth(method='glm', method.args=list(family="binomial"), se=F)+
#   labs(title=paste(sub('^..','', sc), '-', sp),
#        subtitle=bquote(paste( 'N = ',.(Sc$N),'  ', 
#                               lambda, ' = ', .(Sc$lmda),
#                               '  ESA = ', .(Sc$ESA),
#                               '  Runs = ', .(nRuns) )),
#        colour='n_visits', linetype='detP')+
#   scale_x_continuous(name='Number of grid cells sampled', limits=c(0,ifelse(sp=='Fisher',400,1600)))+
#   scale_y_continuous(name='Power', limits=c(0,1)) 
     
                                                                         