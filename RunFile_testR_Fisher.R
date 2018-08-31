### To run: ####################################################################
##   0.) Download/install Program MARK

##   1.) Set working directory
 #setwd("C:/Users/s91r448/Documents/WeaselSims/")
 #setwd("C:/Users/martha.ellis/Documents/WeaselSims") 

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

###   4.) Set scenario number to run - Only run one at a time!!
   #(a) Normal runs
  nRuns=50                                              ### EDIT
  scenarios_to_test=paste0('./Scenario', 1:6)           ### EDIT
      PList<-list(n_yrs=11,                             ### EDIT
                n_visits=10,                            ### EDIT
                n_visit_test=c(3,5,10),                 ### EDIT   
                detP_test = c(0.2, 0.7),                ### EDIT              
                grid_sample=c(0.05,0.25,0.5,0.75,0.95), ### EDIT
                alt_model=c(0))                         # <- DO NOT EDIT
  FUN<-list("WeaselFun", T)
  addtext='' 
    
#   #(c) AltM runs :  Use this section to create a results file comparing:
#   #       - Sampling every year                    alt_model=0
#   #   vs  - Every other year                       alt_model=2
#   #   vs  - Every third year                       alt_model=3
#   #   vs  - 20\% missing observations              alt_model=1
#  nRuns=50                                                     ### EDIT
#  scenarios_to_test=paste0('./Scenario', c(8,11,14,17))        ### EDIT
#     PList<-list(n_yrs=11,                                     ### EDIT
#                n_visits=10,                                   ### EDIT
#                n_visit_test=c(5),                             ### EDIT
#                detP_test = c(0.2, 0.7),                       ### EDIT       
#                grid_sample=c(0.05,0.25,0.5,0.75,0.95),        ### EDIT
#                alt_model=c(0:3))                              # <- DO NOT EDIT
#  FUN<-list("WeaselFun2", T)
#  addtext='_altM' 
  
#   #(d) AltM=5 runs : Use this section to create a results file comparing:
#   #          - Visiting the SAME sites over time, with full variance-cov matrix
#   #  vs      - Visiting the SAME sites over time, with independent occupancy estimates each year
#   #  vs      - Visiting a random selection of sites each year, with independent occupancy estimates each year 
#   #  alt_model=5 for all, altM = 0 (Same sites, VCV); altM = 4 (Same sites, no VCV); altM = 5 (Rand, no VCV)
#  nRuns=50 
#  scenarios_to_test=paste0('./Scenario', c(8,14)) 
#     PList<-list(n_yrs=11,                                ### EDIT
#                n_visits=10,                              ### EDIT
#                n_visit_test=c(5),                        ### EDIT
#                detP_test = c(0.7),                       ### EDIT       
#                grid_sample=c(1.0),                       # <- DO NOT EDIT 
#                alt_model=c(5))                           # <- DO NOT EDIT 
#  FUN<-list("WeaselFun3", F) 
#  addtext='_altM5'   


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
  SPP<-'Fisher'
            
  for(j in 1:length(SPP)){
    sp<-SPP[j]
    
    iFile<-paste0('rSPACE_sc',scN,'_',sp,'_x')
    oFile<-paste0(sp,'_Scenario', scN,'_results', addtext, '.txt')
    
    
    # 2.) Analyze collated data files
    testReplicates(sc, PList, 
                     base.name=iFile, 
                     function_name=FUN[[1]], jttr=FUN[[2]],
                     results.file=oFile,
                     skipConfirm=T, overwrite=T, 
                     sample_matrix=MX, n_runs=nRuns) 

  }}
                                                                      