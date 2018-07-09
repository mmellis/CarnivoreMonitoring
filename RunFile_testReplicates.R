#combine individual types                                                      x
#randomize each encounter history file --- use set.seed() to start analyses    x
#write analysis function

#scatter number of grid cells selected

setwd("C:/Users/s91r448/Documents/WeaselSims/")

source('./ScriptDir/scr/Test_Samples.R')
source('./ScriptDir/scr/AnalysisFun.R')

sc='./Scenario10'
sp<-c('Fisher','Marten')[2]
nRuns=3
oFile<-paste0(sp,gsub('\\./','_',sc),'_results.txt')

set.seed(1)

PList<-list(n_yrs=11,
            n_visits=5,
            n_visit_test=c(3,5),                    
            detP_test = c(0.7,0.2),                              
            grid_sample=c(0.05,0.15,0.25,0.35,0.45,0.55,0.75,0.95), 
            alt_model=c(0))                                     
             
             
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
                 base.name='rSPACE_sc10_Marten_x', 
                 function_name="WeaselFun", jttr=T,
                 results.file=oFile,
                 skipConfirm=T) 


# 3.) Check out results
 library(dplyr)
 library(ggplot2)
 CI = qnorm(0.9)
 dta<-read.table(paste0(sc,'/output/',oFile), header=T)
  dta<-dta %>% mutate(ind=as.numeric(-trend > CI*trendSE))
  ggplot(dta, aes(x=n_grid, y=ind, 
     colour=factor(n_visits), linetype=factor(detP)))+
     #stat_smooth(se=F, size=0.5)+
     stat_smooth(method='glm', method.args=list(family="binomial"), se=F)
     
                                                                         