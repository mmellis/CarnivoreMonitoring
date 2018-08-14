# Relabel all mistake 0.978 runs to 0.933 values and move to 0.933 folders
sc<-c(4L, 5L, 6L, 10L, 11L, 12L, 16L, 17L, 18L, 22L, 23L, 24L)
new_sc<-c(1L, 2L, 3L, 7L, 8L, 9L, 13L, 14L, 15L, 19L, 20L, 21L)

setwd("C:/Users/martha.ellis/Documents/WeaselSims/PTails")

for(i in 2:length(sc)){
  ## Change replicate file labels
  fnames0<-dir(path=paste0('./Scenario',sc[i]), pattern='rSPACE_sc', full.names=T)
    fnames<-gsub(paste0('sc',sc[i]), paste0('sc',new_sc[i]), fnames0)
    repNo<-as.numeric(gsub('^.*\\_x|\\.txt$','',fnames))+50
    fnames<-gsub('_x.*\\.txt$','',fnames)
    fnames<-paste0(fnames,"_x",repNo, '.txt')
    fnames<-gsub(paste0('Scenario',sc[i]), paste0('Scenario',new_sc[i]), fnames)
  
  # Move replicate files to new scenario folder 
  file.rename(fnames0, fnames)
  
  ## Append n_total counts to previous version
  dta<-read.table(dir(path=paste0('./Scenario',sc[i],'/output/'), pattern='nTotal', full.names=T), header=F)
    V<-do.call('rbind',strsplit(paste(dta$V1), split='_'))
    V[,1]<-paste0('sc',new_sc[i])
    V[,3]<-paste0('x',as.numeric(gsub('x','',V[,3]))+50)
    dta$V1<-apply(V, 1, paste, collapse='_')
    
  # Copy to nTotal in new scenario
  write.table(dta, row.names=F, col.names=F, append=T, quote=F, 
   file=dir(path=paste0('./Scenario',new_sc[i],'/output/'), pattern='nTotal', full.names=T))
   
  ## Relabel scenario and run numbers in results file
  fnames0<-sort(dir(path=paste0('./Scenario',sc[i],'/output/'), pattern='results\\.txt', full.names=T))
  fnames1<-sort(dir(path=paste0('./Scenario',new_sc[i],'/output/'), pattern='results\\.txt', full.names=T))
  stopifnot(length(fnames0)==length(fnames1))
  
  for(ii in 1:length(fnames0)){
  dta<-read.table(fnames0[ii], header=T)
    V<-do.call('rbind',strsplit(paste(dta$rn), split='_'))
    V[,2]<-paste0('sc',new_sc[i])
    V[,4]<-paste0('x',as.numeric(gsub('x','',V[,4]))+50)
    dta$rn<-apply(V, 1, paste, collapse='_')

  # Copy to nTotal in new scenario
  write.table(dta, row.names=F, col.names=F, append=T, quote=F, file=fnames1[ii])  
  }}
  