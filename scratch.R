Reg<-read.table('rSPACE_Marten1_x1.txt', header=T)
Tailed<-read.table('rSPACE_Tails_Marten1_x1.txt', header=T)
Eff50<-read.table('rSPACE_Eff50_Marten1_x1.txt', header=T)
Tail50<-read.table('rSPACE_Tail50_Marten1_x1.txt', header=T)

library(ggplot2)

table(Reg$N_0)
table(Tailed$N_0)
table(Eff50$N_0)
table(Tail50$N_0)

dd<-rbind(rbind(rbind(data.frame(grd=Tailed$grd, Reg[-1,], Tails=F, EffA=F), 
  data.frame(Tailed,  Tails=T, EffA=F)), 
  data.frame(grd=Tailed$grd, Eff50[-1,], Tails=F, EffA=T)), 
  data.frame(Tail50, Tails=T, EffA=T))

ggplot(subset(dd,enc_0>0), aes(x=enc_0))+geom_histogram()+facet_grid(Tails~EffA, labeller=label_both)


library(dplyr)

 files<-dir(pattern='nTotal')
outT<-function(fname, index=1){
  nfn<-read.table(fname, header=F,as.is=T)[index,]
   nfn<-list(nfn[1,1], gather(nfn[,-1],'V','N')$N)
   
  out<-read.table(dir(pattern=paste0('rSPACE_',nfn[[1]])), header=T)
  message(paste0('rSPACE_',nfn[[1]]))
  psi<-gather(select(out, starts_with("enc")), "Year", "enc")  %>%
    mutate(Year = as.numeric(sub('enc_', '', Year))) 
      
  sumPsi<- psi  %>%
    group_by(Year) %>%
    summarise(truePsi_1visit = sum(enc)/length(enc),
              truePsi_5visit = sum(1-(1-enc)^5)/length(enc),
              truePsi_Asymptotic = sum(enc>0)/length(enc) ) 
                         
    return(left_join(data.frame(Year=0:10, N=nfn[[2]]), sumPsi, by="Year"))                                           
    }
    
 
mx<-select(out, starts_with("enc"))
n <- nrow(mx)   # number of sites
Yrs <- ncol(mx)   # number of primary periods
J <- 5          # number of secondary periods

site <- 1:n
years <- data.frame(matrix(rep(c(paste0('0', 0:9), '10'), each=n), n, Yrs))
occasions <- data.frame(matrix(rep(1:(J*Yrs), each=n), n, J*Yrs))

y<-matrix(0, n, J*Yrs)
for(sit in 1:n)
  for(vis in 1:(J*Yrs))
    y[sit,vis] <- rbinom(1,1,prob=mx[sit,rep(1:Yrs,each=J)[vis]])

umf <-unmarkedMultFrame(y=y,
    yearlySiteCovs=list(year=years),
    numPrimary=Yrs)    
 summary(umf) ### 
 
 fm <- colext(psiformula=~1, 
               gammaformula =~1, 
               epsilonformula =~year-1, 
               pformula =~year-1, data=umf)  # fit a model 
  
