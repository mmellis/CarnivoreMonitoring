### Analysis function for each observed encounter history       
WeaselFun3<-function(n_yrs, ch=NULL, n_visit=NULL, sample_yr=0, FPC=1, grdID=NULL, ...){
   require(RMark)
## Run with ch=NULL to set up output  
  sim_results<-data.frame(p_est=0, trend=0, trendSE=0, singular=0, matrix(0,1,n_yrs), nG=0, altM=0)
  if(is.null(ch)) return(sim_results)

## If ch has a value, do the rest of it...

## Subfunctions ################################################################
  tryN<-function(expr) tryCatch(expr, error=function(e) return(NULL))
  tryNA<-function(expr) tryCatch(expr, error=function(e) return(NA))          #
  tryW<-function(expr) suppressWarnings(tryN(expr))                            #
  tryM<-function(expr) suppressMessages(tryCatch(expr,                         #
          error=function(e){                                                   #
            if(grepl("mark.exe", e$message)){                                  #
              stop(e$message, call.=F)                                         #
            } else return(NULL)}))                                             #
                                                                               #
  time_int<-function(n_visit, n_yrs){                                          #
    tmp<-rep(0,n_visit)                                                        #
    tmp[n_visit] = 1                                                           #
    tmp<-rep(tmp,n_yrs)                                                        #
    return(tmp[-n_yrs*n_visit])                                                #
    }                                                                          #
                                                                               #
  FPC_trendSE<-function(Random.effects.model, k, FPC){                         #
    trendSE<-Random.effects.model$beta[2,2]                                    #
    process.variance<-Random.effects.model$sigma^2                             #
    sampling.variance<-k*trendSE^2 - process.variance                          #
    trendSE<-sqrt((process.variance+FPC*sampling.variance)/k)                  #
    return(trendSE)                                                            #
    }                                                                          #
    
drop_years <- function(ch, n_visits, dropvec=rep(c(F,T),                       #
                length.out=nchar(ch[1])/n_visits), samples=NULL){            #
  dropvec<-matrix(dropvec, ncol=length(dropvec),nrow=length(ch), byrow=T)      #
  if(!is.null(samples)){                                                       #
    dropvec<-as.matrix(samples[1:length(ch),]) #
  }                                                                            #

    ch_split<-do.call('rbind', strsplit(ch, split=''))                         #
    ch <- unlist(sapply(1:length(ch), function(x) {                            #
      ch_split[x,rep(!!(dropvec[x,]),each=n_visits)]="."                        #
      paste(ch_split[x,],collapse="")} ))                                     #

    if(!any(grepl('[[:digit::]', ch_split[,nchar(ch[1])-c(0:2)])))
      ch<-substr(ch, 1, nchar(ch)-n_visits)
  return(ch)} 
################################################################################

                                                               
  additional.args<-list(...)
    sample_matrix<-additional.args$sample_matrix

  hold<-expand.grid(nG=c(0.05,0.25,0.5,0.75,0.95), 
             altM=c(0,4,5)) 
     hold$nG<-jitter(hold$nG, amount = 0.1)
     hold$nG<-ifelse(hold$nG < 0.05, 0.05, ifelse(hold$nG > 0.95, 0.95, hold$nG))
     hold$nG<-round(length(ch)*hold$nG)
    chL<-sample(ch)
    chR<-lapply(1:(nchar(ch[1])/3), function(x) { substr(ch, seq(1,33,by=3)[x], seq(3,33, by=3)[x])})
      chR<-lapply(chR, sample)
      chR<-apply(do.call(rbind, chR),2, function(x) paste(x, collapse=''))
  
  for(i in 1:nrow(hold)){
     if(hold$altM[i] == 5){
        ch<-chR
     } else { ch<-chL } 
    ch<-ch[1:hold$nG[i]]
    sample_yr<-hold$altM[i]

  if(sample_yr == 0){
          mark_data<-data.frame(ch=ch,freq=rep(1,length(ch)),stringsAsFactors=F)
          test_processed=process.data(mark_data,model="RDOccupEG",time.intervals=time_int(n_visit,n_yrs))
          test_ddl=make.design.data(test_processed)
            test_ddl$Epsilon$eps=-1
            test_ddl$Gamma$eps=1
  
            p.session=list(formula=~session)
            Epsilon.random.shared=list(formula=~-1+eps:time, share=TRUE)
          model.parameters=list(Epsilon=Epsilon.random.shared,p=p.session)

          RDoccupancy=tryM(mark(test_processed,
                                test_ddl,
                                model.parameters=model.parameters,
                                silent=T, delete=T,output=F))
                             
          derived_psi <- tryN(RDoccupancy$results$derived[[1]][,1])
          derived_psi_vcv <-tryN(RDoccupancy$results$derived.vcv[[1]])
          P_est <-tryN(RDoccupancy$results$real$estimate[which(row.names(RDoccupancy$results$real)=="p g1 s1 t1")])
   
          Trend_DM=cbind(1,1:n_yrs)
          Random.effects.model<-tryW(var.components.reml(theta=derived_psi,design=Trend_DM,vcv=derived_psi_vcv))

   }  else if(sample_yr %in% c(4,5) ) {    #random sample

    ConvertToMatrix<-function(ch){
      ch<-strsplit(ch, split='')
      ch<-lapply(ch, function(x) as.numeric(x))
      ch<-do.call(rbind, ch)
    return(ch)}
         
    loglik <- function(param, h)
      {
        s   <- dim(h)[1] # nr sites
        k   <- dim(h)[2] # nr sampling occasions
        psi <- 1/(1+exp(-param[1]))  # to probabilities
        p   <- 1/(1+exp(-param[2]))
        d  	<- sum(sum(h)) # summary statistics
        Sd  <- sum(rowSums(h)>0)
        loglik <- -(Sd*log(psi)+d*log(p)+(k*Sd-d)*log(1-p)+(s-Sd)*log((1-psi)+psi*(1-p)^k))
        return(loglik)
      }
   
   fit<-function(h1){ 
    fm<-optim(par=runif(2), fn=loglik, h=h1, hessian=T)
      pars <- 1/(1+exp(-fm$par))    # to probabilities
      VC<-tryN(solve(fm$hessian, silent=T))
    if (is.null(VC))
      {
        SEs <- rep(NA,2)
      }else{
        SEs <- c(sqrt(diag(VC))*pars*(1-pars))
      }
    return(list(psi1=list(est=pars[1], se=SEs[1]), p1=pars[2])) }
            

   st<-seq(1,nchar(ch[1]), by=n_visit)
   en<-seq(n_visit, nchar(ch[1]), by=n_visit)
  h<-lapply(1:length(st), function(x) substr(ch, st[x], en[x]) ) 
   
   
   fits<-lapply(h, function(x) tryN(fit(ConvertToMatrix(x))))
   
   derived_psi<-sapply(1:n_yrs, function(x) fits[[x]]$psi1$est)
   derived_psi_vcv<-sapply(1:n_yrs, function(x) fits[[x]]$psi1$se)
    derived_psi_vcv<-diag(derived_psi_vcv, nrow=n_yrs)
   P_est<-fits[[1]]$p1
   
   Trend_DM = cbind(1,1:n_yrs)
   Random.effects.model<-tryW(var.components.reml(theta=derived_psi,design=Trend_DM,vcv=derived_psi_vcv))
 
   } 
   
   sim_results$nG[i]<-hold$nG[i]
   sim_results$altM[i]<-hold$altM[i]
   if(!is.null(Random.effects.model)){
     sim_results$trend[i]          <- Random.effects.model$beta[2,1]
     sim_results$trendSE[i]        <- FPC_trendSE(Random.effects.model, k=nrow(Trend_DM), FPC)
   }
   if(!is.null(derived_psi)){
     sim_results[i,grep('X',names(sim_results))] <- matrix(derived_psi,nrow=1)
     sim_results$p_est[i]          <- P_est
     sim_results$singular[i]       <- tryNA(length(RDoccupancy$results$singular))
  } }
   
   return(sim_results)
}
   