plotOutput<-function(sc=NULL){

  stopifnot({
    require(dplyr)
    require(ggplot2) 
    require(tidyr) })
  
  if(is.null(sc)){
    use<-file.choose()
    print(use)
  } else { 
    use<-dir(pattern='_results\\.txt', path=sc, full.names=T, recursive=T)
    print(use) 
  } 
  
  dev.new(width=10, height=7, record=T); par(ask=T) 
  for(i in use){  
  dta<-read.table(i, header=T)
  CI = qnorm(0.9)
 
 # pdf(width=10, height=7, onefile=T, file=paste0(dirname(use),'/plots.pdf'))
   # 1.) P_est versus number of visits
  print(
   ggplot(dta, aes(x=factor(n_visits), colour=factor(detP), y=p_est))+
     geom_boxplot()+
     labs(x='n_visits', colour='detP', title=basename(i))
   )  
   #2.) Est occupancy vs year by sampling instensity
   print(
   dta %>%
    mutate(SamplingIntensity=case_when(
      n_grid == min(n_grid)  ~ "low",
      n_grid == max(n_grid)  ~ "high",
      between(n_grid, fivenum(dta$n_grid)[2], fivenum(dta$n_grid)[4])  ~ "int" ,
      TRUE ~ "drop")
    ) %>% filter(SamplingIntensity != 'drop', n_visits==5) %>%
        select(SamplingIntensity, detP,n_visits, starts_with('X')) %>%
        gather('Yr','X',4:14) %>%
        mutate(Yr=as.numeric(sub('X','',Yr))) %>%     
 ggplot(aes(x=factor(Yr), y=X, colour=factor(n_visits)))+
  geom_jitter(alpha=0.2)+
  geom_boxplot()+
    ylim(0,1)+
      facet_grid(detP~SamplingIntensity)+
      labs(x='Year', y='Est Occupancy', title=basename(i), colour='n_visits')
      )
   
  #3.) Power curves
  print(
  dta %>% mutate(ind=as.numeric(-trend > CI*trendSE)) %>%
   ggplot(aes(x=n_grid, y=ind, colour=factor(n_visits), linetype=factor(detP)))+
   stat_smooth(se=F, size=0.5)+
   stat_smooth(method='glm', method.args=list(family="binomial"), se=F)+
   labs(colour='n_visits', linetype='detP', title=basename(i))+
   scale_x_continuous(name='Number of grid cells sampled', limits=c(0,ifelse(max(dta$n_grid)>400, 1600,400)))+
   scale_y_continuous(name='Power', limits=c(0,1))
   )
   }
 par(ask=F)
 return(dta)} 
