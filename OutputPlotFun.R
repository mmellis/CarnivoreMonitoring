plotOutput<-function(sc=NULL, opt=0, use=NULL){

  stopifnot({
    require(dplyr)
    require(ggplot2) 
    require(tidyr)})
  
  if(is.null(use)){
  if(is.null(sc)){
    use<-file.choose()
    print(use)
  } else { 
    use<-dir(pattern='_results\\.txt', path=sc, full.names=T, recursive=T)
    print(use) 
  } }
  
  if(opt<2){
    dev.new(width=10, height=7, record=T)
    par(ask=T)
  }
    
  for(i in use){  
  dta<-read.table(i, header=T)
  CI = qnorm(0.9)
 # pdf(width=10, height=7, onefile=T, file=paste0(dirname(use),'/plots.pdf'))
   # 1.) P_est versus number of visits

   P1 <- ggplot(dta, aes(x=factor(n_visits), colour=factor(detP), y=p_est))+
     geom_boxplot()+ ylim(0,0.7)+
     labs(x='n_visits', colour='detP', title=basename(i))+theme_bw()

   #2.) Est occupancy vs year by sampling instensity

   P2 <- dta %>%
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
      facet_grid(detP~SamplingIntensity)+theme_bw()+
      labs(x='Year', y='Est Occupancy', title=basename(i), colour='n_visits')
      
  #3.) Estimated trend
  nG<-max(dta$n_grid)/0.95
  ctr<-median(filter(dta, n_visits==10, detP==0.7, n_grid==max(n_grid))$trend, na.rm=T)
  P3<-dta %>% mutate(n_grid = case_when(
      n_grid == max(n_grid) ~ "0.95", 
      between(n_grid, .15*nG, .35*nG) ~ "0.25",
      between(n_grid, .7*nG, .8*nG) ~ "0.75",      
      between(n_grid, 0.45*nG, 0.55*nG) ~ "0.50",
      n_grid == min(n_grid) ~ "0.05", 
      T ~ "")) %>% filter(n_grid!="") %>%
  ggplot(aes(x=n_grid, y=trend))+theme_bw()+
    geom_hline(yintercept=0, linetype='dashed')+
    geom_boxplot()+
   geom_pointrange(aes(ymin=(trend-CI*trendSE), ymax=(trend+CI*trendSE)),
    position=position_jitter(width = 0.4, height = 0), shape=20, alpha=0.1)+
   facet_grid(detP~n_visits) +
     coord_cartesian(ylim = ctr+c(-0.1,0.1)) + 
   labs(x='Proportion of landscape sampled', y='Est trend', title=basename(i))
   
  #4.) Power curves
  P4<- dta %>% mutate(ind=as.numeric(-trend > CI*trendSE)) %>%
   ggplot(aes(x=n_grid, y=ind, colour=factor(n_visits), linetype=factor(detP)))+
   stat_smooth(se=F, size=0.5)+
   stat_smooth(method='glm', method.args=list(family="binomial"), se=F)+
   labs(colour='n_visits', linetype='detP', title=basename(i))+theme_bw()+
   scale_x_continuous(name='Number of grid cells sampled', limits=c(0,ifelse(max(dta$n_grid)>400, 1600,400)))+
   scale_y_continuous(name='Power', limits=c(0,1))

   if(opt<=2){
     print(P1) 
     print(P2)
     print(P3) 
     print(P4) }
   
   if(opt==3){
     grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 3)))
      vplayout <- function(x, y) grid::viewport(layout.pos.row = x, layout.pos.col = y)
        print(P1, vp = vplayout(1, 1))
        print(P2, vp = vplayout(1, 2:3))
        print(P3, vp = vplayout(2, 2:3))
        print(P4, vp = vplayout(2, 1))
        grid::grid.newpage()
   } }
   
   if(opt<2)  
    par(ask=F)
   
   if(opt==1){ return(dta) } else { return(invisible()) }
  } 

hold<-dir(pattern='._results\\.txt', full.names=T, recursive=T)
 hold<-hold[!grepl('RevGrid|PTails|KmTails',hold)]
 hold<-hold[order(as.numeric((gsub('^.*_Scenario|_results.*$', '', hold))))]
  pdf(width=12, height=7, onefile=T)    
     plotOutput(use=hold,opt=3)
  graphics.off()
    