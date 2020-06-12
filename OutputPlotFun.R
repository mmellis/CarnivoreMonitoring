plotOutput<-function(sc=NULL, opt=0, use=NULL, singular.cutoff=5){

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

   P1 <- dta %>%
          mutate(n_visits = factor(n_visits),
                 detP = factor(detP, labels=c("low","moderate"))) %>%
   ggplot(aes(x=n_visits, colour=detP, y=p_est))+
     geom_boxplot()+ ylim(0,0.7)+
     labs(x='Number of visits', y='Estimated detection probability',
          colour='Detection effort')+theme_bw()  +
          scale_colour_manual(values=c('navy','#41b6c4'))

   #2.) Est occupancy vs year by sampling instensity

   P2 <- dta %>%
    mutate(SamplingIntensity=case_when(
      n_grid == min(n_grid)  ~ "Sampling: Low",
      n_grid == max(n_grid)  ~ "Sampling: High",
      between(n_grid, fivenum(dta$n_grid)[2], fivenum(dta$n_grid)[4])  ~ "Sampling: Med" ,
      TRUE ~ "drop"),
      detP = factor(detP, labels=c('Detection: Low', 'Detection: Moderate'))) %>% 
      filter(SamplingIntensity != 'drop', n_visits==5, singular < singular.cutoff) %>%
        select(SamplingIntensity, detP,n_visits, starts_with('X')) %>%
        gather('Yr','X',4:14) %>%
        mutate(Yr=as.numeric(sub('X','',Yr)),
               SamplingIntensity=factor(SamplingIntensity, 
                   levels=paste('Sampling:', c('Low','Med','High')))) %>%     
 ggplot(aes(x=factor(Yr), y=X))+
  geom_jitter(alpha=0.2)+
  geom_boxplot(outlier.shape=NA)+
    ylim(0,1)+
      facet_grid(detP~SamplingIntensity, drop=F)+
      theme_bw()+
      labs(x='Year', y='Est Occupancy', colour='n_visits')
      
      
  #3.) Estimated trend
  nG<-max(dta$n_grid)/0.95
  ctr<-median(filter(dta, n_visits==10, 
                          detP==0.7, 
                          n_grid==max(n_grid))$trend, na.rm=T)
  P3<-dta %>% mutate(n_grid = case_when(
      n_grid == max(n_grid) ~ "0.95", 
      between(n_grid, .15*nG, .35*nG) ~ "0.25",
      between(n_grid, .7*nG, .8*nG) ~ "0.75",      
      between(n_grid, 0.45*nG, 0.55*nG) ~ "0.50",
      n_grid == min(n_grid) ~ "0.05", 
      T ~ ""),
      detP = factor(detP, labels=c('Detection: Low', 'Detection: Moderate')),
      n_visits =factor(n_visits, labels=c('Visits: 3', 'Visits: 5', 'Visits: 10')) ) %>% 
      filter(n_grid!="", singular < singular.cutoff) %>%
  ggplot(aes(x=n_grid, y=trend))+theme_bw()+
    geom_hline(yintercept=0, linetype='dashed')+
    geom_boxplot()+
   geom_pointrange(aes(ymin=(trend-CI*trendSE), ymax=(trend+CI*trendSE)),
    position=position_jitter(width = 0.4, height = 0), shape=20, alpha=0.1)+
   facet_grid(detP ~ n_visits, drop=F)+
     coord_cartesian(ylim = ctr+c(-0.1,0.1)) + 
   labs(x='Proportion of landscape sampled', y='Est trend')
   
  #4.) Power curves
  P4<- dta %>% mutate(ind=as.numeric(-trend > CI*trendSE)) %>%
    filter(singular<singular.cutoff) %>%
   ggplot(aes(x=n_grid, y=ind, colour=factor(n_visits), 
               linetype=factor(detP, labels=c('low','moderate'))))+
   stat_smooth(method='loess', se=F)+
   #stat_smooth(method='glm', method.args=list(family="binomial"), se=F)+
   labs(colour='Visits', linetype='Detection')+theme_bw()+
   scale_x_continuous(name='Number of grid cells sampled', limits=c(0,ifelse(max(dta$n_grid)>400, 1600,400)))+
   scale_y_continuous(name='Power', limits=c(0,1))+
   scale_color_manual(values=c('#a1dab4','#41b6c4','#253494'))+
   scale_linetype_manual(values=c('dashed','solid'))

   if(opt<=2){
     print(P1) 
     print(P2)
     print(P3) 
     print(P4) }
   
   if(opt>=3){
     grid::pushViewport(grid::viewport(layout = grid::grid.layout(3, 3, heights = unit(c(0.5, 5, 5), "null"))))
      vplayout <- function(x, y) grid::viewport(layout.pos.row = x, layout.pos.col = y)
        grid.text(basename(i), gp=gpar(fontsize=20, fontface='bold'),
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
        print(P1, vp = vplayout(2, 1))
        print(P2, vp = vplayout(2, 2:3))
        print(P3, vp = vplayout(3, 2:3))
        print(P4, vp = vplayout(3, 1))
       if(opt==3) grid::grid.newpage()
   } }
   
   if(opt<2)  
    par(ask=F)
   
   if(opt==1){ return(dta) } else { return(invisible()) }
  } 

hold<-dir(pattern='._results\\.txt', full.names=T, recursive=T)
 hold<-hold[!grepl('CopyDir|RevGrid|PTails|KmTails',hold)]
 hold<-hold[order(as.numeric((gsub('^.*_Scenario|_results.*$', '', hold))))]
  pdf(width=12, height=7, onefile=T)    
     plotOutput(use=hold, opt=3, singular.cutoff=3)
  graphics.off()
  
   pdf(width=12, height=7, onefile=T)    
     plotOutput(sc='./Scenario19/',opt=3)
     plotOutput(sc='./Scenario21/',opt=3)
     plotOutput(sc='./Scenario22/',opt=3)
     plotOutput(sc='./Scenario23/',opt=3)
     plotOutput(sc='./Scenario24/',opt=3) 
  graphics.off()
    