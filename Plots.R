 dta<-rbind(data.frame(read.table('./Scenario7/output/Marten_Scenario7_results.txt', header=T), Tails='1km'),
            data.frame(read.table('./Scenario7_PropTails/output/Marten_Scenario7_results.txt', header=T), Tails='Prop'))
            
            
#### Output graphs
# 1.) Estimated detection probability by n_visits and detP             
dev.new(width=10,height=7)
 ggplot(dta, aes(x=factor(n_visits), colour=factor(detP), y=p_est))+
  geom_boxplot()+
  facet_wrap(~Tails)
  

#2.) Occupancy over time maximum sampling        
dev.new(width=10, height=7)
  filter(dta, n_grid==max(n_grid), n_visits==10, detP==max(detP)) %>%
        select(Tails, detP, starts_with('X')) %>%
        gather('Yr','X',3:13) %>%
        mutate(Yr=as.numeric(sub('X','',Yr))) %>%
  ggplot(aes(x=factor(Yr), y=X, colour=Tails))+geom_jitter(alpha=0.2)+geom_boxplot()+ylim(0,1)+
  labs(title='Scenario 7 - Marten', subtitle='N = 150, lmda = 0.933, EA = 6.25, nRuns=50', 
       x='Year', y='Est Occupancy')

#3.) Occupancy over time medium sampling
dev.new(width=10, height=7)
 filter(dta, (n_grid>=100 & n_grid<=600), n_visits==10) %>%
        select(Tails, detP, starts_with('X')) %>%
        gather('Yr','X',3:13) %>%
        mutate(Yr=as.numeric(sub('X','',Yr))) %>%      
 ggplot(aes(x=factor(Yr), y=X, colour=Tails))+geom_jitter(alpha=0.2)+geom_boxplot()+ylim(0,1)+
  labs(title='Scenario 7 - Marten', subtitle='N = 150, lmda = 0.933, EA = 6.25, nRuns=50', 
       x='Year', y='Est Occupancy')
       
#4.) Power curves
dev.new(width=10, height=7)
dta %>% mutate(ind=as.numeric(-trend > CI*trendSE)) %>%
ggplot(aes(x=n_grid, y=ind, colour=factor(n_visits), linetype=factor(detP)))+
   stat_smooth(se=F, size=0.5)+
   stat_smooth(method='glm', method.args=list(family="binomial"), se=F)+
   labs(title=paste(sub('^..','', sc), '-', sp),
        subtitle=bquote(paste( 'N = ',.(Sc$N),'  ', 
                               lambda, ' = ', .(Sc$lmda),
                               '  ESA = ', .(Sc$ESA),
                               '  Runs = ', .(nRuns) )),
        colour='n_visits', linetype='detP')+
   scale_x_continuous(name='Number of grid cells sampled', limits=c(0,ifelse(sp=='Fisher',400,1600)))+
   scale_y_continuous(name='Power', limits=c(0,1))        
        
