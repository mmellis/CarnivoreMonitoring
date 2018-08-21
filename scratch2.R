ch<-c("000100000000000100001000001000000", "011011110010101110101101111111110", 
"010000000000100001000001010010010", "000000000000000000000000000000000", 
"000000000101001100100011000001010", "000000000000000000000000000000000", 
"011110011010000000111000000000000", "000000000000000000000000000000000", 
"000000000000000000000000000000000", "001000000000000000000000000000000", 
"011100011011000000000100100000000", "111000000100010110000111000001100", 
"000000000000000000101010001010000", "000000000000000000000000000000000", 
"000000000000000000000000000000000", "111110100100000011111000011001100", 
"001010000000000000000000000000000", "010010000101000001011000111101011", 
"000000000000000000000000000000000", "101100001101100110101100011001100"
)

WeaselFun2(n_yrs=11, ch=ch, n_visit=3, sample_yr=1, FPC=1)
n_yrs=11; ch=ch, n_visit=3

dta<-read.table(file.choose(), header=T)
  dta$ind<-as.numeric(-dta$trend>qnorm(0.9)*dta$trendSE)
  dta$hold<-rep(1:nrow(dta), each=2)[1:nrow(dta)]
  
gather(dta, 'year', 'psi', X1:X11) %>% 
mutate(year=as.numeric(gsub('X','', year))) %>%
  ggplot(aes(x=year, y=psi, colour=factor(alt_model), group=interaction(alt_model, n_visits, n_grid, detP)))+geom_line()+geom_jitter() +stat_smooth(method='lm', se=F)
  
ggplot(dta, aes(x=hold, y=trend, ymax=trend+qnorm(0.9)*trendSE, ymin=trend-qnorm(0.9)*trendSE,
                 colour=factor(alt_model)))+
                 geom_pointrange(position=position_jitter(height=0, width=0.1))+
                 ylim(-.1,.1)
                 
### Shuffle ch and replace                 
 ch0<-lapply(1:(nchar(ch[1])/3), function(x) { substr(ch, seq(1,33,by=3)[x], seq(3,33, by=3)[x])})
 ch0<-lapply(ch0, sample)
 ch0<-apply(do.call(rbind, ch0),2, function(x) paste(x, collapse=''))