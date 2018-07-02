### Exploring movement parameters ##############################################
plotMovement<-function(moveDist, moveDistQ, maxDistQ,dist2, relP, name=NULL){
  z<-seq(-4,4,by=0.01)
  dz<-dnorm(z)

  scl <- moveDist / qnorm((1+moveDistQ)/2)
  cutoff <- qnorm((1+maxDistQ)/2)

  dy<-ifelse(abs(z)>=cutoff, 0, dz)
  dy<-dy/sum(dy)
  lw<-(1-relP)*max(dy)
  
  distance<-z[dy>0]*scl
  distance<-c(rep(distance[1]-dist2,2), distance[1], distance, rev(distance)[1],rep(rev(distance)[1]+dist2,2))
  prob<-c(0,lw,lw,dy[dy>0],lw,lw,0)
  plot(distance, prob, type='l',xlab='Distance (km)', 
   main=paste0(name, ':   moveDist = ', moveDist, ',  moveQ = ', moveDistQ,',  maxQ = ', maxDistQ))                                 
}

pdf('FisherMarten_1kmTails.pdf')
  plotMovement(Fisher$moveDist[1], Fisher$moveDistQ[1], Fisher$maxDistQ[1],1,0.95, 'Fisher')
  plotMovement(Fisher$moveDist[2], Fisher$moveDistQ[2], Fisher$maxDistQ[2],1,0.95, 'Fisher')
  plotMovement(Marten$moveDist[1], Marten$moveDistQ[1], Marten$maxDistQ[1],1,0.95, 'Marten')
  plotMovement(Marten$moveDist[2], Marten$moveDistQ[2], Marten$maxDistQ[2],1,0.95, 'Marten')
dev.off()

pdf('FisherMarten_ProportionalTails.pdf')
  plotMovement(Fisher$moveDist[1], Fisher$moveDistQ[1], Fisher$maxDistQ[1],Fisher$moveDist[1],0.95, 'Fisher')
  plotMovement(Fisher$moveDist[2], Fisher$moveDistQ[2], Fisher$maxDistQ[2],Fisher$moveDist[2],0.95, 'Fisher')
  plotMovement(Marten$moveDist[1], Marten$moveDistQ[1], Marten$maxDistQ[1],Marten$moveDist[1],0.95, 'Marten')
  plotMovement(Marten$moveDist[2], Marten$moveDistQ[2], Marten$maxDistQ[2],Marten$moveDist[2],0.95, 'Marten')
dev.off()

# If maxQ and moveQ match, they determine the difference between the 
#    max probability and min probability. 
pC<-function(Q){ 
 C<-(dnorm(0)-dnorm(qnorm((1+Q)/2)))/dnorm(0)*100
 cat('%difference: ', C, '\n')}
