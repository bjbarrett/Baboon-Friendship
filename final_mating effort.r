#Modeling starts line 48. First cycling (flat). Lactation line 185.
setwd("Z:/Vroni/Olive Baboons/analyses/DSI/")
me=read.table(file="mating effort info.csv",header=T,sep=",")
cf=read.table(file="input_parenting effort_cyclflat.csv",header=T,sep=",")
l=read.table(file="input_parenting effort_lac.csv",header=T,sep=",")


me$cf_DSI=cf$DSI[match(paste(me$inf,me$momid,me$partner),paste(cf$infid,cf$momid,cf$partner))]
me$cf_rDSI=cf$rDSI[match(paste(me$inf,me$momid,me$partner),paste(cf$infid,cf$momid,cf$partner))]
me$l_DSI=l$DSI[match(paste(me$inf,me$momid,me$partner),paste(l$inf,l$momid,l$partner))]

nrow(me)
471

#add smooth score_doc fro next infant conception
scores <- read.csv("Z:/Vroni/Olive Baboons/analyses/Elo/Elo scores_males_analysis across groups_Jan13Dec16_standardised.csv", header=T, stringsAsFactors=F, sep=",")

me$doc_next_infant=as.Date(me$doc_next_infant,format="%m/%d/%Y")
me$score_doc30=scores$s_scores_smooth30[match(paste(me$doc_next_infant,me$partner),paste(scores$Date,scores$ID))]s
me$score_doc=scores$s_scores[match(paste(me$doc_next_infant,me$partner),paste(scores$Date,scores$ID))]

#complictaed way of deriving from rank file whether male was present at doc
me$res_status_next=rep(NA,nrow(me))
for (i in 1:nrow(me)){
me$res_status_next[i]=ifelse(identical(scores[scores$Date==me$doc_next_infant[i] & scores$ID==me$partner[i], as.character(me$group_next[i])]==1,logical(0)),0,1)
}

me=me[!is.na(me$score_doc30),]
sum(is.na(me$res_status_next))


setwd("Z:/Vroni/Olive Baboons/analyses/DSI/")  
write.table(me, file = "input_mating effort.csv",sep = ",",row.names=F)

cf=me[!is.na(me$cf_rDSI),]
l=me[!is.na(me$l_DSI),]

cf=cf[cf$realdad!=1,]
l=l[l$realdad!=1,]

cf$dyad=paste(cf$momid,cf$partner)
l$dyad=paste(l$momid,l$partner)

write.table(cf, file = "input_mating effort_cf.csv",sep = ",",row.names=F)

write.table(l, file = "input_mating effort_lac.csv",sep = ",",row.names=F)
####################################################################
library(rethinking)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
setwd("Z:/Vroni/Olive Baboons/analyses/DSI/")
cf=read.table(file="input_mating effort_cf.csv",header=T,sep=",")


cf$partner=droplevels(cf$partner)
cf$male_index <- as.integer(as.factor(cf$partner))

cf$PHG <- ifelse(cf$group_next=="PHG" , 1 , 0 )


cf$s_score_doc=(cf$score_doc - mean(cf$score_doc))/sd(cf$score_doc)
cf$s_rDSI=(cf$cf_rDSI - mean(cf$cf_rDSI))/sd(cf$cf_rDSI)
#cf$s_nmales=(cf$nmales - mean(cf$nmales))/sd(cf$nmales)

##################################################

##################################################
plot(m_effort_cf6)
par(mfrow=c(1,1))

m_effort_cf6 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_rDSI + b_rank*s_score_doc + b_phg*PHG + (b_nmales + b_nmalesXdsi*s_rDSI)*nmales + a_id[male_index],
	c(a,b_rank,b_phg,b_dsi,b_nmales,b_nmalesXdsi) ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)
	
) , data=cf ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) #have only done short runs so far
plot(precis(m_effort_cf6))

plot(m_effort_cf5)
par(mfrow=c(1,1))

m_effort_cf5 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_rDSI + b_rank*s_score_doc + b_phg*PHG + b_nmales*nmales + a_id[male_index],
	c(a,b_rank,b_phg,b_dsi,b_nmales) ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)
	
) , data=cf ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) #have only done short runs so far
plot(precis(m_effort_cf5))

m_effort_cf4 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a +  b_rank*s_score_doc + b_phg*PHG + a_id[male_index],
	c(a,b_rank,b_phg) ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)
	
) , data=cf ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) #have only done short runs so far


m_effort_cf3 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_rDSI + b_phg*PHG + a_id[male_index],
	c(a,b_phg,b_dsi) ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)
	
) , data=cf ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) #have only done short runs so far


m_effort_cf2 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_rDSI + (b_rank + b_rankXdsi*s_rDSI)*s_score_doc + b_phg*PHG + a_id[male_index],
	c(a,b_rank,b_phg,b_dsi,b_rankXdsi) ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)
	
) , data=cf ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) #have only done short runs so far


#check trace plots
plot(m_effort_cf)
par(mfrow=c(1,1))

m_effort_cf <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_rDSI + b_rank*s_score_doc + b_phg*PHG + a_id[male_index],
	c(a,b_rank,b_phg,b_dsi) ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)
	
) , data=cf ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) #have only done short runs so far

output_cf=precis(m_effort_cf , depth=2 , digits=2)@output

#check trace plots
plot(m_effort_cf)
par(mfrow=c(1,1))


#Intercepts-only model
m_effort_cf_null <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,
	logit(p) <- a + a_id[male_index],
	a ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)	
) , data=cf ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE)

comp_cf=compare(m_effort_cf4,m_effort_cf3,m_effort_cf2,m_effort_cf,m_effort_cf_null)


#Check for whether the prior is too constraining

par(mar=c(3.5,2.5,1,1))
post=extract.samples(m_effort_cf)

norm=rnorm(1000,0,1)
dens(post$b_dsi,col="red",xlim=c(-5,5),xlab="posterior vs. prior")
dens(post$b_rank,add=T,col="blue")
dens(post$b_phg,add=T,col="orange")
dens(norm,add=T)
legend("topright",legend = c("dsi", "rank","phg","dnorm(0,1)"),col=c("red","blue","orange","black"),lty=1,lwd=2)

#Plot
par(mar=c(3.5,4,1,1))
a_male_z <- matrix(0,1000,length(unique(cf$male_index)))

dsi.seq <- seq(min(cf$s_rDSI),max(cf$s_rDSI),length.out=30)

d.pred <- data.frame(	
	male_index=rep(1,length(dsi.seq)),
	s_score_doc=rep(mean(cf$s_score_doc),length(dsi.seq)),
	s_rDSI=dsi.seq,
	PHG=rep(mean(cf$PHG),length(dsi.seq))	
	)

pred <- link( m_effort_cf , data=d.pred, n=1000 ,replace=
	    list(a_id=a_male_z))

pred.median <- apply(pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )

par(mar=c(3.5,2.5,1,1))
plot( nextdad ~ s_rDSI , data=cf , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1.2,xlim=c(min(cf$s_rDSI),max(cf$s_rDSI)))

lines( dsi.seq , pred.median ,lwd=2)
shade(pred.HPDI,dsi.seq,col=alpha("black",0.2))

xt=seq(0,10,by=1)
xxt=(xt-mean(cf$cf_rDSI))/sd(cf$cf_rDSI)
mtext(side=1,line=1.8,text="Dyadic sociality index - Cycling (flat)",cex=1)
mtext(side=2,line=1,text="Probability of siring next infant",cex=1)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1, cex.axis=1,at= xxt,labels=xt,line=-0.2,col=NA)
axis(2, tck=-0.015,cex=1.2,labels=NA,at=c(0,1))
axis(2, cex.axis=1,at= c(0,1),labels=c(0,1),line=-.2,col=NA)

#alternatively plot 100 random posterior predictions
plot( nextdad ~ s_rDSI , data=cf , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1.2,xlim=c(min(cf$s_rDSI),max(cf$s_rDSI)))
pred.lines=pred[sample(nrow(pred),100,replace=F),]
for (i in 1:100){
lines( dsi.seq , pred.lines[i,] ,lwd=1,col=alpha("black",0.5))
}
lines( dsi.seq , pred.median ,lwd=3,col="black")

xt=seq(0,10,by=1)
xxt=(xt-mean(cf$cf_rDSI))/sd(cf$cf_rDSI)
mtext(side=1,line=1.8,text="Dyadic sociality index - Cycling (flat)",cex=1)
mtext(side=2,line=1,text="Probability of siring next infant",cex=1)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1, cex.axis=1,at= xxt,labels=xt,line=-0.2,col=NA)
axis(2, tck=-0.015,cex=1.2,labels=NA,at=c(0,1))
axis(2, cex.axis=1,at= c(0,1),labels=c(0,1),line=-.2,col=NA)

#plot influence of rank
rank.seq <- seq(min(cf$s_score_doc),max(cf$s_score_doc),length=100)

d.pred <- data.frame(	
	male_index=rep(1,length(rank.seq)),
	s_score_doc=rank.seq,
	s_rDSI=rep(mean(cf$s_rDSI),length(rank.seq)),
	PHG=rep(mean(cf$PHG),length(rank.seq))
)

pred <- link( m_effort_cf , data=d.pred, n=1000 ,replace=
	list(a_id=a_male_z))

pred.median <- apply(pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )
par(mar=c(3.5,3,1,1))
plot( nextdad ~ s_score_doc , data=cf , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1.2#,xlim=c(-0.6,4.8))
)
lines( rank.seq , pred.median ,lwd=2)
shade(pred.HPDI,rank.seq,col=alpha("black",0.2))

############################################################################
##############################################################################
lac=read.table(file="input_mating effort_lac.csv",header=T,sep=",")

lac$partner=droplevels(lac$partner)
lac$male_index <- as.integer(as.factor(lac$partner))

lac$PHG <- ifelse(lac$group_next=="PHG" , 1 , 0 )

lac$s_score_doc=(lac$score_doc - mean(lac$score_doc))/sd(lac$score_doc)
lac$s_lDSI=(lac$l_DSI - mean(lac$l_DSI))/sd(lac$l_DSI)
lac$s_nmales=(lac$nmales - mean(lac$nmales))/sd(lac$nmales)

#################################
################################
m_effort_lac4 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_rank*s_score_doc + b_phg*PHG + a_id[male_index],
	c(a,b_phg,b_rank) ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)
	
) , data=lac ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE)

plot(precis(m_effort_lac4))

m_effort_lac3 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_lDSI + b_phg*PHG + a_id[male_index],
	c(a,b_phg,b_dsi) ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)
	
) , data=lac ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE)

plot(precis(m_effort_lac3))

m_effort_lac2 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_lDSI + (b_rank + b_rankXdsi*s_lDSI)*s_score_doc + b_phg*PHG + a_id[male_index],
	c(a,b_rank,b_phg,b_dsi,b_rankXdsi) ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)
	
) , data=lac ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE)

plot(precis(m_effort_lac2))

m_effort_lac <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_lDSI + b_rank*s_score_doc + b_phg*PHG + a_id[male_index],
	c(a,b_rank,b_phg,b_dsi) ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)
	
) , data=lac ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE)

plot(precis(m_effort_lac))

output_lac=precis(m_effort_lac , depth=2 , digits=2)@output

#intercept-only comparison
m_effort_lac_null <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,
	logit(p) <- a + a_id[male_index] ,
	a ~ dnorm(0,1),
	a_id[male_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)	
) , data=lac ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE)

comp_lac=compare(m_effort_lac4,m_effort_lac3,m_effort_lac2,m_effort_lac,m_effort_lac_null)

#check if prior is too limiting

par(mar=c(3.5,2,1,1))
post=extract.samples(m_effort_lac)

norm=rnorm(10000,0,1)
dens(post$b_dsi,col="red",xlim=c(-5,5),xlab="posterior vs. prior")
dens(post$b_rank,add=T,col="blue")
dens(post$b_phg,add=T,col="orange")
dens(norm,add=T)
legend("topright",legend = c("dsi", "rank","phg","dnorm(0,1)"),col=c("red","blue","orange","black"),lty=1,lwd=2)

#Plot
par(mar=c(3.5,3,1,1))
a_male_z <- matrix(0,1000,length(unique(lac$male_index)))

dsi.seq <- seq(min(lac$s_lDSI),max(lac$s_lDSI),length.out=30)

d.pred <- data.frame(	
	male_index=rep(1,length(dsi.seq)),
	s_score_doc=rep(mean(lac$s_score_doc),length(dsi.seq)),
	s_lDSI=dsi.seq,
	PHG=rep(mean(lac$PHG),length(dsi.seq))	
)

pred <- link( m_effort_lac3 , data=d.pred, n=1000 ,replace=
	list(a_id=a_male_z))



pred.median <- apply(pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )
par(mar=c(3.5,2,1,1))
plot( nextdad ~ s_lDSI , data=lac , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1.2)

lines( dsi.seq , pred.median ,lwd=2)
shade(pred.HPDI,dsi.seq,col=alpha("black",0.2))

xt=seq(0,9,by=1)
xxt=(xt-mean(lac$l_DSI))/sd(lac$l_DSI)
mtext(side=1,line=1.7,text="Dyadic sociality index - Lactation",cex=1)
mtext(side=2,line=1,text="Probability of siring next infant",cex=1)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1, cex.axis=1,at= xxt,labels=xt,line=-.3,col=NA)
axis(2, tck=-0.015,cex=1.2,labels=NA,at=c(0,1))
axis(2, cex.axis=1,at= c(0,1),labels=c(0,1),line=-.1,col=NA)

#alternatively plot 100 random predictions
plot( nextdad ~ s_lDSI , data=lac , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1.2)
pred.lines=pred[sample(nrow(pred),100,replace=F),]
for (i in 1:100){
lines( dsi.seq , pred.lines[i,] ,lwd=0.5,col=alpha("black",0.5))
}
lines( dsi.seq , pred.median ,lwd=2,col="black")

xt=seq(0,9,by=1)
xxt=(xt-mean(lac$l_DSI))/sd(lac$l_DSI)
mtext(side=1,line=1.7,text="Dyadic sociality index - Lactation",cex=1)
mtext(side=2,line=1,text="Probability of siring next infant",cex=1)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1, cex.axis=1,at= xxt,labels=xt,line=-.3,col=NA)
axis(2, tck=-0.015,cex=1.2,labels=NA,at=c(0,1))
axis(2, cex.axis=1,at= c(0,1),labels=c(0,1),line=-.1,col=NA)

#Plot rank
rank.seq <- seq(min(lac$s_score_doc),max(lac$s_score_doc),length.out=30)

d.pred <- data.frame(	
	male_index=rep(1,length(rank.seq)),
	s_score_doc=rank.seq,
	s_lDSI=rep(mean(lac$s_lDSI),length(rank.seq)),
	PHG=rep(mean(lac$PHG),length(rank.seq))
	)
	
pred <- link( m_effort_lac , data=d.pred, n=1000 ,replace=
	list(a_id=a_male_z))

pred.median <- apply(pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )
par(mar=c(3.5,3,1,1))
plot( nextdad ~ rank.seq, data=lac , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1.2#,xlim=c(-0.6,4.8))
)
lines( rank.seq , pred.median ,lwd=2)
shade(pred.HPDI,rank.seq,col=alpha("black",0.2))
#################################

for i in 1:max(lac$male_index){
	plot(s_lDSI ~ lac$ )
}
