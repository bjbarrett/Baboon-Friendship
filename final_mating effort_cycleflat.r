#Models start line 48.
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
me$score_doc30=scores$s_scores_smooth30[match(paste(me$doc_next_infant,me$partner),paste(scores$Date,scores$ID))]
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
library(Cairo)												
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
setwd("Z:/Vroni/Olive Baboons/analyses/DSI/")
cf=read.table(file="input_mating effort_cf_red.csv",header=T,sep=",")

setdiff(cf$momid,cf$partner)
cf$dyad=droplevels(cf$dyad)
cf$partner=droplevels(cf$partner)
cf$momid=droplevels(cf$momid)
c$male_index <- as.integer

cf$PHG <- ifelse(cf$group_next=="PHG" , 1 , 0 )
cf$dyad_index <- as.integer(as.factor(cf$dyad))
cf$prt1_index <- as.integer(as.factor(cf$momid))
cf$prt2_index <- as.integer(as.factor(cf$partner))
cf$male_index <- as.integer(as.factor(cf$dadid))

#cf$s_score_doc30=(cf$score_doc30 - mean(cf$score_doc30))/sd(cf$score_doc30)
cf$s_score_doc=(cf$score_doc - mean(cf$score_doc))/sd(cf$score_doc)
cf$s_rDSI=(cf$cf_rDSI - mean(cf$cf_rDSI))/sd(cf$cf_rDSI)
cf$t_nmales=1/cf$nmales
#cf$s_age_doc=(cf$age_doc - mean(cf$age_doc))/sd(cf$age_doc)
#################################
m_effort2 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- t_nmales + a + b_dsi*s_rDSI + b_rank*s_score_doc + b_phg*PHG + a_id[prt1_index] + a_id[prt2_index] + a_dyad[dyad_index],

	c(a,b_rank,b_phg,b_dsi) ~ dnorm(0,2),
	a_id[prt1_index] ~ dnorm(0,sigma_id),
	a_dyad[dyad_index] ~ dnorm(0, sigma_dyad),
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2)
	
) , data=cf ,warmup=3000,iter=10000,chains=2,cores=2, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE)

output2=precis(m_effort2 , depth=2 , digits=2)@output #N_effs too low
##################################

m_effort1 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- t_nmales + a + b_dsi*s_rDSI + b_rank*s_score_doc + b_phg*PHG + a_id[prt1_index] + a_id[prt2_index] + a_dyad[dyad_index],

	c(a,b_rank,b_phg,b_dsi) ~ dnorm(0,1),
	a_id[prt1_index] ~ dnorm(0,sigma_id),
	a_dyad[dyad_index] ~ dnorm(0, sigma_dyad),
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2)
	
) , data=cf ,warmup=3000,iter=10000,chains=2,cores=2, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE)

###this  is my modified code of m_effort1_BJB####
m_effort1_BJB <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_nmales*t_nmales + b_dsi*s_rDSI + b_rank*s_score_doc + b_phg*PHG + a_male[male_index],
	c(a,b_rank,b_phg,b_dsi,b_nmales) ~ dnorm(0,1),
	a_male[prt1_index] ~ dnorm(0,sigma_id),
	sigma_id ~ dcauchy(0,2)
	
) , data=cf ,warmup=1000,iter=2000,chains=2,cores=2, control=list(adapt_delta=0.95), WAIC=TRUE)


output1=precis(m_effort1 , depth=2 , digits=2)@output
plot(m_effort1)
par(mfrow=c(1,1))

#Plot DSI and 
par(mfrow=c(1,2))
par(mar=c(3.5,2,1,1))
post=extract.samples(m_effort1)

norm=rnorm(10000,0,1)
dens(post$b_dsi,col="red",xlim=c(-5,5),xlab="posterior vs. prior")
dens(post$b_rank,add=T,col="blue")
dens(post$b_phg,add=T,col="orange")
dens(norm,add=T)
legend("topright",legend = c("dsi", "rank","phg","dnorm(0,1)"),col=c("red","blue","orange","black"),lty=1,lwd=2)

a_prt1_z <- matrix(0,1000,length(unique(cf$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(cf$dyad_index)))


dsi.seq <- seq(min(cf$s_rDSI),max(cf$s_rDSI),length.out=100)

d.pred <- data.frame(	
	prt1_index=rep(1,length(dsi.seq)),
	prt2_index=rep(1,length(dsi.seq)),
	dyad_index=rep(1,length(dsi.seq)),
	s_score_doc=rep(mean(cf$s_score_doc),length(dsi.seq)),
	s_rDSI=dsi.seq,
	PHG=rep(mean(cf$PHG),length(dsi.seq)),
	t_nmales=rep(mean(cf$t_nmales),length(dsi.seq))	
)

pred <- link( m_effort1 , data=d.pred, n=1000 ,replace=
	list(a_id=a_prt1_z , a_dyad=a_dyad_z))

pred.mean <- apply(pred , 2 , mean )
pred.HPDI <- apply( pred , 2 , HPDI )
par(mar=c(3.5,2,1,1))
plot( nextdad ~ s_rDSI , data=cf , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1.2,xlim=c(-0.7,4.4))

lines( dsi.seq , pred.mean ,lwd=2)
shade(pred.HPDI,dsi.seq,col=alpha("black",0.2))


xt=seq(0,9,by=1)
xxt=(xt-mean(cf$cf_rDSI))/sd(cf$cf_rDSI)
mtext(side=1,line=1.7,text="Dyadic sociality index - Cycling (flat)",cex=1)
mtext(side=2,line=1,text="Probability of siring next infant",cex=1)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1, cex.axis=1,at= xxt,labels=xt,line=.2,col=NA)
axis(2, tck=-0.015,cex=1.2,labels=NA,at=c(0,1))
axis(2, cex.axis=1,at= c(0,1),labels=c(0,1),line=-.1,col=NA)
##############
rank.seq <- seq(min(cf$s_score_doc),max(cf$s_score_doc),length=100)

d.pred <- data.frame(
	
	prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	s_score_doc=rank.seq,
	s_rDSI=rep(mean(cf$s_rDSI),length(rank.seq)),
	PHG=rep(mean(cf$PHG),length(rank.seq)),
	t_nmales=rep(mean(cf$t_nmales),length(rank.seq))	
)

pred <- link( m_effort1 , data=d.pred, n=1000 ,replace=
	list(a_id=a_prt1_z , a_dyad=a_dyad_z))

pred.mean <- apply(pred , 2 , mean )
pred.HPDI <- apply( pred , 2 , HPDI )
par(mar=c(3.5,3,1,1))
plot( nextdad ~ s_score_doc , data=cf , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1.2#,xlim=c(-0.6,4.8))
)
lines( rank.seq , pred.mean ,lwd=2)
shade(pred.HPDI,rank.seq,col=alpha("black",0.2))

#############################################################################
m_effort3 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- t_nmales + a + b_dsi*s_rDSI + b_rank*s_score_doc + b_phg*PHG + a_id[prt1_index] + a_id[prt2_index] + a_dyad[dyad_index],

	c(a,b_rank,b_phg,b_dsi) ~ dnorm(0,0.8),
	a_id[prt1_index] ~ dnorm(0,sigma_id),
	a_dyad[dyad_index] ~ dnorm(0, sigma_dyad),
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2)
	
) , data=cf ,warmup=3000,iter=10000,chains=2,cores=2, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE)

output3=precis(m_effort3 , depth=2 , digits=2)@output
plot(m_effort1)
par(mfrow=c(1,1))

post=extract.samples(m_effort3)

norm=rnorm(10000,0,0.8)
dens(post$b_dsi,col="red",xlim=c(-5,5))
dens(post$b_rank,add=T,col="blue")
dens(post$b_phg,add=T,col="orange")
dens(norm,add=T)
legend("topright",legend = c("dsi", "rank","phg","dnorm(0,0.8)"),col=c("red","blue","orange","black"),lty=1,lwd=2)

#####################################################################



