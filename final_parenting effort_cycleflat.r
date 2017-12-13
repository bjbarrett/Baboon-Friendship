library(rethinking)
library(Cairo)												
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
########################################################################
setwd("Z:/Vroni/Olive Baboons/analyses/rDSI/")
d=read.table(file="input_parenting effort_cyclflat_red.csv",header=T,sep=",")

nrow(d[d$realdad==1 & d$nextdad==1,]) #how many excluded because current and next dad
d=d[!(d$realdad==1 & d$nextdad==1),]		#exclude current and next dads

d$dyad=droplevels(d$dyad)
d$partner=droplevels(d$partner)
d$momid=droplevels(d$momid)

#check if if each id in each column
setdiff(unique(d$partner),unique(d$momid))

d$PHG <- ifelse(d$group=="PHG" , 1 , 0 )
d$dyad_index <- as.integer(as.factor(d$dyad))
d$prt1_index <- as.integer(as.factor(d$momid))
d$prt2_index <- as.integer(as.factor(d$partner))

d$s_rank=(d$score_cyclf - mean(d$score_cyclf ))/sd(d$score_cyclf )

d$father=d$realdad
d$nodad=ifelse(d$father==1 | d$nextdad==1,0,1 )

#DSI is here reduced to three correlated parameters hence rDSI
###############
p_cycf1 <- map2stan(
alist(

rDSI ~ dzagamma2( p, mu , scale ),
  logit(p) ~ 	ap + bp_pat*father + bp_next*nextdad  + bp_group*PHG + bp_rank*s_rank +
                ap_id[prt1_index] + ap_id[prt2_index] + ap_dyad[dyad_index],
				
  log(mu)  ~	am + bm_pat*father + bm_next*nextdad  + bm_group*PHG + bm_rank*s_rank +
				am_id[prt1_index] + am_id[prt2_index] + am_dyad[dyad_index],
  
    c(ap,am,bm_pat,bp_pat,bp_group,bm_group,bp_next,bm_next,bp_rank,bm_rank) ~ dnorm(0,2),
	c(ap_id,am_id)[prt1_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(ap_dyad,am_dyad)[dyad_index] ~ dmvnormNC(sigma_dyad , Rho_dyad ),
	c(Rho_id,Rho_dyad) ~ dlkjcorr(3),	
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),
data=d, cores=2 , chains=2 , warmup=3500, iter=7000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)

plot(p_cycf1)
par(mfrow = c(1, 1))

output_cycf1=precis(p_cycf1, depth=2 , digits=2)@output
plot(precis(p_cycf1, pars=c("bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group"),depth=2))

#Plot prep paternity
a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

d.pred_father <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=1,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nodad <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nextdad <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=1,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

link_father <- link(p_cycf1, n=1000 , data=d.pred_father,replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z), WAIC=TRUE)
pred_father <- (1-link_father$p)*link_father$mu
median(pred_father)
HPDI(pred_father)

link_nodad <- link(p_cycf1, n=1000 , data=d.pred_nodad, replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z), WAIC=TRUE)
pred_nodad <- (1-link_nodad$p)*link_nodad$mu
median(pred_nodad)
HPDI(pred_nodad)

link_nextdad <- link(p_cycf1, n=1000 , data=d.pred_nextdad, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z), WAIC=TRUE)
pred_nextdad <- (1-link_nextdad$p)*link_nextdad$mu
median(pred_nextdad)
HPDI(pred_nextdad)

#Plot
par(mfrow=c(1,1) , mar=c(0,0,0,0) , oma=c(0,0,0,0))
dens(pred_nodad, xlim=c(0,10.5) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,ylim=c(-0.12,3))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0,10,by=1))
axis(1, cex.axis=0.8,at= c(0,10),labels=c(0,10),line=0.1,col=NA)
ll <- d$rDSI[d$father==1]
points(ll, rep(-.05,length(ll)), pch=15 , col=col.alpha("blue", alpha=0.4) , cex=0.75 )

shade( density(pred_father) , lim= as.vector(HPDI(pred_father, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad) , lim= as.vector(HPDI(pred_nextdad, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad) , lim= as.vector(HPDI(pred_nodad, prob=0.9999)) , col = col.alpha("orange", 0.5))

ll <- d$rDSI[d$nextdad==1]
points(ll, rep(-.13,length(ll)), pch=17 , col=col.alpha("red", alpha=0.4) , cex=0.75 )
ll <- d$rDSI[d$nodad==1]
points(ll, rep(-.20,length(ll)), pch=16 , col=col.alpha("orange", alpha=0.4) , cex=0.75 )

abline(v=median(pred_father) , lty=1,lwd=1)
abline(v=median(pred_nodad) , lty=2,lwd=1)
abline(v=median(pred_nextdad) , lty=4,lwd=1)

mtext(side=1,line=0.5,cex=0.8,text="Dyadic sociality index - Cycling flat")
legend(6,2, legend = c("", "",""),
       col=c(1,1)  , lty=c(1,4,2),
       lw=1 , cex=1, bty="n",y.intersp=0.9)

legend(6.25,2, legend = c("sires of current infant","sires of next infant","other males"),
       col=c(col.alpha("blue", 0.5) ,col.alpha("red", 0.5), col.alpha("orange", 0.5) ) , pch=c(15,17,16),
       pt.cex=1.5 , cex=0.8, bty="n",y.intersp=1.2,x.intersp=2)
	   
#Plot prep for rank
rank.seq=seq(min(d$s_rank),max(d$s_rank),length=30)

d.pred_rank <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(mean(d$father),length(rank.seq)),
	nextdad=rep(mean(d$nextdad),length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))	
)

link_rank <- link(p_cycf1, n=1000 , data=d.pred_rank, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z), WAIC=TRUE)
pred_rank <- (1-link_rank$p)*link_rank$mu
pred.median=apply(pred_rank , 2 , median )
pred.HPDI=apply(pred_rank , 2 , HPDI )

#Plot rank
par(mar=c(3,3,0.5,1.5))
#colour=ifelse(d$father==1,"blue",ifelse(d$nextdad==1,"red","orange"))
#plot( rDSI ~ s_rank , data=d , col=alpha(colour,0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='')
plot( rDSI ~ s_rank , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='')
lines( rank.seq , pred.median )
shade(pred.HPDI,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=1.5)
mtext("rDSI", side=2, line=1.5)

#############################P_LAC2################################################
p_cycf2 <- map2stan(
alist(

rDSI ~ dzagamma2( p, mu , scale ),
  logit(p) ~ 	ap  + bp_pat*father + bp_next*nextdad + bp_group*PHG + bp_rank*s_rank+
				AP + BPf*father + BPf*nextdad + BPr*s_rank,
				
                AP ~ ap_id[prt1_index] + ap_id[prt2_index] + ap_dyad[dyad_index],
				BPf ~ bpf_id[prt1_index] + bpf_id[prt2_index],
				BPr ~ bpr_id[prt1_index] + bpr_id[prt2_index],
								
  log(mu)  ~	am  + bm_pat*father + bm_next*nextdad + bm_group*PHG + bm_rank*s_rank +
				AM + BMf*father + BMf*nextdad + BMr*s_rank,
				
                AM ~ am_id[prt1_index] + am_id[prt2_index] + am_dyad[dyad_index],
				BMf ~ bmf_id[prt1_index] + bmf_id[prt2_index],
				BMr ~ bmr_id[prt1_index] + bmr_id[prt2_index],
  
    c(ap,am,bm_pat,bp_pat,bp_group,bm_group,bp_next,bm_next,bp_rank,bm_rank) ~ dnorm(0,2),
	c(ap_id,am_id,bpf_id,bmf_id,bpr_id,bmr_id)[prt1_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(ap_dyad,am_dyad)[dyad_index] ~ dmvnormNC(sigma_dyad , Rho_dyad ),
	c(Rho_id,Rho_dyad) ~ dlkjcorr(3),	
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),
data=d, cores=2 , chains=2 , warmup=3500, iter=7000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)

plot(p_cycf2)
par(mfrow = c(1, 1))

output_cycf2=precis(p_cycf2, depth=2 , digits=2)@output
plot(precis(p_cycf2, pars=c("bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group"),depth=2))

a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

d.pred_father <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=1,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nodad <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nextdad <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=1,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

link_father <- link(p_cycf2, n=1000 , data=d.pred_father,replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z), WAIC=TRUE)
pred_father <- (1-link_father$p)*link_father$mu
median(pred_father)
HPDI(pred_father)

link_nodad <- link(p_cycf2, n=1000 , data=d.pred_nodad, replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z), WAIC=TRUE)
pred_nodad <- (1-link_nodad$p)*link_nodad$mu
median(pred_nodad)
HPDI(pred_nodad)

link_nextdad <- link(p_cycf2, n=1000 , data=d.pred_nextdad, replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z), WAIC=TRUE)
pred_nextdad <- (1-link_nextdad$p)*link_nextdad$mu
median(pred_nextdad)
HPDI(pred_nextdad)


#Plot
par(mfrow=c(1,1) , mar=c(0,0,0,0) , oma=c(0,0,0,0))
dens(pred_nodad, xlim=c(0,10.5) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,ylim=c(-0.12,3.6))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0,10,by=1))
axis(1, cex.axis=0.8,at= c(0,10),labels=c(0,10),line=0.1,col=NA)
ll <- d$rDSI[d$father==1]
points(ll, rep(-.05,length(ll)), pch=15 , col=col.alpha("blue", alpha=0.4) , cex=0.75 )

shade( density(pred_father) , lim= as.vector(HPDI(pred_father, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad) , lim= as.vector(HPDI(pred_nextdad, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad) , lim= as.vector(HPDI(pred_nodad, prob=0.9999)) , col = col.alpha("orange", 0.5))

ll <- d$rDSI[d$nextdad==1]
points(ll, rep(-.13,length(ll)), pch=17 , col=col.alpha("red", alpha=0.4) , cex=0.75 )
ll <- d$rDSI[d$nodad==1]
points(ll, rep(-.20,length(ll)), pch=16 , col=col.alpha("orange", alpha=0.4) , cex=0.75 )

abline(v=median(pred_father) , lty=1,lwd=1)
abline(v=median(pred_nodad) , lty=2,lwd=1)
abline(v=median(pred_nextdad) , lty=4,lwd=1)

mtext(side=1,line=0.5,cex=0.8,text="Dyadic sociality index - Cycling flat")
legend(6,2, legend = c("", "",""),
       col=c(1,1)  , lty=c(1,4,2),
       lw=1 , cex=1, bty="n",y.intersp=0.9)

legend(6.25,2, legend = c("sires of current infant","sires of next infant","other males"),
       col=c(col.alpha("blue", 0.5) ,col.alpha("red", 0.5), col.alpha("orange", 0.5) ) , pch=c(15,17,16),
       pt.cex=1.5 , cex=0.8, bty="n",y.intersp=1.2,x.intersp=2)

#Plot prep for rank
rank.seq=seq(min(d$s_rank),max(d$s_rank),length=30)

d.pred_rank <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(mean(d$father),length(rank.seq)),
	nextdad=rep(mean(d$nextdad),length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))	
)

link_rank <- link(p_cycf2, n=1000 , data=d.pred_rank, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z), WAIC=TRUE)
pred_rank <- (1-link_rank$p)*link_rank$mu
pred.median=apply(pred_rank , 2 , median )
pred.HPDI=apply(pred_rank , 2 , HPDI )

#Plot rank
par(mar=c(3,3,0.5,1.5))
plot( rDSI ~ s_rank , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='')
lines( rank.seq , pred.median )
shade(pred.HPDI,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=1.5)
mtext("rDSI", side=2, line=1.5)

#########################################P_LAC3#####################################################################
p_cycf3 <- map2stan(
alist(

rDSI ~ dzagamma2( p, mu , scale ),
  logit(p) ~ 	ap + bp_rank*s_rank + bp_pat*father + bp_next*nextdad + bp_group*PHG + 
				bp_rank_pat*s_rank*father + bp_rank_next*s_rank*nextdad + AP,
				
                AP ~ ap_id[prt1_index] + ap_id[prt2_index] + ap_dyad[dyad_index],
				
								
  log(mu)  ~	am + bm_rank*s_rank + bm_pat*father + bm_next*nextdad + bm_group*PHG +
				bm_rank_pat*s_rank*father + bm_rank_next*s_rank*nextdad + AM,
				
                AM ~ am_id[prt1_index] + am_id[prt2_index] + am_dyad[dyad_index],
			
  
    c(ap,am,bm_pat,bp_pat,bp_group,bm_group,bm_rank,bp_rank,bp_next,bm_next,bp_rank_pat,bm_rank_pat,bp_rank_next, bm_rank_next) ~ dnorm(0,2),
	c(ap_id,am_id)[prt1_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(ap_dyad,am_dyad)[dyad_index] ~ dmvnormNC(sigma_dyad , Rho_dyad ),
	c(Rho_id,Rho_dyad) ~ dlkjcorr(3),	
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),
data=d, cores=2 , chains=2 , warmup=3500, iter=7000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)

plot(p_cycf3)
par(mfrow = c(1, 1))

output_cycf3=precis(p_cycf3, depth=2 , digits=2)@output
plot(precis(p_cycf3, pars=c("bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group","bp_rank_pat","bm_rank_pat","bp_rank_next","bm_rank_next"),depth=2))

#prep plot paternity
a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

d.pred_father <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=1,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nodad <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nextdad <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=1,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

link_father <- link(p_cycf3, n=1000 , data=d.pred_father,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred_father <- (1-link_father$p)*link_father$mu
median(pred_father)
HPDI(pred_father)

link_nodad <- link(p_cycf3, n=1000 , data=d.pred_nodad, replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred_nodad <- (1-link_nodad$p)*link_nodad$mu
median(pred_nodad)
HPDI(pred_nodad)

link_nextdad <- link(p_cycf3, n=1000 , data=d.pred_nextdad, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred_nextdad <- (1-link_nextdad$p)*link_nextdad$mu
median(pred_nextdad)
HPDI(pred_nextdad)

#Plot
par(mfrow=c(1,1) , mar=c(0,0,0,0) , oma=c(0,0,0,0))
dens(pred_nodad, xlim=c(0,10.5) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,ylim=c(-0.12,3.6))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0,10,by=1))
axis(1, cex.axis=0.8,at= c(0,10),labels=c(0,10),line=0.1,col=NA)
ll <- d$rDSI[d$father==1]
points(ll, rep(-.05,length(ll)), pch=15 , col=col.alpha("blue", alpha=0.4) , cex=0.75 )

shade( density(pred_father) , lim= as.vector(HPDI(pred_father, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad) , lim= as.vector(HPDI(pred_nextdad, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad) , lim= as.vector(HPDI(pred_nodad, prob=0.9999)) , col = col.alpha("orange", 0.5))

ll <- d$rDSI[d$nextdad==1]
points(ll, rep(-.13,length(ll)), pch=17 , col=col.alpha("red", alpha=0.4) , cex=0.75 )
ll <- d$rDSI[d$nodad==1]
points(ll, rep(-.20,length(ll)), pch=16 , col=col.alpha("orange", alpha=0.4) , cex=0.75 )

abline(v=median(pred_father) , lty=1,lwd=1)
abline(v=median(pred_nodad) , lty=2,lwd=1)
abline(v=median(pred_nextdad) , lty=4,lwd=1)

mtext(side=1,line=0.5,cex=0.8,text="Dyadic sociality index - Lactation")
legend(6,2, legend = c("", "",""),
       col=c(1,1)  , lty=c(1,4,2),
       lw=1 , cex=1, bty="n",y.intersp=0.9)

legend(6.25,2, legend = c("sires of current infant","sires of next infant","other males"),
       col=c(col.alpha("blue", 0.5) ,col.alpha("red", 0.5), col.alpha("orange", 0.5) ) , pch=c(15,17,16),
       pt.cex=1.5 , cex=0.8, bty="n",y.intersp=1.2,x.intersp=2)
	   
#prep plot rank
a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))
rank.seq=seq(min(d$s_rank),max(d$s_rank),length=30)

d.pred_father <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(1,length(rank.seq)),
	nextdad=rep(0,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)

d.pred_nodad <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(0,length(rank.seq)),
	nextdad=rep(0,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)

d.pred_nextdad <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(0,length(rank.seq)),
	nextdad=rep(1,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)

link2 <- link( p_cycf3 , data=d.pred_father , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred1 <- (1-link2$p)*link2$mu
pred.median1 <- apply(pred1 , 2 , median )
pred.HPDI1 <- apply( pred1 , 2 , HPDI )

link2 <- link( p_cycf3 , data=d.pred_nodad , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred2 <- (1-link2$p)*link2$mu
pred.median2 <- apply(pred2 , 2 , median )
pred.HPDI2 <- apply( pred2 , 2 , HPDI )

link2 <- link( p_cycf3 , data=d.pred_nextdad , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred3 <- (1-link2$p)*link2$mu
pred.median3 <- apply(pred3 , 2 , median )
pred.HPDI3 <- apply( pred3 , 2 , HPDI )

#plot rank
par(mfrow = c(1, 3),oma = c( 2.5, 1, 0, 0 ))

par(mar=c(1,3,0.5,0.5))
plot( rDSI[d$father==1] ~ s_rank[d$father==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median1 )
shade(pred.HPDI1,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("rDSI", side=2, line=1.9,cex=0.8)
text(-1,11,"sires current infant")

par(mar=c(1,3,0.5,0.5))
plot( rDSI[d$nodad==1] ~ s_rank[d$nodad==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median2 )
shade(pred.HPDI2,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("rDSI", side=2, line=1.9,cex=0.8)
text(-1.1,11,"other males")

par(mar=c(1,3,0.5,0.5))
plot( rDSI[d$nextdad==1] ~ s_rank[nextdad==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median3 )
shade(pred.HPDI3,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("rDSI", side=2, line=1.9,cex=0.8)
text(-1.1,11,"sires next infant")
par(mfrow = c(1, 1),oma=c(0,0,0,0))

###########################################P-LAC4#################################################################
p_cycf4 <- map2stan(
alist(

rDSI ~ dzagamma2( p, mu , scale ),

    logit(p) ~ 	ap + bp_rank_pat*s_rank*father + bp_rank_next*s_rank*nextdad  + bp_pat*father + bp_next*nextdad + bp_rank*s_rank + bp_group*PHG + 
				#AP + (BPf + BPfr*s_rank)*father + (BPf + BPfr*s_rank)*nextdad + BPr*s_rank,									
				AP + (BPf + BPfr*s_rank)*father + (BPn + BPnr*s_rank)*nextdad + BPr*s_rank,									

				AP ~ ap_id[prt1_index] + ap_id[prt2_index] + ap_dyad[dyad_index], 
				BPf ~ bpf_id[prt1_index] + bpf_id[prt2_index],
				BPr ~ bpr_id[prt1_index] + bpr_id[prt2_index],
				BPn ~ bpn_id[prt1_index] + bpn_id[prt2_index],
				BPfr ~ bpfr_id[prt1_index] + bpfr_id[prt2_index],
				BPnr ~ bpnr_id[prt1_index] + bpnr_id[prt2_index],

    log(mu) ~ 	am + bm_rank_pat*s_rank*father + bm_rank_next*s_rank*nextdad + bm_pat*father + bm_next*nextdad + bm_rank*s_rank + bm_group*PHG +   
				#AM + (BMf + BMfr*s_rank)*father + (BMf + BMfr*s_rank)*nextdad + BMr*s_rank,										 
				AM + (BMf + BMfr*s_rank)*father + (BMn + BMnr*s_rank)*nextdad + BMr*s_rank,										 

				AM ~ am_id[prt1_index] + am_id[prt2_index] + am_dyad[dyad_index], 
				BMf ~ bmf_id[prt1_index] + bmf_id[prt2_index]  ,
				BMr ~ bmr_id[prt1_index] + bmr_id[prt2_index] ,
				BMn ~ bmn_id[prt1_index] + bmn_id[prt2_index],
				BMfr ~ bmfr_id[prt1_index] + bmfr_id[prt2_index]  ,		
				BMnr ~ bmnr_id[prt1_index] + bmnr_id[prt2_index],		

  
    c(ap,am,bp_rank_pat,bm_rank_pat,bp_rank_next,bm_rank_next,bp_pat,bm_pat,bp_next,bm_next,bp_rank,bm_rank,bp_group,bm_group) ~ dnorm(0,2),
	c(ap_id,am_id,bpf_id,bmf_id,bpfr_id,bmfr_id,bpr_id,bmr_id,bpnr_id,bmnr_id,bpn_id,bmn_id)[prt1_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(ap_dyad,am_dyad)[dyad_index] ~ dmvnormNC( sigma_dyad , Rho_dyad ),
	c(Rho_id,Rho_dyad) ~ dlkjcorr(3),	
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
	),
data=d, cores=2 , chains=2, warmup=5000, iter=10000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99,max_treedepth = 15), WAIC=TRUE)

output_cycf4=precis(p_cycf4, depth=2 , digits=2)@output
plot(precis(p_cycf4, pars=c("bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group","bp_rank_pat","bm_rank_pat","bp_rank_next","bm_rank_next"),depth=2))

#model comparisons
comp4=compare(p_cycf1,p_cycf2,p_cycf3,p_cycf4)
coefs4=coeftab(p_cycf1,p_cycf2,p_cycf3,p_cycf4)
plot(coefs4,pars=c("bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group","bp_rank_pat","bm_rank_pat","bp_rank_next","bm_rank_next"),cex=0.7) 

#prep plot paternity
a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

d.pred_father <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=1,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nodad <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nextdad <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=1,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

link_father <- link(p_cycf4, n=1000 , data=d.pred_father,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred_father <- (1-link_father$p)*link_father$mu
median(pred_father)
HPDI(pred_father)

link_nodad <- link(p_cycf4, n=1000 , data=d.pred_nodad, replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred_nodad <- (1-link_nodad$p)*link_nodad$mu
median(pred_nodad)
HPDI(pred_nodad)

link_nextdad <- link(p_cycf4, n=1000 , data=d.pred_nextdad, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred_nextdad <- (1-link_nextdad$p)*link_nextdad$mu
median(pred_nextdad)
HPDI(pred_nextdad)

#Plot
par(mfrow=c(1,1) , mar=c(0,0,0,0) , oma=c(0,0,0,0))
dens(pred_nodad, xlim=c(0,10.5) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,ylim=c(-0.12,2))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0,10,by=1))
axis(1, cex.axis=0.8,at= c(0,10),labels=c(0,10),line=0.1,col=NA)
ll <- d$rDSI[d$father==1]
points(ll, rep(-.05,length(ll)), pch=15 , col=col.alpha("blue", alpha=0.4) , cex=0.75 )

shade( density(pred_father) , lim= as.vector(HPDI(pred_father, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad) , lim= as.vector(HPDI(pred_nextdad, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad) , lim= as.vector(HPDI(pred_nodad, prob=0.9999)) , col = col.alpha("orange", 0.5))

ll <- d$rDSI[d$nextdad==1]
points(ll, rep(-.13,length(ll)), pch=17 , col=col.alpha("red", alpha=0.4) , cex=0.75 )
ll <- d$rDSI[d$nodad==1]
points(ll, rep(-.20,length(ll)), pch=16 , col=col.alpha("orange", alpha=0.4) , cex=0.75 )

abline(v=median(pred_father) , lty=1,lwd=1)
abline(v=median(pred_nodad) , lty=2,lwd=1)
abline(v=median(pred_nextdad) , lty=4,lwd=1)

mtext(side=1,line=0.5,cex=0.8,text="Dyadic sociality index - Cycling (flat)")
legend(6,2, legend = c("", "",""),
       col=c(1,1)  , lty=c(1,4,2),
       lw=1 , cex=1, bty="n",y.intersp=0.9)

legend(6.25,2, legend = c("sires of current infant","sires of next infant","other males"),
       col=c(col.alpha("blue", 0.5) ,col.alpha("red", 0.5), col.alpha("orange", 0.5) ) , pch=c(15,17,16),
       pt.cex=1.5 , cex=0.8, bty="n",y.intersp=1.2,x.intersp=2)
	   
#prep plot rank
a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))
rank.seq=seq(min(d$s_rank),max(d$s_rank),length=30)

d.pred_father <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(1,length(rank.seq)),
		nextdad=rep(0,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)

d.pred_nodad <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(0,length(rank.seq)),
	nextdad=rep(0,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)

d.pred_nextdad <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(0,length(rank.seq)),
	nextdad=rep(1,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)

link2 <- link( p_cycf4 , data=d.pred_father , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred1 <- (1-link2$p)*link2$mu
pred.median1 <- apply(pred1 , 2 , median )
pred.HPDI1 <- apply( pred1 , 2 , HPDI )

link2 <- link( p_cycf4, data=d.pred_nodad , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred2 <- (1-link2$p)*link2$mu
pred.median2 <- apply(pred2 , 2 , median )
pred.HPDI2 <- apply( pred2 , 2 , HPDI )

link2 <- link( p_cycf4 , data=d.pred_nextdad , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred3 <- (1-link2$p)*link2$mu
pred.median3 <- apply(pred3 , 2 , median )
pred.HPDI3 <- apply( pred3 , 2 , HPDI )

#plot rank
par(mfrow = c(1, 3),oma = c( 2.5, 0.5, 0, 0 ))

par(mar=c(1,3,0.5,0.5))
plot( rDSI[d$father==1] ~ s_rank[d$father==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median1 )
shade(pred.HPDI1,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("rDSI", side=2, line=1.9,cex=0.8)
text(-1,11,"sires current infant")

par(mar=c(1,3,0.5,0.5))
plot( rDSI[d$nodad==1] ~ s_rank[d$nodad==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median2 )
shade(pred.HPDI2,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("rDSI", side=2, line=1.9,cex=0.8)
text(-1.1,11,"other males")

par(mar=c(1,3,0.5,0.5))
plot( rDSI[d$nextdad==1] ~ s_rank[nextdad==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median3 )
shade(pred.HPDI3,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("rDSI", side=2, line=1.9,cex=0.8)
text(-1.1,11,"sires next infant")
par(mfrow = c(1, 1),oma=c(0,0,0,0))
################################################
dyadlist <- (unique(d$dyad_index))
idlist <- (unique(d$prt1_index))

d.pred_father.s <- list(
    prt1_index=idlist,
	prt2_index=idlist,
	dyad_index=rep(1,length(idlist)),
	father=rep(1,length(idlist)),
	nextdad=rep(0,length(idlist)),
	s_rank=rep(mean(d$s_rank),length(idlist)),
	PHG=rep(mean(d$PHG),length(idlist))
)

d.pred_nodad.s <- list(
    prt1_index=idlist,
	prt2_index=idlist,
	dyad_index=rep(1,length(idlist)),
	father=rep(0,length(idlist)),
	nextdad=rep(0,length(idlist)),
	s_rank=rep(mean(d$s_rank),length(idlist)),
	PHG=rep(mean(d$PHG),length(idlist))
)

d.pred_nextdad.s <- list(
    prt1_index=idlist,
	prt2_index=idlist,
	dyad_index=rep(1,length(idlist)),
	father=rep(0,length(idlist)),
	nextdad=rep(1,length(idlist)),
	s_rank=rep(mean(d$s_rank),length(idlist)),
	PHG=rep(mean(d$PHG),length(idlist))
)

link_father <- link(p_cycf4, n=1000 , data=d.pred_father.s,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred_father <- (1-link_father$p)*link_father$mu

link_nodad <- link(p_cycf4, n=1000 , data=d.pred_nodad.s, replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred_nodad <- (1-link_nodad$p)*link_nodad$mu

link_nextdad <- link(p_cycf4, n=1000 , data=d.pred_nextdad.s, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred_nextdad <- (1-link_nextdad$p)*link_nextdad$mu

par(cex=0.75, mar=c(3,4,.5,.5))
plot(1,1 , ylim=c(0,4) , xlim=c(0,length(idlist)) , col="white" , cex.axis=1 , xaxt='n' , xlab='' , ylab='')
for (i in 1:length(idlist)){
	lines( y=c(HPDI(pred_father[,i])[1] ,HPDI(pred_father[,i])[2] ) , x=c(i,i) , col="black" )
	points( i , median(pred_father[,i]) , pch=15 , cex=0.5 )
	lines( y=c(HPDI(pred_nodad[,i])[1] ,HPDI(pred_nodad[,i])[2] ) , x=c(i+0.2,i+0.2) , col="red" )
	points( i+0.2 , median(pred_nodad[,i]) , pch=15 , cex=0.5 , col="red"  )
	lines( y=c(HPDI(pred_nextdad[,i])[1] ,HPDI(pred_nextdad[,i])[2] ) , x=c(i+0.4,i+0.4) , col="blue" )
	points( i+0.4 , median(pred_nextdad[,i]) , pch=15 , cex=0.5 , col="blue"  )
}
