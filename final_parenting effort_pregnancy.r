library(rethinking)
library(Cairo)												
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
########################################################################
#setwd("Z:/Vroni/Olive Baboons/analyses/DSI/")
#d=read.table(file="input_parenting effort_pregnancy_red.csv",header=T,sep=",")
d <- read.csv()
nrow(d[d$realdad==1 & d$nextdad==1,]) #how many excluded because current and next dad
d=d[!(d$realdad==1 & d$nextdad==1),]		#exclude current and next dads

d$dyad=droplevels(d$dyad)
d$partner=droplevels(d$partner)
d$momid=droplevels(d$mom)

#check if if each id in each column
setdiff(unique(d$partner),unique(d$mom)) ###this is not working-- indexes for partners are off
##for example individual s yy an yr have different indices when they are the partner and  mom (36,37) columns (41,40)-- this can't be the case

####BJB suggestions

#below will show you instances where entires are missing from other columns

#STEP 1- create 2 new coumns stripped from relationship to start to make manipulations on outside of R

d_old <- d #save orig dataset
d$prt1 <- d$mom
d$prt2 <- d$partner

#STEP 2
#with these unique columns see which ones are present in on columns but not the other

#unique(d$mom[!(d$mom %in% d$partner)]) ##running this will show that eb um xh ut are in mom but missing from partner column
#unique(d$partner[!(d$mom %in% d$partner)]) #running this will show that dz fx jj s1 tp yn yr yy are in partner but missing from mom column
unique(d$prt1[!(d$prt1 %in% d$prt2)]) 
unique(d$prt2[!(d$prt1 %in% d$prt2)]) 
unique(d$prt1[!(d$prt2 %in% d$prt1)]) 
unique(d$prt2[!(d$prt2 %in% d$prt1)]) 
##check to make sure the lengths of each columns are equal

#length(unique(d$mom)) #this should equal below #
#length(unique(d$partner)) #this should equal above #
length(unique(d$prt1)) #this should equal below #
length(unique(d$prt2)) #this should equal above #
sort(unique(d$prt1))
sort(unique(d$prt2))

##if the above returns integer(0) (all entires represented in both columns) or are not equal lengths preceed to STEP 3, otherwise continue below
##write current csv to harddrive
write.csv(d,"parentpreg.csv")

#OUTSIDE OF R, open that csv, make edits manually in prt1 and prt2 columns
##go back into R
d <- read.csv(file="/Users/BJB/Dropbox/Veronika Joan DSI male/parentpreg.csv" , header=TRUE) #reload outside of R edited file
#d <- d[,2:max(length(d))]##drop column imported into csv land if needed

#run above tests again till same # and names are duplicated in prt1 and prt2

##STEP 3
#assign new individual and dyad indices
d$prt1_index <- as.integer(as.factor(d$mom))
d$prt2_index <- as.integer(as.factor(d$partner))

##below is code to automatically creat dyad names and dyad indices-- i recommned this to avoid andy manual errors
library(dplyr)
d$dyad2 <- apply(d[,1:2], 1, function(s) paste0(sort(s), collapse='')) #sorts and created dyad names based off entries in coumns 1 and 2--currently is mom and partner but might need to change later
d$dyad_index <- as.integer(as.factor(d$dyad2))

d$PHG <- ifelse(d$group=="PHG" , 1 , 0 )
#d$dyad_index <- as.integer(as.factor(d$dyad)) #change to d$dyad2 if needed
d$prt1_index <- as.integer(as.factor(d$mom))
d$prt2_index <- as.integer(as.factor(d$partner))

#####END BJB Suggestions########

d$s_rank=(d$score_preg - mean(d$score_preg ))/sd(d$score_preg )
d$father=d$realdad
d$nodad=ifelse(d$father==1 | d$nextdad==1,0,1 )
###########################################################################################################################

p_preg1 <- map2stan(
alist(

DSI ~ dzagamma2( p, mu , scale ),
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
data=d, cores=2 , chains=2 , warmup=3000, iter=6000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)

plot(p_preg1)
par(mfrow = c(1, 1))

output_preg1=precis(p_preg1, depth=2 , digits=2)@output
plot(precis(p_preg1, pars=c("ap","am","bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group"),depth=2))

#Plot prep paternity
a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

d.pred_father1 <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=1,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nodad1 <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nextdad1 <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=1,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

link_father <- link(p_preg1, n=1000 , data=d.pred_father1,replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z), WAIC=TRUE)
pred_father <- (1-link_father$p)*link_father$mu
median(pred_father)
HPDI(pred_father)

link_nodad <- link(p_preg1, n=1000 , data=d.pred_nodad1, replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z), WAIC=TRUE)
pred_nodad <- (1-link_nodad$p)*link_nodad$mu
median(pred_nodad)
HPDI(pred_nodad)

link_nextdad <- link(p_preg1, n=1000 , data=d.pred_nextdad1, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z), WAIC=TRUE)
pred_nextdad <- (1-link_nextdad$p)*link_nextdad$mu
median(pred_nextdad)
HPDI(pred_nextdad)

#Plot
par(mfrow=c(1,1) , mar=c(0,0,0,0) , oma=c(0,0,0,0))
dens(pred_nodad, xlim=c(0,10.5) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,ylim=c(-0.12,5))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0,10,by=1))
axis(1, cex.axis=0.8,at= c(0,10),labels=c(0,10),line=0.1,col=NA)
ll <- d$DSI[d$father==1]
points(ll, rep(-.05,length(ll)), pch=15 , col=col.alpha("blue", alpha=0.4) , cex=0.75 )

shade( density(pred_father) , lim= as.vector(HPDI(pred_father, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad) , lim= as.vector(HPDI(pred_nextdad, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad) , lim= as.vector(HPDI(pred_nodad, prob=0.9999)) , col = col.alpha("orange", 0.5))

ll <- d$DSI[d$nextdad==1]
points(ll, rep(-.13,length(ll)), pch=17 , col=col.alpha("red", alpha=0.4) , cex=0.75 )
ll <- d$DSI[d$nodad==1]
points(ll, rep(-.20,length(ll)), pch=16 , col=col.alpha("orange", alpha=0.4) , cex=0.75 )

abline(v=median(pred_father) , lty=1,lwd=1)
abline(v=median(pred_nodad) , lty=2,lwd=1)
abline(v=median(pred_nextdad) , lty=4,lwd=1)

mtext(side=1,line=0.5,cex=0.8,text="Dyadic sociality index - Pregnancy")
legend(6,2, legend = c("", "",""),
       col=c(1,1)  , lty=c(1,4,2),
       lw=1 , cex=1, bty="n",y.intersp=0.9)

legend(6.25,2, legend = c("sires of current infant","sires of next infant","other males"),
       col=c(col.alpha("blue", 0.5) ,col.alpha("red", 0.5), col.alpha("orange", 0.5) ) , pch=c(15,17,16),
       pt.cex=1.5 , cex=0.8, bty="n",y.intersp=1.2,x.intersp=2)
	   
#Plot prep for rank
rank.seq=seq(min(d$s_rank),max(d$s_rank),length=30)

d.pred_rank1 <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(mean(d$father),length(rank.seq)),
	nextdad=rep(mean(d$nextdad),length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))	
)

link_rank <- link(p_preg1, n=1000 , data=d.pred_rank1, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z), WAIC=TRUE)
pred_rank <- (1-link_rank$p)*link_rank$mu
pred.median=apply(pred_rank , 2 , median )
pred.HPDI=apply(pred_rank , 2 , HPDI )

#Plot rank
par(mar=c(3,3,0.5,1.5))
#colour=ifelse(d$father==1,"blue",ifelse(d$nextdad==1,"red","orange"))
#plot( DSI ~ s_rank , data=d , col=alpha(colour,0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='')
plot( DSI ~ s_rank , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='')
lines( rank.seq , pred.median )
shade(pred.HPDI,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=1.5)
mtext("DSI", side=2, line=1.5)

#############################P_preg2################################################
p_preg2 <- map2stan(
alist(

DSI ~ dzagamma2( p, mu , scale ),
  logit(p) ~ 	ap + bp_rank*s_rank + bp_pat*father + bp_next*nextdad + bp_group*PHG +
				AP + BPr*s_rank + BPf*father + BPf*nextdad,
				
                AP ~ ap_id[prt1_index] + ap_id[prt2_index] + ap_dyad[dyad_index],
				BPf ~ bpf_id[prt1_index] + bpf_id[prt2_index],
				BPr ~ bpr_id[prt1_index] + bpr_id[prt2_index],
								
  log(mu)  ~	am + bm_rank*s_rank + bm_pat*father + bm_next*nextdad + bm_group*PHG +
				AM + BMr*s_rank + BMf*father + BMf*nextdad,
				
                AM ~ am_id[prt1_index] + am_id[prt2_index] + am_dyad[dyad_index],
				BMf ~ bmf_id[prt1_index] + bmf_id[prt2_index],
				BMr ~ bmr_id[prt1_index] + bmr_id[prt2_index],
  
    c(ap,am,bm_pat,bp_pat,bp_group,bm_group,bm_rank,bp_rank,bp_next,bm_next) ~ dnorm(0,2),
	c(ap_id,am_id,bpf_id,bmf_id,bpr_id,bmr_id)[prt1_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(ap_dyad,am_dyad)[dyad_index] ~ dmvnormNC(sigma_dyad , Rho_dyad ),
	c(Rho_id,Rho_dyad) ~ dlkjcorr(3),	
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),
data=d, cores=2 , chains=2 , warmup=3000, iter=6000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE)

plot(p_preg2)
par(mfrow = c(1, 1))

output_preg2=precis(p_preg2, depth=2 , digits=2)@output
plot(precis(p_preg2, pars=c("ap","am","bp_rank","bm_rank","bp_pat","bm_pat","bp_no","bm_no","bp_next","bm_next","bp_group","bm_group"),depth=2))

a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))


link_father <- link(p_preg2, n=1000 , data=d.pred_father1,replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z), WAIC=TRUE)
pred_father <- (1-link_father$p)*link_father$mu
median(pred_father)
HPDI(pred_father)

link_nodad <- link(p_preg2, n=1000 , data=d.pred_nodad1, replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z), WAIC=TRUE)
pred_nodad <- (1-link_nodad$p)*link_nodad$mu
median(pred_nodad)
HPDI(pred_nodad)

link_nextdad <- link(p_preg2, n=1000 , data=d.pred_nextdad1, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z), WAIC=TRUE)
pred_nextdad <- (1-link_nextdad$p)*link_nextdad$mu
median(pred_nextdad)
HPDI(pred_nextdad)

#Plot
par(mfrow=c(1,1) , mar=c(0,0,0,0) , oma=c(0,0,0,0))
dens(pred_nodad, xlim=c(0,10.5) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,ylim=c(-0.12,6))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0,10,by=1))
axis(1, cex.axis=0.8,at= c(0,10),labels=c(0,10),line=0.1,col=NA)
ll <- d$DSI[d$father==1]
points(ll, rep(-.05,length(ll)), pch=15 , col=col.alpha("blue", alpha=0.4) , cex=0.75 )

shade( density(pred_father) , lim= as.vector(HPDI(pred_father, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad) , lim= as.vector(HPDI(pred_nextdad, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad) , lim= as.vector(HPDI(pred_nodad, prob=0.9999)) , col = col.alpha("orange", 0.5))

ll <- d$DSI[d$nextdad==1]
points(ll, rep(-.13,length(ll)), pch=17 , col=col.alpha("red", alpha=0.4) , cex=0.75 )
ll <- d$DSI[d$nodad==1]
points(ll, rep(-.20,length(ll)), pch=16 , col=col.alpha("orange", alpha=0.4) , cex=0.75 )

abline(v=median(pred_father) , lty=1,lwd=1)
abline(v=median(pred_nodad) , lty=2,lwd=1)
abline(v=median(pred_nextdad) , lty=4,lwd=1)

mtext(side=1,line=0.5,cex=0.8,text="Dyadic sociality index - Pregnancy")
legend(6,2, legend = c("", "",""),
       col=c(1,1)  , lty=c(1,4,2),
       lw=1 , cex=1, bty="n",y.intersp=0.9)

legend(6.25,2, legend = c("sires of current infant","sires of next infant","other males"),
       col=c(col.alpha("blue", 0.5) ,col.alpha("red", 0.5), col.alpha("orange", 0.5) ) , pch=c(15,17,16),
       pt.cex=1.5 , cex=0.8, bty="n",y.intersp=1.2,x.intersp=2)

#Plot prep for rank

link_rank <- link(p_preg2, n=1000 , data=d.pred_rank1, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z), WAIC=TRUE)
pred_rank <- (1-link_rank$p)*link_rank$mu
pred.median=apply(pred_rank , 2 , median )
pred.HPDI=apply(pred_rank , 2 , HPDI )

#Plot rank
par(mar=c(3,3,0.5,1.5))
plot( DSI ~ s_rank , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='')
lines( rank.seq , pred.median )
shade(pred.HPDI,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=1.5)
mtext("DSI", side=2, line=1.5)

#########################################P_preg3#####################################################################
p_preg3 <- map2stan(
alist(

DSI ~ dzagamma2( p, mu , scale ),
  logit(p) ~ 	ap + bp_rank*s_rank + bp_pat*father + bp_next*nextdad + bp_group*PHG + 
				bp_rank_pat*s_rank*father + bp_rank_next*s_rank*nextdad  + AP,
				
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
data=d, cores=2 , chains=2 , warmup=3000, iter=6000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)

plot(p_preg3)
par(mfrow = c(1, 1))

output_preg3=precis(p_preg3, depth=2 , digits=2)@output
plot(precis(p_preg3, pars=c("ap","am","bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group","bp_rank_pat","bm_rank_pat","bp_rank_next","bm_rank_next"),depth=2))

#prep plot paternity
a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

link_father <- link(p_preg3, n=1000 , data=d.pred_father1,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred_father <- (1-link_father$p)*link_father$mu
median(pred_father)
HPDI(pred_father)

link_nodad <- link(p_preg3, n=1000 , data=d.pred_nodad1, replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred_nodad <- (1-link_nodad$p)*link_nodad$mu
median(pred_nodad)
HPDI(pred_nodad)

link_nextdad <- link(p_preg3 n=1000 , data=d.pred_nextdad1, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred_nextdad <- (1-link_nextdad$p)*link_nextdad$mu
median(pred_nextdad)
HPDI(pred_nextdad)

#Plot
par(mfrow=c(1,1) , mar=c(0,0,0,0) , oma=c(0,0,0,0))
dens(pred_nodad, xlim=c(0,10.5) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,ylim=c(-0.12,3.6))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0,10,by=1))
axis(1, cex.axis=0.8,at= c(0,10),labels=c(0,10),line=0.1,col=NA)
ll <- d$DSI[d$father==1]
points(ll, rep(-.05,length(ll)), pch=15 , col=col.alpha("blue", alpha=0.4) , cex=0.75 )

shade( density(pred_father) , lim= as.vector(HPDI(pred_father, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad) , lim= as.vector(HPDI(pred_nextdad, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad) , lim= as.vector(HPDI(pred_nodad, prob=0.9999)) , col = col.alpha("orange", 0.5))

ll <- d$DSI[d$nextdad==1]
points(ll, rep(-.13,length(ll)), pch=17 , col=col.alpha("red", alpha=0.4) , cex=0.75 )
ll <- d$DSI[d$nodad==1]
points(ll, rep(-.20,length(ll)), pch=16 , col=col.alpha("orange", alpha=0.4) , cex=0.75 )

abline(v=median(pred_father) , lty=1,lwd=1)
abline(v=median(pred_nodad) , lty=2,lwd=1)
abline(v=median(pred_nextdad) , lty=4,lwd=1)

mtext(side=1,line=0.5,cex=0.8,text="Dyadic sociality index - Pregnancy")
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

d.pred_father3 <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(1,length(rank.seq)),
	nextdad=rep(0,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)

d.pred_nodad3 <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(0,length(rank.seq)),
	nextdad=rep(0,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)

d.pred_nextdad3 <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(0,length(rank.seq)),
	nextdad=rep(1,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)

link2 <- link( p_preg3 , data=d.pred_father3 , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred1 <- (1-link2$p)*link2$mu
pred.median1 <- apply(pred1 , 2 , median )
pred.HPDI1 <- apply( pred1 , 2 , HPDI )

link2 <- link( p_preg3 , data=d.pred_nodad3 , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred2 <- (1-link2$p)*link2$mu
pred.median2 <- apply(pred2 , 2 , median )
pred.HPDI2 <- apply( pred2 , 2 , HPDI )

link2 <- link( p_preg3 , data=d.pred_nextdad3 , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z), WAIC=TRUE)
pred3 <- (1-link2$p)*link2$mu
pred.median3 <- apply(pred3 , 2 , median )
pred.HPDI3 <- apply( pred3 , 2 , HPDI )

#plot rank
par(mfrow = c(1, 3),oma = c( 2.5, 1, 0, 0 ))

par(mar=c(1,3,0.5,0.5))
plot( DSI[d$father==1] ~ s_rank[d$father==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median1 )
shade(pred.HPDI1,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)
text(-1,11,"sires current infant")

par(mar=c(1,3,0.5,0.5))
plot( DSI[d$nodad==1] ~ s_rank[d$nodad==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median2 )
shade(pred.HPDI2,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)
text(-1.1,11,"other males")

par(mar=c(1,3,0.5,0.5))
plot( DSI[d$nextdad==1] ~ s_rank[nextdad==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median3 )
shade(pred.HPDI3,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)
text(-1.1,11,"sires next infant")
par(mfrow = c(1, 1),oma=c(0,0,0,0))

###########################################P-preg4#################################################################
p_preg4 <- map2stan(
alist(

DSI ~ dzagamma2( p, mu , scale ),

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
data=d, cores=2 , chains=2, warmup=3000, iter=6000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99,max_treedepth = 15), WAIC=TRUE)

output_preg4=precis(p_preg4, depth=2 , digits=2)@output
plot(precis(p_preg4, pars=c("ap","am","bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group","bp_rank_pat","bm_rank_pat","bp_rank_next","bm_rank_next"),depth=2))

#model comparisons
comp4=compare(p_preg1,p_preg2,p_preg3,p_preg4)
coefs4=coeftab(p_preg1,p_preg2,p_preg3,p_preg4)
plot(coefs4,pars=c("ap","am","bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group","bp_rank_pat","bm_rank_pat","bp_rank_next","bm_rank_next"),cex=0.7) 


#prep plot paternity
link_father <- link(p_preg4, n=1000 , data=d.pred_father1,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred_father <- (1-link_father$p)*link_father$mu
median(pred_father)
HPDI(pred_father)

link_nodad <- link(p_preg4, n=1000 , data=d.pred_nodad1, replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred_nodad <- (1-link_nodad$p)*link_nodad$mu
median(pred_nodad)
HPDI(pred_nodad)

link_nextdad <- link(p_preg4, n=1000 , data=d.pred_nextdad1, WAIC=TRUE,replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred_nextdad <- (1-link_nextdad$p)*link_nextdad$mu
median(pred_nextdad)
HPDI(pred_nextdad)

#Plot
par(mfrow=c(1,1) , mar=c(0,0,0,0) , oma=c(0,0,0,0))
dens(pred_nodad, xlim=c(0,10.5) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,ylim=c(-0.12,5))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0,10,by=1))
axis(1, cex.axis=0.8,at= c(0,10),labels=c(0,10),line=0.1,col=NA)
ll <- d$DSI[d$father==1]
points(ll, rep(-.05,length(ll)), pch=15 , col=col.alpha("blue", alpha=0.4) , cex=0.75 )

shade( density(pred_father) , lim= as.vector(HPDI(pred_father, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad) , lim= as.vector(HPDI(pred_nextdad, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad) , lim= as.vector(HPDI(pred_nodad, prob=0.9999)) , col = col.alpha("orange", 0.5))

ll <- d$DSI[d$nextdad==1]
points(ll, rep(-.13,length(ll)), pch=17 , col=col.alpha("red", alpha=0.4) , cex=0.75 )
ll <- d$DSI[d$nodad==1]
points(ll, rep(-.20,length(ll)), pch=16 , col=col.alpha("orange", alpha=0.4) , cex=0.75 )

abline(v=median(pred_father) , lty=1,lwd=1)
abline(v=median(pred_nodad) , lty=2,lwd=1)
abline(v=median(pred_nextdad) , lty=4,lwd=1)

mtext(side=1,line=0.5,cex=0.8,text="Dyadic sociality index - Pregnancy")
legend(6,2, legend = c("", "",""),
       col=c(1,1)  , lty=c(1,4,2),
       lw=1 , cex=1, bty="n",y.intersp=0.9)

legend(6.25,2, legend = c("sires of current infant","sires of next infant","other males"),
       col=c(col.alpha("blue", 0.5) ,col.alpha("red", 0.5), col.alpha("orange", 0.5) ) , pch=c(15,17,16),
       pt.cex=1.5 , cex=0.8, bty="n",y.intersp=1.2,x.intersp=2)
	   
#prep plot rank
link2 <- link( p_preg4 , data=d.pred_father3 , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred1 <- (1-link2$p)*link2$mu
pred.median1 <- apply(pred1 , 2 , median )
pred.HPDI1 <- apply( pred1 , 2 , HPDI )

link2 <- link( p_preg4 , data=d.pred_nodad3 , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred2 <- (1-link2$p)*link2$mu
pred.median2 <- apply(pred2 , 2 , median )
pred.HPDI2 <- apply( pred2 , 2 , HPDI )

link2 <- link( p_preg4 , data=d.pred_nextdad3 , n=1000 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
pred3 <- (1-link2$p)*link2$mu
pred.median3 <- apply(pred3 , 2 , median )
pred.HPDI3 <- apply( pred3 , 2 , HPDI )

#plot rank
par(mfrow = c(1, 3),oma = c( 2.5, 0.5, 0, 0 ))

par(mar=c(1,3,0.5,0.5))
plot( DSI[d$father==1] ~ s_rank[d$father==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median1 )
shade(pred.HPDI1,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)
text(-1,11,"sires current infant")

par(mar=c(1,3,0.5,0.5))
plot( DSI[d$nodad==1] ~ s_rank[d$nodad==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median2 )
shade(pred.HPDI2,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)
text(-1.1,11,"other males")

par(mar=c(1,3,0.5,0.5))
plot( DSI[d$nextdad==1] ~ s_rank[nextdad==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median3 )
shade(pred.HPDI3,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)
text(-1.1,11,"sires next infant")
par(mfrow = c(1, 1),oma=c(0,0,0,0))

############################################################################
