#Data prepr. Models from line 96.

full_data=read.csv("Z:/Vroni/Olive Baboons/analyses/DSI/BehavTabs_conception_25.10.2017.csv", header=T, check.names=F, stringsAsFactors=F, sep=",")
str(full data)
nrow(full_data)

xdata=full_data[full_data$res_status==1,]
nrow(xdata)

unique(xdata$partner)
unique(xdata$mom)
unique(xdata$inf)

f=c("bf","bf","bf","un","ut","fd","dz","dz")
m=c("yn","yr","yy","ub","um","fh","xh","xc")
rdyad=paste(f,m)
dyad=paste(xdata$mom,xdata$partner)
red_data=xdata[!dyad%in%rdyad,]

xdata=red_data

correlation=cor(cbind(xdata$prxprop,xdata$a1rate,xdata$a2rate,xdata$vg2rate,xdata$g1rate,xdata$g2rate,xdata$c2m2rate))
          [,1]      [,2]      [,3]       [,4]        [,5]       [,6]        [,7]
[1,] 1.0000000 0.4202626 0.5445440 0.85215748 0.286247033 0.83667907 0.585170293
[2,] 0.4202626 1.0000000 0.4553370 0.31658408 0.318932871 0.34194150 0.115214552
[3,] 0.5445440 0.4553370 1.0000000 0.32394434 0.388500240 0.35250969 0.374645151
[4,] 0.8521575 0.3165841 0.3239443 1.00000000 0.051835968 0.96349925 0.444493277
[5,] 0.2862470 0.3189329 0.3885002 0.05183597 1.000000000 0.03784741 0.007933366
[6,] 0.8366791 0.3419415 0.3525097 0.96349925 0.037847407 1.00000000 0.303730132
[7,] 0.5851703 0.1152146 0.3746452 0.44449328 0.007933366 0.30373013 1.000000000


behav_means_PHG=apply(cbind(xdata$prxprop[xdata$group=="PHG"],xdata$a1rate[xdata$group=="PHG"],xdata$a2rate[xdata$group=="PHG"],xdata$vg2rate[xdata$group=="PHG"],xdata$g1rate[xdata$group=="PHG"],xdata$g2rate[xdata$group=="PHG"],xdata$c2m2rate[xdata$group=="PHG"]),2,mean) #get mean rate for each behavioural variable

behav_means_ENK=apply(cbind(xdata$prxprop[xdata$group=="ENK"],xdata$a1rate[xdata$group=="ENK"],xdata$a2rate[xdata$group=="ENK"],xdata$vg2rate[xdata$group=="ENK"],xdata$g1rate[xdata$group=="ENK"],xdata$g2rate[xdata$group=="ENK"],xdata$c2m2rate[xdata$group=="ENK"]),2,mean)

xdata$DSI[xdata$group=="PHG"]=(xdata$prxprop[xdata$group=="PHG"]/behav_means_PHG[1] + xdata$a1rate[xdata$group=="PHG"]/behav_means_PHG[2] + xdata$a2rate[xdata$group=="PHG"]/behav_means_PHG[3] + xdata$vg2rate[xdata$group=="PHG"]/behav_means_PHG[4] + xdata$g1rate[xdata$group=="PHG"]/behav_means_PHG[5] + xdata$g2rate[xdata$group=="PHG"]/behav_means_PHG[6]+ xdata$c2m2rate[xdata$group=="PHG"]/behav_means_PHG[7])/6 

xdata$DSI[xdata$group=="ENK"]=(xdata$prxprop[xdata$group=="ENK"]/behav_means_ENK[1] + xdata$a1rate[xdata$group=="ENK"]/behav_means_ENK[2] + xdata$a2rate[xdata$group=="ENK"]/behav_means_ENK[3] + xdata$vg2rate[xdata$group=="ENK"]/behav_means_ENK[4] + xdata$g1rate[xdata$group=="ENK"]/behav_means_ENK[5] + xdata$g2rate[xdata$group=="ENK"]/behav_means_ENK[6]+ xdata$c2m2rate[xdata$group=="ENK"]/behav_means_ENK[7])/6

xdata$rDSI[xdata$group=="PHG"]=(xdata$prxprop[xdata$group=="PHG"]/behav_means_PHG[1] + xdata$a2rate[xdata$group=="PHG"]/behav_means_PHG[3] + xdata$vg2rate[xdata$group=="PHG"]/behav_means_PHG[4] + xdata$g2rate[xdata$group=="PHG"]/behav_means_PHG[6])/4 

xdata$rDSI[xdata$group=="ENK"]=(xdata$prxprop[xdata$group=="ENK"]/behav_means_ENK[1] + xdata$a2rate[xdata$group=="ENK"]/behav_means_ENK[3] + xdata$vg2rate[xdata$group=="ENK"]/behav_means_ENK[4] + xdata$g2rate[xdata$group=="ENK"]/behav_means_ENK[6])/4


cor(cbind(xdata$DSI,xdata$rDSI))

hist(xdata$DSI)
hist(xdata$rDSI)

#make random effect for dyad
xdata$dyad=paste(xdata$mom,xdata$partner, sep = "_")



##############################################
#This is for relating consort DSI to lactation/preg DSI but I think both measures are too inaccurate in combination while paternity is fine
lac=read.csv("Z:/Vroni/Olive Baboons/analyses/DSI/input_parenting effort.csv", header=T, check.names=F, stringsAsFactors=F, sep=",")

preg=read.csv("Z:/Vroni/Olive Baboons/analyses/DSI/input_parenting effort_pregnancy.csv", header=T, check.names=F, stringsAsFactors=F, sep=",")
#too few dyads where we have consortship and full lac info, take preg instead?

#where=match(paste(xdata$mom,xdata$partner,xdata$inf),paste(lac$momid,lac$partner,lac$inf))

where=match(paste(xdata$mom,xdata$partner,xdata$inf),paste(preg$mom,preg$partner,preg$infid))
sum(is.na(where))
#xdata$lac_DSI=lac$DSI[where]

xdata$preg_DSI=preg$DSI[where]
xdata=xdata[!is.na(xdata$preg_DSI),]

setwd("Z:/Vroni/Olive Baboons/analyses/DSI/")  
write.table(xdata, file = "input_consort effort DSI.csv",sep = ",",row.names=F)

##############################################
#Adding in conception elos
scores <- read.csv("Z:/Vroni/Olive Baboons/analyses/Elo/Elo scores_males_analysis across groups_Jan13Dec16_standardised.csv", header=T, stringsAsFactors=F, sep=",")
behav =read.csv("Z:/Vroni/Olive Baboons/analyses/DSI/input_consort effort DSI_red.csv", header=T, check.names=F, stringsAsFactors=F, sep=",")

scores$Date=as.Date(scores$Date,format="%Y-%m-%d")

behav$doc_infant=as.Date(behav$doc_infant,format="%m/%d/%Y")

#add average score during lactation
behav$score_doc=rep(NA,nrow(behav))

for(i in 1:nrow(behav)){
behav$score_doc[i]=scores$s_scores[(scores$ID==behav$partner[i]|scores$ID==behav$mom[i]) & scores$Date==behav$doc_infant[i]]
behav$score_doc30[i]=scores$s_scores_smooth30[(scores$ID==behav$partner[i]|scores$ID==behav$mom[i]) & scores$Date==behav$doc_infant[i]]
}

setwd("Z:/Vroni/Olive Baboons/analyses/DSI/")  
write.table(behav, file = "input_consort effort DSI_red.csv",sep = ",",row.names=F)

####################################################
library(rethinking)
library(Cairo)											
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("Z:/Vroni/Olive Baboons/analyses/DSI/")  
d=read.csv("Z:/Vroni/Olive Baboons/analyses/DSI/input_consort effort DSI_red.csv", header=T, check.names=F, stringsAsFactors=F, sep=",")

d$PHG <- ifelse(d$group=="PHG" , 1 , 0 )


d2=d[d$hrsobs>0.25,]
d3=d[d$rDSI<60,]

d$dyad_index <- as.integer(as.factor(d$dyad))
d$prt1_index <- as.integer(as.factor(d$mom))
d$prt2_index <- as.integer(as.factor(d$partner))

d2$dyad_index <- as.integer(as.factor(d2$dyad))
d2$prt1_index <- as.integer(as.factor(d2$mom))
d2$prt2_index <- as.integer(as.factor(d2$partner))


d3$dyad_index <- as.integer(as.factor(d3$dyad))
d3$prt1_index <- as.integer(as.factor(d3$mom))
d3$prt2_index <- as.integer(as.factor(d3$partner))

d$s_rDSI=(d$rDSI - mean(d$rDSI))/sd(d$rDSI)
d2$s_rDSI=(d2$rDSI - mean(d2$rDSI))/sd(d2$rDSI)
d3$s_rDSI=(d3$rDSI - mean(d3$rDSI))/sd(d3$rDSI)
d$s_score_doc=(d$score_doc - mean(d$score_doc))/sd(d$score_doc)
d2$s_score_doc=(d2$score_doc - mean(d2$score_doc))/sd(d2$score_doc)
d3$s_score_doc=(d3$score_doc - mean(d3$score_doc))/sd(d3$score_doc)

table(d$PHG,d$realdad)

############################################################
consort <- map2stan( 
alist(
realdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_rank*s_score_doc + b_dsi*s_rDSI + b_group*PHG +
				a_id[prt1_index] + a_id[prt2_index] + a_dyad[dyad_index],

	c(a,b_rank,b_dsi,b_group) ~ dnorm(0,2),
	a_id[prt1_index] ~ dnorm(0,sigma_id),
	a_dyad[dyad_index] ~ dnorm(0,sigma_dyad ),
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2)
		
) , data=d2 ,warmup=4000,iter=10000,chains=2,cores=2, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE)


output=precis(consort , depth=2 , digits=2)@output #n_eff extremely low despite many iterations

################################################################################
consort2 <- map2stan(  #same as 1 but observations =0.25h removed
alist(
realdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_rank*s_score_doc + b_dsi*s_rDSI + b_group*PHG +
				a_id[prt1_index] + a_id[prt2_index] + a_dyad[dyad_index],

	c(a,b_rank,b_dsi,b_group) ~ dnorm(0,1),
	a_id[prt1_index] ~ dnorm(0,sigma_id),
	a_dyad[dyad_index] ~ dnorm(0,sigma_dyad ),
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2)
		
) , data=d2 ,warmup=4000,iter=10000,chains=2,cores=2, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE)

output2=precis(consort2 , depth=2 , digits=2)@output


par(mar=c(2.5,2.5,1,1))
post=extract.samples(consort2)

norm=rnorm(10000,0,1) #compare prior and posteriors, seesm fine
dens(post$b_dsi,col="red",xlim=c(-5,5),xlab="posterior vs. prior")
dens(post$b_rank,add=T,col="blue")
dens(post$b_group,add=T,col="orange")
dens(norm,add=T)
legend("topright",legend = c("dsi", "rank","phg","dnorm(0,1)"),col=c("red","blue","orange","black"),lty=1,lwd=2)

a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

dsi.seq=seq(min(d2$s_rDSI),max(d2$s_rDSI),length=30)

d.pred_dsi <- list(
	prt1_index=	rep(1,length(dsi.seq)),
	prt2_index=	rep(1,length(dsi.seq)),
	dyad_index=	rep(1,length(dsi.seq)),
	s_rDSI=dsi.seq,
	s_score_doc=rep(mean(d2$s_score_doc),length(dsi.seq)),
	PHG=rep(mean(d2$PHG),length(dsi.seq))
)

pred <- link(consort2, n=1000 , data=d.pred_dsi, replace= list(a_id=a_prt1_z , a_dyad=a_dyad_z ), WAIC=TRUE)

pred.median<- apply(pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )

par(mar=c(2.5,2.5,1,1))
plot( realdad ~ s_rDSI , data=d2 , col=alpha("blue",0.5),pch=19 ,yaxt='n',ylab=NA,xlab=NA,xaxt='n',cex=1.2)
lines( dsi.seq , pred.median ,lwd=2)
shade(pred.HPDI,dsi.seq,col=alpha("black",0.2))
mtext(side=1,line=1,text="sDSI - conception",cex=1)
mtext(side=2,line=1,text="Likelihood of paternity",cex=1)
####
a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

rank.seq=seq(min(d$s_score_doc),max(d$s_score_doc),length=30)

d.pred_rank <- list(
	prt1_index=	rep(1,length(rank.seq)),
	prt2_index=	rep(1,length(rank.seq)),
	dyad_index=	rep(1,length(rank.seq)),
	s_score_doc=rank.seq,
	s_rDSI=rep(mean(d$s_rDSI),length(rank.seq)),
	PHG=rep(mean(d$PHG),length(rank.seq))
)

pred <- link(consort2, n=1000 , data=d.pred_rank, replace= list(a_id=a_prt1_z , a_dyad=a_dyad_z ), WAIC=TRUE)

pred.median<- apply(pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )

par(mar=c(3,2.5,1,1))
plot( realdad ~ s_score_doc , data=d , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1.2)
lines( rank.seq , pred.median ,lwd=2)
shade(pred.HPDI,rank.seq,col=alpha("black",0.2))
mtext(side=1,line=0.6,at=c(-1.5,1),text=c("low","high"),cex=1)
mtext(side=1,line=1,text="rank",cex=1)
mtext(side=2,line=1,text="Probability of paternity",cex=1)
####


#############################################
consort3 <- map2stan(  #same as 2 but whole dataset with outlier removed
alist(
realdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_rank*s_score_doc + b_dsi*s_rDSI + b_group*PHG +
				a_id[prt1_index] + a_id[prt2_index] + a_dyad[dyad_index],

	c(a,b_rank,b_dsi,b_group) ~ dnorm(0,1),
	a_id[prt1_index] ~ dnorm(0,sigma_id),
	a_dyad[dyad_index] ~ dnorm(0,sigma_dyad ),
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2)
		
) , data=d3 ,warmup=4000,iter=10000,chains=2,cores=2, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE)

output3=precis(consort3 , depth=2 , digits=2)@output

plot(consort3)
par(mfrow = c(1, 1))

a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

dsi.seq=seq(min(d3$s_rDSI),max(d3$s_rDSI),length=30)

d.pred_dsi <- list(
	prt1_index=	rep(1,length(dsi.seq)),
	prt2_index=	rep(1,length(dsi.seq)),
	dyad_index=	rep(1,length(dsi.seq)),
	s_rDSI=dsi.seq,
	s_score_doc=rep(mean(d3$s_score_doc),length(dsi.seq)),
	PHG=rep(mean(d3$PHG),length(dsi.seq))
)

pred <- link(consort3, n=1000 , data=d.pred_dsi, replace= list(a_id=a_prt1_z , a_dyad=a_dyad_z ), WAIC=TRUE)

pred.median<- apply(pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )


par(mar=c(2.5,2.5,1,1))
plot( realdad ~ s_rDSI , data=d3 , col=alpha("blue",0.5),pch=19 ,yaxt='n',ylab=NA,xlab=NA,cex=1.2,xaxt='n')
lines( dsi.seq , pred.median ,lwd=2)
shade(pred.HPDI,dsi.seq,col=alpha("black",0.2))
mtext(side=1,line=1,text="sDSI - conception",cex=1)
mtext(side=2,line=1,text="Likelihood of paternity",cex=1)
########################################################################
consort4 <- map2stan(  #same as 2 and 3 but whole dataset
alist(
realdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_rank*s_score_doc + b_dsi*s_rDSI + b_group*PHG +
				a_id[prt1_index] + a_id[prt2_index] + a_dyad[dyad_index],

	c(a,b_rank,b_dsi,b_group) ~ dnorm(0,1),
	a_id[prt1_index] ~ dnorm(0,sigma_id),
	a_dyad[dyad_index] ~ dnorm(0,sigma_dyad ),
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2)
		
) , data=d ,warmup=4000,iter=10000,chains=2,cores=2, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE)

output4=precis(consort4 , depth=2 , digits=2)@output

a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

dsi.seq=seq(min(d$s_rDSI),max(d$s_rDSI),length=30)

d.pred_dsi <- list(
	prt1_index=	rep(1,length(dsi.seq)),
	prt2_index=	rep(1,length(dsi.seq)),
	dyad_index=	rep(1,length(dsi.seq)),
	s_rDSI=dsi.seq,
	s_score_doc=rep(mean(d$s_score_doc),length(dsi.seq)),
	PHG=rep(mean(d$PHG),length(dsi.seq))
)

pred <- link(consort4, n=1000 , data=d.pred_dsi, replace= list(a_id=a_prt1_z , a_dyad=a_dyad_z ), WAIC=TRUE)

pred.median<- apply(pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )


par(mar=c(2.5,2.5,1,1))
plot( realdad ~ s_rDSI , data=d , col=alpha("blue",0.5),pch=19 ,yaxt='n',ylab=NA,xlab=NA,cex=1.2,xaxt='n')
lines( dsi.seq , pred.median ,lwd=2)
shade(pred.HPDI,dsi.seq,col=alpha("black",0.2))
mtext(side=1,line=1,text="sDSI - conception",cex=1)
mtext(side=2,line=1,text="Likelihood of paternity",cex=1)

