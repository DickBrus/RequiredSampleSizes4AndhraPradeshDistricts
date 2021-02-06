library(ggplot2)

dat <- read.csv("c:/Users/brus003/OneDrive - Wageningen University & Research/SISIndia/AP_SHC_c1_processed_modelling_v201911.csv")
mandals <- sort(unique(dat$block))
Zn_crit <- 0.9
nlegacy <- mu <- sigma2 <- f <- numeric(length=length(mandals))
for (i in 1:length(mandals)) {
  ids <- which(dat$block==mandals[i])  
  subdat <- dat[ids,]                  
  nlegacy[i] <- sum(!is.na(subdat$Zn)) 
  mu[i] <- mean(log(subdat$Zn),na.rm=TRUE)
  sigma2[i] <- var(log(subdat$Zn),na.rm=TRUE)
  f[i] <- sum(subdat$Zn < Zn_crit, na.rm=TRUE)/nlegacy[i]  # the district's fraction
}

df <- data.frame(district=mandals,n=nlegacy,mu=mu,var=sigma2,f=f)
write.csv(df,file="SummaryStatistics_mandals.csv",row.names=FALSE)


library(SampleSizeMeans)

wmax <- 0.2 #maximum width of confidence interval for mean ln(Zn) of district
conflevel <- 0.95
worstlevel <- 0.80
lambda <- 1/sigma2 #prior estimate of precision parameter
cv <- 0.25 #coefficient of variation of gamma distribution for lambda

nreq.freq <- nreq.alc.bayes <- nreq.alc.mbl <- nreq.acc.bayes <- nreq.acc.mbl <- nreq.woc.bayes <- nreq.woc.mbl <-  numeric(length=length(mandals))
cv <- 0.25
lambda <- 1/sigma2
for (i in 1:length(mandals)) {
  a <- 1/cv^2
  b <- a/lambda[i]
  nreq.freq[i] <- mu.freq(len=wmax, lambda= lambda[i], level=conflevel)
  nreq.alc.bayes[i] <- mu.alc(len=wmax, alpha=a, beta=b, n0=0, level=conflevel)
  nreq.alc.mbl[i] <- mu.mblalc(len=wmax, alpha=a, beta=b, level=conflevel)
  nreq.acc.bayes[i] <- mu.acc(len=wmax, alpha=a, beta=b, n0=0, level=conflevel)
  nreq.acc.mbl[i] <- mu.mblacc(len=wmax, alpha=a, beta=b, level=conflevel)
  nreq.woc.bayes[i] <- mu.modwoc(len=wmax, alpha=a, beta=b, n0=0, level=conflevel, worst.level=worstlevel)
  nreq.woc.mbl[i] <- mu.mblmodwoc(len=wmax, alpha=a, beta=b, level=conflevel, worst.level=worstlevel)
}
df <- data.frame(mandal=mandals, lambda = round(lambda,2),
                 freq = nreq.freq,
                 alc = nreq.alc.bayes, 
                 alc.mbl = nreq.alc.mbl, 
                 acc = nreq.acc.bayes,
                 acc.mbl = nreq.acc.mbl,
                 woc = nreq.woc.bayes,
                 woc.mbl = nreq.woc.mbl)
write.csv(df,file="RequiredSampleSizes_Mean_mandals.csv",row.names=FALSE)

nreq_means <- read.csv(file="RequiredSampleSizes_Mean_mandals.csv",header=TRUE)
dn <- (nreq_means[,c(4,5,6,7,8,9)] - nlegacy)*-1

pdf(file="SampleSurplusMandals_Mean.pdf",width=6,height=4)
hist(dn$alc.mbl,main="Mean of ln(Zn)",xlab="Sample surplus",breaks=seq(from=-750,to=3000,by=250))
dev.off()

# Required sample sizes for estimating areal fraction with Zn deficiency

library(binomSamSize)
library(SampleSizeBinomial) 

wmax <- 0.1
conflevel <- 0.95
alpha <- 1-conflevel
worstlevel <- 0.80
n0 <- nlegacy/100

nreq.wald <- nreq.alc.bayes <- nreq.alc.mbl <- nreq.acc.bayes <- nreq.acc.mbl <- nreq.woc.bayes <- nreq.woc.mbl <-  numeric(length=length(mandals))
for (i in 1:length(mandals)) {
  #see Eq. 28 in Sambucini
  a <- n0[i]*f[i]+1
  b <- n0[i]*(1-f[i])+1
#  nreq.wald[i] <- binomSamSize::ciss.wald(p0=f[i],d=wmax/2,alpha=1-conflevel)
  nreq.wald[i] <- ceiling((qnorm(1-alpha/2)*sqrt(f[i]*(1-f[i]))/(wmax/2))^2+1)
  out <- prop.alc(len=wmax,alpha=a,beta=b,level=conflevel,exact=TRUE)
  nreq.alc.bayes[i] <- out$n
  out <- prop.mblalc(len=wmax,alpha=a,beta=b,level=conflevel,exact=TRUE)
  nreq.alc.mbl[i] <- out$n
  out <- prop.acc(len=wmax,alpha=a,beta=b,level=conflevel,exact=TRUE)
  nreq.acc.bayes[i] <- out$n
  out <- prop.mblacc(len=wmax,alpha=a,beta=b,level=conflevel,exact=TRUE)
  nreq.acc.mbl[i] <- out$n
  out <- prop.modwoc(len=wmax,alpha=a,beta=b,level=conflevel,exact=TRUE,worst.level=worstlevel)
  nreq.woc.bayes[i] <- out$n
  out <- prop.mblmodwoc(len=wmax,alpha=a,beta=b,level=conflevel,exact=TRUE,worst.level=worstlevel)
  nreq.woc.mbl[i] <- out$n
}
df <- data.frame(mandal=mandals, f = round(f,3),
                 wald = nreq.wald,
                 alc = nreq.alc.bayes, 
                 alc.mbl = nreq.alc.mbl, 
                 acc = nreq.acc.bayes,
                 acc.mbl = nreq.acc.mbl,
                 woc = nreq.woc.bayes,
                 woc.mbl = nreq.woc.mbl)
write.csv(df,file="RequiredSampleSizes_ArealFraction_mandals.csv",row.names=FALSE)

nreq_fractions <- read.csv(file="RequiredSampleSizes_ArealFraction_mandals.csv",header=TRUE)
dn <- (nreq_fractions[,c(4,5,6,7,8,9)] - nlegacy)*-1

pdf(file="SampleSurplusMandals_ArealFraction.pdf",width=6,height=4)
hist(dn$alc.mbl,main="Areal fraction with Zn-deficiency",xlab="Sample surplus",,breaks=seq(from=-250,to=3000,by=250))
dev.off()


