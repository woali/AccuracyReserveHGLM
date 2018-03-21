##!!MUST HAVE INTERNET CONNECTION
##authors: Alicja Wolny-Dominiak, Tomasz Żądło

# install.packages("hglm")
# install.packages("statmod")
# install.packages("tweedie")
# install.packages("StatMatch") ##dummy
# install.packages("ggplot2")

library(hglm)
library(statmod)
library(tweedie)
library(StatMatch) ##dummy
library(ggplot2)
set.seed(456)

options(warn = -1)

lboot <- 2
probs <- c(0.5,0.75, 0.9, 0.95, 0.99)

triangle_upper <- read.csv2(file="http://web.ue.katowice.pl/woali/triangle_upper.csv")
triangle_lower <- read.csv2(file="http://web.ue.katowice.pl/woali/triangle_lower.csv")
X.upper <- model.matrix(~as.factor(dev),data=triangle_upper) 
Z.upper <- fact2dummy(as.factor(triangle_upper$origin))
X.lower <- as.matrix(cbind(rep(1,length(triangle_lower$dev)),fact2dummy(as.factor(triangle_lower$dev)))) 
Z.lower <- fact2dummy(as.factor(triangle_lower$origin))

#upper triangle
claim <- triangle_upper$claim

#Tweedie's power (ODP distribution)
p <- 1 

###model
start.upper <-function(y.upper, p){
beta1=hglm(fixed=claim~as.factor(dev),random=~1|as.factor(origin), family=tweedie(var.power=p,link.power=0), 
rand.family=Gamma(log), data=triangle_upper)

beta=beta1$fixef
u=beta1$ranef
v=log(beta1$ranef)
hphi <- beta1$varFix
hmu_ij <- beta1$fv
r_ij=(claim-hmu_ij)/sqrt(hmu_ij^p) 
phi.u <-beta1$varRanef
e <-beta1$varFix
out <-list(beta=beta, u=u, v=v, hphi=hphi, r_ij=r_ij, hmu_ij=hmu_ij, phi.u=phi.u, e=e)
}

beta=start.upper(claim, 1)$beta
u=start.upper(claim, 1)$u
v=start.upper(claim, 1)$v
hphi=start.upper(claim, 1)$hphi

#total reserve R
hy_ij_p <- exp(X.lower%*%beta+Z.lower%*%v[2:10])
Ri <- tapply(hy_ij_p, triangle_lower$origin, sum)
R <- sum(Ri)
R

mean.dev <- tapply(hy_ij_p, triangle_lower$dev, mean)

######### residual boot function
bootResid <- function(y.upper, B, probs){
beta=start.upper(y.upper, 1)$beta
u=start.upper(y.upper, 1)$u
v=start.upper(y.upper, 1)$v
hphi=start.upper(y.upper, 1)$hphi
hmu_ij <- start.upper(y.upper, 1)$hmu_ij
r_ij <-start.upper(y.upper, 1)$r_ij

RB <- R_pB <-RiB <- Ri_pB <- NULL

DjB <-Dj_pB <-RiB <-RDi_pB <-NULL

for (i in 1:B) {
r_ijB <- sample(r_ij, nrow(triangle_upper), replace = TRUE)
hy_ijB <- abs(r_ijB * sqrt(hmu_ij^p) + hmu_ij) 
                
betaB1 <- try(hglm(fixed=hy_ijB~as.factor(dev),random=~1|as.factor(origin), family=tweedie(var.power=p,link.power=0), 
rand.family=Gamma(log), data=triangle_upper))

hphiB <- betaB1$varFix
betaB=betaB1$fixef
uB=betaB1$ranef
vB=log(betaB1$ranef)
   
hy_ij_pB <- exp(X.lower%*%betaB+Z.lower%*%v[2:10]) 
hy_ijB <- rtweedie(length(hy_ij_pB), mu = c(hy_ij_pB), phi = hphiB, power = p) 

RB=rbind(RB, sum(hy_ijB)) 
R_pB=rbind(R_pB, sum(hy_ij_pB)) 

DjB <- rbind(DjB, tapply(hy_ijB, triangle_lower$dev, sum)) 
Dj_pB <- rbind(Dj_pB, tapply(hy_ij_pB, triangle_lower$dev, sum)) 

RiB <- rbind(RiB, tapply(hy_ijB, triangle_lower$origin, sum)) 
Ri_pB <- rbind(Ri_pB, tapply(hy_ij_pB, triangle_lower$origin, sum)) 
}

RMSEP <- sqrt(sum((R_pB - RB)^2)/B)
QAPE <-quantile(abs(R_pB - RB), probs=probs)
QAPE_j <-quantile(abs(DjB-Dj_pB), probs=probs)
QAPE_i <-quantile(abs(RiB-Ri_pB), probs=probs)
QAPE_ij <-quantile(abs(hy_ijB-hy_ij_pB), probs=probs)

out <-list(RB=RB, RpB=R_pB, rmsep=RMSEP, qape=QAPE, qapej=QAPE_j, qapei=QAPE_i, qapeij=QAPE_ij)
}

###To obtain application results:
br <- bootResid(claim, lboot, probs)
br

######### parametric bootstrap function
bootPar <- function(y.upper, B, probs) {
beta=start.upper(y.upper, 1)$beta
u=start.upper(y.upper, 1)$u
v=start.upper(y.upper, 1)$v
hphi=start.upper(y.upper, 1)$hphi

Ri1B <- Ri1BB <- R1B <- R1BB <-Dj1B <- Dj1BB <- ymodel <-NULL

for (i in 1:B) {
u1 = rtweedie(length(u), mu = c(u), phi =  start.upper(claim, 1)$phi.u, power = 2)
v1=log(u1)
hy_ij_p1B <- rtweedie(nrow(triangle_lower), mu = c(exp(X.lower%*%beta+Z.lower%*%v1[2:10])), phi = hphi, power = p)

hy_ij1B <- rtweedie(nrow(triangle_upper), mu = c(exp(X.upper%*%beta+Z.upper%*%v1)), phi = hphi, power = p)

beta1B <- hglm(fixed=hy_ij1B~as.factor(dev),random=~1|as.factor(origin), family=tweedie(var.power=p,link.power=0), 
rand.family=Gamma(log), data=triangle_upper, maxit=5)

beta1BB=beta1B$fixef
u1BB=beta1B$ranef
v1BB=log(u1BB)   
hy_ij_p1BB <- exp(X.lower%*%beta1BB+Z.lower%*%v1BB[2:10]) 

Ri1B <- rbind(Ri1B, tapply(hy_ij_p1B, triangle_lower$origin, sum))
Ri1BB <- rbind(Ri1BB, tapply(hy_ij_p1BB, triangle_lower$origin, sum))

R1B=rbind(R1B, sum(Ri1B[i,]))
R1BB=rbind(R1BB, sum(Ri1BB[i,]))

Dj1B <- rbind(Dj1B, tapply(hy_ij_p1B, triangle_lower$dev, sum))
Dj1BB <- rbind(Dj1BB, tapply(hy_ij_p1BB, triangle_lower$dev, sum)) 
}

RMSEP1 = sqrt(sum((R1B - R1BB)^2)/B)
QAPE1 = quantile(abs(R1B - R1BB), probs=probs)
QAPE1_j <-quantile(abs(Dj1B-Dj1BB), probs=probs)
QAPE1_i <-quantile(abs(Ri1B-Ri1BB), probs=probs)
QAPE1_ij <-quantile(abs(hy_ij_p1B-hy_ij_p1BB), probs=probs)


out <-list( Rtrue=R1B, Rp=R1BB, rmsep1=RMSEP1, qape1=QAPE1, qapej1=QAPE1_j, qapei1=QAPE1_i, qapeij1=QAPE1_ij)
}

###To obtain application results:
bp <- bootPar(claim, lboot, probs)
bp

######################RESULTS
#Results - residual bootstrap
br$rmsep
br$qape
br$qapej
br$qapei
br$qapeij

#Results - parametric bootstrap
bp$rmsep1
bp$qape1
bp$qapej1
bp$qapei1
bp$qapeij1

ciag=c(br$rmsep, br$qape[-5], br$qapej[-5], br$qapei[-5], br$qapeij[-5], bp$rmsep1[-5], bp$qape1[-5], bp$qapej1[-5], bp$qapei1[-5], bp$qapeij1[-5])
MX=max(ciag)
MI=min(ciag)

###Figure
options(scipen=10)
plot(c(2:5),br$qape[-5],type="b",pch=1,ylim=c(MI,1.15*MX),xlim=c(0.92,5.4),cex=1.6,lty=3,ylab="estimates of accuracy measures",xlab="",xaxt="n") ###
axis(1, at=c(1:5),labels=c(
expression(paste("R", widehat(MS),"Es")),
expression(paste("Q", widehat(AP), "Es for ", alpha,"=0.5")),
expression(paste("Q", widehat(AP), "Es for ", alpha,"=0.75")),
expression(paste("Q", widehat(AP), "Es for ", alpha,"=0.9")),
expression(paste("Q", widehat(AP), "Es for ", alpha,"=0.95"))
))
lines(c(2:5),bp$qape1[-5],type="b",pch=16,cex=1.6)###
lines(c(2:5),br$qapej[-5],type="b",pch=2,cex=1.6,lty=3)###
lines(c(2:5),bp$qapej1[-5],type="b",pch=17,cex=1.6)### 
lines(c(2:5),br$qapei[-5],type="b",pch=0,cex=1.6,lty=3)### 
lines(c(2:5),bp$qapei1[-5],type="b",pch=15,cex=1.6)### 
lines(c(2:5),br$qapeij[-5],type="b",pch=6,cex=1.6,lty=3)
lines(c(2:5),bp$qapeij1[-5],type="b",pch=25,bg=1,cex=1.6)
points(1,br$rmsep,pch=5,cex=1.6) ###
points(1,bp$rmsep1,pch=23,bg=1,cex=1.6) ###
legend(0.78,MX*1.18,c(
expression(paste("resid. boot. ", "R", widehat(MS),"E")),
expression(paste("param. boot. ", "R", widehat(MS),"E")),
expression(paste("resid. boot. ", "Q", widehat(AP), E[alpha])),
expression(paste("param. boot. ","Q", widehat(AP),E[alpha])),
expression(paste("resid. boot. ","Q", widehat(AP),E[alpha]^j)),
expression(paste("param. boot. ", "Q", widehat(AP),E[alpha]^j)),
expression(paste("resid. boot. ","Q", widehat(AP),E[alpha]^i)),
expression(paste("param. boot. ", "Q", widehat(AP),E[alpha]^i)),
expression(paste("resid. boot. ","Q", widehat(AP),E[alpha]^ij)),
expression(paste("param. boot. ", "Q", widehat(AP),E[alpha]^ij))
),
pch=c(5,23,1,16,2,17,0,15,6,25),
pt.bg=c(0,1,0,0,0,0,0,0,0,1),ncol=2)
