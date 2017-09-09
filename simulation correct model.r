##FIRST RUN THE CODE FROM THE FILE "procedure.r"

s.true <- 2

## MONTE CARLO TRUE RMSE and QAPE 
RMSEP1.boot <-QAPE1.boot <-QAPE1_j.boot <-QAPE1_i.boot <-QAPE1_ij.boot <-NULL
RMSEP2.boot <-QAPE2.boot <-QAPE2_j.boot <-QAPE2_i.boot <-QAPE2_ij.boot <-NULL

Ri.true <- Ri11 <- Dj.true <- Dj11 <- matrix(NA,s.true,9)
R.true <- R11 <- matrix(NA,s.true,1)
y <- yp <- matrix(NA,45,s.true)
hy_ij_trueM <- hy_ij_p11M <- matrix(NA,45,s.true)

time1 <- proc.time()

for (i in 1:s.true){
u.sim1 <- rtweedie(length(u), mu = c(u), phi =  start.upper(claim, 1)$phi.u, power = 2) 
v.sim1 <- log(u.sim1)

hy_ij_true <- rtweedie(nrow(triangle_lower), mu = c(exp(X.lower%*%beta+Z.lower%*%v.sim1[2:10])), phi = hphi, power = p)
hy_ij_trueM[,i]=hy_ij_true

hy_ij <- rtweedie(nrow(triangle_upper), mu = c(exp(X.upper%*%beta+Z.upper%*%v.sim1)), phi = hphi, power = p)

beta10 <- try(hglm(fixed=hy_ij~as.factor(triangle_upper$dev),random=~1|as.factor(triangle_upper$origin), 
family=tweedie(var.power=p,link.power=0), rand.family=Gamma(log), maxit=5))

beta11 <- beta10$fixef
u11 <- beta10$ranef
v11 <- log(u11)   
hy_ij_p11 <- exp(X.lower%*%beta11+Z.lower%*%v11[2:10])
hy_ij_p11M[,i] <- hy_ij_p11

y[,i] <- hy_ij_p11
yp[,i] <- hy_ij_p11

Ri.true[i,] <- tapply(hy_ij_true, triangle_lower$origin, sum)
Ri11[i,] <- tapply(hy_ij_p11, triangle_lower$origin, sum)

R.true[i] <- sum(Ri.true[i,])
R11[i] <- sum(Ri11[i,])

Dj.true [i,] <- tapply(hy_ij_true, triangle_lower$dev, sum)
Dj11[i,] <- tapply(hy_ij_p11, triangle_lower$dev, sum)

##bootsrap errors
errorRes <- bootResid(hy_ij, lboot,probs)

RMSEP1.boot <- cbind(RMSEP1.boot,errorRes$rmsep)
QAPE1.boot <- cbind(QAPE1.boot,errorRes$qape)
QAPE1_j.boot <- cbind(QAPE1_j.boot,errorRes$qapej)
QAPE1_i.boot <- cbind(QAPE1_i.boot,errorRes$qapei)
QAPE1_ij.boot <- cbind(QAPE1_ij.boot,errorRes$qapeij)

errorPar <- bootPar(hy_ij, lboot, probs)

RMSEP2.boot <- cbind(RMSEP2.boot,errorPar$rmsep1)
QAPE2.boot <- cbind(QAPE2.boot,errorPar$qape1)
QAPE2_j.boot <- cbind(QAPE2_j.boot,errorPar$qapej1)
QAPE2_i.boot <- cbind(QAPE2_i.boot,errorPar$qapei1)
QAPE2_ij.boot <- cbind(QAPE2_ij.boot,errorPar$qapeij1)
}

time2 <- proc.time()

###results
RMSEP.true <- sqrt(mean((R.true - R11)^2))
QAPE.true <- quantile(abs(R.true - R11), probs=probs)
QAPE_j.true <-quantile(abs(Dj.true-Dj11), probs=probs)
QAPE_i.true <-quantile(abs(Ri.true-Ri11), probs=probs)
QAPE_ij.true <-quantile(abs(c(hy_ij_trueM)-c(hy_ij_p11M)), probs=probs)

#############################
#relative bias of bootstrap estimator RES
RMSEP.bootr <- mean(RMSEP1.boot)
QAPE.bootr <- rowMeans(QAPE1.boot)
QAPE_j.bootr<- rowMeans(QAPE1_j.boot) 
QAPE_i.bootr<- rowMeans(QAPE1_i.boot)
QAPE_ij.bootr<- rowMeans(QAPE1_ij.boot)

bias1r <- 100*((RMSEP.bootr-RMSEP.true)/RMSEP.true)
bias2r <- 100*((QAPE.bootr-QAPE.true)/QAPE.true)
bias3r <- 100*((QAPE_j.bootr-QAPE_j.true)/QAPE_j.true)
bias4r <- 100*((QAPE_i.bootr-QAPE_i.true)/QAPE_i.true)
bias5r <- 100*((QAPE_ij.bootr-QAPE_ij.true)/QAPE_ij.true)

#relative RMSE of bootstrap estimator RES
RMSE1r <- 100*sqrt(mean((RMSEP1.boot-RMSEP.true)^2))/RMSEP.true
RMSE2r <- 100*sqrt(rowMeans((QAPE1.boot-QAPE.true)^2))/QAPE.true
RMSE3r <- 100*sqrt(rowMeans((QAPE1_j.boot-QAPE_j.true )^2))/ QAPE_j.true
RMSE4r <- 100*sqrt(rowMeans((QAPE1_i.boot-QAPE_i.true)^2))/QAPE_i.true
RMSE5r <- 100*sqrt(rowMeans((QAPE1_ij.boot-QAPE_ij.true)^2))/QAPE_ij.true

#relative bias of bootstrap estimator PAR
RMSEP.bootp <- mean(RMSEP2.boot)
QAPE.bootp <- rowMeans(QAPE2.boot)
QAPE_j.bootp <- rowMeans(QAPE2_j.boot) 
QAPE_i.bootp <- rowMeans(QAPE2_i.boot)
QAPE_ij.bootp <- rowMeans(QAPE2_ij.boot)

bias1p <- 100*((RMSEP.bootp-RMSEP.true)/RMSEP.true)
bias2p <- 100*((QAPE.bootp-QAPE.true)/QAPE.true)
bias3p <- 100*((QAPE_j.bootp-QAPE_j.true)/QAPE_j.true)
bias4p <- 100*((QAPE_i.bootp-QAPE_i.true)/QAPE_i.true)
bias5p <- 100*((QAPE_ij.bootp-QAPE_ij.true)/QAPE_ij.true)

#relative RMSE of bootstrap estimator PAR
RMSE1p <- 100*sqrt(mean((RMSEP2.boot-RMSEP.true)^2))/RMSEP.true
RMSE2p <- 100*sqrt(rowMeans((QAPE2.boot-QAPE.true)^2))/QAPE.true
RMSE3p <- 100*sqrt(rowMeans((QAPE2_j.boot-QAPE_j.true )^2))/ QAPE_j.true
RMSE4p <- 100*sqrt(rowMeans((QAPE2_i.boot-QAPE_i.true)^2))/QAPE_i.true
RMSE5p <- 100*sqrt(rowMeans((QAPE2_ij.boot-QAPE_ij.true)^2))/QAPE_ij.true

summaryALL1 <- c(RMSEP.true, QAPE.true, QAPE_j.true, QAPE_i.true, QAPE_ij.true)
summaryALL2 <- c(bias1r,bias2r,bias3r,bias4r,bias5r)
summaryALL3 <- c(RMSE1r,RMSE2r,RMSE3r,RMSE4r,RMSE5r)
summaryALL4 <- c(bias1p,bias2p,bias3p,bias4p,bias5p)
summaryALL5 <- c(RMSE1p,RMSE2p,RMSE3p,RMSE4p,RMSE5p)

#relative biases in %
rel_bias <- cbind(summaryALL2,summaryALL4)
rownames(rel_bias)=c("eRMSE","eQAPE50","eQAPE75","eQAPE90","eQAPE95","eQAPE99",
"eQAPE50_j","eQAPE75_j","eQAPE90_j","eQAPE95_j","eQAPE99_j",
"eQAPE50_i","eQAPE75_i","eQAPE90_i","eQAPE95_i","eQAPE99_i",
"eQAPE50_ij","eQAPE75_ij","eQAPE90_ij","eQAPE95_ij","eQAPE99_ij")
colnames(rel_bias)=c("residual bootstrap estimator","parametric bootstrap estimator")

#relative RMSEs in %
rel_RMSE <- cbind(summaryALL3,summaryALL5)
rownames(rel_RMSE)=c("eRMSE","eQAPE50","eQAPE75","eQAPE90","eQAPE95","eQAPE99",
"eQAPE50_j","eQAPE75_j","eQAPE90_j","eQAPE95_j","eQAPE99_j",
"eQAPE50_i","eQAPE75_i","eQAPE90_i","eQAPE95_i","eQAPE99_i",
"eQAPE50_ij","eQAPE75_ij","eQAPE90_ij","eQAPE95_ij","eQAPE99_ij")
colnames(rel_RMSE)=c("residual bootstrap estimator","parametric bootstrap estimator")

#Monte Carlo time-consuming
time2-time1 #Monte Carlo loop
(time2-time1)/s.true #one iteration of Monte Carlo

###summary
summaryT <- rbind(cbind(RMSEP.true, QAPE.true, QAPE_j.true, QAPE_i.true, QAPE_ij.true),
cbind(RMSEP.bootr, QAPE.bootr, QAPE_j.bootr, QAPE_i.bootr, QAPE_ij.bootr),
cbind(bias1r, bias2r, bias3r, bias4r, bias5r),
cbind(RMSEP.bootp, QAPE.bootp, QAPE_j.bootp, QAPE_i.bootp, QAPE_ij.bootp),
cbind(bias1p, bias2p, bias3p, bias4p, bias5p))

r1 <- as.numeric(abs(R.true - R11))
r2 <- as.numeric(abs(Dj.true-Dj11))
r3 <- as.numeric(abs(Ri.true-Ri11))
r4 <- as.numeric(abs(hy_ij_trueM-hy_ij_p11M))
length(r1)
length(r2)
length(r3)
length(r4)

w <- c(r1,r2,r3,r4)
gr <- c(rep(1,length(r1)),rep(2,length(r2)),rep(3,length(r3)),rep(4,length(r4)))
nazw <- c("QAPE","QAPE_j","QAPE_i","QAPE_ij")

#
options(scipen=10)
plot(0:1,0:1,type="n",xlim=c(0.5,4.5),ylim=c(min(w),max(w)), axes=FALSE,ann=FALSE)
vioplot(r1,r2,r3,r4,col="LIGHTGREY",drawRect=FALSE,names=nazw,axes=FALSE,ann=FALSE,add=T)
title(ylab="values of absolute prediction errors")
axis(side=1,at=1:4,labels=c(
expression(paste("| ", hat(R)," - ","R |")),
expression(paste("| ",hat(R)[j]," - ",R[j]," |")),
expression(paste("| ",hat(R)[i]," - ",R[i]," |")),
expression(paste("| ",hat(Y)[ij]," - ",Y[ij]," |"))
))
axis(side=2,at=c(0,500000,1000000,1500000),labels=c(0,500000,1000000,1500000))
points(x=rep(1,4),y=QAPE.true[-5],cex=1.3,pch=21,bg="grey35")
points(x=rep(2,4),y=QAPE_j.true[-5],cex=1.3,pch=24,bg="grey35")
points(x=rep(3,4),y=QAPE_i.true[-5],cex=1.3,pch=22,bg="grey35")
points(x=rep(4,4),y=QAPE_ij.true[-5],cex=1.3,pch=25,bg="grey35")
points(x=1,y=RMSEP.true,cex=2,pch=23,bg="grey35")
legend(2.3,max(w),c(
expression(paste("RMSE")),
expression(paste("QAP", E[alpha],"  for ",alpha, "={0.5, 0.75, 0.9, 0.95}")),
expression(paste("QAP",E[alpha]^j,"  for ",alpha, "={0.5, 0.75, 0.9, 0.95}")),
expression(paste("QAP",E[alpha]^i,"  for ",alpha, "={0.5, 0.75, 0.9, 0.95}")),
expression(paste("QAP",E[alpha]^ij,"  for ",alpha, "={0.5, 0.75, 0.9, 0.95}"))
),
pch=c(23,21,24,22,25),
pt.bg=c("grey35","grey35","grey35","grey35","grey35"),
pt.cex=c(1.5,1,1,1,1))

#Tables
round(rel_bias,2)
round(rel_RMSE,2)

#Relative bias of predictor
relative.bias.true <- 100*mean(R11-R.true )/R.true


