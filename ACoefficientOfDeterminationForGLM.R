####CREATION OF R SUQARED OF MAGEE Rlr####
#l(Y,mu) = logverosimiglianza del modello E[Y|X]=mu(X) per i dati osservati (Y;X)
#l(Y,mu1) = logverosimiglianza del modello con sola intercetta

Rlr = function(fit,distr){
  fam = distr
  dati = fit$model
  n = fit$df.null+1
  fit0 = update(fit,.~1,family=fam,data=dati)
  value = 1 - exp((2/n)*(as.numeric(logLik(fit0))-as.numeric(logLik(fit))))
  return(value)
}

####CREATION OF R SUQARED OF NAGELKERKE####

Rn = function(fit,distr){
  fam = distr
  dati = fit$model
  n = fit$df.null+1
  fit0 = update(fit,.~1, family=fam, data=dati)
  value = (1 - exp((2/n)*(logLik(fit0)-logLik(fit))))/(1-exp((2/n)*(logLik(fit0))))
  return(as.numeric(value))
}

####CREATION OF R SUQARED OF KULLBACK-LEIBLER####

Rkl = function(fit){
  {
    sse1 = fit$deviance
    sse0 = fit$null.deviance
    rsq = 1-(sse1/sse0)
    return(rsq)
  }
}

#ADJUSTED

Rkladj = function(fit){
  sse1 = fit$deviance
  sse0 = fit$null.deviance
  rsq.adj = 1-((sse1/sse0)*(fit$df.null/fit$df.residual))
  return(rsq.adj)
}

####CREATION OF R SUQARED OF DABAO ZHANG####
#possibility to select different families 

##funzione per calcolare la varianza V'(mu) che serve per il calcolo della distanza
deriv.fun = function(x,fam){
  if(fam=="normal"){
    v2 = 0
    v1 = 0
    v0 = 1
  }
  else if(fam=="binomial"){
    v2 = -1
    v1 = 1
  }
  else if(fam=="poisson"){
    v2 = 0
    v1 = 1
  }
  else if(fam=="Gamma"){
    v1 = 0
    v2 = 1
  }
  else if(fam=="quasibinomial"){ #V(mu) = mu*(1-mu)
    v2 = -1
    v1 = 1
  }
  else if(fam=="quasipoisson"){ #V(mu) = mu
    v2 = 0
    v1 = 1
  }
  else if(fam=="negative.binomial"){   # V(mu)=mu+mu*mu/theta
    v2 = 1/theta
    v1 = 1
  }
  ris = (2*v2*x)+v1
  return(ris)
}

#calcolo distanza fra due punti
dist.v = function(da,db,fam){
  if(fam=="normal"){
    v1=0
    dvab = (1+v1^2)*((db-da)^2)
  }
  else if(fam=="binomial"){
    v2=-1
    dvab = (1/(16*(v2^2)))*(log((db+sqrt(1+db^2))/(da+sqrt(1+da^2)))+(db*sqrt(1+db^2))-(da*sqrt(1+da^2)))^2
  }
  else if(fam=="poisson"){
    v1=1
    dvab = (1+v1^2)*((db-da)^2)
  }
  else if(fam=="Gamma"){
    v2=1
    dvab = (1/(16*(v2^2)))*(log((db+sqrt(1+db^2))/(da+sqrt(1+da^2)))+(db*sqrt(1+db^2))-(da*sqrt(1+da^2)))^2
  }
  else if(fam=="quasibinomial"){ 
    v2 = -1
    v1 = 1
    dvab = (1/(16*(v2^2)))*(log((db+sqrt(1+db^2))/(da+sqrt(1+da^2)))+(db*sqrt(1+db^2))-(da*sqrt(1+da^2)))^2
  }
  else if(fam=="quasipoisson"){ 
    v2 = 0
    v1 = 1
    dvab = (1+v1^2)*((db-da)^2)
  }
  else if(fam=="negative.binomial"){   # V(mu)=mu+mu*mu/theta
    v2 = 1/theta
    dvab = (1/(16*(v2^2)))*(log((db+sqrt(1+db^2))/(da+sqrt(1+da^2)))+(db*sqrt(1+db^2))-(da*sqrt(1+da^2)))^2
  }
  return(dvab)
}

Rv = function(y,fit,fam){
  fit.v = fitted(fit)
  fit0 = fitted(update(fit,.~1))
  if(fam=="normal"|fam=="poisson"){
    num = dist.v(y,fit.v,fam)
    den = dist.v(y,fit0,fam)
  }
  else{  #per le altre distrib. devo calcolare V'(.) per trovare dv, quindi è più complesso
    da = deriv.fun(y,fam)
    db = deriv.fun(fit.v,fam)
    num = dist.v(da,db,fam)
    dc = deriv.fun(fit0,fam)
    den = dist.v(da,dc,fam)
  }
  value = 1 - (sum(num)/sum(den))
  return(value)
}

#ADJUSTED

Rvadj = function(y,fit,fam){
  n = length(y)
  p = length(coef(fit))
  fit.v = fitted(fit)
  fit0 = fitted(glm(y~1,family=fam))
  if(fam=="normal"|fam=="poisson"){
    num = dist.v(y,fit.v,fam)
    den = dist.v(y,fit0,fam)
  }
  else{
    da = deriv.fun(y,fam)
    db = deriv.fun(fit.v,fam)
    num = dist.v(da,db,fam)
    dc = deriv.fun(fit0,fam)
    den = dist.v(da,dc,fam)
  }
  value = 1 - ((sum(num)/(n-p))/(sum(den)/(n-1)))
  return(value)
}

#EMPIRICAL STUDIES#
####GENERATION OF 100 DATASET####
library(performance)

Rsgen = function(Nsim=100,beta.max=8,fam="binomial",Seed=281196,display=TRUE,compare=FALSE){
  set.seed(Seed)
  X = c(rep(1,25),rep(-1,25))
  beta = seq(0,beta.max,by=0.1) #beta.max=8 -> così vedo che comportamento hanno gli indici a +Inf
  Rlr.final = Rn.final = Rkl.final = Rv.final = Rkladj.final = Rvadj.final = rep(NA,length(beta)) #pre allocazione dei vettori Rquadro finali (per ogni beta)
  for(i in 1:length(beta)){
    Rlr.gen = Rn.gen = Rkl.gen = Rv.gen = Rkladj.gen = Rvadj.gen = rep(NA,length(beta))
    if(i%%5==0) print(i) #a che punto siamo?
    current.beta = beta[i]
    eta = current.beta*X
    if(fam=="binomial")  {p = plogis(eta)}
    else if(fam=="poisson")  {p = exp(eta)}
    else if(fam=="Gamma")  {p = 0.01/(current.beta+1+eta)} 
    Y = matrix(NA,50,Nsim)
    for(j in 1:Nsim){
      if(fam=="binomial")  {Y[,j] = rbinom(50,1,p)}
      else if(fam=="poisson")  {Y[,j] = rpois(50,p)}
      else if(fam=="Gamma") {Y[,j] = rgamma(50,shape=100,scale=p)} 
      current.model = glm(Y[,j]~X, family=fam)
      Rlr.gen[j] = Rlr(current.model,distr=fam)
      #Rn.gen[j] = Rn(current.model,distr=fam)
      Rn.gen[j] = r2_nagelkerke(current.model)
      Rkl.gen[j] = Rkl(current.model)
      Rv.gen[j] = Rv(Y[,j],current.model,fam)
      Rkladj.gen[j] = Rkladj(current.model)
      Rvadj.gen[j] = Rvadj(Y[,j],current.model,fam)
    }
    Rlr.final[i] = mean(Rlr.gen) 
    Rn.final[i] = mean(Rn.gen) 
    Rkl.final[i] = mean(Rkl.gen) 
    Rv.final[i] = mean(Rv.gen) 
    Rkladj.final[i] = mean(Rkladj.gen) 
    Rvadj.final[i] = mean(Rvadj.gen) 
  }
  R.square = list(Rlr = Rlr.final, Rn = Rn.final, Rkl = Rkl.final, Rv = Rv.final, Rkladj = Rkladj.final, Rvadj = Rvadj.final)
  if(display){
    plot(beta,Rlr.final,type="l",lwd=1,lty=1,ylim = c(0,1),xlab=expression(beta),ylab=expression("Coefficienti R "^2))
    lines(beta,Rn.final,lwd=1,lty=5)
    lines(beta,Rkl.final,lwd=2,lty=1)
    lines(beta,Rv.final,lwd=2,lty=5)
    legend((beta.max-1.2),0.2,c("RLR","RN","RKL","RV"),lty=c(1,5,1,5),lwd=c(1,1,2,2),bty="n")
  }
  if(compare){
    par(mfrow=c(1,2))
    plot(beta,Rkl.final,type="l",lwd=1,ylim = c(0,1),xlab=expression(beta),ylab="Coefficiente R quadro")
    lines(beta,Rkladj.final,lwd=1,col=2)
    plot(beta,Rv.final,type="l",lwd=1,ylim = c(0,1))
    lines(beta,Rvadj.final,lwd=1,col=2)
  }
  par(mfrow=c(1,1))
  return(R.square)
}

Rfinal.bin = Rsgen(beta.max=6,fam="binomial",compare=F)
Rfinal.bin
Rfinal.poi = Rsgen(beta.max=6,fam="poisson",compare=F)
Rfinal.poi
Rfinal.gamma = Rsgen(beta.max=6,fam="Gamma",compare=F)
Rfinal.gamma

####ADJUSTED COEFFICIENT SIMULATION####

#per generare le x, usiamo la funzione messa a disposizione dallo stesso Zhang

library(rsq) #libreria per estrazione delle x

Radjtest3 = function(Nsim=100,beta.max=5,Seed=281196,fam="binomial"){
  set.seed(Seed)
  beta = seq(0,beta.max,by=0.1)
  Rkladj.final = Rvadj.final = Rkl.final = Rv.final = rep(NA,length(beta))
  for(i in 1:length(beta)){
    Rkladj.gen = Rvadj.gen = Rkl.gen = Rv.gen = rep(NA,length(beta))
    if(i%%5==0) print(i) #a che punto siamo?
    current.beta = beta[i]
    if(fam=="binomial") {X = simglm(family="binomial",lambda=current.beta,n=50,p=3)}
    else if(fam=="poisson") {X = simglm(family="poisson",lambda=current.beta,n=50,p=3)}
    else if(fam=="Gamma") {X = simglm(family="Gamma",lambda=current.beta,n=50,p=3)}
    X1 = X$yx$x.1
    X2 = X$yx$x.2
    X3 = X$yx$x.3
    X123 = cbind(X1,X2,X3)
    eta = X123%*%as.matrix(X$beta[-1])
    if(fam=="binomial") {p = plogis(eta)}
    else if(fam=="poisson") {p = exp(eta)}
    else if(fam=="Gamma") {p = 0.01/(current.beta+1+eta)}
    Y = matrix(NA,50,Nsim)
    for(j in 1:Nsim){
      if(fam=="binomial") {Y[,j] = rbinom(50,1,p)}
      else if(fam=="poisson") {Y[,j] = rpois(50,p)}
      else if(fam=="Gamma") {Y[,j] = rgamma(50,shape=100,scale=p)}
      complex.model = glm(Y[,j]~X123, family=fam)
      simple.model = glm(Y[,j]~X123[,1], family=fam)
      Rkladjt = Rkladj(complex.model)
      Rvadjt = Rvadj(Y[,j],complex.model,fam)
      Rklt = Rkl(complex.model)
      Rvt = Rv(Y[,j],complex.model,fam)
      Rkladj.gen[j] = Rkladjt - Rkladj(simple.model)
      Rvadj.gen[j] = Rvadjt - Rvadj(Y[,j],simple.model,fam)
      Rkl.gen[j] = Rklt - Rkl(simple.model)
      Rv.gen[j] = Rvt - Rv(Y[,j],simple.model,fam)
    }
    Rkladj.final[i] = mean(na.omit(Rkladj.gen)) 
    Rvadj.final[i] = mean(Rvadj.gen) 
    Rkl.final[i] = mean(Rkl.gen) 
    Rv.final[i] = mean(Rv.gen)
  }
  R.square = list(Rkladj = Rkladj.final, Rvadj = Rvadj.final, Rkl = Rkl.final, Rv = Rv.final)
  return(R.square)
}


layout(matrix(c(1,2,3,4), nrow=4, byrow=TRUE), heights=c(1.2,1.2,1.2,0.3))
par(mai=rep(0.3, 4))

RB = Radjtest3(fam="binomial")
RB$Rkladj
beta = seq(0,5,by=0.1)
lo = smooth.spline(RB$Rkladj~beta,spar=0.35)
plot(beta, RB$Rkladj,xlab=expression(beta),ylab=expression(paste("R"^2,"(X"[1],"X"[2],"X"[3],")",-"R"^2,"(X"[1],")")),type="n",ylim=c(-0.015,0.05),main="Modello Binomiale")
lines(lo,lwd=1,lty=1)
lo = smooth.spline(RB$Rkl~beta,spar=0.35)
lines(lo,lwd=2,lty=1)
lo = smooth.spline(RB$Rvadj~beta,spar=0.35)
lines(lo,lwd=1,lty=5)
lo = smooth.spline(RB$Rv~beta,spar=0.35)
lines(beta, RB$Rv,lwd=2,lty=5)

RP = Radjtest3(fam="poisson")
RP$Rkladj
beta = seq(0,5,by=0.1)
lo = smooth.spline(RP$Rkladj~beta,spar=0.2)
plot(beta, RP$Rkladj,xlab=expression(beta),ylab=expression(paste("R"^2,"(X"[1],"X"[2],"X"[3],")",-"R"^2,"(X"[1],")")),type="n",ylim=c(-0.015,0.055),main="Modello Poisson")
lines(lo,lwd=1,lty=1)
lo = smooth.spline(RP$Rkl~beta,spar=0.2)
lines(lo,lwd=2,lty=1)
lo = smooth.spline(RP$Rvadj~beta,spar=0.2)
lines(lo,lwd=1,lty=5)
lo = smooth.spline(RP$Rv~beta,spar=0.2)
lines(beta, RP$Rv,lwd=2,lty=5)

RG = Radjtest3(fam="Gamma")
RG$Rkladj
beta = seq(0,5,by=0.1)
lo = smooth.spline(RG$Rkladj~beta,spar=0.1)
plot(beta, RG$Rkladj,xlab=expression(beta),ylab=expression(paste("R"^2,"(X"[1],"X"[2],"X"[3],")",-"R"^2,"(X"[1],")")),type="n",ylim=c(-0.015,0.055),main="Modello Gamma")
lines(lo,lwd=1,lty=1)
lo = smooth.spline(RG$Rkl~beta,spar=0.1)
lines(lo,lwd=2,lty=1)
lo = smooth.spline(RG$Rvadj~beta,spar=0.1)
lines(lo,lwd=1,lty=5)
lo = smooth.spline(RG$Rv~beta,spar=0.1)
lines(beta, RG$Rv,lwd=2,lty=5)

par(mai=c(0,0,0,0))
plot.new()
legend(x="center",ncol=4,legend=c("Rkladj","Rkl","Rvadj","Rv"),lty=c(1,1,5,5),
       lwd=c(1,2,1,2))

####TEST BINOMIAL####

Rfinal.bin$Rlr[length(seq(0,6,by=0.1))] #0.744457 -> tende a 0.75
Rfinal.bin$Rlr[1] #dovrebbe essere zero -> 0.01704168

Rfinal.bin$Rn[length(seq(0,6,by=0.1))] #0.9926695
Rfinal.bin$Rn[1] #dovrebbe essere zero -> 0.02285787

Rfinal.bin$Rkl[length(seq(0,6,by=0.1))] #0.9854457
Rfinal.bin$Rkl[1] #dovrebbe essere zero -> 0.01277485

Rfinal.bin$Rv[length(seq(0,6,by=0.1))] #0.9907529
Rfinal.bin$Rv[1] #dovrebbe essere zero -> 0.01711641

####TEST POISSON####

Rfinal.poi$Rlr[length(seq(0,6,by=0.1))] #1
Rfinal.poi$Rlr[1] #dovrebbe essere zero -> 0.01758655

Rfinal.poi$Rn[length(seq(0,6,by=0.1))] #1
Rfinal.poi$Rn[1] #dovrebbe essere zero -> 0.01898508

Rfinal.poi$Rkl[length(seq(0,6,by=0.1))] #0.9982738
Rfinal.poi$Rkl[1] #dovrebbe essere zero -> 0.01566021

Rfinal.poi$Rv[length(seq(0,6,by=0.1))] #0.9952668
Rfinal.poi$Rv[1] #dovrebbe essere zero -> 0.01798161

####APPLICATION TO A REAL DATASET####
#crabs.dat

dati = read.table("crabs.dat",header=T)
head(dati)
dati$color = as.factor(dati$color)
dati$spine = as.factor(dati$spine)
attach(dati)
dati = dati[,-1]
#facciamo una breve analisi esplorativa
table(y)

satyes = as.numeric(y!=0) #per la binomiale uso questa var risposta!
satyes
head(dati)

####TEST REAL BINOMIAL####

fit0 = glm(satyes~weight+width+color+spine+I(weight^2)+I(width^2),family = binomial)
AIC(fit0)
BIC(fit0)
Rlr(satyes,fit0,fam="binomial")
Rn(satyes,fit0,fam="binomial")
Rkl(fit0)
Rkladj(fit0)
Rv(satyes,fit0,fam="binomial")
Rvadj(satyes,fit0,fam="binomial")

fit1 = glm(satyes~weight+width+color+I(weight^2)+I(width^2),family = binomial)
AIC(fit1)
BIC(fit1)
Rlr(satyes,fit1,fam="binomial")
Rn(satyes,fit1,fam="binomial")
Rkl(fit1)
Rkladj(fit1)
Rv(satyes,fit1,fam="binomial")
Rvadj(satyes,fit1,fam="binomial")

fit2 = glm(satyes~width+color+I(width^2),family = binomial)
AIC(fit2)
BIC(fit2)
Rlr(satyes,fit2,fam="binomial")
Rn(satyes,fit2,fam="binomial")
Rkl(fit2)
Rkladj(fit2)
Rv(satyes,fit2,fam="binomial")
Rvadj(satyes,fit2,fam="binomial")

fit3 = glm(satyes~weight+color+I(weight^2),family = binomial)
AIC(fit3)
BIC(fit3)
Rlr(satyes,fit3,fam="binomial")
Rn(satyes,fit3,fam="binomial")
Rkl(fit3)
Rkladj(fit3)
Rv(satyes,fit3,fam="binomial")
Rvadj(satyes,fit3,fam="binomial")

fit4 = glm(satyes~weight+width+I(weight^2)+I(width^2),family = binomial)
AIC(fit4)
BIC(fit4)
Rlr(satyes,fit4,fam="binomial")
Rn(satyes,fit4,fam="binomial")
Rkl(fit4)
Rkladj(fit4)
Rv(satyes,fit4,fam="binomial")
Rvadj(satyes,fit4,fam="binomial")

fit8 = glm(satyes~color,family = binomial)
AIC(fit8)
BIC(fit8)
Rlr(satyes,fit8,fam="binomial")
Rn(satyes,fit8,fam="binomial")
Rkl(fit8)
Rkladj(fit8)
Rv(satyes,fit8,fam="binomial")
Rvadj(satyes,fit8,fam="binomial")

fit7 = glm(satyes~spine,family = binomial)
AIC(fit7)
BIC(fit7)
Rlr(satyes,fit7,fam="binomial")
Rn(satyes,fit7,fam="binomial")
Rkl(fit7)
Rkladj(fit7)
Rv(satyes,fit7,fam="binomial")
Rvadj(satyes,fit7,fam="binomial")

fit5 = glm(satyes~weight,family = binomial)
AIC(fit5)
BIC(fit5)
Rlr(satyes,fit5,fam="binomial")
Rn(satyes,fit5,fam="binomial")
Rkl(fit5)
Rkladj(fit5)
Rv(satyes,fit5,fam="binomial")
Rvadj(satyes,fit5,fam="binomial")

fit6 = glm(satyes~width,family = binomial)
AIC(fit6)
BIC(fit6)
Rlr(satyes,fit6,fam="binomial")
Rn(satyes,fit6,fam="binomial")
Rkl(fit6)
Rkladj(fit6)
Rv(satyes,fit6,fam="binomial")
Rvadj(satyes,fit6,fam="binomial")

####DATASET REALE POISSON####

fit0 = glm(y~weight+width+color+spine+I(weight^2)+I(width^2),family = poisson)
AIC(fit0)
BIC(fit0)
Rlr(y,fit0,fam="poisson")
Rn(y,fit0,fam="poisson")
Rkl(fit0)
Rv(y,fit0,fam="poisson")
Rkladj(fit0)
Rvadj(y,fit0,fam="poisson")

#eccetera: i calcoli sembrano combaciare con quelli dell'articolo

####ANALYSIS OF CREDIT.DAT####
#BINOMIAL MODEL

credit = read.table("credit.dat",header=TRUE)

#ANALISI ESPLORATIVA
head(credit)
str(credit)
credit$Y = as.numeric(credit$Y=="mal")

#mal=1, buen=0
#casi e controlli fissati, si vedono come variano le concomitanti rispetto ai due strati
attach(credit)

#sembra che i clienti insolventi tendano a chiedere prestiti più a lungo termine:
boxplot(Mes~Y)
#ammontare del prestito per categoria di debitori
boxplot(DM~Y)
table(Y,Cuenta)
barplot(t(prop.table(table(Y,Cuenta),1)),beside=T,legend=T,ylim=c(0,0.6))
#eccetera per le altre variabili

#la logistica permette di fare analisi con y fissata e x osservate
#adottiamo un modello di regressione logistica con selezione tipo backward (come in crabs.dat)

#SELEZIONE BACKWARD
#STEP 1: modello completo

glm0 = glm(Y~.,family = "binomial",data=credit)
summary(glm0)
AIC(glm0)
BIC(glm0)
Rlr(glm0,distr="binomial")
Rn(glm0,distr="binomial")
Rkl(glm0)
Rv(Y,glm0,fam="binomial")
Rkladj(glm0)
Rvadj(Y,glm0,fam="binomial")
drop1(glm0,test="Chisq")

#STEP 2: esce DM

glm1 = update(glm0,.~.-DM)
summary(glm1)
AIC(glm1)
BIC(glm1)
Rlr(glm1,distr="binomial")
Rn(glm1,distr="binomial")
Rkl(glm1)
Rv(Y,glm1,fam="binomial")
Rkladj(glm1)
Rvadj(Y,glm1,fam="binomial")
drop1(glm1,test="Chisq")

#STEP 3: esce Sexo

glm2 = update(glm1,.~.-Sexo)
summary(glm2)
AIC(glm2)
BIC(glm2)
Rlr(glm2,distr="binomial")
Rn(glm2,distr="binomial")
Rkl(glm2)
Rv(Y,glm2,fam="binomial")
Rkladj(glm2)
Rvadj(Y,glm2,fam="binomial")
drop1(glm2,test="Chisq")

#STEP 4: faccio uscire Uso per vedere se indici colgono la differenza

glm3 = update(glm2,.~.-Uso)
summary(glm3)
AIC(glm3)
BIC(glm3)
Rlr(glm3,distr="binomial")
Rn(glm3,distr="binomial")
Rkl(glm3)
Rv(Y,glm3,fam="binomial")
Rkladj(glm3)
Rvadj(Y,glm3,fam="binomial")
drop1(glm3,test="Chisq")

glm4 = update(glm3,.~.-Estc)
summary(glm4)
AIC(glm4)
BIC(glm4)
Rlr(glm4,distr="binomial")
Rn(glm4,distr="binomial")
Rkl(glm4)
Rv(Y,glm4,fam="binomial")
Rkladj(glm4)
Rvadj(Y,glm4,fam="binomial")
drop1(glm4,test="Chisq")

glm5 = update(glm1,.~.-Estc)
summary(glm5)
AIC(glm5)
BIC(glm5)
Rlr(glm5,distr="binomial")
Rn(glm5,distr="binomial")
Rkl(glm5)
Rv(Y,glm5,fam="binomial")
Rkladj(glm5)
Rvadj(Y,glm5,fam="binomial")
drop1(glm5,test="Chisq")

glm6 = update(glm1,.~.-Uso)
summary(glm6)
AIC(glm6)
BIC(glm6)
Rlr(glm6,distr="binomial")
Rn(glm6,distr="binomial")
Rkl(glm6)
Rv(Y,glm6,fam="binomial")
Rkladj(glm6)
Rvadj(Y,glm6,fam="binomial")
drop1(glm6,test="Chisq")

glm7 = update(glm2,.~.-Estc)
summary(glm7)
AIC(glm7)
BIC(glm7)
Rlr(glm7,distr="binomial")
Rn(glm7,distr="binomial")
Rkl(glm7)
Rv(Y,glm7,fam="binomial")
Rkladj(glm7)
Rvadj(Y,glm7,fam="binomial")
drop1(glm7,test="Chisq")

####ANALYSIS DATASET CRABS.DAT####
#POISSON MODEL

crabs = read.table("crabs.dat",header=T)
head(crabs)
crabs$color = as.factor(crabs$color)
crabs$spine = as.factor(crabs$spine)
attach(crabs)
crabs = crabs[,-1]
#facciamo una breve analisi esplorativa
table(y)
head(crabs)

#STEP 1: MODELLO COMPLETO
fit0 = glm(y~weight+width+color+spine+I(weight^2)+I(width^2),family = poisson)
AIC(fit0)
BIC(fit0)
Rlr(fit0,distr="poisson")
Rn(fit0,distr="poisson")
Rkl(fit0)
Rv(y,fit0,fam="poisson")
Rkladj(fit0)
Rvadj(y,fit0,fam="poisson")
drop1(fit0,test="Chisq")

#STEP 2: esce spine

fit1 = update(fit0,.~.-spine)
summary(fit1)
AIC(fit1)
BIC(fit1)
Rlr(fit1,distr="poisson")
Rn(fit1,distr="poisson")
Rkl(fit1)
Rv(y,fit1,fam="poisson")
Rkladj(fit1)
Rvadj(y,fit1,fam="poisson")
drop1(fit1,test="Chisq")

#STEP 3: esce I(weight^2)

fit2 = update(fit1,.~.-I(weight^2))
summary(fit2)
AIC(fit2)
BIC(fit2)
Rlr(y,fit2,fam="poisson")
Rn(y,fit2,fam="poisson")
Rkl(fit2)
Rv(y,fit2,fam="poisson")
Rkladj(fit2)
Rvadj(y,fit2,fam="poisson")
drop1(fit2,test="Chisq")

#STEP 4: esce color

fit3 = update(fit2,.~.-color)
summary(fit3)
AIC(fit3)
BIC(fit3)
Rlr(fit3,distr="poisson")
Rn(fit3,distr="poisson")
Rkl(fit3)
Rv(y,fit3,fam="poisson")
Rkladj(fit3)
Rvadj(y,fit3,fam="poisson")
drop1(fit3,test="Chisq")


#STEP 4: esce color

fit4 = update(fit1,.~.-width)
fit4 = update(fit4,.~.-I(width^2))
summary(fit4)
AIC(fit4)
BIC(fit4)
Rlr(fit4,distr="poisson")
Rn(fit4,distr="poisson")
Rkl(fit4)
Rv(y,fit4,fam="poisson")
Rkladj(fit4)
Rvadj(y,fit4,fam="poisson")
drop1(fit4,test="Chisq")


####ANALYSIS DATASET RATS.DAT####
#MODELLO DI QUASI VEROSIMIGLIANZA

rats = read.table("Rats.dat",header=T)
head(rats)

#effetto regime alimentare sullo sviluppo fetale dei ratti
#ratto femmina assegnato casualmente a uno dei quattro gruppi
#s = numero di feti morti nella nidiata
#n = numero totale dei feti nella nidiata
#h = livello omoglobina della madre

str(rats)
attach(rats)

#y = variabile d'interesse = probabilità di morte del feto: s/n

y = s/n

plot(y~as.factor(group))

rats0 = glm(y~group*h, weights=n, family=binomial)
summary(rats0)

pers.res = resid(rats0, type="pearson")
sum(pers.res^2)
rats0$df.residual
phitilde = sum(pers.res^2)/rats0$df.residual
phitilde #2.75 >> 1

#calcolo comunque i coefficienti

AIC(rats0)
BIC(rats0)
Rlr(rats0,distr="binomial")
Rn(rats0,distr="binomial")
Rkl(rats0)
Rv(y,rats0,fam="binomial")
Rkladj(rats0)
Rvadj(y,rats0,fam="binomial")
drop1(rats0,test="Chisq")

#sembra quindi esserci sovradispersione -> stretto legame con parametro di dispersione, formula 2.53

library(aod)
rats0.ql = glm(y~group, family=quasibinomial, weights=n)
summary(rats0.ql)
Rv(y,rats0.ql,fam="quasibinomial")
Rvadj(y,rats0.ql,fam="quasibinomial")

rats1.ql = glm(y~h, family=quasibinomial, weights=n)
summary(rats1.ql)
Rv(y,rats1.ql,fam="quasibinomial")
Rvadj(y,rats1.ql,fam="quasibinomial")

rats2.ql = glm(y~group+h, family=quasibinomial, weights=n)
summary(rats2.ql)
Rv(y,rats2.ql,fam="quasibinomial")
Rvadj(y,rats2.ql,fam="quasibinomial")

rats3.ql = glm(y~group*h, family=quasibinomial, weights=n)
summary(rats3.ql)
Rv(y,rats3.ql,fam="quasibinomial")
Rvadj(y,rats3.ql,fam="quasibinomial")
