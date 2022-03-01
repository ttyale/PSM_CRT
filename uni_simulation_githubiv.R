################################################################################################################################################################################################
##### Simulation code for the joint model with principal stratification model and outcome regression model with a continuous outcome and two study covariates in cluster randomized trials #####
################################################################################################################################################################################################

require("lme4")
require("mvtnorm")
###################################
## 1. Data generating process  ####
###################################

# N is the total number of individuals
# s is the number of clusters. 

gendata<-function(N=1500,s=60){
    
  m=N/s  #mean cluster size
  tau2<-1 #cluster variance
  sigma2<-5 #outcome variance
 
  cl<-rep(1:s,each=m) #cluster index
  x1<-rnorm(N,0,2) #covariate 1: normal continuous 
  x2<-runif(N,-5,5)  #covariate 2: uniform continuous 
  X<-as.matrix(data.frame(rep(1,N),x1,x2))   #design matrix
  
  #strata model parameters
  beta=c(-1,0.3,0.5) #group00 (never survivor) vs group11 (always survivor)
  gamma=c(-0.8,0.6,0.4) #group10 (protected) vs group11 (always survivor)
  
  
  #principal strata probability
  p00=exp(X%*%beta)/(exp(X%*%beta)+exp(X%*%gamma)+rep(1,N))
  p10=exp(X%*%gamma)/(exp(X%*%beta)+exp(X%*%gamma)+rep(1,N))
  p11=1-p00-p10
  
  #principal strata membership matrix
  Gm<-t(apply(cbind(p00,p10,p11),1, function(x) as.numeric(rmultinom(1,1,x))))
  
  #individual treatment indicator
  d<-rep(sample(rep(c(0,1),each=s/2),replace=F),each=m)
  
  #outcome model parameters
  alpha111<-c(1.5,0.5,0.8)  #always survivor in treatment  
  alpha101<-c(0.2,0.3,0.6)  #protected in treatment
  alpha110<-c(-1.5,0.9,0.5) #always survivor in control  
  
  #residual error and cluster-level random effect
  resid<-rnorm(N,0,sqrt(sigma2))
  eta<-rnorm(s,0,sqrt(tau2))
  etan<-rep(eta,each=m)
  
  #continuous outcome
  Y<-rep(NA,N) #outcome variable (NA indicates death)
  Y[intersect(which(Gm[,3]==1),which(d==1))]<-X[intersect(which(Gm[,3]==1),which(d==1)),]%*%alpha111
  Y[intersect(which(Gm[,2]==1),which(d==1))]<-X[intersect(which(Gm[,2]==1),which(d==1)),]%*%alpha101
  Y[intersect(which(Gm[,3]==1),which(d==0))]<-X[intersect(which(Gm[,3]==1),which(d==0)),]%*%alpha110
  Y<-Y+etan+resid

  #principal strata membership (0: always non-survivor,1: protected, 2: always survivor)
  G=Gm[,2]+2*Gm[,3]
  
  #survivor average causal effect
  TE<-colMeans(X[which(G==2),]%*%alpha111-X[which(G==2),]%*%alpha110)
  
  return(data.frame(cbind(Y,x1,x2,d,cl,G,TE)))
}

###################################
## 2. MCMC model estimation  ######
###################################

#df: the simulated data 
#S: length of chain
#dau1: scaling parameter for chain move

uvsampler<-function(df=df,S=10000,dau1=0.05){
  
  #df<-gendata(N=1500,s=60)
  
  # extract data information
  N=nrow(df) # number of individuals
  cl=df$cl   # cluster ID
  cs=as.numeric(table(cl)) # size of cluster
  s<-length(unique(df$cl)) # number of clusters
  
  mbar=N/s #average cluster size
  p=3 # length of regression coefficients
  X<-as.matrix(cbind(rep(1,N),df$x1,df$x2)) # design matrix
  Y=df$Y # outcome
  Nst<-length(which(!is.na(Y))) # non-missing individuals
  D<-df$d # treatment indicator
  

    #initial values for outcome regression model
  lmm<-summary(lmer(Y~x1+x2+(1|cl),data=df))
  alpha111<-alpha101<-alpha110<-as.numeric(lmm$coefficients[,1])
  sigma2<-lmm$sigma^2
  tau2<-as.numeric(lmm$varcor$cl)^2
  # random effect at individual level
  etan<-rep(eta<-rnorm(s,0,sqrt(tau2)),each=mbar)
  G=df$G
  
  #true model parameters
  #alpha111<-c(1.5,0.5,0.8)   
  #alpha101<-c(0.2,0.3,0.6)
  #alpha110<-c(-1.5,0.9,0.5) 
  #sigma2<-5
  #tau2<-1
  
  #initial values for principal strata model
  lmm1<-summary(glmer(G==0~x1+x2+(1|cl),family=binomial,data=df[which(G!=1),]))
  beta<-as.numeric(lmm1$coefficients[,1])
  
  lmm2<-summary(glmer(G==1~x1+x2+(1|cl),family=binomial,data=df[which(G>0),]))
  gamma<-as.numeric(lmm2$coefficients[,1])
  
  #true model parameters
  #beta=c(-1,0.3,0.5) #group00 vs group11
  #gamma=c(-0.8,0.6,0.4) #group10 vs group11
  
  ##weakly informative priors
  #outcome model regression coefficient
  #a111<-a101<-a110<-rep(0,p) set all normal prior means of regression parameters as 0
  Sigma111<-Sigma101<-Sigma110<-diag(p)*1000 # normal variance of regession parameter
  a=0.001  #prior variance parameters for IG
  b=0.001
  cc=0.001
  dd=0.001
  
  #principal strata model regression coefficients
  beta0=gamma0=rep(0,p) #normal prior with mean zero
  Lambda=Gamma=diag(p)*1000
  
  g=0.001  #prior variance parameter
  h=0.001
  
  ##store chains for parameters/estimates with length S
  ALPHA111<-ALPHA101<-ALPHA110<-matrix(NA,p,S)
  SIGMA<-TAU<-TE<-rep(0,S) #sample TE by estimating the potential outcomes of G=2.
  
  BETA<-GAMMA<-matrix(NA,p,S)
  GV<-matrix(NA,3,S)
  ACC<-rep(0,S)
  
  #GT<-matrix(NA, nrow(df), S)   #store G trajectory
  #rownames(GT)<-df$GP
  
  for(i in 1:S){
    
    #update alpha
    tgt<-intersect(which(D==1),which(G==2)) 
    v=solve(crossprod(X[tgt,])/sigma2+solve(Sigma111))
    m=v%*%(t(X[tgt,])%*%(Y[tgt]-etan[tgt])/sigma2)
    alpha111=rmvnorm(1,m,v)
    ALPHA111[,i]<-alpha111
    
    tgt<-intersect(which(D==1),which(G==1)) 
    v=solve(crossprod(X[tgt,])/sigma2+solve(Sigma101))
    m=v%*%(t(X[tgt,])%*%(Y[tgt]-etan[tgt])/sigma2)
    alpha101=rmvnorm(1,m,v)
    ALPHA101[,i]<-alpha101
    
    tgt<-intersect(which(D==0),which(G==2)) 
    v=solve(crossprod(X[tgt,])/sigma2+solve(Sigma110))
    m=v%*%(t(X[tgt,])%*%(Y[tgt]-etan[tgt])/sigma2)
    alpha110=rmvnorm(1,m,v)
    ALPHA110[,i]<-alpha110
    
    #update eta
    for(j in 1:s){
      if (D[which(cl==j)[1]]==0){
        tgt3<-intersect(which(G==2),intersect(which(D==0),which(cl==j)))
        v=(length(tgt3)/sigma2+1/tau2)^(-1)
        m=v*(sum(Y[tgt3]-X[tgt3,]%*%t(alpha110))/sigma2)
        
      }else{
        tgt1<-intersect(which(G==2),intersect(which(D==1),which(cl==j)))
        tgt2<-intersect(which(G==1),intersect(which(D==1),which(cl==j)))
        v=(length(c(tgt1,tgt2))/sigma2+1/tau2)^(-1)
        m=v*(sum(Y[tgt1]-X[tgt1,]%*%t(alpha111))+sum(Y[tgt2]-X[tgt2,]%*%t(alpha101)))/sigma2
      }
      eta[j]<-rnorm(1,m,sqrt(v))
    }
    etan<-rep(eta,cs)
    
    #update tau2
    tau2<-rgamma(1,shape=(a+s/2),rate=(b+sum(eta^2)/2))^(-1)
    TAU[i]<-tau2
    
    #update sigma2
    tgt1<-intersect(which(G==2),which(D==1))
    tgt2<-intersect(which(G==1),which(D==1))
    tgt3<-intersect(which(G==2),which(D==0))
    rate<-dd+0.5*(sum((Y[tgt1]-X[tgt1,]%*%t(alpha111)-etan[tgt1])^2))+
      0.5*(sum((Y[tgt2]-X[tgt2,]%*%t(alpha101)-etan[tgt2])^2))+
      0.5*(sum((Y[tgt3]-X[tgt3,]%*%t(alpha110)-etan[tgt3])^2))
    
    sigma2<-rgamma(1,shape=(cc+Nst/2),rate=rate)^(-1)
    SIGMA[i]<-sigma2
    
    #potential outcomes
    dff<-X[c(tgt1,tgt3),]
    #X1<-as.matrix(cbind(rep(1,nrow(df1)),df1$x1,df1$x2))
    
    Y1<-dff%*%t(alpha111)
    Y0<-dff%*%t(alpha110)
    te<-mean(Y1-Y0)
    TE[i]<-te 
    
    ######################################
    ####MH 
    
    #update beta/gamma
    beta0=gamma0=rep(0,p)
    Lambda=Gamma=diag(p)*1000
    burn=S/5
    
    #mu0<-rep(0,length(beta))			# Prior mean of gamma (for all classes)	
    #S0<-s*diag(length(beta)) 			# Assume diffuse prior variance
    
    covkbeta<-covkgamma<-diag(p)			# Proposal covariance (updated after 2*burn)
    if (i>(2*burn)){
      covkbeta<-cov(t(BETA[1:p,(burn+1):(2*burn)]))
      covkgamma<-cov(t(GAMMA[1:p,(burn+1):(2*burn)]))
    } 
    
    gnew_beta<-beta + rmvt(1,sigma=dau1*covkbeta,11)  	# Draw from symmetric MV t-dist   #dau is the tuning parameter 
    gnew_gamma<-gamma + rmvt(1,sigma=dau1*covkgamma,11)  	# Draw from symmetric MV t-dist
    
    lpold_beta<-dmvnorm(beta,beta0,Lambda,log=T)		# Old Log Prior Beta
    lpnew_beta<-dmvnorm(gnew_beta,beta0,Lambda,log=T)	# New Log Prior
    
    lpold_gamma<-dmvnorm(gamma,gamma0,Lambda,log=T)		# Old Log Prior Gamma
    lpnew_gamma<-dmvnorm(gnew_gamma,beta0,Lambda,log=T)	# New Log Prior
    
    # Old and new class probabilities (from G-Logit Multinomial)
    etaold_00<-X%*%beta
    etaold_10<-X%*%gamma
    
    pold_00<-exp(etaold_00)/(1+exp(etaold_00)+exp(etaold_10))
    pold_10<-exp(etaold_10)/(1+exp(etaold_00)+exp(etaold_10))
    pold_11<-1/(1+exp(etaold_00)+exp(etaold_10))
    
    etanew_00<-X%*%t(gnew_beta)
    etanew_10<-X%*%t(gnew_gamma)
    
    pnew_00<-exp(etanew_00)/(1+exp(etanew_00)+exp(etanew_10))
    pnew_10<-exp(etanew_10)/(1+exp(etanew_00)+exp(etanew_10))
    pnew_11<-1/(1+exp(etanew_00)+exp(etanew_10))
    
    # Old and new multinomial likelihoods
    lold_00<-sum(log(pold_00)[which(G==0)])	# Class k loglike = sum[I(c_i=k)*log(p_ik)]
    lold_10<-sum(log(pold_10)[which(G==1)])
    lold_11<-sum(log(pold_11)[which(G==2)])
    
    lnew_00<-sum(log(pnew_00)[which(G==0)])
    lnew_10<-sum(log(pnew_10)[which(G==1)])
    lnew_11<-sum(log(pnew_11)[which(G==2)])
    
    # MH acceptance ratio on log scale
    ratio<-sum(lnew_00)+sum(lnew_10)+sum(lnew_11)+sum(lpnew_beta)+sum(lpnew_gamma)-(sum(lold_00)+sum(lold_10)+sum(lold_11)+sum(lpold_beta)+sum(lpold_gamma))
    
    if(log(runif(1))<ratio) {
      beta=c(gnew_beta)
      gamma=c(gnew_gamma)
      ACC[i]<-1
    }else{
      ACC[i]<-0
    }
    
    BETA[,i]<-beta
    GAMMA[,i]<-gamma
    
    #Update membership
    #membership for death in control
    tgt<-intersect(which(is.na(Y)),which(D==0))
    A=exp(X[tgt,]%*%beta)/(exp(X[tgt,]%*%beta)+exp(X[tgt,]%*%gamma)+1) #00
    B=exp(X[tgt,]%*%gamma)/(exp(X[tgt,]%*%beta)+exp(X[tgt,]%*%gamma)+1) #10
    
    GP<-rep(1,length(tgt))-rbinom(length(tgt),1,I(A/(A+B)))
    G[tgt]<-GP
    
    #membership for alive in treatment
    tgt<-intersect(which(!is.na(Y)),which(D==1))
    D1<-dnorm(as.numeric(Y[tgt]-etan[tgt]-X[tgt,]%*%t(alpha101))/sqrt(sigma2))
    D2<-dnorm(as.numeric(Y[tgt]-etan[tgt]-X[tgt,]%*%t(alpha111))/sqrt(sigma2))
    
    A=exp(X[tgt,]%*%gamma)/(exp(X[tgt,]%*%beta)+exp(X[tgt,]%*%gamma)+1)
    B=1/(exp(X[tgt,]%*%beta)+exp(X[tgt,]%*%gamma)+1)
    
    P=A*D1/(A*D1+B*D2)     
    GP<-rep(2,length(tgt))-rbinom(length(tgt),1,P)
    G[tgt]<-GP
    
    #store G percentage
    GV[,i]<-c(length(which(G==0)),length(which(G==1)),length(which(G==2)))/N
    if (i%%500==0) print(i) #print every 500 steps
  
  }  
  
  return(list(ALPHA111=ALPHA111,ALPHA101=ALPHA101,ALPHA110=ALPHA110,BETA=BETA,GAMMA=GAMMA,TAU=TAU,SIGMA=SIGMA,GV=GV,TE=TE,ACC=ACC))
}

###################################
## 3. Simulation Studies ##########
###################################

sim0<-function(S=10000,nsim=50, N=1500, s=60,seed=1234567){
  #p=3
  set.seed(seed)
  ALPHA111<-ALPHA101<-ALPHA110<-BETA<-GAMMA<-GV<-array(NA,c(nsim,3,S))
  TAU<-SIGMA<-TE<-ACC<-matrix(NA,nrow=nsim,ncol=S)
  
  for(i in 1:nsim){
    print(paste("replicate:",i))
    df<-gendata(N=N,s=s)
    
    otp<-try(uvsampler(df=df,S=S, dau1=0.05)) #avoid potential cases with extreme data/no intial value; not observed in simulations.
    if (class(otp) =="try-error"){
      next
    }
   
    ALPHA111[i,,]<-otp$ALPHA111
    ALPHA101[i,,]<-otp$ALPHA101
    ALPHA110[i,,]<-otp$ALPHA110
    BETA[i,,]<-otp$BETA
    GAMMA[i,,]<-otp$GAMMA
    
    TAU[i,]<-otp$TAU
    SIGMA[i,]<-otp$SIGMA
    TE[i,]<-otp$TE
    GV[i,,]<-otp$GV
    ACC[i,]<-otp$ACC
  }
  
  output<-list(ALPHA111=ALPHA111,ALPHA101=ALPHA101,ALPHA110=ALPHA110,BETA=BETA,GAMMA=GAMMA,
               GV=GV,TAU=TAU,SIGMA=SIGMA,TE=TE,ACC=ACC)
  
  return(output)
}



######################################################
## 4. Output summary and ICC estimation ##############
######################################################

## simulation result summary
sumstat<-function(sim1,afterburn=2501:10000){
  pestimate<-c(
    mean(rowMeans(sim1$ALPHA111[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA111[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA111[,3,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$ALPHA101[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA101[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA101[,3,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$ALPHA110[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA110[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA110[,3,afterburn]),na.rm=T),
    
    
    mean(rowMeans(sim1$BETA[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$BETA[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$BETA[,3,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$GAMMA[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$GAMMA[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$GAMMA[,3,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$TAU[, afterburn]),na.rm=T),
    mean(rowMeans(sim1$SIGMA[,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$TE[,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$GV[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$GV[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$GV[,3,afterburn]),na.rm=T))
  
  param<-c("alpha111","alpha111","alpha111",
           "alpha101","alpha101","alpha101",
           "alpha110","alpha110","alpha110",
           
           "beta","beta","beta",
           "gamma","gamma","gamma",
           "tau","sigma", "te",
           "G0","G1","G2")
  
  true_val<-c(1.5, 0.5, 0.8,
              0.2, 0.3, 0.6,
              -1.5, 0.9, 0.5,
              
              -1,0.3,0.5,
              -0.8,0.6,0.4,
              1, 5, 2.846862,
              
              0.2112513, 0.2648087, 0.5239400)
  
  rbias<-(pestimate-true_val)/true_val*100
  
  coverage<-c()
  a.true1=c(1.5, 0.5, 0.8)
  
  cla1=apply(sim1$ALPHA111[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$ALPHA111[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$ALPHA111[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  a1<-mean(a.true1[1]>=cla1[1,] & a.true1[1]<=cla1[2,],na.rm=T)
  a2<-mean(a.true1[2]>=cla2[1,] & a.true1[2]<=cla2[2,],na.rm=T)
  a3<-mean(a.true1[3]>=cla3[1,] & a.true1[3]<=cla3[2,],na.rm=T)
  
  coverage<-c(coverage,a1,a2,a3)
  
  a.true2=c(0.2, 0.3, 0.6)
  
  cla1=apply(sim1$ALPHA101[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$ALPHA101[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$ALPHA101[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  a1<-mean(a.true2[1]>=cla1[1,] & a.true2[1]<=cla1[2,],na.rm=T)
  a2<-mean(a.true2[2]>=cla2[1,] & a.true2[2]<=cla2[2,],na.rm=T)
  a3<-mean(a.true2[3]>=cla3[1,] & a.true2[3]<=cla3[2,],na.rm=T)
  
  coverage<-c(coverage,a1,a2,a3)
  
  a.true3=c(-1.5, 0.9, 0.5)
  
  cla1=apply(sim1$ALPHA110[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$ALPHA110[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$ALPHA110[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  a1<-mean(a.true3[1]>=cla1[1,] & a.true3[1]<=cla1[2,],na.rm=T)
  a2<-mean(a.true3[2]>=cla2[1,] & a.true3[2]<=cla2[2,],na.rm=T)
  a3<-mean(a.true3[3]>=cla3[1,] & a.true3[3]<=cla3[2,],na.rm=T)
  
  coverage<-c(coverage,a1,a2,a3)
  
  #beta
  b.true=c(-1,0.3,0.5)
  cla1=apply(sim1$BETA[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$BETA[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$BETA[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  b1<-mean(b.true[1]>=cla1[1,] & b.true[1]<=cla1[2,],na.rm=T)
  b2<-mean(b.true[2]>=cla2[1,] & b.true[2]<=cla2[2,],na.rm=T)
  b3<-mean(b.true[3]>=cla3[1,] & b.true[3]<=cla3[2,],na.rm=T)
  
  coverage<-c(coverage,b1,b2,b3)
  
  #gamma
  g.true=c(-0.8,0.6,0.4)
  cla1=apply(sim1$GAMMA[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$GAMMA[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$GAMMA[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  g1<-mean(g.true[1]>=cla1[1,] & g.true[1]<=cla1[2,],na.rm=T)
  g2<-mean(g.true[2]>=cla2[1,] & g.true[2]<=cla2[2,],na.rm=T)
  g3<-mean(g.true[3]>=cla3[1,] & g.true[3]<=cla3[2,],na.rm=T)
  
  coverage<-c(coverage,g1,g2,g3)
  
  #tau.true=c(1,0.3*sqrt(2),2)
  
  tau.true=1
  cla1=apply(sim1$TAU[,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  g1<-mean(tau.true[1]>=cla1[1,] & tau.true[1]<=cla1[2,],na.rm=T)
  
  coverage<-c(coverage,g1) 
  
  #sigma.true=c(5,0.3*sqrt(50),10)
  
  sigma.true=5
  
  cla1=apply(sim1$SIGMA[,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  g1<-mean(sigma.true[1]>=cla1[1,] & sigma.true[1]<=cla1[2,],na.rm=T)
  
  coverage<-c(coverage,g1)   
  
  #te.true=c(2.846862,3.642532)
  te.true=2.846862
  
  cla1=apply(sim1$TE[,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  g1<-mean(te.true[1]>=cla1[1,] & te.true[1]<=cla1[2,],na.rm=T)
  coverage<-c(coverage,g1) 
  
  g.true=c(0.2112513, 0.2648087, 0.5239400)
  
  cla1=apply(sim1$GV[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$GV[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$GV[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  g1<-mean(g.true[1]>=cla1[1,] & g.true[1]<=cla1[2,],na.rm=T)
  g2<-mean(g.true[2]>=cla2[1,] & g.true[2]<=cla2[2,],na.rm=T)
  g3<-mean(g.true[3]>=cla3[1,] & g.true[3]<=cla3[2,],na.rm=T)
  
  coverage<-c(coverage,g1,g2,g3) 
  
  return(data.frame(param,true_val,pestimate,rbias,coverage))
}


## ICC estimation
icc_calculation<-function(sim1,afterburn=2501:10000){
  
  RHO<-sim1$TAU[, afterburn]/(sim1$TAU[, afterburn]+sim1$SIGMA[, afterburn])
  
  #mean ICC
  icc<-c(rho<-mean(rowMeans(RHO),na.rm=T))
  
  #coverage
  val.true=c(1/6)
  rbias<-(icc-val.true)/val.true*100
  cla1=apply(RHO,1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  g1<-mean(val.true[1]>=cla1[1,] & val.true[1]<=cla1[2,],na.rm=T)
  
  coverage<-g1
  
  return(cbind( val.true,icc, rbias, coverage))
  
}

####################################################################
## 5. Example simulation execution and result summary ##############
####################################################################
output<-sim0(S=5000,nsim=50, N=1500, s=60)
sumstat(output,afterburn=2501:5000)
icc_calculation(output,afterburn=2501:5000)


