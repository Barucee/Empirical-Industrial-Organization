library(tidyverse)
library(mvtnorm)
library(nloptr)
library(polyMatrix)
library(rootSolve)
library(stargazer)
library(MASS)
library(kableExtra)
library(descr)
library(xtable)
##EIO project

  ######    PART 1  ######




##Creation a dummy that flags if carrier is legacy or not
table(Data_US.2017_2$AIRLINE)
df1<-Data_US.2017_2%>%mutate(legacy=case_when(AIRLINE=="AA"~"1",
                                              AIRLINE=="UA"~"1",
                                              AIRLINE=="DL"~"1",
                                              AIRLINE=="UA"~"1",
                                              AIRLINE=="AS"~"0",
                                              AIRLINE=="B6"~"0",
                                              AIRLINE=="F9"~"0",
                                              AIRLINE=="G4"~"0",
                                              AIRLINE=="NK"~"0",
                                              AIRLINE=="OO"~"0",
                                              AIRLINE=="SY"~"0",
                                              AIRLINE=="VX"~"0",
                                              AIRLINE=="WN"~"0",
                                              ))

df1$legacy<-as.numeric(df1$legacy)
df1<-df1%>%mutate(LCC=case_when(legacy==1~0,
                                legacy==0~1))

    #Creating the controls


    ###Airport presence

#Get count of how many times each airline served the markets present in our data
df2<-df1%>%group_by(AERODEP,AIRLINE, AEROARR )%>%summarise(count=n())
df3<-df2%>%group_by(AERODEP, AIRLINE)%>%summarise(count=n())
#Total count of markets served by each airport
df4<-df1%>%group_by(AERODEP, AEROARR)%>%summarise(count=n())

df4<-df4%>%group_by(AERODEP)%>%summarise(count=n())

#Final presence of each airline at each airport
test1<-left_join(df3, df4, by="AERODEP")
test1<-test1%>%mutate(AERO_PRESENCE=count.x/count.y)

    ##Creating final database

final<-df1%>%group_by(AERODEP,AEROARR, MARKET)%>%summarise(N_leg=sum(legacy), N_LCC=sum(LCC), dist=mean(Distance), pop_dep=mean(POPDEP), pop_arr=mean(POPARR), pax=mean(PAX), mdep=mean(MAINDEP), marr=mean(MAINARR), marketsize=mean(SQRT.PRODPOP.))

    ##Adding two dummy variables on if there are a strictly positive number of legacy/LCC carriers in the market
final<-final%>%mutate(d_leg=ifelse(N_leg>0,1,0),
                      d_LCC=ifelse(N_LCC>0,1,0))
final$pop_mean2<-(log(final$pop_dep/1E6)+log(final$pop_arr/1E6))/2

#Frequency table
table(final$N_leg, final$N_LCC)%>%kable(format="latex", caption="Counts of market structures")

  ######  PART 2 ######

#nloptr only using market size and not including pop mean anymore

n=nrow(final)


###### Expl variables #########
## name the data  ##
distance=final$dist/100
dist2=distance^2
market_size=log(final$marketsize)


#### N1, N2 #############
N_FSC=as.matrix(final$N_leg)
N_LCC=as.matrix(final$N_LCC)

##Indicators that N_FSC>0 and N_LCC>0
D_FSC=as.matrix(final$d_leg)
D_LCC=as.matrix(final$d_LCC)
####Matrix
Xm=as.matrix(cbind(rep(1,n),distance,dist2, market_size))

########## profit functions ##########
#theta=c(theta1,delta_FF,delta_FS,delta_FL,thetaS,delta_SF,delta_SL,theta3,delta_LF,delta_LS,delta_LL)
d1=10
d2=11


######
Krho=0.999
#####

# Profit functions: FSC
func_f <- function(Xm,D_FSC, Nb_FSC, Nb_LCC, para1){
  beta_FSC <- para1[1:4]
  delta_FF <- para1[5]
  delta_FL <- para1[6]
  f <- Xm %*% as.matrix(beta_FSC) - delta_FF * (log(Nb_FSC)*D_FSC) - delta_FL * (log(Nb_LCC + 1))/Nb_FSC
  return(-f)
}




# Profit function: LCC
func_g <- function(Xm,D_LCC, Nb_FSC, Nb_LCC, para2){
  beta_LCC <- para2[1:4]
  delta_LF <- para2[5]
  delta_LL <- para2[6]
  g <- Xm %*% as.matrix(beta_LCC)  - delta_LF * (log(Nb_FSC + 1))/Nb_LCC - delta_LL * (log(Nb_LCC)*D_LCC)
  return(-g)
}



#### Log likelihood of obs number i
LL_one2 <- function(theta,i){
  rho <- theta[13]
  corr <- diag(2)
  corr[lower.tri(corr)] <- 0.99
  corr[upper.tri(corr)] <- 0.99
  
  if(N_FSC[i] > 0 & N_LCC[i] > 0){
    prob <- pmvnorm(lower = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i], theta[1:6]), 
                              func_g(Xm[i,],D_LCC[i], Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i], theta[7:14])),
                    upper = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i], theta[1:6]), 
                              func_g(Xm[i,],D_LCC[i], Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i]+1, theta[7:12])), mean = c(0,0), corr) -
      pmvnorm(lower = c(func_f(Xm[i,],D_FSC[i],  Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i]-1,  theta[1:6]), 
                        func_g(Xm[i,],D_LCC[i],  Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i], theta[7:12])), 
              upper = c(func_f(Xm[i,],D_FSC[i],  Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i], theta[1:6]), 
                        func_g(Xm[i,],D_LCC[i],  Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i], theta[7:12])), mean = c(0,0), corr)
    if ( prob <= 0| is.na(prob)==T) {
      prob <- 1e-40
    }
    ln_prob <- log(prob)
    
  }else if(N_FSC[i] > 0 & N_LCC[i] == 0){
    prob <- pmvnorm(lower = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i], Nb_LCC = 0, theta[1:6]),-Inf),
                    upper = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i]+1, Nb_LCC = 0, theta[1:6]), 
                              func_g(Xm[i,],D_LCC[i], Nb_FSC = N_FSC[i], Nb_LCC = 1, theta[7:12])), mean = c(0,0), corr) 
    if (prob <= 0| is.na(prob)==T) {
      prob <- 1e-40
    }
    ln_prob <- log(prob)
    
  }else if(N_FSC[i] == 0 & N_LCC[i] >0){
    prob <- pmvnorm(lower = c(-Inf, 
                              func_g(Xm[i,],D_LCC[i], Nb_FSC = 0, Nb_LCC = N_LCC[i], theta[7:12])),
                    upper = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = 1, Nb_LCC = N_LCC[i], theta[1:6]), 
                              func_g(Xm[i,],D_LCC[i], Nb_FSC = 0, Nb_LCC = N_LCC[i]+1, theta[7:12])), mean = c(0,0), corr) -
      pmvnorm(lower = c(func_f(Xm[i,],D_FSC[i],  Nb_FSC = 1, Nb_LCC = N_LCC[i]-1,  theta[1:6]), 
                        func_g(Xm[i,],D_LCC[i],  Nb_FSC = 0, Nb_LCC = N_LCC[i], theta[7:12])), 
              upper = c(func_f(Xm[i,],D_FSC[i],  Nb_FSC = 1, Nb_LCC = N_LCC[i], theta[1:6]), 
                        func_g(Xm[i,],D_LCC[i],  Nb_FSC = 1, Nb_LCC = N_LCC[i], theta[7:12])), mean = c(0,0), corr)
    if( prob <= 0 | is.na(prob)==T){
      prob <- 1e-40
    }
    ln_prob <- log(prob)
    
  }else{
    prob <- pmvnorm(lower = c(-Inf ,-Inf),upper = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = 1, Nb_LCC = 0, theta[1:6]), 
                                                    func_g(Xm[i,],D_LCC[i], Nb_FSC = 0, Nb_LCC = 1, theta[7:12])), mean = c(0,0), corr) 
    
    prob <- pmvnorm(lower = c(func_f(Xm[i,],D_FSC[i],  Nb_FSC = 1, Nb_LCC = 0,  theta[1:6]), 
                              func_g(Xm[i,],D_LCC[i],  Nb_FSC = 0, Nb_LCC = 1,  theta[7:12])), 
                    upper = rep(Inf, 2), mean = c(0,0), corr)
    if ( prob <= 0| is.na(prob)==T) {
      prob <- 1e-40
    }
    ln_prob <- log(prob)
  }
  (ln_prob)
}


neg_log_likelihood2 <- function(theta) {
  
  neg_ll2 <- 0
  
  
  for (i in 1:3983) {
    
    neg_ll2 <- neg_ll2 - LL_one2(theta,i)
  }
  
  
  return(neg_ll2)
}


#Estimating theta_init3
p8<-polr(factor(N_leg)~dist +dist_2+ log(market_size), data=final, method="probit")
p9<-polr(factor(N_LCC)~dist + dist_2+ log(market_size), data=final, method="probit")


##
theta_init3= as.matrix(c(2.5,0.00079,-0.000000055,-0.829,0.3,0.15,-5,0.00213,-0.000000521,-0.328,0.07,0.3,2))

##
#Creating inequality constraints for parameters

v<-function(theta){
  g1<-theta[6]-theta[5]
  g2<-theta[11] - theta[12]
  g3<- -theta[6]
  g4<- -theta[5]
  g5<- -theta[11]
  g6<- -theta[12]
  return(c(g1,g2, g3, g4, g5, g6))
}

#Optimization

opts <-list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-4,maxeval=10000,"print_level"=2)
res3<-nloptr(x0=theta_init3,eval_f=neg_log_likelihood2, eval_g_ineq = v,opts=opts)
print(res3)
res3$solution


#Computing the standard errors

  #Deriving hessian of log likelihood funcion
h<-hessian(neg_log_likelihood2, res3$solution)

  #Taking away the correlation parameter, and computing inverse of hessian to obtain the std errors
h2<-h[-13,-13]
h3<-solve(h2)
diag(h3)

#Building stargazer table for results

  #fake regressions
d <- data.frame(y = runif(30),
                x = runif(30),
                a = runif(30),
                b = runif(30),
                c = runif(30),
                e = runif(30),
                f = runif(30),
                g = runif(30),
                h = runif(30),
                i = runif(30),
                j = runif(30),
                k = runif(30),
                l = runif(30),
                m = runif(30),
                n = runif(30),
                o = runif(30),
                p = runif(30),
                q = runif(30),
                r = runif(30),
                s = runif(30),
                t = runif(30),
                u = runif(30),
                v = runif(30),
                w = runif(30))
j<-lm(y ~ a+b+c+e+f+g+h+i+j+k+l+m-1, data=d)
i<-lm(x ~ a+b+c+e+f+g+h+i+j+k+l+m-1, data=d)
o<-lm(y ~ a+b+c+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t+u+v+w-1, data=d)
p<-lm(x ~ a+b+c+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t+u+v+w-1, data=d)


  #General results
stargazer(j,covariate.labels=c("cst_{Leg}", "distance_{Leg}","distance^{2}_{Leg}","marketsize_{Leg}","delta_{Leg, Leg}","delta_{Leg, LCC}",
                             "cst_{LCC}", "distance_{LCC}","distance^{2}_{LCC}","marketsize_{LCC}","delta_{LCC, Leg}","delta_{LCC, LCC}"),
          dep.var.labels = "",
          coef=list(res3$solution[-13]),
          se=list(diag(h3)), type="latex", title="NL Optimization estimates",
          omit.stat = "all"
          )

#Different competitive effects depending on number of entrants
  
  #Intra FSC
Comp_eff_INTR_FSC<-function(x){
  m<-c()
  for (i in 1:x){

    y<-res3$solution[5]*log(i)
    m<-c(m,y)
 
}
  return(m)
}
Delta_FF<-Comp_eff_INTR_FSC(6)
  #Inter FSC
Comp_eff_INTE_FSC<-function(x){
  m<-c()
  for (i in 1:x){
    
    y<-res3$solution[6]*log(i+1)
    m<-c(m,y)
    
  }
  return(m)
}
Delta_Fl<-Comp_eff_INTE_FSC(6)

  #Intra LCC
Comp_eff_INTR_LCC<-function(x){
  m<-c()
  for (i in 1:x){
    
    y<-res3$solution[12]*log(i)
    m<-c(m,y)
    
  }
  return(m)
}
Delta_LL<-Comp_eff_INTR_LCC(6)
  #Inter LCC
Comp_eff_INTE_LCC<-function(x){
  m<-c()
  for (i in 1:x){
    
    y<-res3$solution[11]*log(i+1)
    m<-c(m,y)
    
  }
  return(m)
}
Delta_LF<-Comp_eff_INTE_LCC(6)

#Combinig parameters for table
Eff_FSC<-round(as.vector(c(Delta_FF, Delta_Fl)),3)
Eff_LCC<-round(as.vector(c(Delta_LL,Delta_LF)),3)

#Table
stargazer(j,i,covariate.labels=c("delta^{1}_{f,f}", "delta^{2}_{f,f}","delta^{3}_{f,f}","delta^{4}_{f,f}","delta^{5}_{f,f}","delta^{6}_{f,f}",
                                 "delta^{1}_{f,-f}", "delta^{2}_{f,-f}","delta^{3}_{f,-f}","delta^{4}_{f,-f}","delta^{5}_{f,-f}","delta^{6}_{f,-f}"),
          dep.var.labels = c("FSC", "LCC"),
          coef=c(list(Eff_FSC), list(Eff_LCC)),
          se=c(list(rep(NA,12)), list(rep(NA,12))),
          omit.stat = "all",
          type="latex", title="Competitive effects")




  ######  PART 3 ######

#####Estimating the ratios 

  #Intra_FSC
INTRA_FCC <- function(x) {
  m <- c()
  for(i in 1:x) {
    y <- exp( ( res3$solution[5] *( log(i) -log(i-1) ) /res3$solution[1] ) ) * ((i-1)/i)
    m <- c(m, y)
  }
  return(m)
}

#Values
INTRA_FSC<-INTRA_FCC(11)
INTRA_FSC


  #Intra_LCC

INTRA_LCC <- function(x) {
  m <- c()
  for(i in 1:x) {
    y <- exp( ( res3$solution[12] *( log(i) -log(i-1) ) /res3$solution[7] ) ) * ((i-1)/i)
    m <- c(m, y)
  }
  return(m)
}

#Values
INTRA_LCC<-INTRA_LCC(11)
INTRA_LCC

  #Inter_FSC

INTER_FCC <- function(x) {
  m <- c()
  for(i in 1:x) {
    y <- exp( ( res3$solution[6] *( log(i) -log(i-1) ) /res3$solution[1] ) )
    m <- c(m, y)
  }
  return(m)
}

#Values
INTER_FSC<-INTER_FCC(11)
INTER_FSC
  #Inter_LCC

INTER_LCC <- function(x) {
  m <- c()
  for(i in 1:x) {
    y <- exp( ( res3$solution[11] *( log(i) -log(i-1) ) /res3$solution[7] ) )
    m <- c(m, y)
  }
  return(m)
}

#Values
INTER_LCC<-INTER_LCC(11)
INTER_LCC
#Making the table
INT_LCC<-c(INTRA_LCC, INTER_LCC)
INT_FSC<-c(INTRA_FSC, INTER_FSC)
stargazer(o,p, covariate.labels=c("INTRA^{1}", "INTRA^{2}","INTRA^{3}","INTRA^{4}","INTRA^{5}","INTRA^{6}","INTRA^{7}","INTRA^{8}","INTRA^{9}","INTRA^{10}","INTRA^{11}",
                                  "INTER^{1}", "INTER^{2}","INTER^{3}","INTER^{4}","INTER^{5}","INTER^{6}","INTER^{7}","INTER^{8}","INTER^{9}","INTER^{10}","INTER^{11}"),
          dep.var.labels = c("FSC", "LCC"),
          coef=c(list(INT_FSC), list(INT_LCC)),
          se=c(list(rep(NA,12)), list(rep(NA,12))),
          omit.stat = "all",
          type="latex", title="Competitive effects")


#Final nloptr using population mean rather than market size (counterfactual)

n=nrow(final)


###### Expl variables #########
## name the data  ##
distance=final$dist/100
dist2=distance^2
population_mean=(log(final$pop_dep/1E6)+log(final$pop_arr/1E6))/2


#### N1, N2 #############
N_FSC=as.matrix(final$N_leg)
N_LCC=as.matrix(final$N_LCC)

##Indicators that N_FSC>0 and N_LCC>0
D_FSC=as.matrix(final$d_leg)
D_LCC=as.matrix(final$d_LCC)
####Matrix
Xm=as.matrix(cbind(rep(1,n),distance,dist2, population_mean))

########## profit functions ##########
#theta=c(theta1,delta_FF,delta_FS,delta_FL,thetaS,delta_SF,delta_SL,theta3,delta_LF,delta_LS,delta_LL)
d1=10
d2=11


######
Krho=0.999
#####

# Profit functions: FSC
func_f <- function(Xm,D_FSC, Nb_FSC, Nb_LCC, para1){
  beta_FSC <- para1[1:4]
  delta_FF <- para1[5]
  delta_FL <- para1[6]
  f <- Xm %*% as.matrix(beta_FSC) - delta_FF * (log(Nb_FSC)*D_FSC) - delta_FL * (log(Nb_LCC + 1))/Nb_FSC
  return(-f)
}




# Profit function: LCC
func_g <- function(Xm,D_LCC, Nb_FSC, Nb_LCC, para2){
  beta_LCC <- para2[1:4]
  delta_LF <- para2[5]
  delta_LL <- para2[6]
  g <- Xm %*% as.matrix(beta_LCC)  - delta_LF * (log(Nb_FSC + 1))/Nb_LCC - delta_LL * (log(Nb_LCC)*D_LCC)
  return(-g)
}



#### Log likelihood of obs number i
LL_one3 <- function(theta,i){
  rho <- theta[13]
  corr <- diag(2)
  corr[lower.tri(corr)] <- 0.99
  corr[upper.tri(corr)] <- 0.99
  
  if(N_FSC[i] > 0 & N_LCC[i] > 0){
    prob <- pmvnorm(lower = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i], theta[1:6]), 
                              func_g(Xm[i,],D_LCC[i], Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i], theta[7:14])),
                    upper = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i], theta[1:6]), 
                              func_g(Xm[i,],D_LCC[i], Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i]+1, theta[7:12])), mean = c(0,0), corr) -
      pmvnorm(lower = c(func_f(Xm[i,],D_FSC[i],  Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i]-1,  theta[1:6]), 
                        func_g(Xm[i,],D_LCC[i],  Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i], theta[7:12])), 
              upper = c(func_f(Xm[i,],D_FSC[i],  Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i], theta[1:6]), 
                        func_g(Xm[i,],D_LCC[i],  Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i], theta[7:12])), mean = c(0,0), corr)
    if ( prob <= 0| is.na(prob)==T) {
      prob <- 1e-40
    }
    ln_prob <- log(prob)
    
  }else if(N_FSC[i] > 0 & N_LCC[i] == 0){
    prob <- pmvnorm(lower = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i], Nb_LCC = 0, theta[1:6]),-Inf),
                    upper = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i]+1, Nb_LCC = 0, theta[1:6]), 
                              func_g(Xm[i,],D_LCC[i], Nb_FSC = N_FSC[i], Nb_LCC = 1, theta[7:12])), mean = c(0,0), corr) 
    if (prob <= 0| is.na(prob)==T) {
      prob <- 1e-40
    }
    ln_prob <- log(prob)
    
  }else if(N_FSC[i] == 0 & N_LCC[i] >0){
    prob <- pmvnorm(lower = c(-Inf, 
                              func_g(Xm[i,],D_LCC[i], Nb_FSC = 0, Nb_LCC = N_LCC[i], theta[7:12])),
                    upper = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = 1, Nb_LCC = N_LCC[i], theta[1:6]), 
                              func_g(Xm[i,],D_LCC[i], Nb_FSC = 0, Nb_LCC = N_LCC[i]+1, theta[7:12])), mean = c(0,0), corr) -
      pmvnorm(lower = c(func_f(Xm[i,],D_FSC[i],  Nb_FSC = 1, Nb_LCC = N_LCC[i]-1,  theta[1:6]), 
                        func_g(Xm[i,],D_LCC[i],  Nb_FSC = 0, Nb_LCC = N_LCC[i], theta[7:12])), 
              upper = c(func_f(Xm[i,],D_FSC[i],  Nb_FSC = 1, Nb_LCC = N_LCC[i], theta[1:6]), 
                        func_g(Xm[i,],D_LCC[i],  Nb_FSC = 1, Nb_LCC = N_LCC[i], theta[7:12])), mean = c(0,0), corr)
    if( prob <= 0 | is.na(prob)==T){
      prob <- 1e-40
    }
    ln_prob <- log(prob)
    
  }else{
    prob <- pmvnorm(lower = c(-Inf ,-Inf),upper = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = 1, Nb_LCC = 0, theta[1:6]), 
                                                    func_g(Xm[i,],D_LCC[i], Nb_FSC = 0, Nb_LCC = 1, theta[7:12])), mean = c(0,0), corr) 
    
    prob <- pmvnorm(lower = c(func_f(Xm[i,],D_FSC[i],  Nb_FSC = 1, Nb_LCC = 0,  theta[1:6]), 
                              func_g(Xm[i,],D_LCC[i],  Nb_FSC = 0, Nb_LCC = 1,  theta[7:12])), 
                    upper = rep(Inf, 2), mean = c(0,0), corr)
    if ( prob <= 0| is.na(prob)==T) {
      prob <- 1e-40
    }
    ln_prob <- log(prob)
  }
  (ln_prob)
}


neg_log_likelihood3 <- function(theta) {
  
  neg_ll2 <- 0
  
  
  for (i in 1:3983) {
    
    neg_ll2 <- neg_ll2 - LL_one3(theta,i)
  }
  
  
  return(neg_ll2)
}


#Estimating theta_init4
final$pop_mean<-(log(final$pop_dep/1E6)+log(final$pop_arr/1E6))/2
p8<-polr(factor(N_leg)~dist +dist_2+ pop_mean, data=final, method="probit")
p9<-polr(factor(N_LCC)~dist + dist_2+ pop_mean, data=final, method="probit")



##
theta_init4= as.matrix(c(2.5,0.00079,-0.000000055,-1,0.3,0.15,-5,0.00213,-0.000000521,-0.51,0.07,0.3,2))
theta_init4= as.matrix(c(2.5,0.00079,-0.000000055,-1,0.3,0.15,-5,0.00213,-0.000000521,-0.51,0.07,0.3,2))

##
#Creating inequality constraints for parameters

v<-function(theta){
  g1<-theta[6]-theta[5]
  g2<-theta[11] - theta[12]
  g3<- -theta[6]
  g4<- -theta[5]
  g5<- -theta[11]
  g6<- -theta[12]
  return(c(g1,g2, g3, g4, g5, g6))
}

#Optimization

opts <-list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-4,maxeval=10000,"print_level"=2)
res4<-nloptr(x0=theta_init4,eval_f=neg_log_likelihood2, eval_g_ineq = v,opts=opts)
#Hessian
hess<-hessian(neg_log_likelihood3,res4$solution)
hess2<-hess[-13,-13]
hess3<-solve(hess2)
#Table
stargazer(j,covariate.labels=c("cst_{Leg}", "distance_{Leg}","distance^{2}_{Leg}","popmean_{Leg}","delta_{Leg, Leg}","delta_{Leg, LCC}",
                               "cst_{LCC}", "distance_{LCC}","distance^{2}_{LCC}","popmean_{LCC}","delta_{LCC, Leg}","delta_{LCC, LCC}"),
          dep.var.labels = "",
          coef=list(res4$solution[-13]),
          se=list(diag(hess3)), type="latex", title="NL Optimization estimates",
          omit.stat = "all"
)


#Changing order of entry
n=nrow(final)


###### Expl variables #########
## name the data  ##
distance=final$dist/100
dist2=distance^2
market_size=log(final$marketsize)


#### N1, N2 #############
N_FSC=as.matrix(final$N_leg)
N_LCC=as.matrix(final$N_LCC)

##Indicators that N_FSC>0 and N_LCC>0
D_FSC=as.matrix(final$d_leg)
D_LCC=as.matrix(final$d_LCC)
####Matrix
Xm=as.matrix(cbind(rep(1,n),distance,dist2, market_size))

########## profit functions ##########
#theta=c(theta1,delta_FF,delta_FS,delta_FL,thetaS,delta_SF,delta_SL,theta3,delta_LF,delta_LS,delta_LL)
d1=10
d2=11


######
Krho=0.999
#####

# Profit functions: FSC
func_f <- function(Xm,D_FSC, Nb_FSC, Nb_LCC, para1){
  beta_FSC <- para1[1:4]
  delta_FF <- para1[5]
  delta_FL <- para1[6]
  f <- Xm %*% as.matrix(beta_FSC) - delta_FF * (log(Nb_FSC)*D_FSC) - delta_FL * (log(Nb_LCC + 1))/Nb_FSC
  return(-f)
}




# Profit function: LCC
func_g <- function(Xm,D_LCC, Nb_FSC, Nb_LCC, para2){
  beta_LCC <- para2[1:4]
  delta_LF <- para2[5]
  delta_LL <- para2[6]
  g <- Xm %*% as.matrix(beta_LCC)  - delta_LF * (log(Nb_FSC + 1))/Nb_LCC - delta_LL * (log(Nb_LCC)*D_LCC)
  return(-g)
}



#### Log likelihood of obs number i
LL_one4 <- function(theta,i){
  rho <- theta[13]
  corr <- diag(2)
  corr[lower.tri(corr)] <- 0.99
  corr[upper.tri(corr)] <- 0.99
  
  if(N_FSC[i] > 0 & N_LCC[i] > 0){
    prob <- pmvnorm(lower = c(func_g(Xm[i,],D_LCC[i], Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i], theta[7:14]),
                              func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i], theta[1:6]) 
                              ),
                    upper = c(func_g(Xm[i,],D_LCC[i], Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i]+1, theta[7:12]),
                              func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i], theta[1:6]) 
                              ), mean = c(0,0), corr) -
      pmvnorm(lower = c(func_g(Xm[i,],D_LCC[i],  Nb_FSC = N_FSC[i], Nb_LCC = N_LCC[i], theta[7:12]),
                        func_f(Xm[i,],D_FSC[i],  Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i]-1,  theta[1:6]) 
                        ), 
              upper = c(func_g(Xm[i,],D_LCC[i],  Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i], theta[7:12]),
                        func_f(Xm[i,],D_FSC[i],  Nb_FSC = N_FSC[i]+1, Nb_LCC = N_LCC[i], theta[1:6]) 
                        ), mean = c(0,0), corr)
    if ( prob <= 0| is.na(prob)==T) {
      prob <- 1e-40
    }
    ln_prob <- log(prob)
    
  }else if(N_FSC[i] > 0 & N_LCC[i] == 0){
    prob <- pmvnorm(lower = c(func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i], Nb_LCC = 0, theta[1:6]),-Inf),
                    upper = c(func_g(Xm[i,],D_LCC[i], Nb_FSC = N_FSC[i], Nb_LCC = 1, theta[7:12]),
                              func_f(Xm[i,],D_FSC[i], Nb_FSC = N_FSC[i]+1, Nb_LCC = 0, theta[1:6])
                              ), mean = c(0,0), corr) 
    if (prob <= 0| is.na(prob)==T) {
      prob <- 1e-40
    }
    ln_prob <- log(prob)
    
  }else if(N_FSC[i] == 0 & N_LCC[i] >0){
    prob <- pmvnorm(lower = c(-Inf, 
                              func_g(Xm[i,],D_LCC[i], Nb_FSC = 0, Nb_LCC = N_LCC[i], theta[7:12])),
                    upper = c(func_g(Xm[i,],D_LCC[i], Nb_FSC = 0, Nb_LCC = N_LCC[i]+1, theta[7:12]),
                              func_f(Xm[i,],D_FSC[i], Nb_FSC = 1, Nb_LCC = N_LCC[i], theta[1:6])
                              ), mean = c(0,0), corr) -
      pmvnorm(lower = c(func_g(Xm[i,],D_LCC[i],  Nb_FSC = 0, Nb_LCC = N_LCC[i], theta[7:12]),
                        func_f(Xm[i,],D_FSC[i],  Nb_FSC = 1, Nb_LCC = N_LCC[i]-1,  theta[1:6]), 
                        ), 
              upper = c(func_g(Xm[i,],D_LCC[i],  Nb_FSC = 1, Nb_LCC = N_LCC[i], theta[7:12]),
                        func_f(Xm[i,],D_FSC[i],  Nb_FSC = 1, Nb_LCC = N_LCC[i], theta[1:6]) 
                        ), mean = c(0,0), corr)
    if( prob <= 0 | is.na(prob)==T){
      prob <- 1e-40
    }
    ln_prob <- log(prob)
    
  }else{
    prob <- pmvnorm(lower = c(-Inf ,-Inf),upper = c(func_g(Xm[i,],D_LCC[i], Nb_FSC = 0, Nb_LCC = 1, theta[7:12]),
                                                    func_f(Xm[i,],D_FSC[i], Nb_FSC = 1, Nb_LCC = 0, theta[1:6]), 
                                                    ), mean = c(0,0), corr) 
    
    prob <- pmvnorm(lower = c(func_g(Xm[i,],D_LCC[i],  Nb_FSC = 0, Nb_LCC = 1,  theta[7:12]),
                              func_f(Xm[i,],D_FSC[i],  Nb_FSC = 1, Nb_LCC = 0,  theta[1:6]), 
                              ), 
                    upper = rep(Inf, 2), mean = c(0,0), corr)
    if ( prob <= 0| is.na(prob)==T) {
      prob <- 1e-40
    }
    ln_prob <- log(prob)
  }
  (ln_prob)
}


neg_log_likelihood4 <- function(theta) {
  
  neg_ll2 <- 0
  
  
  for (i in 1:3983) {
    
    neg_ll2 <- neg_ll2 - LL_one4(theta,i)
  }
  
  
  return(neg_ll2)
}

opts <-list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-4,maxeval=10000,"print_level"=2)
res5<-nloptr(x0=theta_init3,eval_f=neg_log_likelihood4, eval_g_ineq = v,opts=opts)

