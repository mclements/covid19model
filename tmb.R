library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(EnvStats)
library(bayesplot)

countries <- c(
  "Denmark",
  "Italy",
  "Germany",
  "Spain",
  "United_Kingdom",
  "France",
  "Norway",
  "Belgium",
  "Austria", 
  "Sweden",
  "Switzerland"
)

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0) {
  args = 'base'
} 
StanModel = args[1]

print(sprintf("Running %s",StanModel))

## Reading all data
d=readRDS('data/COVID-19-up-to-date.rds')

## get CFR
cfr.by.country = read.csv("data/weighted_fatality.csv")
cfr.by.country$country = as.character(cfr.by.country[,2])
cfr.by.country$country[cfr.by.country$country == "United Kingdom"] = "United_Kingdom"

serial.interval = read.csv("data/serial_interval.csv")
covariates = read.csv('data/interventions.csv', stringsAsFactors = FALSE)
covariates <- covariates[1:11, c(1,2,3,4,5,6, 7, 8)]

## making all covariates that happen after lockdown to have same date as lockdown
covariates$schools_universities[covariates$schools_universities > covariates$lockdown] <- covariates$lockdown[covariates$schools_universities > covariates$lockdown]
covariates$travel_restrictions[covariates$travel_restrictions > covariates$lockdown] <- covariates$lockdown[covariates$travel_restrictions > covariates$lockdown] 
covariates$public_events[covariates$public_events > covariates$lockdown] <- covariates$lockdown[covariates$public_events > covariates$lockdown]
covariates$sport[covariates$sport > covariates$lockdown] <- covariates$lockdown[covariates$sport > covariates$lockdown]
covariates$social_distancing_encouraged[covariates$social_distancing_encouraged > covariates$lockdown] <- covariates$lockdown[covariates$social_distancing_encouraged > covariates$lockdown]
covariates$self_isolating_if_ill[covariates$self_isolating_if_ill > covariates$lockdown] <- covariates$lockdown[covariates$self_isolating_if_ill > covariates$lockdown]

p <- ncol(covariates) - 1
forecast = 0

DEBUG = FALSE
if(DEBUG == FALSE) {
  N2 = 75 # Increase this for a further forecast
}  else  {
  ### For faster runs:
  # countries = c("Austria","Belgium") #,Spain")
  N2 = 75
}
# countries = c("Italy","United_Kingdom","Spain","Norway","Austria","Switzerland")

dates = list()
reported_cases = list()
stan_data = list(M=length(countries),N=NULL,p=p,x1=poly(1:N2,2)[,1],x2=poly(1:N2,2)[,2],
                 y=NULL,covariate1=NULL,covariate2=NULL,covariate3=NULL,covariate4=NULL,covariate5=NULL,covariate6=NULL,covariate7=NULL,deaths=NULL,f=NULL,
                 N0=6,cases=NULL,LENGTHSCALE=7,SI=serial.interval$fit[1:N2],
                 EpidemicStart = NULL) # N0 = 6 to make it consistent with Rayleigh
deaths_by_country = list()


for(Country in countries) {
  CFR=cfr.by.country$weighted_fatality[cfr.by.country$country == Country]
  
  covariates1 <- covariates[covariates$Country == Country, 2:8]
  
  d1=d[d$Countries.and.territories==Country,]
  d1$date = as.Date(d1$DateRep,format='%d/%m/%Y')
  d1$t = decimal_date(d1$date) 
  d1=d1[order(d1$t),]
  index = which(d1$Cases>0)[1]
  index1 = which(cumsum(d1$Deaths)>=10)[1] # also 5
  index2 = index1-30
  
  print(sprintf("First non-zero cases is on day %d, and 30 days before 5 days is day %d",index,index2))
  d1=d1[index2:nrow(d1),]
  stan_data$EpidemicStart = c(stan_data$EpidemicStart,index1+1-index2)
  
  
  for (ii in 1:ncol(covariates1)) {
    covariate = names(covariates1)[ii]
    d1[covariate] <- (as.Date(d1$DateRep, format='%d/%m/%Y') >= as.Date(covariates1[1,covariate]))*1  # should this be > or >=?
  }
  
  dates[[Country]] = d1$date
  # hazard estimation
  N = length(d1$Cases)
  print(sprintf("%s has %d days of data",Country,N))
  forecast = N2 - N
  if(forecast < 0) {
    print(sprintf("%s: %d", Country, N))
    print("ERROR!!!! increasing N2")
    N2 = N
    forecast = N2 - N
  }
  
  h = rep(0,forecast+N) # discrete hazard rate from time t = 1, ..., 100
  if(DEBUG) { # OLD -- but faster for testing this part of the code
    mean = 18.8
    cv = 0.45
    
    for(i in 1:length(h))
      h[i] = (CFR*pgammaAlt(i,mean = mean,cv=cv) - CFR*pgammaAlt(i-1,mean = mean,cv=cv)) / (1-CFR*pgammaAlt(i-1,mean = mean,cv=cv))
  } else { # NEW
    mean1 = 5.1; cv1 = 0.86; # infection to onset
    mean2 = 18.8; cv2 = 0.45 # onset to death
    ## assume that CFR is probability of dying given infection
    x1 = rgammaAlt(5e6,mean1,cv1) # infection-to-onset ----> do all people who are infected get to onset?
    x2 = rgammaAlt(5e6,mean2,cv2) # onset-to-death
    f = ecdf(x1+x2)
    convolution = function(u) (CFR * f(u))
    
    h[1] = (convolution(1.5) - convolution(0)) 
    for(i in 2:length(h)) {
      h[i] = (convolution(i+.5) - convolution(i-.5)) / (1-convolution(i-.5))
    }
  }
  s = rep(0,N2)
  s[1] = 1 
  for(i in 2:N2) {
    s[i] = s[i-1]*(1-h[i-1])
  }
  f = s * h
  
  y=c(as.vector(as.numeric(d1$Cases)),rep(-1,forecast))
  reported_cases[[Country]] = as.vector(as.numeric(d1$Cases))
  deaths=c(as.vector(as.numeric(d1$Deaths)),rep(-1,forecast))
  cases=c(as.vector(as.numeric(d1$Cases)),rep(-1,forecast))
  deaths_by_country[[Country]] = as.vector(as.numeric(d1$Deaths))
  covariates2 <- as.data.frame(d1[, colnames(covariates1)])
  # x=1:(N+forecast)
  covariates2[N:(N+forecast),] <- covariates2[N,]
  
  ## append data
  stan_data$N = c(stan_data$N,N)
  stan_data$y = c(stan_data$y,y[1]) # just the index case!
  # stan_data$x = cbind(stan_data$x,x)
  stan_data$covariate1 = cbind(stan_data$covariate1,covariates2[,1])
  stan_data$covariate2 = cbind(stan_data$covariate2,covariates2[,2])
  stan_data$covariate3 = cbind(stan_data$covariate3,covariates2[,3])
  stan_data$covariate4 = cbind(stan_data$covariate4,covariates2[,4])
  stan_data$covariate5 = cbind(stan_data$covariate5,covariates2[,5])
  stan_data$covariate6 = cbind(stan_data$covariate6,covariates2[,6])
  stan_data$covariate7 = cbind(stan_data$covariate7,covariates2[,7]) 
  stan_data$f = cbind(stan_data$f,f)
  stan_data$deaths = cbind(stan_data$deaths,deaths)
  stan_data$cases = cbind(stan_data$cases,cases)
  
  stan_data$N2=N2
  stan_data$x=1:N2
  if(length(stan_data$N) == 1) {
    stan_data$N = as.array(stan_data$N)
  }
}

stan_data$covariate2 = 0 * stan_data$covariate2 # remove travel bans
stan_data$covariate4 = 0 * stan_data$covariate5 # remove sport

#stan_data$covariate1 = stan_data$covariate1 # school closure
stan_data$covariate2 = stan_data$covariate7 # self-isolating if ill
#stan_data$covariate3 = stan_data$covariate3 # public events
# create the `any intervention` covariate
stan_data$covariate4 = 1*as.data.frame((stan_data$covariate1+
                                          stan_data$covariate3+
                                          stan_data$covariate5+
                                          stan_data$covariate6+
                                          stan_data$covariate7) >= 1)
stan_data$covariate5 = stan_data$covariate5 # lockdown
stan_data$covariate6 = stan_data$covariate6 # social distancing encouraged
stan_data$covariate7 = 0 # models should only take 6 covariates

if(DEBUG) {
  for(i in 1:length(countries)) {
    write.csv(
      data.frame(date=dates[[i]],
                 `school closure`=stan_data$covariate1[1:stan_data$N[i],i],
                 `self isolating if ill`=stan_data$covariate2[1:stan_data$N[i],i],
                 `public events`=stan_data$covariate3[1:stan_data$N[i],i],
                 `government makes any intervention`=stan_data$covariate4[1:stan_data$N[i],i],
                 `lockdown`=stan_data$covariate5[1:stan_data$N[i],i],
                 `social distancing encouraged`=stan_data$covariate6[1:stan_data$N[i],i]),
      file=sprintf("results/%s-check-dates.csv",countries[i]),row.names=F)
  }
}

## NEW CODE: ignore u
## library(rstan)
## library(data.table)
## library(lubridate)
## library(gdata)
## library(EnvStats)
## library(bayesplot)
## load("results/base-1296366.Rdata")
library(TMB)
covariateNames <- c("school closure",
                    "self isolating if ill",
                    "public events",
                    "government makes any intervention",
                    "lockdown",
                    "social distancing encouraged")
src <- "#include <TMB.hpp>
// template<class T, size_t N>
// constexpr size_t size(T (&)[N]) { return N; }
// #define READ_VECTOR(v,inv) vector<Type> v(size(inv)); \\
//  {for (int i = 0; i < size(inv); i ++) \\
//    v(i) = inv[i];} 
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(M); // number of countries
  DATA_INTEGER(N0); // number of days for which to impute infections
  DATA_IVECTOR(N); // length M; days of observed data for country m. each entry must be <= N2
  DATA_INTEGER(N2); // days of observed data + # of days to forecast
  DATA_MATRIX(deaths); // dim N2, M; reported deaths -- the rows with i > N contain -1 and should be ignored
  DATA_MATRIX(f); // dim N2, M; h * s; improper density from incidence to death
  DATA_MATRIX(x0); // school closure; dim N2, M
  DATA_MATRIX(x1); // self-isolating if ill; dim N2, M
  DATA_MATRIX(x2); // restrict public events; dim N2, M
  DATA_MATRIX(x3); // government makes any intervention; dim N2, M
  DATA_MATRIX(x4); // lockdown; dim N2, M
  DATA_MATRIX(x5); // social distancing encouraged; dim N2, M
  DATA_IVECTOR(ix); // vector of indices for included covariates
  DATA_IVECTOR(EpidemicStart); // index for first death from the epidemic by country; length M
  DATA_VECTOR(SI); // length N2; fixed pre-calculated SI using emprical data from Neil Fergusson
  PARAMETER(alpha); // log(R0) fixed effect
  // PARAMETER_VECTOR(u); // random intercept for log(Rt); length M
  // PARAMETER(log_sigma_u); // sd for u
  PARAMETER_VECTOR(beta); // the coefficients for effects of interventions; length 6
  PARAMETER(logtau); // mean for the daily infections
  PARAMETER_VECTOR(v); // random effect for log(mean daily infections); length M
  PARAMETER(logpsi); // negative-binomial variance parameter for the number of deaths
  Type convolution;
  matrix<Type> prediction(N2,M); // initialised to 0
  matrix<Type> E_deaths(N2,M);
  matrix<Type> Rt(N2,M);
  vector<Type> x(ix.size());
  for (int m=0; m<M; m++) { // by country
    for (int i=0; i<N0; i++) { // by time for 0..N0-1
      prediction(i,m) = exp(logtau+v(m));
    }
    for (int i=N0; i<N2; i++) {// by time for N0..N2-1
      for (int j=0; j<ix.size(); j++)
          x(j) = (ix(j)==0) ? x0(i,m) :
                 (ix(j)==1) ? x1(i,m) :
                 (ix(j)==2) ? x2(i,m) :
                 (ix(j)==3) ? x3(i,m) :
                 (ix(j)==4) ? x4(i,m) : x5(i,m) ;
      Rt(i,m) = exp(alpha-(x*beta).sum());
      // Rt(i,m) = exp(alpha+u(m)-(x*beta).sum());
      convolution=Type(0);
      for(int j=0; j<i-1; j++) {
	convolution += prediction(j, m)*SI(i-j);
      }
      prediction(i, m) = Rt(i,m) * convolution;
    }
    E_deaths(0, m)= Type(1.0e-9);
    for (int i=1; i<N2; i++) {
      E_deaths(i,m)= Type(0);
      for(int j=0; j<i-1; j++) {
	E_deaths(i,m) += prediction(j,m)*f(i-j,m);
      }
    }
  }
  Type logl = Type(0);
  for(int m=0; m<M; m++) {
    for(int i=EpidemicStart[m]-1; i<N(m); i++) {
      logl += dnbinom2(deaths(i,m),E_deaths(i,m),E_deaths(i,m)+E_deaths(i,m)*E_deaths(i,m)/exp(logpsi),1); 
    }
    // logl += dnorm(u(m), Type(0), exp(log_sigma_u), true);
    logl += dnorm(v(m), Type(0), sqrt(log(Type(2))), true);
  }
  return -logl;
}"
write(src, file="tmb.cpp")
compile("tmb.cpp") # slow compilation
dyn.load(dynlib("tmb"))
## dyn.unload(dynlib("tmb"))
##

stan_data2 <- with(stan_data, list(M=M, N0=N0, N=N, N2=N2, deaths=deaths, f=f,
                                   x0=covariate1, x1=covariate2, x2=covariate3,
                                   x3=as.matrix(covariate4), x4=covariate5, x5=covariate6,
                                   EpidemicStart=EpidemicStart, SI=SI, ix=0:5))
flaxman <- function(ix=0:5, drop=NULL) {
    stopifnot(!any(ix<0 | ix>5))
    if(!is.null(drop)) {
        stopifnot(!any(drop<0 | drop>5))
        ix <- setdiff(ix,drop)
    }
    stan_data2$ix <- ix
    init <- list(alpha=log(2.4),
                 ## u=rep(0,11), log_sigma_u=0,
                 beta=rep(0,length(ix)), logtau=log(30),
                 v=rep(0,11), 
                 logpsi=1)
    f <- MakeADFun(data=stan_data2,
                   parameters=init,
                   method="nlminb",
                   ## random=c("u","v"),                
                   random="v",
                   DLL="tmb", silent=TRUE)
    fit <- nlminb(f$par,f$fn,f$gr)
    fit
    h <- optimHess(fit$par, f$fn, f$gr)
    covariateNames <- c("school closure",
                        "self isolating if ill",
                        "public events",
                        "government makes any intervention",
                        "lockdown",
                        "social distancing encouraged")
    nms <- names(fit$par)
    nms[nms=="beta"] <- covariateNames[ix+1]
    printCoefmat(transform("rownames<-"(data.frame(Estimate=fit$par,Std.Err=sqrt(diag(solve(h)))),
                                        nms),
                           "Z value"=Estimate/Std.Err,
                           "Pr(>z)"=2*(1-pnorm(abs(Estimate/Std.Err)))),
                 has.Pvalue=TRUE, P.values=TRUE)
    cat("log(l) =",-fit$objective,"\n")
    cat("df =",length(fit$par),"\n")
    invisible(list(fit=fit,h=h))
} # a formula interface would be more pleasant

fit <- flaxman() # schools+, self-isolating(+), public(-), govt(-), lockdown++, social-
sum(fit$fit$par[c(3,4,5,7)]) # Sweden
##
flaxman(drop=3) # schools+, self-isolating+, public(-), lockdown++, social--
flaxman(drop=2:3) # schools+, self-isolating+, lockdown++, social--
flaxman(drop=c(2:3,5)) # schools+, self-isolating(+), lockdown+
flaxman(c(0,4)) # schools+, lockdown+
##
flaxman(drop=5) # schools++, self-isolating+, public(-), govt-, lockdown++
flaxman(drop=c(2,5)) # schools++, self-isolating+, govt-, lockdown++
flaxman(drop=c(2,3,5)) # schools+, self-isolating(+), lockdown(+)
##
flaxman(3) # govt+
flaxman(c(3,5)) # govt(+), social(+)
flaxman(c(3,5)) # govt(+), social(+)
flaxman(c(0,3,4)) # schools+, govt(-), lockdown(+)
##
flaxman(c(1,2,3,5)) # self-isolating(-), public+, govt(+), social(-)
flaxman(c(1,2,3)) # self-isolating(-), public+, govt(+)
flaxman(c(2,3)) # public+, govt(+)

## investigate covariates
with(stan_data2,
     lapply(1:ncol(x0), function(i) cor(cbind(x0[,i],x1[,i],x2[,i],x3[,i],x4[,i],x5[,i]))))
stan_data3 <- with(stan_data2, data.frame(x0=as.vector(x0), x1=as.vector(x1), x2=as.vector(x2), x3=as.vector(x3), x4=as.vector(x4), x5=as.vector(x5))) 
summary(glm(x5~x0+x1+x2+x3+x4, data=stan_data3, family=binomial))
## x5
xtabs(~x0+x5,data=stan_data3) # concordance between social distancing and school closures (discordant in Sweden?)
xtabs(~x1+x5,data=stan_data3) # social distancing => self isolating if ill
xtabs(~x4+x5,data=stan_data3) # lockdown => social distancing
## x1
xtabs(~x1+x0,data=stan_data3) # good concordance between school closures and self-isolating
xtabs(~x4+x0,data=stan_data3) # lockdown => school closures
xtabs(~x5+x0,data=stan_data3) # good concordance between school closures and social distancing


