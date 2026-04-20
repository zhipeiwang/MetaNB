# return JAGS model text for weak realistic priors
get_tri_model_weak <- function(include_EVPI = FALSE, J = 1000) {

  model_string <- "
# Weak realistic priors

model{
for (i in 1:n_study){

n_event[i] ~ dbin(prev[i],n[i])
tp[i] ~ dbin(sens[i],n_event[i])
tn[i] ~ dbin(spec[i],n_nonevent[i])

logit( prev[i] )<- logitp[i]
logit( sens[i] )<-logitsens[i]
logit( spec[i] )<-logitspec[i]

logitp[i] ~ dnorm(etap,precp)
logitsens[i] ~ dnorm(etasens[i], precsens)
etasens[i]<-lambdasens0+lambdasens1*logitp[i]
logitspec[i] ~dnorm(etaspec[i], precspec)
etaspec[i]<-lambdaspec0+lambdaspec1*logitp[i]+lambdaspec2*logitsens[i]

NB[i] <- sens[i] * prev[i] - (1 - spec[i]) * (1 - prev[i]) * t / (1-t)
NB_TA[i] <- prev[i] - (1 - prev[i]) * t / (1-t)
RU[i] <- (NB[i] - max(NB_TA[i], 0) ) / (prev[i] - max(NB_TA[i], 0))
}

# compute summary statistics

mean.p<-etap
mean.sens<-lambdasens0+lambdasens1*etap
mean.spec<-lambdaspec0+lambdaspec1*etap+lambdaspec2*mean.sens
pooledprev <- exp(mean.p)/(1+exp(mean.p))
pooledsens <- exp(mean.sens)/(1+exp(mean.sens))
pooledspec <- exp(mean.spec)/(1+exp(mean.spec))


pooledNB<-pooledsens*pooledprev-(1-pooledspec)*(1-pooledprev)*t/(1-t)
pooledNB_TA<-pooledprev-(1-pooledprev)*t/(1-t)
pooledRU<-(pooledNB - max(pooledNB_TA, 0) ) / (pooledprev - max(pooledNB_TA, 0))

# pooled NB / NB_TA / RU at known prevalence
pooledNB_known<-pooledsens*prev_known-(1-pooledspec)*(1-prev_known)*t/(1-t)
pooledNB_TA_known<-prev_known-(1-prev_known)*t/(1-t)
pooledRU_known <- (pooledNB_known - max(pooledNB_TA_known, 0)) / (prev_known - max(pooledNB_TA_known, 0))


# priors
etap~dnorm(mu_etap, tau_etap)
lambdasens0~dnorm(mu_lambdasens0, tau_lambdasens0)
lambdaspec0~dnorm(mu_lambdaspec0, tau_lambdaspec0)

# Fisher or Uniform priors for correlations
zss~dnorm(mu_zss, tau_zss)
corr.sens.spec<-(exp(2*zss)-1)/(exp(2*zss)+1)
zsp~dnorm(mu_zsp, tau_zsp)
corr.spec.prev<-(exp(2*zsp)-1)/(exp(2*zsp)+1)
corr.sens.prev~ dunif(a_csp, b_csp)

# Halfnormal priors for between-cluster heterogeneity
varprev~dnorm(0, tau_varprev)I(0, )
varsens~dnorm(0, tau_varsens)I(0, )
varspec~dnorm(0, tau_varspec)I(0, )

# to monitor
logvarsens<-log(varsens)
logvarspec<-log(varspec)
logvarprev<-log(varprev)

# Implied relations
precp<-1/varprev
###
varsens.c<- varsens-pow(lambdasens1,2)*varprev
precsens<-1/varsens.c
###
varspec.c<- varspec-pow(lambdaspec1,2)*varprev-pow(lambdaspec2,2)*varsens-2*lambdaspec1*lambdaspec2*lambdasens1*varprev
precspec<-1/varspec.c

lambdasens1<-corr.sens.prev*sqrt(varsens)/sqrt(varprev)
###
lambdaspec1<-((corr.spec.prev-corr.sens.prev*corr.sens.spec)*sqrt(varspec))/(sqrt(varprev)*(1-pow(corr.sens.prev,2)))
###
lambdaspec2<-((corr.sens.spec-corr.sens.prev*corr.spec.prev)*sqrt(varspec))/(sqrt(varsens)*(1-pow(corr.sens.prev,2)))


# Predict new triade of sens and spec and prev
logitp.new~ dnorm(etap,precp)
logitsens.new~ dnorm(etasens.new, precsens)
etasens.new<-lambdasens0+lambdasens1*logitp.new
logitspec.new~dnorm(etaspec.new, precspec)
etaspec.new<-lambdaspec0+lambdaspec1*logitp.new+lambdaspec2*logitsens.new

prevnew<-exp(logitp.new)/(1+exp(logitp.new))
sensnew<-exp(logitsens.new)/(1+exp(logitsens.new))
specnew<-exp(logitspec.new)/(1+exp(logitspec.new))

NBnew<-sensnew*prevnew-(1-specnew)*(1-prevnew)*t/(1-t)
NBnew_TA<-prevnew-(1-prevnew)*t/(1-t)

probuseful<-equals(max(max(NBnew,NBnew_TA),0), NBnew)

RUnew<-(NBnew - max(NBnew_TA, 0) ) / (prevnew - max(NBnew_TA, 0))

# Predict new triad at known prevalence
logitp.new.known<-logit(prev_known)

logitsens.new.known ~ dnorm(etasens.new.known,precsens)
etasens.new.known<-lambdasens0+lambdasens1*logitp.new.known

logitspec.new.known~dnorm(etaspec.new.known,precspec)
etaspec.new.known<-lambdaspec0+lambdaspec1*logitp.new.known+lambdaspec2*logitsens.new.known

sensnew.known<-exp(logitsens.new.known)/(1+exp(logitsens.new.known))
specnew.known<-exp(logitspec.new.known)/(1+exp(logitspec.new.known))

NBnew_known<- sensnew.known*prev_known-(1-specnew.known)*(1-prev_known)*t/(1-t)
NBnew_TA_known<-prev_known-(1-prev_known)*t/(1-t)

probuseful_known<-equals(max(max(NBnew_known,NBnew_TA_known),0),NBnew_known)

RUnew_known<-(NBnew_known - max(NBnew_TA_known, 0) ) / (prev_known - max(NBnew_TA_known, 0))


  "
if (include_EVPI) {
model_string <- paste0(model_string, "
# --- EVPI population-level block ---
for (j in 1:J) {
  logitp.new2[j] ~ dnorm(etap, precp)
  etasens.new2[j] <- lambdasens0 + lambdasens1 * logitp.new2[j]
  logitsens.new2[j] ~ dnorm(etasens.new2[j], precsens)
  etaspec.new2[j] <- lambdaspec0 + lambdaspec1 * logitp.new2[j]
                     + lambdaspec2 * logitsens.new2[j]
  logitspec.new2[j] ~ dnorm(etaspec.new2[j], precspec)

  prevnew2[j] <- exp(logitp.new2[j]) / (1 + exp(logitp.new2[j]))
  sensnew2[j] <- exp(logitsens.new2[j]) / (1 + exp(logitsens.new2[j]))
  specnew2[j] <- exp(logitspec.new2[j]) / (1 + exp(logitspec.new2[j]))

  NBnew2[j]    <- sensnew2[j] * prevnew2[j]
                  - (1 - specnew2[j]) * (1 - prevnew2[j]) * t / (1 - t)
  NBnew_TA2[j] <- prevnew2[j]
                  - (1 - prevnew2[j]) * t / (1 - t)
                  }
  ENBnew    <- mean(NBnew2[])
  ENBnew_TA <- mean(NBnew_TA2[])
    ")
}

model_string <- paste0(model_string, "}")

return(model_string)
}

