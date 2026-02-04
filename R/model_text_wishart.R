# return JAGS model text for Wishart prior
get_tri_model_wishart <- function(include_EVPI = compute_EVPI, J = 1000) {

  model_string <- "
# Inverse Wishart prior for the variance-covariance matrix

model{
# ---------- Likelihood ----------
for (i in 1:n_study){

# n[i]<-n_event[i]+n_nonevent[i]
n_event[i] ~ dbin(prev[i],n[i])
tp[i] ~ dbin(sens[i],n_event[i])
tn[i] ~ dbin(spec[i],n_nonevent[i])

# Multivariate normal on logits with precision matrix T
logits[i,1:3] ~ dmnorm(mu[1:3],T[1:3,1:3])

# Back-transform to probabilities
prev[i] <- exp(logits[i,1])/(1+exp(logits[i,1]))
sens[i] <- exp(logits[i,2])/(1+exp(logits[i,2]))
spec[i] <- exp(logits[i,3])/(1+exp(logits[i,3]))

# Per study
NB[i] <- sens[i] * prev[i] - (1 - spec[i]) * (1 - prev[i]) * t / (1-t)
NB_TA[i] <- prev[i] - (1 - prev[i]) * t / (1-t)
RU[i] <- (NB[i] - max(NB_TA[i], 0) ) / (prev[i] - max(NB_TA[i], 0))
}

# ---------- Pooled probabilities ----------
pooledprev <- exp(mu[1])/(1+exp(mu[1]))
pooledsens <- exp(mu[2])/(1+exp(mu[2]))
pooledspec <- exp(mu[3])/(1+exp(mu[3]))

# ---------- Derive covariance & correlations from precision ----------
# (T is precision; tau is covariance = inverse(T))
tau[1:3,1:3] <- inverse(T[1:3,1:3])

tau2prev<- tau[1,1]
tau2sens<- tau[2,2]
tau2spec<- tau[3,3]

cov12 <- tau[1,2]
corr12 <- cov12/sqrt(tau2sens*tau2prev) # corr(prev, sens)
cov13 <- tau[1,3]
corr13 <- cov13/sqrt(tau2spec*tau2prev) # corr(prev, spec)
cov23 <- tau[2,3]
corr23 <- cov23/sqrt(tau2spec*tau2sens) # corr(sens, spec)

# compute summary statistics

pooledNB<-pooledsens*pooledprev-(1-pooledspec)*(1-pooledprev)*t/(1-t)
pooledNB_TA<-pooledprev-(1-pooledprev)*t/(1-t)
pooledRU<-(pooledNB - max(pooledNB_TA, 0) ) / (pooledprev - max(pooledNB_TA, 0))


# pooled NB / NB_TA / RU at known prevalence
pooledNB_ref<-pooledsens*prev_ref-(1-pooledspec)*(1-prev_ref)*t/(1-t)
pooledNB_TA_ref<-prev_ref-(1-prev_ref)*t/(1-t)
pooledRU_ref <- (pooledNB_ref - max(pooledNB_TA_ref, 0)) / (prev_ref - max(pooledNB_TA_ref, 0))


# ---------- Priors ----------
mu[1:3] ~ dmnorm(mn[1:3],prec[1:3,1:3])
T[1:3,1:3] ~ dwish(R[1:3,1:3],df)

# ---------- Predict new triade of sens and spec and prev ----------

logitsnew[1:3] ~ dmnorm(mu[1:3],T[1:3,1:3])
prevnew<-exp(logitsnew[1])/(1+exp(logitsnew[1]))
sensnew<-exp(logitsnew[2])/(1+exp(logitsnew[2]))
specnew<-exp(logitsnew[3])/(1+exp(logitsnew[3]))
NBnew<-sensnew*prevnew-(1-specnew)*(1-prevnew)*t/(1-t)
NBnew_TA<-prevnew-(1-prevnew)*t/(1-t)
NBnew_ref<-sensnew*prev_ref-(1-specnew)*(1-prev_ref)*t/(1-t)
NBnew_TA_ref<-prev_ref-(1-prev_ref)*t/(1-t)

probuseful<-equals(max(max(NBnew,NBnew_TA),0), NBnew)
probuseful_ref<-equals(max(max(NBnew_ref,NBnew_TA_ref),0), NBnew_ref)

RUnew<-(NBnew - max(NBnew_TA, 0) ) / (prevnew - max(NBnew_TA, 0))
RUnew_ref<-(NBnew_ref - max(NBnew_TA_ref, 0) ) / (prev_ref - max(NBnew_TA_ref, 0))

  "
if (include_EVPI) {
  model_string <- paste0(model_string, "
# --- EVPI population-level block ---
for (j in 1:J) {
  logitsnew2[j,1:3] ~ dmnorm(mu[1:3], T[1:3,1:3])

  prevnew2[j] <- exp(logitsnew2[j,1]) / (1 + exp(logitsnew2[j,1]))
  sensnew2[j] <- exp(logitsnew2[j,2]) / (1 + exp(logitsnew2[j,2]))
  specnew2[j] <- exp(logitsnew2[j,3]) / (1 + exp(logitsnew2[j,3]))

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
