library(tidyverse)
library(ggplot2)
# TRANSITION MATRIX FOR THE problem
tmat = matrix(FALSE, 3,3)
colnames(tmat) = rownames(tmat) = c("1","2","3")
tmat[1,c(2,3)] = TRUE
tmat[2,c(1,3)] = TRUE

#tmat

# define the hazard functions (one for each of the three events) ---------
set.seed(1)
alpha1 <- 0.8 # weibull parameter 1
alpha2 <- 1.2 #weibull parameter 2
k1 <- 1 # constant 1
k2 <- 1 #constant 2

#intercepts
int12 <- 1.3
int13 <- 0.5
int21 <- 1.1
int23 <- 0.9

f2 =  c("weibullsurv","weibullsurv","exponentialsurv","exponentialsurv")
fam_models = f2

# for plotting baseline haz
h12_0p <- function(t) alpha1 * t^{alpha1 - 1} * exp(int12)
h13_0p <- function(t) alpha2 * t^{alpha2 - 1} * exp(int13)
h21_0p <- function(t) rep(k1, length(t)) * exp(int21)
h23_0p <- function(t) rep(k2, length(t)) * exp(int23)

# plot baseline hazards
tt = seq(0,5,0.01)
data.frame(t = tt,
           h12_0  = h12_0p(tt),
           h13_0  = h13_0p(tt),
           h21_0  = h21_0p(tt),
           h23_0  = h23_0p(tt)
) %>%
  pivot_longer(-t) %>%
  ggplot() + geom_line(aes(t,value, group= name, color = name, linetype=name)) + 
  scale_linetype_manual(values=c("solid","solid","solid" , "dashed"))

h12_0 <- function(t) alpha1 * t^{alpha1 - 1}
h13_0 <- function(t) alpha2 * t^{alpha2 - 1} 
h21_0 <- function(t) rep(k1, length(t)) 
h23_0 <- function(t) rep(k2, length(t)) 

# compute the cumulative hazards ------------------------------------------
## These are the integrals from t0 to t of the functions above
H12_0p = function(t,t0) t^{alpha1} - t0^{alpha1} * exp(int12)
H13_0p = function(t,t0) t^{alpha2} - t0^{alpha2} * exp(int13)
H21_0p = function(t,t0) k1*(t - t0) * exp(int21)
H23_0p = function(t,t0) k2*(t - t0) * exp(int23)

#tt = seq(0,2,0.01)
data.frame(t = tt,
           H12  = H12_0p(tt, 0),
           H13  = H13_0p(tt, 0),
           H21  = H21_0p(tt, 0),
           H23  = H23_0p(tt, 0)
) %>%
  pivot_longer(-t) %>%
  ggplot() + geom_line(aes(t,value, group= name, color = name, linetype=name)) +
  scale_linetype_manual(values=c("solid","solid","solid" , "dashed"))

H12_0 = function(t,t0) t^{alpha1} - t0^{alpha1}
H13_0 = function(t,t0) t^{alpha2} - t0^{alpha2}
H21_0 = function(t,t0) k1*(t - t0)
H23_0 = function(t,t0) k2*(t - t0)


# simulation start
# covariates
N = 300
X = cbind(scale(runif(N)), sample(c(0,1), N, replace = T))

beta12 = c(0.5,0.5)
beta13 = c(-.5,0.3)
beta21 = c(0.2,-0.2)
beta23 = c(0.4,-0.1)

# linear predictor
eta12 <- int12+X %*% beta12 
eta13 <- int13+X %*% beta13
eta21 <- int21+X %*% beta21
eta23 <- int23+X %*% beta23

# function to generate event times using uniroot
generate.my.times <- function(n, max.int, temp,t0) { # gives n independent replicates of T
  s_time <- NULL
  i <- 1
  
  while(length(s_time) < n) {
    u <- runif(1)
    ## If endpoints are of opposite sign, call uniroot:
    if (temp(0,t0,log(1-u))*temp(max.int,t0,log(1 - u)) < 0) {
      res = uniroot(temp, c(0,max.int),
                    tol = 0.0001, t0 = t0, y = log(1-u))
      if(res$root>0)
      {
        s_time[i] = res$root
        i = i+1
      }
    }
    # else 
    #   {
    #     cat("values at endpoint not of opposite sign\n")
    #   }
  }
  return(s_time)
}


# looping through all patients, simulating the states that happen to them
data = data.frame()
for(j in 1:N){
  TT = 7   ## longest time 
  state = sample(c(1,2),1,prob = c(0.6,0.4)) ## 60% chance of state 1, 40% of state 2
  new_state = state
  times = 0
  ctime = 0
  from_state = c()
  to_state = c()
  all_times = c() 
  duration  = c()
  i = 1
  while(ctime<TT & new_state!=3) 
  {
    i = i+1
    if(state[i-1]==1) # checks which state got sampled above
    {
      possible_tras <- as.vector(which(tmat[1,])) # gives possible states to jump to
      
      H12 <- function(t,t0) H12_0(t,t0) * exp(eta12[j])
      H13 <- function(t,t0) H13_0(t,t0) * exp(eta13[j])
      temp = function(t,t0,y) H12(t,t0) + H13(t,t0) + y 
      
      t = generate.my.times(1, # generate 1 time
                            max.int =  TT + 2, ## interval searches for uniroot from 0 to TT + 2 
                            temp  =  temp, # temp is made above
                            t0 = times[i-1])
      #print(t)
      prob = h12_0(t) * exp(eta12[j])/(h12_0(t)* exp(eta12[j]) + 
                                         h13_0(t)* exp(eta13[j]))
      new_state = sample(x = possible_tras, size = 1, prob = c(prob, 1-prob)) # samples either state 2 or 3 with calculated probability
      state = c(state, new_state)
      times = c(times,t)
      duration = c(duration,t-times[i-1])
      from_state = c(from_state, rep(state[i-1],2)) 
      to_state = c(to_state, c(new_state, setdiff(possible_tras, new_state))) 
      all_times = c(all_times, c(t, NA)) # assigns time NA to event that didnt happen
      ctime =  t
      #print(paste(state[i-1], " to ", state[i]))
      
    }
    else if (state[i-1]==2)
    {
      possible_tras <- as.vector(which(tmat[2,])) # gives possible states to jump to
      
      H21 <- function(t,t0) H21_0(t,t0) * exp(eta21[j])
      H23 <- function(t,t0) H23_0(t,t0) * exp(eta23[j])
      temp = function(t,t0,y) H21(t,t0) + H23(t,t0) + y 
      
      t = generate.my.times(1, # generate 1 time
                            max.int =  TT + 2, ## interval searches for uniroot from 0 to TT + 2 
                            temp  =  temp, # temp is made above
                            t0 = times[i-1])
      #print(t)
      prob = h21_0(t) * exp(eta21[j])/(h21_0(t)* exp(eta21[j]) + 
                                         h23_0(t)* exp(eta23[j]))
      #print(prob)
      new_state = sample(x = possible_tras, size = 1, prob = c(prob, 1-prob)) # samples either state 2 or 3 with calculated probability
      state = c(state, new_state)
      times = c(times,t)
      duration = c(duration,t-times[i-1])
      from_state = c(from_state, rep(state[i-1],2)) 
      to_state = c(to_state, c(new_state, setdiff(possible_tras, new_state))) 
      all_times = c(all_times, c(t, NA))
      ctime =  t
      #print(paste(state[i-1], " to ", state[i]))
      
    }
  }
  df = data.frame(id = j, #temporary dataframe to store info for each patient j
                  Tstart = rep(times[-length(state)], each = 2), # last observed time is removed
                  Tstop  = rep(times[-1], each = 2), # first time is removed
                  duration = rep(times[-1], each = 2) - 
                    rep(times[-length(state)], each = 2),
                  from = from_state,
                  to = to_state,
                  status = as.numeric(!is.na(all_times))) #makes TRUE/FALSE, then converts to 1/0
  
  df = df %>% mutate(status = ifelse(Tstop > TT, 0, status), #make status=0 if Tstop is over the max value TT
                     Tstop = ifelse(Tstop > TT, TT,Tstop))  #make Tstop value be TT if the recorded Tstop value was larger than TT
  
  data = rbind(data, df) #puts all the small temporary data frames into one big one
}
data = left_join(data, data.frame(id = 1:N, cov1 = X[,1], cov2 = X[,2])) #add cov1 & cov2 to the end of data
#plot(data$Tstop)
head(data)

ggplot(data, aes(x=data$Tstop)) + geom_histogram(bins = 50) + xlab("Time") + ggtitle("Histogram of simulated event times")

# look at history of some patients
np = 12 # makes a nice grid
ii = sample(1:N,np) %>% sort() # takes a random sample of 12 patients, sorts in increasing order of ID
data %>% dplyr::filter(id %in% ii, status == 1) %>%
  dplyr::select(id,Tstart,Tstop, from,to) %>%
  rename(T0_start = Tstart,
         T0_stop = Tstop,
         T1_start = from,
         T1_stop = to) %>%
  pivot_longer(cols = -id,
               names_to = c('.value', 'brand'),
               names_sep = "_") %>%
  ggplot()+geom_step(aes(T0,T1, group = id)) + facet_wrap(.~id) +
  xlim(c(0,2)) + ylim(c(1,3))+ xlab("time")+ylab("states")

# percentage of each transition
per_trans <- data %>% group_by(from, to) %>% summarise(n = n(), m = sum(status)) %>%
  mutate(percent = m/n * 100) %>%
  dplyr::select(from, to, percent)  

data12 = data %>% dplyr::filter(from  == 1, to == 2) 
data13 = data %>% dplyr::filter(from  == 1, to == 3) 
data21 = data %>% dplyr::filter(from  == 2, to == 1) 
data23 = data %>% dplyr::filter(from  == 2, to == 3) 

fixed.eff <- data.frame(
  int12 = c(rep(1,dim(data12)[1]), rep(NA,dim(data13)[1]), rep(NA,dim(data21)[1]), rep(NA,dim(data23)[1])),
  int13 = c(rep(NA,dim(data12)[1]), rep(1,dim(data13)[1]), rep(NA,dim(data21)[1]), rep(NA,dim(data23)[1])),
  int21 = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]), rep(1,dim(data21)[1]), rep(NA,dim(data23)[1])),
  int23 = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]), rep(NA,dim(data21)[1]), rep(1,dim(data23)[1])),
  
  cov1_12 = c(data12$cov1, rep(NA,dim(data13)[1]), rep(NA,dim(data21)[1]), rep(NA,dim(data23)[1])),
  cov1_13 = c(rep(NA,dim(data12)[1]), data13$cov1, rep(NA,dim(data21)[1]), rep(NA,dim(data23)[1])),
  cov1_21 = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]), data21$cov1, rep(NA,dim(data23)[1])),
  cov1_23 = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]), rep(NA,dim(data21)[1]), data23$cov1),
  
  cov2_12 = c(data12$cov2, rep(NA,dim(data13)[1]), rep(NA,dim(data21)[1]), rep(NA,dim(data23)[1])),
  cov2_13 = c(rep(NA,dim(data12)[1]), data13$cov2, rep(NA,dim(data21)[1]), rep(NA,dim(data23)[1])),
  cov2_21 = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]), data21$cov2, rep(NA,dim(data23)[1])),
  cov2_23 = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]), rep(NA,dim(data21)[1]), data23$cov2)
)

# ------------------------------------------------------------------------------------------
# how y data should look like to fit current dataset
library(INLA)


y.surv12 <- inla.surv(time = c(data12$Tstop, rep(NA,dim(data13)[1]),rep(NA,dim(data21)[1]), rep(NA,dim(data23)[1])),
                      event = c(data12$status, rep(NA,dim(data13)[1]), rep(NA,dim(data21)[1]),rep(NA,dim(data23)[1])),
                      truncation = c(data12$Tstart, rep(NA,dim(data13)[1]), rep(NA,dim(data21)[1]),rep(NA,dim(data23)[1])))

y.surv13 <- inla.surv(time = c(rep(NA,dim(data12)[1]), data13$Tstop,rep(NA,dim(data21)[1]), rep(NA,dim(data23)[1])),
                      event = c(rep(NA,dim(data12)[1]), data13$status,rep(NA,dim(data21)[1]), rep(NA,dim(data23)[1])),
                      truncation = c(rep(NA,dim(data12)[1]), data13$Tstart,rep(NA,dim(data21)[1]), rep(NA,dim(data23)[1])))

y.surv21 <- inla.surv(time = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]),data21$Tstop, rep(NA,dim(data23)[1])),
                      event = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]),data21$status, rep(NA,dim(data23)[1])),
                      truncation = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]),data21$Tstart, rep(NA,dim(data23)[1])))

y.surv23 <- inla.surv(time = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]),rep(NA,dim(data21)[1]), data23$Tstop),
                      event = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]),rep(NA,dim(data21)[1]), data23$status),
                      truncation = c(rep(NA,dim(data12)[1]), rep(NA,dim(data13)[1]),rep(NA,dim(data21)[1]), data23$Tstart))

y.joint <- list(y.surv12, y.surv13, y.surv21, y.surv23)

# ------------------------------------------------------------------------------------------
jointdata <- c(fixed.eff) # julie: should I include ID here??
jointdata$Y <- y.joint

#Model fit
formula <- Y ~ -1 + int12+int13+int21+int23 + 
  cov1_12+
  cov2_12+
  cov1_13+
  cov2_13+
  cov1_21+
  cov2_21+
  cov1_23 + 
  cov2_23

variant = 0
# fam_models: h12~weibull, h13~weibull, h21~exponential, h23~exponential
Jointmodel = inla(formula, family = fam_models,
                  control.family = list(list(variant = variant),list(variant = variant),list(),list()),
                  #control.predictor = list(list(prior = list()),list(),list(),list()),
                  #verbose=TRUE, 
                  data = jointdata,
                  control.compute=list(config = TRUE, dic=TRUE, waic=TRUE, mlik=TRUE))

#summary(Jointmodel)
s.f <- Jointmodel$summary.fixed
s.h <- Jointmodel$summary.hyperpar

# Jointmodel$dic$dic
# Jointmodel$waic$waic


library(ggplot2)
library('latex2exp')
# plotting estimation and conf int of fixed effects
nnn = c("beta_12(0)", "beta_13(0)", "beta_21(0)", "beta_23(0)")
nn =  paste(paste("beta", rep(c("12","13","21","23"),each = 2),sep="_"),
            paste("(",rep(c(1,2),3),")",sep = ""), sep = "")
hyp = c("alpha", "alpha_2")
names <- c(nnn, nn, hyp)
estimates <- rbind(Jointmodel$summary.fixed[c(1:12),c(1,3,5)], Jointmodel$summary.hyperpar[c(1:2), c(1,3,5)])
data.frame(true = c(int12,int13,int21,int23,beta12[1], beta12[2],beta13[1],beta13[2],beta21[1],beta21[2],beta23[1],beta23[2], alpha1, alpha2),
           estimates,
           i =  names)%>%
  ggplot() + geom_errorbar(aes(x  = i, ymin = X0.025quant, ymax = X0.975quant)) +
  geom_point(aes(i,true), color = "red") +
  xlab("") + ylab("") +
  scale_x_discrete(labels=c('alpha'=parse(text = TeX('$\\alpha_1$')),
                            'alpha_2'=parse(text = TeX('$\\alpha_2$')),
                            'beta_12(0)'=parse(text = TeX('$\\beta^{12}_0$')),
                            'beta_13(0)'=parse(text = TeX('$\\beta^{13}_0$')),
                            'beta_21(0)'=parse(text = TeX('$\\beta^{21}_0$')),
                            'beta_23(0)'=parse(text = TeX('$\\beta^{23}_0$')),
                            'beta_12(1)'=parse(text = TeX('$\\beta^{12}_1$')),
                            'beta_12(2)'=parse(text = TeX('$\\beta^{12}_2$')),
                            'beta_13(1)'=parse(text = TeX('$\\beta^{13}_1$')),
                            'beta_13(2)'=parse(text = TeX('$\\beta^{13}_2$')),
                            'beta_21(1)'=parse(text = TeX('$\\beta^{21}_1$')),
                            'beta_21(2)'=parse(text = TeX('$\\beta^{21}_2$')),
                            'beta_23(1)'=parse(text = TeX('$\\beta^{23}_1$')),
                            'beta_23(2)'=parse(text = TeX('$\\beta^{23}_2$'))))

# plotting estimated hazard functions using INLA
library(patchwork)
#only need to sample once
samples <- inla.posterior.sample(Jointmodel, n = 1000) # samples from "joint posterior"
#mm <- max(sdata$time) # not using this i think

# transition 1 --> 2 weibull
func <- function(...){ 
  t = seq(0,3.5,0.001)
  alpha = theta[1]
  lambda = exp(int12)
  alpha * t^(alpha-1) * lambda
}
aa <- inla.posterior.sample.eval(func, samples) 
t = seq(0,3.5,0.001)
#time = seq(0,mm,0.001)

# plotting with ggplot
df12 <- data.frame(time = t, 
                   post_mean = apply(aa,1,mean),
                   q1 = apply(aa,1,quantile, 0.025),
                   q2 = apply(aa,1,quantile, 0.975))
p12 <- ggplot() + geom_line(data = df12, aes(time, post_mean)) +
  geom_ribbon(data=df12, aes(time, ymin= q1, ymax = q2), alpha = 0.2) +
  geom_line(data=data.frame(time = t) %>%
              mutate(haz = alpha1 *t^(alpha1-1) * exp(int12)), aes(time, haz), color = "red") +
  ggtitle("state 1 - state 2") + ylab("Posterior mean")

# transition 1 --> 3 weibull
func <- function(...){ 
  t = seq(0,3.5,0.001)
  alpha = theta[2]
  lambda = exp(int13)
  alpha * t^(alpha-1) * lambda
}
aa <- inla.posterior.sample.eval(func, samples) 
aa[1,] <- 0
t = seq(0,3.5,0.001)
#time = seq(0,mm,0.001)
func0 <- function(...){ 
  alpha = theta[2]
  alpha
}
alpha=inla.posterior.sample.eval(func0, samples) 
# plotting with ggplot
df13 <- data.frame(time = t, 
                   post_mean = apply(aa,1,mean),
                   q1 = apply(aa,1,quantile, 0.025),
                   q2 = apply(aa,1,quantile, 0.975))
p13 <- ggplot() + geom_line(data = df13, aes(time, post_mean)) +
  geom_ribbon(data=df13, aes(time, ymin= q1, ymax = q2), alpha = 0.2) +
  geom_line(data=data.frame(time = t) %>%
              mutate(haz = alpha2 *t^(alpha2-1) * exp(int13)), aes(time, haz), color = "red") +
  ggtitle("state 1 - state 3") + ylab("Posterior mean")


# transition 2 --> 1 exponential
func <- function(...)
{
  exp(int21)
}
aa <- inla.posterior.sample.eval(func, samples) 
t = seq(0,1.5,0.001)

# plotting with ggplot:
df21 <- data.frame(time = t, 
                   post_mean = apply(aa,1,mean),
                   q1 = apply(aa,1,quantile, 0.025),
                   q2 = apply(aa,1,quantile, 0.975))
p21 <- ggplot() + geom_line(data = df21, aes(time, post_mean)) +
  geom_ribbon(data=df21, aes(time, ymin= q1, ymax = q2), alpha = 0.2) +
  geom_line(data=data.frame(time = t) %>%
              mutate(haz = exp(int21)), aes(time, haz), color = "red") +
  ggtitle("state 2 - state 1") + coord_cartesian(ylim = c(0,10)) + ylab("Posterior mean")


# transition 2 --> 3 exponential
func <- function(...)
{
  exp(int23)
}
aa <- inla.posterior.sample.eval(func, samples) 
t = seq(0,1.5,0.001)

# plotting with ggplot:
df23 <- data.frame(time = t, 
                   post_mean = apply(aa,1,mean),
                   q1 = apply(aa,1,quantile, 0.025),
                   q2 = apply(aa,1,quantile, 0.975))
p23 <- ggplot() + geom_line(data = df23, aes(time, post_mean)) +
  geom_ribbon(data=df23, aes(time, ymin= q1, ymax = q2), alpha = 0.2) +
  geom_line(data=data.frame(time = t) %>%
              mutate(haz = exp(int23)), aes(time, haz), color = "red") +
  ggtitle("state 2 - state 3") + coord_cartesian(ylim = c(0,10)) + ylab("Posterior mean")


p12 + p13 + p21 + p23 







