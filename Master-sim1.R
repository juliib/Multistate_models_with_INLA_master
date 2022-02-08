set.seed(5)
# function to simulate event times ----------------------------------------
#' generate.my.times
#'
#'
#' @param n The number of patients to simulate
#' @param max.int The longest observable time
#' @param temp The cumulative all-cause hazard + y
#' @param t0 The time one wants to start simulating from
#'
#' @return Vector of event times
#' @export
#'
#' @examples
generate.my.times <- function(n, max.int, temp) { # gives n independent replicates of T
  s_time <- NULL
  i <- 1
  t = 0
  while(length(s_time) < n) {
    while(t==0)
    {
      u <- runif(1)
      ## If endpoints are of opposite sign, call uniroot:
      if (temp(0,log(1-u))*temp(max.int,log(1-u)) < 0) 
      {
        res = uniroot(temp, c(0,max.int), tol = 0.0001, y = log(1-u))
        t = res$root
      }
    }
    s_time[i] = t
    i = i+1
  }
  
  return(s_time)
}


# set up param ------------------------------------------------------------

N=300 #number of patients

tt = seq(0,5,0.01)

# covariates
X = cbind(rnorm(N), sample(c(0,1), N, replace = T))
beta01 = c(2,0.4)
beta02 = c(0.6,1)

# intercepts
int1 <- 1
int2 <- 2

# linear predictor
eta1 <- int1 + X %*% beta01
eta2 <- int2 + X %*% beta02

# hyperparameters
alpha1  = 1 # makes baseline haz = 1 -> exponential distributed
alpha2 = 1.2


# baseline hazards + plot
h01_0 <- function(t) alpha1 * t^{alpha1 - 1} * exp(int1)
h02_0 <- function(t) alpha2 * t^{alpha2 - 1} * exp(int2)

data.frame(t = tt,
           h01_0  = h01_0 (tt),
           h02_0  = h02_0 (tt)
) %>%
  pivot_longer(-t) %>%
  ggplot() + geom_line(aes(t,value, group= name, color = name)) 

# cumulative baseline hazards + plot
H01_0 <- function(t) t^alpha1 * exp(int1)
H02_0 <- function(t) t^alpha2 * exp(int2)

data.frame(t = tt,
           H01_0  = H01_0 (tt),
           H02_0  = H02_0 (tt)
) %>%
  pivot_longer(-t) %>%
  ggplot() + geom_line(aes(t,value, group= name, color = name))

# TRANSITION VECTOR FOR THE MODEL
tvec <- matrix(FALSE, 1,3)
colnames(tvec) <- c("0","1","2")
rownames(tvec) <- "1" 
dimnames(tvec) <- list(from = "0", 
                       to = c("0","1", "2"))
tvec[1,c(2,3)] = TRUE

#tvec

sdata = data.frame(ID = 1:N, x1 = X[,1], x2 = X[,2], time = NA, event = NA)
for(i in 1:N){
  alpha01 = function(t)
    alpha1 * t^(alpha1-1) *exp(eta1[i])
  alpha02 = function(t)
    alpha2 * t^(alpha2-1) * exp(eta2[i])
  alpha0 = function(t)
    alpha01(t) + alpha02(t)
  temp = function(t, y)
    exp(eta1[i]) * t^(alpha1) + exp(eta2[i]) *t^(alpha2) + y
  
  sim_t = generate.my.times(1, 5, temp) # sim 1 time on max interval 5
  prob_1 = alpha01(sim_t)/(alpha01(sim_t) +  alpha02(sim_t))
  event = sample(c(1,2), 1, prob = c(prob_1, 1-prob_1))
  sdata$time[i] = sim_t
  sdata$event[i] = event
}
#table(sdata$event)
#head(sdata)
#max_time <- max(sdata$time)
#sdata$time <- sdata$time/max_time

ggplot(sdata, aes(x=sdata$time)) + geom_histogram(bins = 50) + xlab("Time") + ggtitle("Histogram of simulated event times")

sdata1 = sdata %>% 
  mutate(status = ifelse(event==1,1,0))

sdata2 = sdata %>% 
  mutate(status = ifelse(event==2,1,0))


fixed.eff <- data.frame(
  int1 = rep(c(1,0),each = N),
  int2 = rep(c(0,1),each = N),
  X1_01=c(sdata$x1, rep(0,N)),
  X1_02=c(rep(0,N),sdata$x1),
  X2_01=c(sdata$x2, rep(0,N)),
  X2_02=c(rep(0,N),sdata$x2))



library(INLA)
y.surv1 <- inla.surv(time = c(sdata1$time,rep(NA,N)), 
                     event = c(sdata1$status,rep(NA,N)))
y.surv2 <- inla.surv(time = c(rep(NA,N),sdata2$time), 
                     event = c(rep(NA,N),sdata2$status))
y.joint<-list(y.surv1,y.surv2)
jointdata <- c(fixed.eff)
jointdata$Y <- y.joint 

#Model fit
formula=Y~-1 + int1 + int2 + X1_01 +  X2_01 + X1_02 + X2_02 

Jointmodel= inla(formula, 
                 family = c("exponentialsurv","weibullsurv"),
                 data = jointdata,
                 control.compute=list(config = TRUE, dic=TRUE, waic=TRUE, mlik=TRUE))
#summary(Jointmodel)
s.f <- Jointmodel$summary.fixed
s.h <- Jointmodel$summary.hyperpar

# plotting estimation and conf int 
library(ggplot2)
library('latex2exp')
nnn = c("beta_01(0)", "beta_02(0)")
nn =  c("beta_01(1)", "beta_01(2)", "beta_02(1)", "beta_02(2)")
hyp = c("alpha")
names <- c(nnn, nn, hyp)
estimates <- rbind(Jointmodel$summary.fixed[c(1:6),c(1,3,5)], Jointmodel$summary.hyperpar[c(1:1), c(1,3,5)])
data.frame(true = c(int1,int2,beta01[1],beta01[2],beta02[1],beta02[2],alpha2),
           estimates,
           i =  names)%>%
  ggplot() + geom_errorbar(aes(x  = i, ymin = X0.025quant, ymax = X0.975quant)) +
  geom_point(aes(i,true), color = "red") +
  xlab("") + ylab("") +
  scale_x_discrete(labels=c('alpha'=parse(text = TeX('$\\alpha$')),
                            'beta_01(0)'=parse(text = TeX('$\\beta^{01}_0$')),
                            'beta_02(0)'=parse(text = TeX('$\\beta^{02}_0$')),
                            'beta_01(1)'=parse(text = TeX('$\\beta^{01}_1$')),
                            'beta_01(2)'=parse(text = TeX('$\\beta^{01}_2$')),
                            'beta_02(1)'=parse(text = TeX('$\\beta^{02}_1$')),
                            'beta_02(2)'=parse(text = TeX('$\\beta^{02}_2$'))))

library(patchwork)
# we only sample once
samples <- inla.posterior.sample(Jointmodel, n = 1000) # samples from "joint posterior"
# transition 0 --> 2
func <- function(...)
{
  t = seq(0,1.5,0.001) # time needs to be defined here
  alpha = theta[1]
  lambda= exp(int2)
  alpha * t^(alpha-1) * lambda
}
aa <- inla.posterior.sample.eval(func, samples) 
aa[1,] <- 0
t = seq(0,1.5,0.001) # time also needs to be defined outside the function

# plotting with ggplot
df02 <- data.frame(time = t, 
                   post_mean = apply(aa,1,mean),
                   q1 = apply(aa,1,quantile, 0.025),
                   q2 = apply(aa,1,quantile, 0.975))
p02 <- ggplot() + geom_line(data = df02, aes(time, post_mean)) +
  geom_ribbon(data=df02, aes(time, ymin= q1, ymax = q2), alpha = 0.2) +
  geom_line(data=data.frame(time = t) %>%
              mutate(haz = alpha2 *t^(alpha2-1) * exp(int2)), aes(time, haz), color = "red") +
  ggtitle("state 0 - state 2") + ylab("Posterior mean")


# transition 0 --> 1
func <- function(...)
{
  exp(int1)
}
aa <- inla.posterior.sample.eval(func, samples) 
t = seq(0,1.5,0.001)

# plotting with ggplot:
df01 <- data.frame(time = t, 
                   post_mean = apply(aa,1,mean),
                   q1 = apply(aa,1,quantile, 0.025),
                   q2 = apply(aa,1,quantile, 0.975))
p01 <- ggplot() + geom_line(data = df01, aes(time, post_mean)) +
  geom_ribbon(data=df01, aes(time, ymin= q1, ymax = q2), alpha = 0.2) +
  geom_line(data=data.frame(time = t) %>%
              mutate(haz = exp(int1)), aes(time, haz), color = "red") +
  ggtitle("state 0 - state 1") + coord_cartesian(ylim = c(0,10)) + ylab("Posterior mean")


p01 + p02





