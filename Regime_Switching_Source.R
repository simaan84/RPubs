## ----warning=F,message=F-------------------------------------------------

# This code was extracted from the Rmd file published at http://rpubs.com/simaan84/regime_switching
library(knitr)
library(kableExtra)
library(dplyr)
theta_v <- data.frame(t(c(2.00,-2.00,1.00,2.00,0.95,0.85)))
names(theta_v) <- c("$\\mu_1$","$\\mu_2$","$\\sigma_1$","$\\sigma_2$","$p_{11}$","$p_{22}$")
kable(theta_v, "html", booktabs = F,escape = T) %>% 
        kable_styling(position = "center")

## ------------------------------------------------------------------------
p11 <- theta_v[1,5]	
p22 <- theta_v[1,6]
P <- matrix(c(p11,1-p22,1-p11,p22),2,2)
P[1,]

## ----fig.align="center"--------------------------------------------------
set.seed(13)
T_end <- 10^2

s0 <- 1 
st <- function(i) sample(1:2,1,prob = P[i,])

s <- st(s0)
for(t in 2:T_end) {
  s <- c(s,st(s[t-1]))
}
plot(s, pch = 20,cex = 0.5)

## ------------------------------------------------------------------------
P_stat <- Reduce(function(M1,M2) M1%*%M2 ,lapply(1:100,function(x) P ))
P_stat[1,]

## ------------------------------------------------------------------------
mean(s==1)

## ----fig.align="center"--------------------------------------------------
set.seed(11)
x1 <- rnorm(T_end,theta_v[1,1],theta_v[1,3])
set.seed(17)
x2 <- rnorm(T_end,theta_v[1,2],theta_v[1,4])

t_index <- 1:T_end
x <- rep(0,T_end)
x[s==1] <- x1[s==1]
x[s==2] <- x2[s==2]
plot(x~t_index, pch = 20)
points(x[s == 2]~t_index[s==2],col = 2)

## ------------------------------------------------------------------------
HMM_Lik <- function(theta) {
  mu <- theta[1:2] # means
  sig <- theta[3:4] # volatilities
  p <- t(matrix(c(theta[5],1-theta[6],1-theta[5],theta[6]),2,2)) # transition matrix
  xi <- c(0.5,0.5) # prior about the filter
  
  Xi <- numeric()
  L <- numeric()
  for(t in 1:T_end) {
    phi1 <- dnorm((x[t]-mu[1])/sig[1])/sig[1]
    phi2 <- dnorm((x[t]-mu[2])/sig[2])/sig[2]	
    phi <- c(phi1,phi2)
    xi_next <- t(p)%*%xi
    l_t <- t(phi)%*%xi_next
    L <- c(L,l_t)
    # update the filter for next period
    xi <- c((1/l_t))*(phi*xi_next) # inference prob
    Xi <- rbind(Xi,t(xi))
  }
  LL <- sum(log(L))
  list(Xi = Xi,LL = LL)
}

## ------------------------------------------------------------------------
theta_known <- unlist(theta_v[1,])
Filter <- HMM_Lik(theta_known)$Xi
Filter <- cbind(t_index,Filter)
colnames(Filter) <- c("$t$","$\\xi_{t \\mid t, 1}$","$\\xi_{t \\mid t, 2}$")
rownames(Filter) <- NULL
kable(round(head(Filter),3), "html", booktabs = F,escape = F) %>% 
        kable_styling(position = "center")


## ------------------------------------------------------------------------
all(round(apply(Filter[,-1],1,sum),9) == 1)

## ----fig.align="center"--------------------------------------------------
plot(Filter[,3]~t_index, type = "l", ylab = expression(xi[2]))
points(Filter[s==2,3]~t_index[s==2],pch = 20, col = 2)

## ------------------------------------------------------------------------
m1 <- m2 <- mean(x,na.rm = T)
s1 <- s2 <- sd(x,na.rm = T)
p1 <- p2 <- 0.5
theta0 <- c(m1,m2,s1,s2,p1,p2)

## ------------------------------------------------------------------------
# set the constraints
ui <- matrix(0,4,6)
ui[,5] <- c(1,0,-1,0)
ui[,6] <- c(0,1,0,-1)

# stack in matrix and vector
A <- t(matrix(ui,ncol= 4))
A <- ui
B <- c(0.01,0.01,-0.99,-0.99)

## ------------------------------------------------------------------------
all(A%*%theta0 >= B)

## ------------------------------------------------------------------------
L.lik <- function(theta) -HMM_Lik(theta)$LL
opt <- constrOptim(theta0,L.lik, NULL,ui = A, ci = B )
opt

## ----fig.align="center"--------------------------------------------------
plot(opt$par ~ theta_known,pch = 20,cex=2,ylab="MLE",xlab = "True")
abline(a=0,b=1,lty=2)

## ------------------------------------------------------------------------
mod <- lm(x~1)
mod_est <- data.frame(mu = mod$coefficients, sig = summary(mod)$sig)
names(mod_est) <- c("$\\hat{\\mu}$","$\\hat{\\sigma}$")
rownames(mod_est) <- NULL
kable(mod_est, "html", booktabs = F,escape = F) %>% 
        kable_styling(position = "center")


## ------------------------------------------------------------------------
EX <- 0.75*2 + 0.25*-2
EX

## ------------------------------------------------------------------------
EX2 <- (2^2 + 1^2)*0.75 + ((-2)^2 + 2^2)*0.25
VX <- EX2 - EX^2
sqrt(VX)

## ----warning=F,message=F-------------------------------------------------
library(MSwM)
mod.mswm <- msmFit(mod,k=2,p=0,sw=c(TRUE,TRUE),control=list(parallel=TRUE))

## ------------------------------------------------------------------------
mod.mswm

## ----fig.align="center"--------------------------------------------------
theta_mswm <- c(unlist(mod.mswm@Coef),mod.mswm@std,diag(mod.mswm@transMat))
plot(opt$par ~ theta_known,pch = 20,cex=2,ylab="MLE",xlab = "True")
points(theta_mswm~theta_known,pch = 1,col = 2, cex = 2,lwd = 2)
abline(a=0,b=1,lty=2)
legend("topleft",c("Manual","MSwM"), pch = c(20,1), col = 1:2)

## ----fig.align="center"--------------------------------------------------
par(mar = 2*c(1,1,1,1),mfrow = c(2,1))
plotProb(mod.mswm,2)

## ----fig.align="center"--------------------------------------------------
Filter <- HMM_Lik(opt$par)$Xi
Filter <- cbind(t_index,Filter)
Filter <- data.frame(Filter)
Filter$Regime_1 <- (Filter[,2]>=0.5)*1
xx <- t_index[Filter$Regime_1 == 1] # plot regimes
plot(x~t_index,type ="l",col = 0,xlim=c(1,100))
rect(xx-1,-10,xx,10,col = "lightgray",lty = 0)
lines(x~t_index)
points(x[s==2] ~ t_index[s==2],col = 1,pch = 20)

## ------------------------------------------------------------------------
mean(Filter$Regime_1 == (s==1)*1)

