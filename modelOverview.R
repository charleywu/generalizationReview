#Model overview
#Charley Wu September 2021
#house keeping
rm(list=ls())

packages <- c('plyr', 'jsonlite', 'reshape', 'DEoptim', "matrixcalc", "fields", 'RColorBrewer', 'ggplot2', 'gridExtra', 'cowplot', 'stats', 'stringr')
lapply(packages, require, character.only = TRUE)
source("models.R")

#Load envs
envs <- lapply(fromJSON("data/kernel2.json", flatten=TRUE),  FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,121), c('x2', 'y', 'x1'))))
choices <- expand.grid('x1'=0:10, 'x2'=0:10) #create index for each unique location to some integer from 1:121

#number formating function
numformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.0f", val)) }
##############################################################################################################
#Observations
##############################################################################################################
envNum <- 3 #fix
trial <- 20 #Which trial is currently being predicted?
scalingFactor <- 85
#Create a set of randomn choices. In the future, replace with either handpicked observations, participant data, or sampled choices under a policy
chosen <- sample(1:121, trial) 
x1 <- choices[chosen,'x1']
x2 <- choices[chosen,'x2']
#create observation matrix
X<-as.matrix(cbind(x1,x2))
#initialize Xtest
Xnew<-as.matrix(expand.grid(0:10,0:10))
#make sure X is a matrix
X<-as.matrix(X)
y <- envs[[envNum]][chosen,'y'] + rnorm(trial, 0, sd = sqrt(0.001)) #TODO: check that this observation variance is consistent with the noise added in the experiment w.r.t. to reward scaling
y[y<.01] <- .01
y[y>1] <- 1
yrescaled <- y*scalingFactor #Think of this as canonical scaling, where in the experiment rewards are scaled to randomly sampled max values


#Plot 1. Subject's search space, which has been revealed until the specific trial number defined by "trial"
#construct dataframe
d1 <- data.frame(cbind(Xnew,rep(NA, 121))) # fill all 121 spots with NA
colnames(d1) <- c("x1", "x2", "y")
#replace NA with observations
for (row in seq(1:length(y))){ 
  obs <- X[row,]
  d1$y[d1$x1==obs[1] & d1$x2==obs[2]] <- yrescaled[row]
}
d1$ytext <- numformat(d1$y)
d1$ytext[d1$ytext=='NA'] <- ""
p1 <- ggplot(d1, aes(x = x1+1, y = x2+1, fill = y)) + 
  geom_tile(color='black', width=1, height=1) +
  theme_bw() +
  coord_equal() +
  xlim(0.5,11.5) +
  ylim(0.5,11.5) + 
  geom_text(aes(x = x1+1, y = x2+1, label = ytext))+
  #ggtitle(paste0('Revealed (t=',trial,")")) +
  #scale_fill_gradientn(name='Payoff', colours = hm.palette(100), values = seq(0, 100, length=9),  rescaler = function(x, ...) x, oob = identity) +
  scale_fill_distiller(palette = "Spectral",limits=c(-.05*scalingFactor,1.05*scalingFactor), na.value = 'white', breaks=c(0,.25*scalingFactor,.50*scalingFactor,.75*scalingFactor,1*scalingFactor))+
  labs(fill="Payoff")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = 'none')
p1

##############################################################################################################
#GP predictions
##############################################################################################################

#Gaussian Process function
#X.test: matrix for predcitions
#theta: vector of hyper-parameters (lambda, Sf, Sn)
#X; matrix of observations
#y: vector of observed outcomes
#kernel: used kernel function, can be "rbf", "oru", or "mat"
gpr <- function(X.test, theta, X,Y, k){
  #make it a matrix
  Xstar <- as.matrix(X.test)
  #dimensions
  d <- ncol(as.matrix(X))
  #calculate capital K
  K <- k(X,X,theta) 
  #Check if matrix is positive semi-definite
  if (is.positive.definite(K)){
    #KK <- cov.inverse.chol(K) #use Cholesky
    KK.inv <- chol2inv(chol(K)) #MASS implementation of Cholesky
  } else {
    KK.inv <- cov.inverse.svd(K) #use SVD
  }
  #times y
  Ky <- KK.inv %*% Y
  #apply the kernel
  result <- apply(Xstar, 1, function(x){
    XX <- matrix(x,nrow=1) 
    Kstar <- k(X, XX, theta)
    Kstarstar <- k(XX,XX,theta)
    #get mean vector
    mu <- t(Kstar) %*% Ky
    #get covariance
    cv <- Kstarstar - (t(Kstar) %*% KK.inv %*% Kstar) #BUG: sometimes cv<0, leading to NaN when we return sqrt(cv)
    #DEBUG
    if (cv<0){ cv <- abs(cv)} #TEMPORARY SOLUTION: MANUALLY SET CV TO BE POSITIVE IF NEGATIVE
    #return means and variance
    return(c(mu, cv))
  })
  #as a data frame with names mu and sig
  prediction <- as.data.frame(t(result))
  prediction[is.na(prediction)] <- 0.01 #remove NaN items with the noise variance 
  colnames(prediction) <- c("mu", "sig")
  return(prediction)
}


#Radial Basis Kernel
rbf <- function(X1,X2,theta){
  #transfer to matrices
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  #check dimensions
  if(ncol(X1) != ncol(X2)){
    stop("X1 and X2 must contain input values of the same dimension.")
  } else if(!all(theta>=0)){
    stop("All parameters must be >= 0.")
  }
  #get dimensions
  N1 <- nrow(X1)
  N2 <- nrow(X2)
  d <- ncol(X1)
  #initialize sigma
  sigma <-  matrix(rep(0, N1*N2),nrow=N1)
  #observational variance
  sf <- theta[d+1]
  #noise variance
  sn <- theta[d+2]
  #loop through
  for(i in 1:d){
    #length scale
    l <- theta[i] #Note: assumes a unique length scale for each dimension
    #x-diff
    xdiff <- (outer(X1[,i],X2[,i],function(x,y) x - y)/l)^2
    sigma <- sigma + xdiff
  }
  #RBF function
  if(identical(X1,X2)){
    id <- diag(rep(1,N1))
    sigma.final <- sf*exp(-0.5*sigma) + sn*id
  } else {
    sigma.final <- sf*exp(-0.5*sigma)
  }
  #return final covariance matrix
  return(sigma.final)
}
class(rbf)<- c(class(rbf), "GP") #identify the rbf kernel as a gp model

#Model specifications
kernel=rbf
acq = ucb
lambda <- .78
signalVariance <- 1
errorVariance <- 0.001
beta <- .5
tau <- 0.09

#Compute posterior predictions
GPpost <- gpr(X.test=Xnew, theta=c(lambda, lambda, signalVariance, errorVariance), X=X, Y=y-.5, k=kernel)
GPpost$mu <- GPpost$mu + 0.5

#Plot 2. Posterior mean
d2 <- melt(matrix((GPpost[,1]*scalingFactor), nrow=11, ncol=11)) #*100 is to preserve same scaling factor
names(d2) <- c('X1', 'X2', 'value')
p2<-ggplot(d2, aes(x = X1, y = X2, fill = value)) + 
  geom_tile(color='black', width=1, height=1) +
  theme_bw() +
  xlim(0.5,11.5) +
  ylim(0.5,11.5) + 
  coord_equal() +
  #scale_fill_gradientn(name = "Exp. Payoff", colours = hm.palette(100),values = seq(0, 100, length=9)) +
  scale_fill_distiller(palette = "Spectral",limits=c(-.05*scalingFactor,1.05*scalingFactor), na.value = 'white', breaks=c(0,25,50,75,100))+
  labs(fill="Payoff")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
p2



#Plot 3. Posterior variance
d3 <- melt(matrix((GPpost[,2]*scalingFactor), nrow=11, ncol=11))  #*100 is to preserve same scaling factor
names(d3) <- c('X1', 'X2', 'value')
p3<-ggplot(d3, aes(x = X1, y = X2, fill = value)) + 
  geom_tile(color='black', width=1, height=1) +
  theme_bw() +
  xlim(0.5,11.5) +
  ylim(0.5,11.5) + 
  coord_equal() +
  #ggtitle('Uncertainty') +
  #scale_fill_gradientn(name = "Exp. Payoff", colours = hm.palette(100),values = seq(0, 100, length=9)) +
  scale_fill_distiller(palette = "Spectral",limits=c(-.05*scalingFactor,1.05*scalingFactor), na.value = 'white', breaks=c(0,25,50,75,100))+
  labs(fill="Variance")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
p3

##############################################################################################################
#UCB
##############################################################################################################

#Compute UCB values
utilityVec <- acq(GPpost$mu, GPpost$sig, beta) #

#Plot 4. UCB value
d4 <- melt(matrix((utilityVec  * scalingFactor), nrow=11, ncol=11))  
names(d4) <- c('X1', 'X2', 'value')
p4<-ggplot(d4, aes(x = X1, y = X2, fill = value)) + 
  geom_tile(color='black', width=1, height=1) +
  theme_bw() +
  xlim(0.5,11.5) +
  ylim(0.5,11.5) + 
  coord_equal() +
  #ggtitle('UCB') +
  #scale_fill_gradientn(name = "Exp. Payoff", colours = hm.palette(100),values = seq(0, 100, length=9)) +
  scale_fill_distiller(palette = "Spectral")+
  labs(fill="UCB")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = 'none')
p4


##############################################################################################################
#Softmax
##############################################################################################################

p <- utilityVec - max(utilityVec) #subtract max value to prevent overflow
p <- exp(p/tau)
p <- p/sum(p)
#avoid underflow by setting a floor and a ceiling
p <- (pmax(p, 0.00001))
p <- (pmin(p, 0.99999))
#Next choice at t = trial + 1
print(sum(p))


#Plot 3. Softmax surface of acquisition function
d5<- melt(matrix(p, nrow=11, ncol=11))
names(d5) <- c('X1', 'X2', 'value')
d5<- round(d5, 3)
p5<-ggplot(d5, aes(x = X1, y = X2, fill = value)) + 
  geom_tile(color='black', width=1, height=1) +
  theme_bw() +
  xlim(0.5,11.5) +
  ylim(0.5,11.5) + 
  coord_equal() +
  #ggtitle('Softmax') +
  #scale_fill_gradientn(name='P(choice)',colours = hm.palette(100),values = seq(0, 1, length=9)) +
  scale_fill_distiller(palette = "Spectral", na.value = 'white' )+
  labs(fill="P(Choice)")+
  #annotate("text", x = nextChoice[1] + 1, y = nextChoice[2] + 1, label = "X") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
p5

plots <- plot_grid(p1,p2,p3,p5, nrow = 1)
plots

ggsave('plots/modelOverview.pdf', width = 10, height = 4, units = 'in')
