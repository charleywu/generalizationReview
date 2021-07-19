#Compare Performance of Diffusion Kernel to RBF kernel in a bandit context

rm(list=ls())
packages <- c('dplyr', "ComplexHeatmap", "ggplot2", "RColorBrewer", "igraph", 'matrixcalc', 'Matrix', 'lattice', 'jsonlite', 'viridis')
lapply(packages, require, character.only = TRUE)

source("models.R")
source('exportImportGraph.R')

#Aesthetics
palf <- colorRampPalette(c("white", "black")) 
######################################################################################################################################################
# Create some graph stimuli
######################################################################################################################################################
#Generate some lattice network with broken edges
numGraphs <- 40
dimensions <- c(8,8)
pruningFactor <- 0.4 #proportion of edges to delete

graphs <- list()
found <- 0
#rejection sampling to make sure that all generated networks are connected (i.e., without orphan nodes)
while (found < numGraphs){
  lattice <- make_lattice(dimvector =dimensions) #candidate
  delete_edges <- length(E(lattice)) * pruningFactor #how mny edges to delete
  pruned <- igraph::delete.edges(lattice, sample(seq(1:length(E(lattice))), delete_edges)) #randomly delete some edges
  #plot(lattice, vertex.size=5, vertex.label=NA) 
  if (igraph::components(pruned)$no==1){#Check that it only has a single connected component
    found <- found + 1
    graphs[[found]] <- pruned
  }
}
#plot example
g <- graphs[[1]]
plot(g, vertex.size=10, vertex.label=NA, margin=0)


#############################################################################
# Sample rewards from a GP-Diffusion Kernel prior
#############################################################################
library(ggnet)
#generating parameters
betaValue <- 2

#define node positions
gridLoc <- expand.grid(1:dimensions[1],1:dimensions[2])
positionScale <- 80
nodePositions <- data.frame(id = 1:prod(dimensions), x =gridLoc$Var1*positionScale, y =gridLoc$Var2*positionScale ) 

#output
graphJSON <- list()

#loop through graphs
for (i in 1:numGraphs){
  g <- graphs[[i]]
  #Sample environment from GP prior
  k_lattice <- diffusionKernel(g,beta = betaValue)
  #heatmap(k_lattice, col = palf(100)) #kernel matrix
  #Sample a reward distribution from the prior
  prior <- mvrnorm(1, rep(0, length(V(g))), k_lattice)
  V(g)$value = 1 #used for scaling
  #V(g)$reward = normalizeMinMax(prior, 1,rewardMax)
  V(g)$reward = prior #unormalized!
  #V(g)$label <- round(normalizeMinMax(prior, 1,rewardMax))
  V(g)$label <- round(prior, digits = 2)
  V(g)$color = cut(V(g)$reward, seq(-.5,0.5,length.out = 9)) #only used in R for discretization of color values
  V(g)$id <- seq(1:length(V(g)))
  V(g)$x <-  nodePositions$x-1
  V(g)$y <-  -nodePositions$y-1
  #p <- ggnet2(g,   mode = c("x", "y"), color = "color", palette='OrRd', label = 'label', label.size = 4, legend.position = "None", edge.color = 'black')
  #ggsave(paste0('environmentPlots/env', i, '.beta', betaValue,'.pdf'), p, width = 6, height = 6,  useDingbat=F)
  #out back into graphs
  graphs[[i]] <- g
  #Save graph as json
  object_json <- fromJSON(exportGraph(g))
  object_json$example <- c(i)
  graphJSON[[i]] <-object_json
}


#############################################################################
# Run model Comparison
#############################################################################
#Parameters
ucbBeta <- 0.5 #UCB beta
tau <-  0.1 #softmax tau
sigma <- 0.1 #SD of gaussian noise added to observations
X_star <- seq(1,64)
replications <- 100

lattice <- make_lattice(dimvector =dimensions) #make proxi

# Start the clock!
ptm <- proc.time()

#loop through graphs
for (i in 1:numGraphs){
  diffPerformance <- rep(0,25)
  rbfPerformance <- rep(0,25)
  g <- graphs[[i]]
  #ggnet2(g,   mode = c("x", "y"), color = "color", palette='OrRd', label = 'label', label.size = 4, legend.position = "None", edge.color = 'black')
  #Compute kernels
  k_diff <- diffusionKernel(g,beta = betaValue)
  k_rbf <-  diffusionKernel(lattice,beta = betaValue)
  for (rep in 1:replications){
    #sample initial value
    init <- sample(V(g),1)
    Xdiff <- Xrbf <- rep(0,25)
    Xdiff[1] <- Xrbf[1] <-  init
    Ydiff <- Yrbf <-  rep(0,25)
    Ydiff[1] <- Yrbf[1] <- init$reward+ rnorm(1,0,sigma)
    for (j in 2:25){ #loop through trials
      #Compute posteriors
      postDiff <- gpr(X.test = X_star, X = Xdiff[1:(j-1)], Y = Ydiff[1:(j-1)], k = k_diff, mu_0 = 0)
      postRBF <- gpr(X.test = X_star, X = Xrbf[1:(j-1)], Y = Yrbf[1:(j-1)], k = k_rbf, mu_0 = 0)
      #Utility Vecs
      Qdiff <- ucb(postDiff$mu, postDiff$var, ucbBeta)
      Qdiff <- Qdiff - max(Qdiff) #prevent overflow
      Qrbf <- ucb(postRBF$mu, postRBF$var, ucbBeta)
      Qrbf <- Qrbf - max(Qrbf) #prevent overflow
      #Choice probabilities
      pDiff <- exp(Qdiff/tau)
      pDiff <- pDiff/sum(pDiff)
      pRBF <- exp(Qrbf/tau)
      pRBF <- pRBF / sum(pRBF)
      #Make choice
      diffChoice <- sample(V(g),1, prob = pDiff, replace = TRUE)
      Xdiff[j] <-  diffChoice
      Ydiff[j] <- diffChoice$reward + rnorm(1,0,sigma)
      rbfChoice <- sample(V(g), 1, prob = pRBF, replace=TRUE)
      Xrbf[j] <-  rbfChoice
      Yrbf[j] <- rbfChoice$reward+ rnorm(1,0,sigma)
    }
    diffPerformance <- diffPerformance + Ydiff
    rbfPerformance <- rbfPerformance + Yrbf
  }
}

# Stop the clock
proc.time() - ptm

#average performance
diffPerformance <- diffPerformance / (numGraphs*replications)
rbfPerformance <- rbfPerformance / (numGraphs*replications)

simDF <- data.frame(trial = seq(1:25), performance = diffPerformance, model = 'Diffusion Kernel')
simDF <-rbind(simDF, data.frame(trial = seq(1:25), performance = rbfPerformance, model = 'RBF Kernel'))

ggplot(simDF, aes(trial, performance, color = model, fill = model))+
  geom_smooth(se = F)+
  theme_classic()
  
install.packages("beepr")
library(beepr)
beep()
  