#Graph regularization operators
# Smola & Kondor 2003

rm(list=ls())

packages <- c('ggplot2', 'tidyr','dplyr', "MASS",'jsonlite', 'matrixcalc', 'Matrix', 'igraph', 'expm','lattice', 'cowplot', 'viridis')
lapply(packages, require, character.only = TRUE)
lattice.options(default.theme = standard.theme(color = FALSE))
source('utilities.R')
source('models.R')


######################################################################
#Define some graphs
######################################################################
experimentGraphFiles <- 'data/exp1/graphs.json'
networks <- fromJSON(experimentGraphFiles)
edges <- networks[[1]]
nodes <- networks[[2]]

graphList = list()

for (i in 1:length(edges)){
  g <- graph_from_edgelist(as.matrix(edges[[i]]), directed=F) 
  V(g)$reward  <- nodes[[i]]$reward
  #plot(g, vertex.size = V(g)$reward)
  graphList[[i]] <-g
}

#alternative: Only look at lattice graphs
#g <- make_lattice(dimvector =c(8,8)) #Lattice graph

#Define a single environment to use for analyses
g <- graphList[[1]]
n_nodes <- length(V(g))
plot(g)
######################################################################
#Define graph primitives
######################################################################
#degree matrix
D <- diag(strength(g)) 
#Weight matrix
W <- as_adjacency_matrix(g, sparse=F) #unweighted
#W <- as_adjacency_matrix(g, attr="weight", sparse=F) #for weighted graphs

#Construct laplacian
L <- D - W
L <- sqrt(solve(D))%*%L%*%sqrt(solve(D)) #if degree normalized Laplacian
levelplot(L,col.regions = gray(0:100/100))
######################################################################
#Construct kernels
######################################################################

#1. Regularized Laplacian (used in Kemp & Tenenbaum, 2008, 2009)
sigma <- 2
K_rl <- solve(diag(n_nodes) + sigma^2 * L)

levelplot(K_rl,col.regions = gray(0:100/100))

#2. Diffusion kernel (defined in terms of Kondor & Smola, 2003)
alpha <- sigma^2/2
eig <- eigen(-alpha*L) #Eigenvalue decomposition
K_df <- eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors) #+ diag(rep(noise, length(V(network))))
K_df <- round(K_df, 10) #necessary since loss of precision can make other functions not recognize K as positive semi-definite 
levelplot(K_df,col.regions = gray(0:100/100))
#check that this corresponds to the alternative implementation in Kondor & Lafferty (2002)
K_alt <- diffusionKernel(g, alpha = alpha, normalizedLaplacian=T)
levelplot(K_alt,col.regions = gray(0:100/100))
all.equal(K_df, K_alt) #Holds even for non-lattice graphs

#3. p-step random walk
a <- 2
p <- 1
K_rw <- matrix.power(a*diag(n_nodes) - L, p)
levelplot(K_rw,col.regions = gray(0:100/100))


######################################################################
#Compute correlation plots
######################################################################

replications <- 10000
sigmaList <- c(1,2,4,8)
simDF <- data.frame()

#compute correlation plots for regularized laplacian kernel and diffusion kernel
for (i in length(graphList)){
  g <- graphList[[i]]
  #Group pairs of nodes by distance
  distanceTable <- distances(g) #distances between each node
  distanceTable[lower.tri(distanceTable)] <- NA #set lower symmetric triangle to NA
  uniqueDistances <- unique(c(distanceTable[!is.na(distanceTable)]))
  #Graph primatives
  D <- diag(strength(g)) #degree matrix
  W <- as_adjacency_matrix(g, sparse=F) #unweighted adjacency matrix
  #W <- as_adjacency_matrix(g, attr="weight", sparse=F) #for weighted graphs
  L <- D - W #Construct laplacian
  L <- sqrt(solve(D))%*%L%*%sqrt(solve(D)) #if degree normalized Laplacian
  for (sigmaSquared in sigmaList){
    alpha <- sigmaSquared/2 #define alpha in terms of sigma
    K_rl <- solve(diag(n_nodes) + sigmaSquared * L)
    K_df <-diffusionKernel(g,alpha = alpha, normalizedLaplacian = T) 
    prior_rl <- mvrnorm(replications, rep(0, length(V(g))), K_rl) #sample from the prior
    prior_df <- mvrnorm(replications, rep(0, length(V(g))), K_df) #sample from the prior
    for (d in uniqueDistances){
      pairs <- which(distanceTable==d, arr.ind=T)
      #Regularized Laplacian
      valuesA <- prior_rl[,pairs[,1]]#rewards
      valuesB <- prior_rl[,pairs[,2]]
      dat <- data.frame(graph =  i, kernel = 'RL', sigmaSquared = sigmaSquared, alpha = NA, A = NA, p = NA, distance = d, correlation= cor(c(valuesA),c(valuesB)))
      simDF<- rbind(simDF,dat)
      #Diffusion Kernel
      valuesA <- prior_df[,pairs[,1]]#rewards
      valuesB <- prior_df[,pairs[,2]]
      dat <- data.frame(graph =  i, kernel = 'DF', sigmaSquared = sigmaSquared, alpha = alpha,A = NA, p = NA, distance = d, correlation= cor(c(valuesA),c(valuesB)))
      simDF<- rbind(simDF,dat)
    }
  }
}

plotsimDF <- simDF %>% group_by(kernel,sigmaSquared, distance) %>% summarize(Correlation = mean(correlation))
plotsimDF$sigmaSquared <- factor(plotsimDF$sigmaSquared)

#Regularized Laplacian
p1 <- ggplot(subset(plotsimDF, kernel == 'RL'), aes(x=distance, y = Correlation, color = sigmaSquared, shape = sigmaSquared))+
  geom_line(size = 1) +
  geom_point(size = 2) +
  ylab("Pearson Correlation")+xlab("Distance")+
  theme_classic()+
  xlim(c(0,5))+
  #scale_x_continuous(breaks = round(seq(min(plotsimDF$distance), max(plotsimDF$distance), by = 1),1))+
  scale_color_brewer(palette="Set1", direction=-1, name='sigma^2')+
  scale_shape_manual(values = c(15,16,17,18),name = 'sigma^2')+
  #scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(text = element_text(size=14,  family="sans"), legend.justification = c(1, 1), legend.position = c(1, 1))+
  theme(strip.background=element_blank(), legend.key=element_rect(color=NA)) +
  ggtitle("Regularized Laplacian")
p1

p2 <- ggplot(subset(plotsimDF, kernel == 'DF'), aes(x=distance, y = Correlation, color = sigmaSquared, shape = sigmaSquared))+
  geom_line(size = 1) +
  geom_point(size = 2) +
  ylab("Pearson Correlation")+xlab("Distance")+
  theme_classic()+
  xlim(c(0,5))+
  #scale_x_continuous(breaks = round(seq(min(plotsimDF$distance), max(plotsimDF$distance), by = 1),1))+
  scale_color_brewer(palette="Set1", direction=-1, name='alpha', labels = sigmaList/2)+
  scale_shape_manual(values = c(15,16,17,18),name = 'alpha', labels = sigmaList/2)+
  #scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(text = element_text(size=14,  family="sans"),  legend.justification = c(1, 1),legend.position = c(1, 1))+
  theme(strip.background=element_blank(), legend.key=element_rect(color=NA)) +
  ggtitle("Diffusion Kernel")
p2


#Now compute p-step random walk kernel
aList <- c(2)
pList <- c(1,2,3,4,5,6)

#compute correlation plots for random walk kernel
for (i in length(graphList)){
  g <- graphList[[i]]
  #Group pairs of nodes by distance
  distanceTable <- distances(g) #distances between each node
  distanceTable[lower.tri(distanceTable)] <- NA #set lower symmetric triangle to NA
  uniqueDistances <- unique(c(distanceTable[!is.na(distanceTable)]))
  #Graph primatives
  D <- diag(strength(g)) #degree matrix
  W <- as_adjacency_matrix(g, sparse=F) #unweighted adjacency matrix
  #W <- as_adjacency_matrix(g, attr="weight", sparse=F) #for weighted graphs
  L <- D - W #Construct laplacian
  L <- sqrt(solve(D))%*%L%*%sqrt(solve(D)) #if degree normalized Laplacian
  for (A in aList){
    for (p in pList){
      K_rw <- matrix.power(A*diag(n_nodes) - L, p)
      prior_rw <- mvrnorm(replications, rep(0, length(V(g))), K_rw)
      #Group pairs of nodes by distance
      distanceTable <- distances(g) #distances between each node
      distanceTable[lower.tri(distanceTable)] <- NA #set lower symmetric triangle to NA
      uniqueDistances <- unique(c(distanceTable[!is.na(distanceTable)]))
      for (d in uniqueDistances){
        pairs <- which(distanceTable==d, arr.ind=T)
        #Regularized Laplacian
        valuesA <- prior_rw[,pairs[,1]]#rewards
        valuesB <- prior_rw[,pairs[,2]]
        dat <- data.frame(graph=i, kernel = 'RW', sigmaSquared = NA, alpha = NA, A = A, p = p, distance = d, correlation= cor(c(valuesA),c(valuesB)))
        simDF<- rbind(simDF,dat)
      }
    }
  }
}

plotsimDF <- subset(simDF, kernel=='RW') %>% group_by(kernel,A,p, distance) %>% summarize(Correlation = mean(correlation))
plotsimDF$A <- factor(plotsimDF$A)
plotsimDF$p <- factor(plotsimDF$p)

p3 <- ggplot(plotsimDF, aes(x=distance, y = Correlation, color = p, shape = p))+
  geom_line(size = 1) +
  geom_point(size = 2) +
  ylab("Pearson Correlation")+xlab("Distance")+
  theme_classic()+
  xlim(c(0,5))+
  #scale_x_continuous(breaks = round(seq(min(plotsimDF$distance), max(plotsimDF$distance), by = 1),1))+
  scale_color_viridis_d( direction=-1, name = 'Steps')+
  scale_shape_manual(values = c(0,1,2,3,4,5,6), name = 'Steps')+
  #scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(text = element_text(size=14,  family="sans"), legend.justification = c(1, 1), legend.position = c(1, 1))+
  theme(strip.background=element_blank(), legend.key=element_rect(color=NA))+
  ggtitle('Random walk')
p3


p <- plot_grid(p1,p2,p3, nrow=1)
p
ggsave(p, filename='plots/correlationDecayNormalizedLaplacian.pdf', width = 10, height = 3, units = 'in', useDingbats=F)



