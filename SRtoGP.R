#Connections between the Successor Representation and diffusion kernel
#Charley Wu 2018

rm(list=ls())
packages <- c('dplyr', "ggplot2", 'GGally', 'network', 'sna', "RColorBrewer", "igraph", 'matrixcalc', 'Matrix', 'lattice', 'jsonlite', 'viridis', 'cowplot', 'ggnetwork')
lapply(packages, require, character.only = TRUE)

source("models.R")

############################################################
#Comparing SR and Diffusion kernel eigenvectors for different example graphs
############################################################
#Make graphs
#1. Standard lattice
g1 <-  make_lattice(dimvector =c(8,8)) #8x8 lattice
#2. prune some edges
found <- F
while (found != T){
  lattice <- make_lattice(dimvector =c(8,8)) #candidate
  delete_edges <- length(E(lattice)) * .2 #how mny edges to delete
  lattice <- igraph::delete.edges(lattice, sample(seq(1:length(E(lattice))), delete_edges)) #randomly delete some edges
  #plot(lattice, vertex.size=5, vertex.label=NA) 
  if (igraph::components(lattice)$no==1){#Check that it only has a single connected component
    found <- found + 1
    g2 <- lattice
  }
}
#TODO Make these look nicer
ggnet(g1)
ggnet(g2)

#Compute SR and diffusion kernel
gamma <- .99
alpha = 1
#SR
SR1 <- SuccRep(g1,gamma)
SR2 <- SuccRep(g2,gamma)
#diffusionkernel
DF1 <- diffusionKernel(g1, alpha)
DF2 <- diffusionKernel(g2, alpha)

#Compute eigenvectors
#SR
eigSR1 <- eigen(SR1)
eigSR2 <- eigen(SR2)
#diffusionkernel
eigDF1 <- eigen(DF1)
eigDF2 <- eigen(DF2)

#Plot eigenvectors
grid <- expand.grid(x1 = seq(1,8), x2 = seq(1,8))
#iterate over each column, convert each vector to a matrix, and 
SReigs <- data.frame(eigenVector=c(as.vector(eigSR1$vectors),as.vector(eigSR2$vectors)), eigenValue = c(rep(eigSR1$value, each = 64),rep(eigSR1$value, each = 64)), 
                     eigenValueRank = c(rep(seq(1,64), each = 64),rep(seq(1,64), each = 64)), x1 = c(rep(grid$x1, 64),rep(grid$x1, 64)), x2 = c(rep(grid$x2, 64),rep(grid$x2, 64)),
                    graph = rep(c(1,2),each=64*64), model = 'SR')
DFeigs <- data.frame(eigenVector=c(as.vector(eigDF1$vectors),as.vector(eigDF2$vectors)), eigenValue = c(rep(eigDF1$value, each = 64),rep(eigDF2$value, each = 64)), 
                     eigenValueRank = c(rep(seq(1,64), each = 64),rep(seq(1,64), each = 64)), x1 = c(rep(grid$x1, 64),rep(grid$x1, 64)), x2 = c(rep(grid$x2, 64),rep(grid$x2, 64)),
                     graph = rep(c(1,2),each=64*64), model='DF')
eigDF <- rbind(SReigs, DFeigs)
eigDF$eigenValueRank <- factor(eigDF$eigenValueRank)

pSR1 <- ggplot(subset(eigDF, graph==1 & model=='SR'), aes(x = x1, y = x2, fill = eigenVector))+
  geom_tile()+
  scale_fill_viridis_c(name='')+
  facet_wrap(~eigenValueRank)+
  theme_void()+
  ggtitle('Successor Representation')
pSR1

pSR2 <- ggplot(subset(eigDF, graph==2 & model=='SR'), aes(x = x1, y = x2, fill = eigenVector))+
  geom_tile()+
  scale_fill_viridis_c(name='')+
  facet_wrap(~eigenValueRank)+
  theme_void()+
  ggtitle('Successor Representation')
pSR2

pDF1 <- ggplot(subset(eigDF, graph==1 & model=='DF'), aes(x = x1, y = x2, fill = eigenVector))+
  geom_tile()+
  scale_fill_viridis_c(name='')+
  facet_wrap(~eigenValueRank)+
  theme_void()+
  ggtitle('Diffusion Kernel')
pDF1

pDF2 <- ggplot(subset(eigDF, graph==2 & model=='DF'), aes(x = x1, y = x2, fill = eigenVector))+
  geom_tile()+
  scale_fill_viridis_c(name='')+
  facet_wrap(~eigenValueRank)+
  theme_void()+
  ggtitle('Diffusion Kernel')
pDF2

p <- plot_grid(pSR1,pDF1, pSR2, pDF2, nrow = 2 )
ggsave('plots/SRGPeigenvectors.pdf', p, height = 12 , width = 12, units = 'in')


###########################################################################################################################################
# What does exponentiation do to eigenvectors?
############################################################################################################################################
#Extracted from diffusionKernel()
D <- diag(strength(g1)) #degree matrix
W <- as_adjacency_matrix(g1, sparse=F)
L <- D - W
L <- sqrt(solve(D))%*%L%*%sqrt(solve(D)) #if degree normalized Laplacian

#Eigenvectors
eigL <- eigen(L) #use negative to reverse the ordering of the 

laplacianDF <- data.frame(eigenVector=as.vector(eigL$vectors), eigenValue = rep(eigL$value, each = 64),
           eigenValueRank = rep(seq(1,64), each = 64), x1 = rep(grid$x1, 64), x2 = rep(grid$x2, 64), model = 'NormalizedLaplacian')


pLa <- ggplot(laplacianDF, aes(x = x1, y = x2, fill = eigenVector))+
  geom_tile()+
  scale_fill_viridis_c(name='')+
  facet_wrap(~eigenValueRank)+
  theme_void()+
  ggtitle('Normalized Laplacian')
pLa

#Compute matrix exponentiation
expL <- eigL$vectors %*% diag(exp(alpha * eigL$values)) %*% t(eigL$vectors)
eigExpL <- eigen(expL)

exponentiateddF <- data.frame(eigenVector=as.vector(eigExpL$vectors), eigenValue = rep(eigExpL$value, each = 64),
                          eigenValueRank = rep(seq(1,64), each = 64), x1 = rep(grid$x1, 64), x2 = rep(grid$x2, 64), model = 'Diffusion Kernel (alpha=1)')

pDF <- ggplot(exponentiateddF, aes(x = x1, y = x2, fill = eigenVector))+
  geom_tile()+
  scale_fill_viridis_c(name='')+
  facet_wrap(~eigenValueRank)+
  theme_void()+
  ggtitle('Diffusion Kernel')
pDF

p <- plot_grid(pLa,pDF,  rows=1  )
ggsave('plots/LaplacianEigenvectors.pdf', p, height = 6 , width = 12, units = 'in')
###########################################################################################################################################
#Circa 2018
# Create some graph stimuli
############################################################################################################################################
#Aesthetics
palf <- colorRampPalette(c("white", "black")) 
###########
#Generate some lattice network with broken edges
numGraphs <- 40
dimensions <- c(3,3)
pruningFactor <- 0.2 #proportion of edges to delete

graphs <- list()
found <- 0
#rejection sampling to make sure that all generated networks are connected (i.e., without orphan nodes)
while (found < numGraphs){
  lattice <- make_lattice(dimvector =dimensions) #candidate
  delete_edges <- length(E(lattice)) * pruningFactor #how mny edges to delete
  lattice <- igraph::delete.edges(lattice, sample(seq(1:length(E(lattice))), delete_edges)) #randomly delete some edges
  #plot(lattice, vertex.size=5, vertex.label=NA) 
  if (igraph::components(lattice)$no==1){#Check that it only has a single connected component
    found <- found + 1
    graphs[[found]] <- lattice
  }
}
#plot example
g <- graphs[[1]]
plot(g, vertex.size=10, vertex.label=NA, margin=0)


#media network as another example
nodes <- read.csv("data/Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("data/Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net <- simplify(net, remove.multiple = F, remove.loops = T) #remove loops
net <- as.undirected(net, ) #convert to undirected graph

#plot(net, vertex.size=10, vertex.label=NA, margin=0)
set.edge.attribute(net, "color", ifelse(bip %e% "weights" > 1, "black", "grey75"))
ggnet(net, edge.size = "weight")
E(net)$weight

net.df <- ggnetwork::ggnetwork(net,layout="fruchtermanreingold",weights="weight",niter=50000)
ggplot2::ggplot(net,df,aes(x=x,y=y,xend=xend,yend=yend))+
  ggnetwork::geom_edges(color="gray",size=aes(size=weight))+
  ggnetwork::geom_nodes(color="black")+
  ggnetwork::geom_nodelabel_repel(aes(label=vertex.names,color=vertex.names),fontface="bold",box.padding=unit(1,"lines"))+
  ggplot2::theme_minimal()+ggplot2::theme(axis.text=ggplot2::element_blank(),axis.title=ggplot2::element_blank(),legend.position="none")
######################################################################################################################################################
# Graph primatives
######################################################################################################################################################

I <- diag(rep(1,length(V(g)))) #identity matrix
#Adjacency matrix (weighted if graph is also weighted)
if(is_weighted(g)){
  A <- as_adjacency_matrix(g,  attr="weight", sparse=F)  
}else{
  A <- as_adjacency_matrix(g, sparse=F) 
}
D <- diag(strength(g)) #degree matrix
n <-length(V(g)) #number of nodes

######################################################################################################################################################
# Construct different types of laplacians
######################################################################################################################################################

#Graph Laplacian: it is defined as L = D-A in Kondor & Lafferty (2002) and Smola & Kondor (2003), but as L = A-D in Kondor & Vert (2003)
L_df <- D-A
levelplot(L_df)

#Normalized Laplacian defined the conventional way as I - D^(-1/2)AD^(-1/2)
L_normed <- I - sqrt(solve(D))%*%A%*%sqrt(solve(D))
levelplot(L_normed)

#Normalized laplacian defined in terms of the diffusion Laplacian as \tilde{L} = D^-1/2 L_df D^-1/2 (from Smola & Condor, 2003) and which is also used in Machado et al (2018)
L_normed2  <- sqrt(solve(D))%*%L_df%*%sqrt(solve(D))
levelplot(L_normed2)
all.equal(L_normed, L_normed2) #proof holds

#Are the normalized eigenvectors of L_df the same as the eigenvectors of the normalized Laplacian?
ev1 <- eigen(L_df)
eigenValues1 <- ev1$values
eigenVectors1 <- ev1$vectors
levelplot(eigenVectors1)

ev2 <- eigen(L_normed)
eigenValues2 <- ev2$values
eigenVectors2 <- ev2$vectors
levelplot(eigenVectors2)

all.equal(eigenValues1, eigenValues2) #eigen vectors and values are not equal (but maybe equivalent to some factor)
all.equal(eigenVectors1, (sqrt(solve(D))%*%eigenVectors2%*%sqrt(solve(D)))) #


#Row norm Laplacian aka random walk laplacian (for unweighted undirected graphs)
L_rownorm<- solve(D)%*%L_df
levelplot(L_rownorm)
all.equal(L_rownorm, I-solve(D)%*%A) #Stachen feld thesis: the row normalized laplacian is equal to I - T_rw (the random walk transition matrix of the SR)


######################################################################################################################################################
# Building the Sucessor Representation
######################################################################################################################################################
library(reshape2)
#Random walk transition probabilities from Stachenfeld et al. (2014)
T_rw <-  solve(D)%*%A
gamma <- 0.9
M <- solve(I - gamma* T_rw) #SR representation in closed form

ggplot(melt(M), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  theme_minimal()

ev3 <- eigen(M)#compute the eigen vectors
levelplot(ev3$vectors)

#SOLVED! EQUIVALENCY BETWEEN M(S,S') AND L
all.equal(T_rw, I-(solve(D)%*%L_df)) #equivalency between the T of the SR and the graph laplacian L used in the diffusion kernel
all.equal(M, solve(I - gamma*(I-(solve(D)%*%L_df)))) #!EQUIVALENCY BETWEEN M AND L

#Some other redefinitions
all.equal(T_rw, I - (L_rownorm)) #T in terms of L_RW
all.equal(L_df, D%*%(I- T_rw)) #L in terms of T

left <-D %*% (I -(I - solve(M))/gamma )
right <- L_df
all.equal(left, right) #L in terms of M


######################################################################################################################################################
# Compute eigenvectors use v for eigenvectors and lambda for eigenvalues
######################################################################################################################################################
#Graph laplacian (used in the diffusion kernel)
ev_df <- eigen(L_df)
lambda_df <- ev_df$values
v_df <- ev_df$vectors
levelplot(v_df)
#normalized laplacian
ev_norm <- eigen(L_normed)
lambda_norm <- ev_norm$values
v_norm <- ev_norm$vectors
levelplot(v_norm)
#Row normalized laplacian
ev_rw <- eigen(L_rownorm)
lambda_rw <- ev_rw$values
v_rw <- ev_rw$vectors
levelplot(v_rw)
#SR
ev_sr <- eigen(M)
lambda_sr <- ev_sr$values
v_sr <- ev_sr$vectors
levelplot(v_sr)

######################################################################################################################################################
# Try to verify equivalencies stated in various papers
######################################################################################################################################################


#Stachenfeld (2017) says that the eigenvectors of M are the same as the eigenVectors of T_rw (or conversely in her thesis, with the eigenvectors of L_rw)
levelplot(v_norm) 
levelplot(v_rw)
levelplot(v_sr[,n:1]) #reversed order
all.equal(v_norm, v_sr[,n:1]) #They are similar but not the same. NOTE THAT THE EQUIVALENCE MAY ONLY HOLD IN THE LIMIT OF CONTINUITY

#look at first and last eigenvectors
levelplot(matrix(v_norm[,1], nrow=sqrt(n))) #both grided
levelplot(matrix(v_sr[,n], nrow=sqrt(n)))
#levelplot(matrix(v_rw[,1], nrow=sqrt(n)))

levelplot(ifelse(v_sr<0,0,v_sr)) #Thresholded at 0
levelplot(matrix(M[,n], nrow=sqrt(n))) #The receptive field of the nth node
levelplot(matrix(v_sr[,4], nrow=sqrt(n))) #Simulated grid cell using the ith eigenvector of the SR


#Machado et al. (ICLR 2018) theorem of the equivalency between the eigenvalues and eigenvectors of the normalized laplacian and SR
#This is based on the Stachenfeld et al. (2014) claim that the eigenvectors of L_rw are equal to L_normed by the scaling factor D^1/2, but broadens it by including the discount factor gamma
#EIGEN VALUES ARE RELATED, BUT THE EQUIVALENCY FOR EIGENVECTORS MAY ONLY HOLD IN THE LIMIT OF CONTINUITY
#Check first that the eigen values are related
all.equal(lambda_norm, 1 - (1 - rev(lambda_sr)^-1)*gamma^-1) #Eigen values are related

i = 1 #choose some arbitrary i
j = n-i + 1 #j indexes from the opposite direction where i + j = n + 1
all.equal(v_norm[,j], c(gamma^-1 * sqrt(D)%*%v_sr[,i])) #doesn't hold
plot(v_norm[,j])
plot(c((gamma^-1 * sqrt(D))%*%v_sr[,i]))
levelplot(v_norm)
levelplot((gamma^-1 *sqrt(D))%*%v_sr[,n:1])


plot(v_norm[,j], main='Normalized Laplacian')
plot((gamma^-1 *sqrt(D))%*%v_sr[,i], main='SR-rescaled')

#Where are the large differences?
scaleFactor <- v_norm / (gamma^-1 *sqrt(D))%*%v_sr[,n:1]
levelplot(scaleFactor) #this may be showing the connecting to the 


sqrt(D)%*%eigenVectors2 - eigenVectors3
#Mahadevan (2007) claims that the eigenvectors of L_rw are equivalent to the eigenvectors of the normalized laplacian by the relation:
# lambda_rw = (1 - lambda_norm)*D^-1/2
eigen(T_rw)$values - (1-eigen(L_normed)$values)%*%sqrt(solve(D))
eigen(T_rw)$vectors - (1-eigen(L_normed)$vectors)%*%sqrt(solve(D)) #Doesn't work


############################################################################################################################################################################################################################
#Generalization Curves
################################################################################################################################################################################


####################################################################################
# Load pre-generated graphs from exp1
####################################################################################
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

#save plot
png(file = 'plots/exampleGraph.png', height = 3, width = 3, units = 'in', res = 600)
plot(graphList[[4]], vertex.size=15, vertex.label=NA, margin=0, layout=layout_with_fr)
dev.off()


####################################################################################
# Compute correlations
####################################################################################

#Simulate GP correlations
betaList <- c(0.5,1,2 ,4)
gpSimDF <- data.frame()
replications <- 10000

for (i in 1:length(graphList)){ #loop through graphs
  g <- graphList[[i]] #select graph
  #Compute GP corelations
  for (betaValue in betaList){
    k <- diffusionKernel(g,beta = betaValue) #compute kernel matrix
    prior <- mvrnorm(replications, rep(0, length(V(g))), k) #sample from the prior
    #Group pairs of nodes by distance
    distanceTable <- distances(g) #distances between each node
    distanceTable[lower.tri(distanceTable)] <- NA #set lower symmetric triangle to NA
    uniqueDistances <- unique(c(distanceTable[!is.na(distanceTable)]))
    for (d in uniqueDistances){
      pairs <- which(distanceTable==d, arr.ind=T)
      valuesA <- prior[,pairs[,1]]#rewards
      valuesB <- prior[,pairs[,2]]
      dat <- data.frame(graphNum = i, Beta = betaValue, distance = d, correlation= cor(c(valuesA),c(valuesB)), model = 'GP')
      gpSimDF<- rbind(gpSimDF,dat)
    }
  }
}

plotgpSimDF <- gpSimDF %>% group_by(distance, Beta) %>% summarize(Correlation = mean(correlation))
plotgpSimDF$Beta <- factor(plotgpSimDF$Beta)

myPalette <- colorRampPalette(rev(brewer.pal(4, "Spectral")), space="Lab")
p1 <- ggplot(plotgpSimDF, aes(x=distance, y = Correlation, color = Beta))+
  geom_line(size = 1) +
  #geom_point(size = 2) +
  ylab("Pearson Correlation")+xlab("Distance")+
  theme_classic()+
  #xlim(c(0,10))+
  scale_x_continuous(breaks = round(seq(0, max(plotgpSimDF$distance), by = 1),1), limits = c(0,6))+
  scale_color_manual(values = myPalette(4), name = expression(beta))+
  #scale_color_viridis_d(name = expression(beta))+
  #scale_color_brewer(palette="Set1", direction=-1)+
  scale_linetype_manual(values = c( 'twodash', 'longdash','solid'))+
  #scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(text = element_text(size=14,  family="sans"), legend.justification = c(1, 1), legend.position = c(1, 1))+
  theme(strip.background=element_blank(), legend.key=element_rect(color=NA), legend.title=element_text(size=14))
#ggtitle("Decay of Correlation")
p1

ggsave(p1, filename='plots/GPcorrelationDecay.pdf', width = 3, height = 3, units = 'in')


#Compute SR correlations
gammaList <- c(0.5, 0.6, 0.7, 0.8)
srSimDF <- data.frame()
for (gammaValue in gammaList){
  dat <- data.frame(Gamma = gammaValue, distance = seq(0,6), correlation= gammaValue^seq(0,6), model = 'SR')
  srSimDF<- rbind(srSimDF,dat)
}

plotsrSimDF <- srSimDF %>% group_by(distance, Gamma) %>% summarize(Correlation = mean(correlation))
plotsrSimDF$Gamma <- factor(plotsrSimDF$Gamma)


p2 <- ggplot(plotsrSimDF, aes(x=distance, y = Correlation, color = Gamma))+
  geom_line(size = 1) +
  #geom_point(size = 2) +
  ylab("Reward Generalization")+xlab("Distance")+
  theme_classic()+
  #xlim(c(0,10))+
  scale_x_continuous(breaks = round(seq(0, max(plotsrSimDF$distance), by = 1),1), limits = c(0,6))+
  scale_color_viridis_d(name = expression(gamma), direction = -1)+
  #scale_color_brewer(palette="Set1", direction=-1)+
  scale_linetype_manual(values = c( 'twodash', 'longdash','solid'))+
  #scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(text = element_text(size=14,  family="sans"), legend.justification = c(1, 1), legend.position = c(1, 1))+
  theme(strip.background=element_blank(), legend.key=element_rect(color=NA), legend.title=element_text(size=14))
#ggtitle("Decay of Correlation")
p2

ggsave(p2, filename='plots/SRcorrelationDecay.pdf', width = 3, height = 3, units = 'in')



################################################################################################################
#Comparison to RBF prior on n x n lattice graph
################################################################################################################

reps <- 1000 #number of environments to generate for each length scale
kernelFun <- rbf #choose kernel function

#pre-compute all pairwise distances between points in a 8 x 8 matrix
grid <- as.matrix(expand.grid(1:8,1:8))
pairwiseDistances <- apply(grid, 1, function(i) apply(grid, 1, function(j) dist(rbind(i, j),method='manhattan')))
uniqueDistances <- unique(c(pairwiseDistances))

#create data frame
corData <- data.frame(lambda = numeric(), correlation = numeric(), distance = numeric())

#Loop over different length scales
for (lambda in c(0.5, 1, 2)){
  #generate samples from GP-rbf priors
  sampledEnvs <- mvrnorm(n = reps, mu = rep(0,64), kernelFun(as.matrix(expand.grid(1:8,1:8)),as.matrix(expand.grid(1:8,1:8)),c(lambda,lambda,1,.0001)))
  #calculate the strength of correlations for each distance
  correlations <- c() #placeholder
  for (d in uniqueDistances){
    pairs <- which(pairwiseDistances==d, arr.ind=T) # select all pairs of points where the pairwise distance is equal to d
    valuesA <- sampledEnvs[,c(pairs[,1])] #rewards
    valuesB <- sampledEnvs[,c(pairs[,2])]
    correlations <- rbind(correlations, cor(c(valuesA),c(valuesB)))
  }
  newDat <- data.frame(lambda = rep(lambda,length(uniqueDistances)), correlation= correlations, distance = uniqueDistances)
  corData <- rbind(corData, newDat)
}

corData$Lambda <- factor(corData$lambda)
# Plot


p <- ggplot(corData, aes(x=distance, y = correlation, color = Lambda, shape = Lambda, linetype = Lambda)) + 
  geom_line(size = 1) +
  #geom_point(size = 2) +
  ylab("Pearson Correlation")+xlab("Distance")+
  theme_classic()+
  #xlim(c(0,10))+
  scale_x_continuous(breaks = round(seq(1, 10, by = 1),1), limits = c(0,10))+
  scale_color_brewer(palette="Set1", direction=-1)+
  scale_linetype_manual(values = c( 'twodash', 'longdash','solid'))+
  #scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(text = element_text(size=14,  family="sans"), legend.justification = c(1, 1), legend.position = c(1, 1))+
  theme(strip.background=element_blank(), legend.key=element_rect(color=NA)) +
  ggtitle("RBF")
p

ggsave(filename = 'plots/kernelCorrelationRBF.pdf', p, height = 4, width = 5)



########################################################################################################################################
#GP Inference
#########################################################################################################################################
library(GGally)
library(ggnet)
library(network)
library(sna)
library(intergraph)
library(igraph)
library(RColorBrewer)
library(cowplot)


colfunc <- colorRampPalette(brewer.pal(3, "YlOrRd"))


####################################################################################
# Load pre-generated graphs from exp1
####################################################################################
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

g <- graphList[[sample(1:40,1)]]

#Define fixed node positions
nodePositions  = gplot.layout.fruchtermanreingold(asNetwork(g), NULL)
#nodePositions <- read.csv('plots/nodepositionsFig1.csv')
V(g)$x = nodePositions[, 1]
V(g)$y = nodePositions[, 2]

#save plot
exampleGraph <- ggnet2(g, mode = c("x", "y"), color = '#E69F00', edge.color = 'black')
exampleGraph
ggsave(filename = 'plots/exampleGraph.pdf', exampleGraph, height = 3, width = 3, units = 'in', useDingbats = F)


#observations
#X <- sample(seq(1,length(V(g))), size = 5, replace = F)
X <- c(7,3,5,4,8)
Y <- round(V(g)$reward[X])


#Let's plot what this looks like
V(g)$label <- ""
V(g)$label[X] <- round(V(g)$reward[X])
V(g)$observed <- 'grey'
V(g)$observed[X] <- colfunc(50)[Y]

obsGraph <- ggnet2(g,  mode = c("x", "y"), color = "observed", label = 'label', label.size = 5, legend.position = "None", edge.color = 'black')
obsGraph
#ggsave(filename = 'plots/obsGraph.pdf', obsGraph, height = 3, width = 3, units = 'in', useDingbats = F)

#GP Inference
K <- diffusionKernel(g, beta = 2) #diffusion kernel
X.test <- seq(1:length(V(g))) #node list
posterior <- gpr(X.test, X = X, Y=Y, k = K, mu_0 =25) #zero mean for rewards

#Mean reward
V(g)$gpMu <- V(g)$label
V(g)$gpMu[-X] <- round(posterior$mu)[-X]
#Uncertainty
V(g)$gpSig <-  normalizeMinMax(posterior$sig, 11, 20)
V(g)$gpSig[X] <- 0


g.df <- as.data.frame(list(Vertex=X.test, gpMu=V(g)$gpMu, gpSig = V(g)$gpSig, color = V(g)$observed, x = V(g)$x, y = V(g)$y ))

gpPredGraph <- ggnet2(g,  mode = c("x", "y"), color = "observed", label = 'gpMu', label.size = 5, legend.position = "None", edge.color = 'black')+
  geom_point(size = V(g)$gpSig, alpha = 0.15) 
gpPredGraph
#ggsave(filename = 'plots/gpPredGraph.pdf', gpPredGraph, height = 3, width = 3, units = 'in', useDingbats = F)


#Successor Representation
M <- SuccRep(g, gamma = 0.7)
Rewards <- rep(0, length(V(g))) #default unobserved nodes to mean?
Rewards[X]<- Y - 25 #put it observations
SRpred <- (M %*% Rewards) + 25

V(g)$SR <- V(g)$label
V(g)$SR[-X] <- round(SRpred[,1])[-X]

SRpredGraph <- ggnet2(g,  mode = c("x", "y"), color = "observed", label = 'SR', label.size = 5, legend.position = "None", edge.color = 'black')
SRpredGraph
#ggsave(filename = 'plots/SRpredGraph.pdf', SRpredGraph, height = 3, width = 3, units = 'in', useDingbats = F)

together <- plot_grid(obsGraph,gpPredGraph,  SRpredGraph,ncol = 1 )
together

ggsave(filename = 'plots/examplePreds.pdf', together, height = 9, width = 3, units = 'in', useDingbats = F)


