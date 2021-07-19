#Charley wu
rm(list=ls())
packages <- c('dplyr', "ggplot2", 'GGally', 'network', 'sna', "RColorBrewer", "igraph", 'matrixcalc', 'Matrix', 'lattice', 'jsonlite', 'viridis', 'cowplot', 'ggnetwork', 'reshape2', 'mosaic')
lapply(packages, require, character.only = TRUE)

source("models.R")


######################################################################################################################################################
# Variational Auto Encoder (cartoon example)
######################################################################################################################################################
set.seed(123)

#create input data
xin.df <- data.frame(x = 1:10, xin =  runif(10))


pXin <- ggplot(xin.df, aes(x=1, y=x, fill=xin)) + 
  geom_tile()+
  theme_minimal()+
  scale_fill_viridis()+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
pXin
 
#low dimension
z <- runif(5)


pZ <- ggplot(data.frame(x = 1:5, z = z), aes(x=1, y=x, fill=z)) + 
  geom_tile()+
  theme_minimal()+
  scale_fill_viridis()+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
pZ


xin.df$xout <- xin.df$xin + rnorm(10, mean = 0, sd = 0.1)
pXout <- ggplot(xin.df, aes(x=1, y=x, fill=xout)) + 
  geom_tile()+
  theme_minimal()+
  scale_fill_viridis()+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
pXout


pVAE <- plot_grid(pXin, pZ, pXout, nrow = 1)
pVAE

ggsave('plots/vae.pdf', width = 8 , height = 6, units = 'in')
######################################################################################################################################################
# Example graph based on Wu et al., 2021, to illustrate SR
######################################################################################################################################################

experimentGraphFiles <- 'data/graphs.json'
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



# assign edge's width as a function of weights.
E(g)$width <- E(g)$weight + min(E(g)$weight) -5 # offset=- 1

pGraph <- ggnet2(g, label = TRUE)
pGraph

ggsave('plots/exampleGraph.pdf', pGraph, width = 3, height = 3, units = 'in')


#singular reward observations
r.df <- data.frame(s = 1:length(V(g)), r =  V(g)$reward)
pReward <- ggplot(r.df, aes(x=1, y=s, fill=r)) + 
  geom_tile()+
  theme_minimal()+
  scale_fill_distiller(palette = "OrRd")+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
pReward
ggsave('plots/SRreward.pdf', pReward, width = 3, height = 3, units = 'in')


#Create graph from random transition matrix 
# a <- matrix(runif(5*5, min=0, max=2), ncol=5)
# diag(a) <- 0 # remove loops.
# a[lower.tri(a)] = t(a)[lower.tri(a)]#make symmetric
# a
# 
# pMatrix <- ggplot(melt(a), aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile()+
#   theme_minimal()+
#   scale_fill_gradient(low = 'white', high = 'black')+
#   theme(axis.line=element_blank(),
#         axis.text.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         legend.position="none",
#         panel.background=element_blank(),
#         panel.border=element_blank(),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         plot.background=element_blank())
# pMatrix
# ggsave('plots/transitionMatrix.pdf', pMatrix, width = 3, height = 3, units = 'in')


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
# Building the Sucessor Representation
######################################################################################################################################################
library(reshape2)
#Random walk transition probabilities from Stachenfeld et al. (2014)
T_rw <-  solve(D)%*%A
gamma <- 0.9
M <- solve(I - gamma* T_rw) #SR representation in closed form

pSR <- ggplot(melt(M), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  theme_minimal()+
  scale_fill_distiller(palette = "PuBu")+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
pSR
ggsave('plots/SR.pdf', pSR, width = 3, height = 3, units = 'in')


#Make value generalizations
v.df <- data.frame(s = 1:length(V(g)), v = rowSums(t(M) * V(g)$reward))
pValue <- ggplot(v.df, aes(x=1, y=s, fill=v)) + 
  geom_tile()+
  theme_minimal()+
  scale_fill_distiller(palette = "OrRd")+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
pValue
ggsave('plots/SRvalue.pdf', pValue, width = 3, height = 3, units = 'in')

ev3 <- eigen(M)#compute the eigen vectors

ggplot(melt(ev3$vectors), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  theme_minimal()+
  scale_fill_gradient(low = 'white', high = 'black')+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

#############################################
# Construct different candidate graphs
#############################################

#Generate some lattice network with some broken edges
set.seed(1234)
graphs <- list()
edges <- c(1,2,3,4)
found <- 0
#rejection sampling to make sure that all generated networks are connected (i.e., without orphan nodes)
while (found < 4){
  lattice <- make_lattice(dimvector = c(3,3)) #candidate
  delete_edges <- edges[found+1]
  lattice <- delete.edges(lattice, sample(seq(1:length(E(lattice))), delete_edges)) #randomly delete some edges
  #plot(lattice, vertex.size=5, vertex.label=NA) 
  if (components(lattice)$no==1){#Check that it only has a single connected component
    found <- found + 1
    graphs[[found]] <- lattice
  }
}

p1 <- ggnet2(graphs[[1]],  node.size = 6, node.color = 'black',  edge.size = 1,vertex.label = NA)
p2 <- ggnet2(graphs[[2]], node.size = 6, node.color = 'black',edge.size = 1,vertex.label = NA)
p3 <- ggnet2(graphs[[3]],  node.size = 6, node.color = 'black', edge.size = 1,vertex.label = NA)
p4 <- ggnet2(graphs[[4]], node.size =6,  node.color = 'black',edge.size = 1, vertex.label = NA)

E(graphs[[1]])
E(graphs[[2]])
E(graphs[[3]])
E(graphs[[4]])

exampleGraphs <- plot_grid(p4,p3,p2,p1, nrow = 1)
exampleGraphs
ggsave('plots/exampleGraphs.pdf', exampleGraphs, width = 12, height = 3, units = 'in')

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