library(igraph)
setwd("/winmounts/walther/AG-Bioinformatics/Dirk/Hanne/Pocketeome/probis_scores")

for(cycle in 1:1) # for repeat runs of randomization experiments
{ 
  
SAVE_OUTPUT=FALSE

#settings for main result set in MS
SHUFFLE=FALSE # randomized similarity matrix => network
SHUFFLE_DISCR=FALSE # shuffle similarity matrix, then discretize (explicit/complicated version of sample_gnm())
DISCRETE=TRUE # unweighted/binary/undirected network
RAND_REWIRE_DEG=FALSE # degree-preserving shuffling of edges, sets DISCRETE TO TRUE and SHUFFLE to FALSE
RAND_NETWORK=FALSE # sample_gnm() generated random network

threshold=0.1 # threshold for setting edges to 1, max_norm ProBIS scores, 0.1 (=paper)/0.75
resol=0.01 # Leiden resolution parameter, 0.01

if(RAND_NETWORK==TRUE) {DISCRETE=TRUE; SHUFFLE=FALSE;}
if(RAND_REWIRE_DEG==TRUE) {DISCRETE==TRUE; SHUFFLE=FALSE;} 
if(SHUFFLE_DISCR==TRUE) {SHUFFLE=TRUE; DISCRETE=FALSE;}


files=c("ARATH_scores.csv","CAEEL_scores.csv","CANAL_scores.csv","DROME_scores.csv","ECOLI_scores.csv","HUMAN_scores.csv","MAIZE_scores.csv","MOUSE_scores.csv","ORYSJ_scores.csv","SOYBN_scores.csv","YEAST_scores.csv")
FS=read.table("foldseek_cluster.csv",header=TRUE,sep=",",row.names=1)
FS <- FS[order(rownames(FS)), ]

# creates a shuffled version of a square, symmetric input matrix, output is also symmetric
shuffle_symmetric_matrix <- function(mat, seed = NULL) {
    
  #if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(mat)
  
  upper_idx <- which(upper.tri(mat))
  upper_values <- mat[upper_idx]
  
  shuffled_values <- sample(upper_values)
  
  # new symmetric matrix
  shuffled_mat <- matrix(0, n, n)
  shuffled_mat[upper_idx] <- shuffled_values
  shuffled_mat <- shuffled_mat + t(shuffled_mat)
  
  # Restore original diagonal
  diag(shuffled_mat) <- diag(mat)
  
  return(shuffled_mat)
}


Nclusters=c()
Ncomm=c()
Npockets=c()
Nsingletons=c()


for(f in files)
{
  d=read.table(f,sep=",",header=T,row.names=1)
  
  sim_matrix=as.matrix(d/max(d))
  
  #h=hist(sim_matrix)
  #plot(h$breaks[2:length(h$breaks)],log(h$counts))
  #readline(prompt="Press [enter] to proceed")
  
  if(DISCRETE==TRUE) {sim_matrix[sim_matrix<=threshold]=0; sim_matrix[sim_matrix>threshold]=1;}
  if(SHUFFLE==TRUE) {
    rN=rownames(sim_matrix)
    cN=colnames(sim_matrix)
    sim_matrix=shuffle_symmetric_matrix(sim_matrix)
    row.names(sim_matrix)=rN
    colnames(sim_matrix)=cN
    if(SHUFFLE_DISCR==TRUE) {sim_matrix[sim_matrix<=threshold]=0; sim_matrix[sim_matrix>threshold]=1;}
    }

  
  #heatmap(sim_matrix, Rowv = NA, Colv = NA,symm = TRUE)
  #check for symmetry
  #upperV=sim_matrix[upper.tri(sim_matrix)]
  #lowerV=t(sim_matrix)[upper.tri(sim_matrix)]
  #cor(upperV,lowerV)
  
  print(paste(f,"N pockets=",ncol(d)))
  
  if(DISCRETE==TRUE | SHUFFLE_DISCR==TRUE) {
      graph = graph_from_adjacency_matrix(sim_matrix, mode = "undirected", weighted = FALSE, diag = FALSE)
    } else {
    graph = graph_from_adjacency_matrix(sim_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  }
  if(RAND_REWIRE_DEG==TRUE) {graph <- rewire(graph, with = keeping_degseq(niter = 10 * gsize(graph)))}
  
  if(RAND_NETWORK==TRUE)
  {
    nNodes=gorder(graph)
    nEdges=gsize(graph)
    V(graph)$name=colnames(sim_matrix)
    graph=sample_gnm(nNodes,nEdges,directed=FALSE,loops=FALSE)
    V(graph)$name=colnames(sim_matrix)
  }
  
  # Leiden clustering
  leiden_clusters = cluster_leiden(graph, resolution_parameter = resol) # resol. controls the size of detected communities: lower values → bigger communities, higher values → smaller ones.
    
  # cluster membership vector
  membership_vec <- membership(leiden_clusters)
  
  #save membership_vector
  fN=paste(f,"_membership_vec.txt",sep="")
  fN=paste(f,"_discretized_",DISCRETE,"_rewired_",RAND_REWIRE_DEG,"_randomShuffle_",SHUFFLE,"_t",threshold,"_r",resol,"_membership_vec.txt",sep="")
  df=data.frame(pocket=names(membership_vec),clusterID=as.vector(membership_vec))
  if(SAVE_OUTPUT==TRUE) {write.table(df, file = fN, row.names = FALSE, col.names = FALSE, quote = FALSE)}
  
  # Tabulate cluster sizes
  cluster_sizes <- table(membership_vec)
  hist(cluster_sizes)
  
  # Count singleton clusters
  num_singletons <- sum(cluster_sizes == 1)
  num_commun <- sum(cluster_sizes >1)
 
  # number of clusters
  Ncl = length(cluster_sizes)
  print(paste("N communities detected=",num_commun,"N singletons=",num_singletons, "N clusters total",Ncl))

  Nclusters = append(Nclusters, Ncl)
  Ncomm = append(Ncomm, num_commun)
  Npockets=append(Npockets,ncol(d))
  Nsingletons=append(Nsingletons,num_singletons)
}

# plot and determine scaling law coefficients
NFS=FS$n_fs_clust
plot(NFS,Nclusters)
m=lm(log(Nclusters)~log(NFS))
plot(log(NFS),log(Nclusters),main=paste("alpha=",m$coefficients[2]))
print(cor.test(log(NFS),log(Nclusters)))
print(paste("alpha=",m$coefficients[2]))
print(paste("intercept=",m$coefficients[1]))

df=data.frame(species=files,n_fs_clust=NFS,n_clusters=Nclusters,n_singletons=Nsingletons,n_communities=Ncomm)
fN=paste("clusterStats","_discretized_",DISCRETE,"_rewired_",RAND_REWIRE_DEG,"_randomShuffle_",SHUFFLE,"_t",threshold,"_r",resol,".txt",sep="")
if(SAVE_OUTPUT==TRUE) {write.table(df, file = fN, row.names = FALSE, col.names = TRUE, quote = FALSE)}
}