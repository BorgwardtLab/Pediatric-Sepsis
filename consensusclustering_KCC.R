library(ConsensusClusterPlus)
library(cluster)
library(dbscan)
library(flexclust)
library(glue)
dir = '/cluster/work/borgw/SPSS/MultiOmicsAnalysis/ConsensusClustering/'
setwd(dir)

mkdirs <- function(fp) {
  if(!file.exists(fp)) {
    mkdirs(dirname(fp))
    dir.create(fp)
  }
} 
views <- c('clinical')
Ks <- c(2, 3, 4, 5, 6)
input = 'raw'
method = 'km'
for (K in Ks) {
  for (view in views){
    input_file = glue('data/ThreeViewsKCC_{view}_K_{K}.csv')
    output_dir = glue("Clustering/{input}_{method}/ThreeViewsKCC_{view}_K_{K}/")
    
    
    mat <- read.csv(input_file, row.names = 1)
    mat_t = t(mat)
    
    if (input=='raw') {
      d=mat_t
      distance="euclidean"
    } else {
      d = sweep(d,1, apply(d,1,median,na.rm=T))
      distancee="pearson"
    }
    
    results =  ConsensusClusterPlus(d, maxK=6, reps=100, pItem=0.8, pFeature=1,
                                    clusterAlg=method, title=output_dir, 
                                    distance=distance,seed=126211,plot="png")
    
    name = glue('{output_dir}/ClusterAssignment.csv',sep="")
    
    assigment <- data.frame(index=row.names(mat),
                            cluster_2=results[[2]][["consensusClass"]],
                            cluster_3=results[[3]][["consensusClass"]],
                            cluster_4=results[[4]][["consensusClass"]],
                            cluster_5=results[[5]][["consensusClass"]],
                            cluster_6=results[[6]][["consensusClass"]])
    write.csv(assigment, file=name)   
  }
}


