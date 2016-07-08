directory="/clusterdata/uqpwains/ibscratch/QC-IDATS"
setwd(directory)
load("mValsSw.RObject")
library("cluster")
library("cluster")
revmval  <- t(mValsSw)

cluster1 = hclust(dist(revmval, method = "euclidean"), method = "average")
png("03_cluster.png")
plot(cluster1,"Clustering")
dev.off()