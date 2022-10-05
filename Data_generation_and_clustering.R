library("MixSim")
library("clustMixType")
library("pdfCluster")
library("kmed")
library("cluster")
library("kamila")
library("FactoMineR")


# data generation via discretising continuous vars
method1 <- function(objs=250, vars=4, clusters=3, cat=0.5, 
                    sph=TRUE, catlevels=2, overlap=0.01) {
  Q <- MixSim(BarOmega=overlap, K=clusters, p=vars, sph=sph, resN=100000)
  contData <- simdataset(objs, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
  df <- data.frame(contData[[1]])
  for (var in 1:(vars*cat)){
    v <- as.numeric(unlist(df[var]))
    df[var] <- cut(v, quantile(v, seq(0, 1, 1/catlevels)), 
                   labels=seq(1, catlevels, 1), include.lowest=TRUE)
  }
  print(df)
  return(list(df, contData[[2]]))
}

# data generation via cluster labels with noise
method2 <- function(objs=250, vars=4, clusters=3, cat=0.5, 
                    sph=TRUE, overlap=0.01,resample=0.2){
  Q <- MixSim(BarOmega=overlap, K=clusters, p=vars*(1-cat), sph=sph, resN=10000)
  contData <- simdataset(objs, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
  df <- data.frame(contData[[1]])
  classes <- contData[[2]]
  for(catVar in 1:(vars*cat)){
    j <- classes
    randomObjects <- sample.int(objs, objs*resample)
    for (i in randomObjects){
      j[i] <- sample.int(clusters, 1)
    }
    df <- cbind(as.factor(j), df)
  }
  return(list(df, classes))
}

algos <- c("Kproto", "PAM_Gower", "PAM_Ahmad Dey", 
                "Hierarchical_Gower", "Hierarchical_Ahmad Dey", 
                "Modha-Spangler", "Tandem", "Kamila")

computeARIScores <- function(repeats=50, inits=30, iters=1000,
                             objs=250, vars=4, clusters=3, cat=0.5, 
                             sph=TRUE, overlap=0.01, catlevels=2){ 
  ARIScores <- data.frame(matrix(0, nrow=repeats, ncol=length(algos)))
  names(ARIScores) <- algos
  
  # folder for saving simulated datasets with these combination of parameters
  paramsName <- paste(as.character(vars),
                      as.character(clusters),
                      as.character(cat*100),
                      as.character(sph),
                      as.character(overlap*100),
                      as.character(catlevels),
                      sep="_")
  fullFolderName <- paste("method1_datasets", paramsName, sep="/")
  dir.create(fullFolderName)
    
  for (i in 1:repeats){
    simulation <- method1(objs=objs, vars=vars, clusters=clusters, cat=cat,
                            sph=sph, overlap=overlap, catlevels=catlevels)
    simFileName <- paste(paramsName, "df", as.character(i), sep="_")
    saveRDS(simulation, file=paste(fullFolderName, simFileName, sep="/"))
    data <- simulation[[1]]
    classes <- simulation[[2]]
  
    # k-prototypes
    kprotoOut <- kproto(data, clusters, verbose=F,
                        iter.max=iters, nstart=inits)
    ARIScores[i, 1] <- adj.rand.index(kprotoOut$cluster, classes)
    
    # PAM_Gower
    gowerDistMatrix <- distmix(data, method="gower",
                               idcat = seq(1, vars*cat, 1),
                               idnum = seq(vars*cat + 1, vars, 1))
    pam_gOut <- pam(gowerDistMatrix, clusters, nstart=inits,
                    cluster.only=T, do.swap=F)
    ARIScores[i, 2] <- adj.rand.index(pam_gOut, classes)
    
    # PAM_Ahmad_Dey
    # excludes combinations with only one categorical variable
    if(!(vars==4 & cat==0.25)){
      adDistMatrix <- distmix(data, method="ahmad",
                              idcat = seq(1, vars*cat, 1),
                              idnum = seq(vars*cat + 1, vars, 1))
      pam_aOut <- pam(adDistMatrix, clusters, nstart=inits,
                      cluster.only=T, do.swap=F)
      ARIScores[i, 3] <- adj.rand.index(pam_aOut, classes)
    }
    else{ARIScores[i, 3] <- NA}
    
    # Hierarchical_Gower
    hc_gOut <- hclust(as.dist(gowerDistMatrix))
    ARIScores[i, 4] <- adj.rand.index(cutree(hc_gOut, k=clusters), classes)
    
    # Hierarchical_Ahmad Dey
    # excludes combinations with only one categorical variable
    if(!(vars==4 & cat==0.25)){
      hc_aOut <- hclust(as.dist(adDistMatrix))
      ARIScores[i, 5] <- adj.rand.index(cutree(hc_aOut, k=clusters),classes)
    }
    else{ARIScores[i, 5] <- NA}
    
    # Modha-Spangler
    # excludes combinations with only one continuous variable
    if(!(vars==4 & cat==0.75)){
    catData <- dummyCodeFactorDf(data[1:(vars*cat)])
    ms <- gmsClust(conData=data[(vars*cat+1):vars],
                   catData=catData, 
                   nclust=clusters, searchDensity=5)
    ARIScores[i, 6] <- adj.rand.index(ms$results$cluster, classes)
    }
    else{ARIScores[i, 6] <- NA}
  
    # FAMD then k-means
    reducedData <- FAMD(data, ncp=clusters-1, graph=F)$ind$coord
    tandemOut <- kmeans(reducedData, clusters, 
                        nstart=inits, iter.max=iters)$cluster
    ARIScores[i, 7] <- adj.rand.index(tandemOut, classes)
    
    # kamila
    kamilaOut <- kamila(conVar=data[(vars*cat+1):vars],
                catFactor=data[1:(vars*cat)],
                numClust=clusters,
                numInit=inits, maxIter=iters)
    ARIScores[i, 8] <- adj.rand.index(kamilaOut$finalMemb, classes)
  }
  saveRDS(ARIScores, file=paste("method1_scores", paramsName, sep="/"))
  return(colMeans(ARIScores))
}


#collecting information
vars <- c(4, 8, 12)
clusters <- c(3, 5, 7)
catlevels <- c(2, 3, 4)
cat <- c(0.25, 0.5, 0.75)
sph <- c(T, F)
overlap <- c(0.01, 0.1, 0.2)
categories <- c("vars", "clusters", "cat", "sph", "overlap", "cat levels")
numParams <- length(categories)
paramsGrid <- expand.grid(vars=vars, clusters=clusters,
                          cat=cat, sph=sph, overlap=overlap,
                          catlevels=catlevels)

meanARImethod1 <- data.frame((matrix(0, nrow=0, ncol=numParams+2)))
meanARIm1 <- data.frame((matrix(0, nrow=0, ncol=numParams+length(algos))))
names(meanARImethod1) <- c(categories, "algorithm", "mean ARI")
names(meanARIm1) <- c(categories, algos)
meanFileNameLong <- paste("meanARImethod1", as.character(start), 
                          as.character(stop), sep="_")
meanFileNameShort <- paste("meanARIm1", as.character(start), 
                           as.character(stop), sep="_")
rows <- 0
alg <- length(algos)

set.seed(80)
for (params in 1:length(paramsGrid[,1])){
  p <- paramsGrid[params, ]
  meanARImethod1[alg*rows+(1:alg), 1:numParams] <- p
  meanARIm1[rows+1, 1:numParams] <- p
  scores <- computeARIScores(objs=10, repeats=3,
                             vars=p$vars, clusters=p$clusters,
                             cat=p$cat, sph=p$sph, overlap=p$overlap,
                             catlevels=p$catlevels)
  meanARIm1[rows+1, numParams+(1:alg)] <- scores
  for (a in 1:alg){
    meanARImethod1[8*rows+a, numParams+1] <- algos[a]
    meanARImethod1[8*rows+a, numParams+2] <- scores[[algos[a]]]
  }
  saveRDS(meanARImethod1, file=meanFileNameLong)
  saveRDS(meanARIm1, file=meanFileNameShort)
  rows <- rows + 1
}
