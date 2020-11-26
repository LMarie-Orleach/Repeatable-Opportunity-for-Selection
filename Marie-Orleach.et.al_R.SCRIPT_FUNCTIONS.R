## Source R script associated to the manuscript entitled 
## "The repeatable opportunity for selection differs between pre- and postcopulatory fitness components" 
## by L. Marie-Orleach, N. Vellnow, and L. Schärer


## Function to relativise a data column
REL <- function (x) { x/mean(x, na.rm=TRUE) }


## Function to compute the fitness components (as explained on Figure 1).
FITCOMP <- function(dataset) {
  dataset$mRS <- dataset$focal_offspring # mRS
  dataset$F   <- dataset$total_offspring # F
  dataset$MS  <- dataset$focal_matings / dataset$total_matings # MS
  dataset$STE <- (dataset$focal_sperm / dataset$total_sperm) / (dataset$focal_matings / dataset$total_matings) # STE
  dataset$SFE <- (dataset$focal_offspring / dataset$total_offspring) / (dataset$focal_sperm / dataset$total_sperm) # SFE
  dataset$SFE[which(!is.finite(dataset$SFE))] <- NaN
  
  return(dataset)
}


## Function to compute the standard variance in a fitness component
STDVAR <- function(fit.comp, dataset) {
  var(s(dataset[,fit.comp])/mean(s(dataset[,fit.comp]), na.rm=TRUE)) 
}


## Function to estimate the binomial sampling error arising from STE and SFE
## see Pelissie et al, Evolution, 2012 (Appendix A) and Marie-Orleach et al, Evolution, 2016 (Supp. Info. E) for more details
BSE <- function(fit.comp, dataset) {
  if(fit.comp=="STE"){
    dataset$V_STS_exp <- ((dataset$focal_sperm/dataset$total_sperm) * (1 - (dataset$focal_sperm/dataset$total_sperm))) / (dataset$MS^2) # expected variance in each STS estimate
    bse_STE_rel <- (sum(dataset$V_STS_exp) / (sum(dataset$total_sperm))) / (mean(dataset$STE)^2) # binomial sampling error in STE
    return (bse_STE_rel)
  }
  
  if(fit.comp=="SFE"){
    dataset$STS       <- dataset$focal_sperm / dataset$total_sperm # sperm transfer success
    dataset$V_PS_exp  <- (dataset$focal_offspring / dataset$total_offspring) * (1 - (dataset$focal_offspring /dataset$total_offspring)) # expected variance in each PS estimate
    dataset$V_SFE_exp <- dataset$V_PS_exp / (dataset$STS^2) # expected variance in each SFE estimate
    dataset_147 <- subset(dataset, focal_sperm!=0) # subset
    bse_SFE_rel <- (sum(dataset_147$V_SFE_exp) / sum(dataset_147$total_offspring)) / (mean(dataset_147$SFE)^2) # binomial sampling error in SFE
    return (bse_SFE_rel)
  }
}


## Function to compute the standard covariance between two fitness components
STDCOV <- function (fit.comp1, fit.comp2, dataset){
  if (fit.comp1=="SFE" | fit.comp2=="SFE") { dataset <- dataset[complete.cases(dataset),] }
  CovArray <- (REL(dataset[,fit.comp1])-mean(REL(dataset[,fit.comp1]), na.rm=TRUE)) * (REL(dataset[,fit.comp2])-mean(REL(dataset[,fit.comp2]), na.rm=TRUE))
  return (mean(CovArray))
}


## Function to compute the total variance explained by our model
TOTVAR <- function(dataset){
  stdvar_F   <- STDVAR ("F"  , dataset)
  stdvar_MS  <- STDVAR ("MS" , dataset)
  stdvar_STE <- STDVAR ("STE", dataset)
  stdvar_SFE <- STDVAR ("SFE", dataset)
  
  bse_STE <- BSE ("STE", dataset)
  bse_SFE <- BSE ("SFE", dataset)
  
  stdcov_F.MS    <- STDCOV ("F",   "MS",  dataset)
  stdcov_F.STE   <- STDCOV ("F",   "STE", dataset)
  stdcov_F.SFE   <- STDCOV ("F",   "SFE", dataset)
  stdcov_MS.STE  <- STDCOV ("MS",  "STE", dataset)
  stdcov_MS.SFE  <- STDCOV ("MS",  "SFE", dataset)
  stdcov_STE.SFE <- STDCOV ("STE", "SFE", dataset)
  
  return(stdvar_F + stdvar_MS + (stdvar_STE-bse_STE) + (stdvar_SFE-bse_SFE) +
    (2*stdcov_F.MS) + (2*stdcov_F.STE) + (2*stdcov_F.SFE) + (2*stdcov_MS.STE) + (2*stdcov_MS.SFE) + (2*stdcov_STE.SFE) +
    bse_STE + bse_SFE)
}


## Function to bootstrapp the variances and covariances to assess their 95 CI
## Outcomes are here expressed as percentages of total variance
BOOT_VARCOVAR <- function (dataset, iteration){
  N <- NROW(dataset)
  
  stor.data <- data.frame(stdvar_mRS=numeric(iteration), stdvar_F=numeric(iteration), stdvar_MS=numeric(iteration), stdvar_STE=numeric(iteration), stdvar_SFE=numeric(iteration),
                        bse_STE=numeric(iteration), bse_SFE=numeric(iteration),
                        stdcov_F.MS=numeric(iteration), stdcov_F.STE=numeric(iteration), stdcov_F.SFE=numeric(iteration), stdcov_MS.STE=numeric(iteration), stdcov_MS.SFE=numeric(iteration), stdcov_STE.SFE=numeric(iteration))
  
  for(i in 1:iteration){
    rand = sample(1:N, N, replace=TRUE)
    dataset$nb <- 1:N
    newdataset = dataset[match(rand, dataset$nb), ]
    
    stdvar_mRS <- STDVAR ("mRS", newdataset)
    stdvar_F   <- STDVAR ("F"  , newdataset)
    stdvar_MS  <- STDVAR ("MS" , newdataset)
    stdvar_STE <- STDVAR ("STE", newdataset)
    stdvar_SFE <- STDVAR ("SFE", newdataset)
    
    bse_STE <- BSE ("STE", newdataset)
    bse_SFE <- BSE ("SFE", newdataset)
    
    stdcov_F.MS    <- STDCOV ("F",   "MS",  newdataset)
    stdcov_F.STE   <- STDCOV ("F",   "STE", newdataset)
    stdcov_F.SFE   <- STDCOV ("F",   "SFE", newdataset)
    stdcov_MS.STE  <- STDCOV ("MS",  "STE", newdataset)
    stdcov_MS.SFE  <- STDCOV ("MS",  "SFE", newdataset)
    stdcov_STE.SFE <- STDCOV ("STE", "SFE", newdataset)
    
    totalvariance <- stdvar_F + stdvar_MS + (stdvar_STE-bse_STE) + (stdvar_SFE-bse_SFE) +
      (2*stdcov_F.MS) + (2*stdcov_F.STE) + (2*stdcov_F.SFE) + (2*stdcov_MS.STE) + (2*stdcov_MS.SFE) + (2*stdcov_STE.SFE) +
      bse_STE + bse_SFE
    
    stor.data[i,1]  = stdvar_mRS / totalvariance
    stor.data[i,2]  = stdvar_F / totalvariance
    stor.data[i,3]  = stdvar_MS / totalvariance
    stor.data[i,4]  = (stdvar_STE-bse_STE) / totalvariance
    stor.data[i,5]  = (stdvar_SFE-bse_SFE) / totalvariance
     
    stor.data[i,6]  = bse_STE / totalvariance
    stor.data[i,7]  = bse_SFE / totalvariance
     
    stor.data[i,8]  = stdcov_F.MS / totalvariance
    stor.data[i,9]  = stdcov_F.STE / totalvariance
    stor.data[i,10] = stdcov_F.SFE / totalvariance
    stor.data[i,11] = stdcov_MS.STE / totalvariance
    stor.data[i,12] = stdcov_MS.SFE / totalvariance
    stor.data[i,13] = stdcov_STE.SFE / totalvariance
  }
  output <- data.frame(source=numeric(13), mean=numeric(13), "Q.025"=numeric(13), "Q.500"=numeric(13), "Q.975"=numeric(13))
  output[,1]  <- names(stor.data)
  
  for (i in 1:13){
    output[i,2] <- mean(stor.data[,i])
    output[i,3] <- quantile(stor.data[,i], .025)
    output[i,4] <- quantile(stor.data[,i], .500)
    output[i,5] <- quantile(stor.data[,i], .975)
  }
  return(output)
}


## Function to bootstrapp the variances to do pairwise comparisons
PWCOMP_VAR <- function (fit.comp1, fit.comp2, dataset, iteration){
  dataset <- FITCOMP (dataset)
  
  ##use the n=150 dataset if SFE is not involved in the pairwise comparison
  if (fit.comp1!="SFE" && fit.comp2!="SFE"){  
        
    N <- NROW(dataset)
    stor.data <- data.frame(var_fit.comp1=numeric(iteration), var_fit.comp2=numeric(iteration))

    for(i in 1:iteration){
      rand = sample(1:N, N, replace=TRUE)
      dataset$nb <- 1:N
      newdataset = dataset[match(rand, dataset$nb), ]

      ##estimate binomial sampling error in STE
      bse_STE <- BSE ("STE", newdataset)
      
      ##estimate standardised variance
      stdvar_fit.comp1 <- STDVAR (fit.comp1, newdataset)
      stdvar_fit.comp2 <- STDVAR (fit.comp2, newdataset)
      
      if (fit.comp1=="STE") { stdvar_fit.comp1 = stdvar_fit.comp1 - bse_STE }
      if (fit.comp2=="STE") { stdvar_fit.comp2 = stdvar_fit.comp2 - bse_STE }
      stor.data[i,1] = stdvar_fit.comp1
      stor.data[i,2] = stdvar_fit.comp2
    }
    
  ##use the n=147 dataset if SFE is involved in the pairwise comparison
  } else if(fit.comp1=="SFE" | fit.comp2=="SFE"){
    
    dataset <- dataset[complete.cases(dataset),]
    N <- NROW(dataset)
    stor.data <- data.frame(var_fit.comp1=numeric(iteration), var_fit.comp2=numeric(iteration))

    for(i in 1:iteration){
      rand = sample(1:N, N, replace=TRUE)
      dataset$nb <- 1:N
      newdataset = dataset[match(rand, dataset$nb), ]
        
      ##estimate binomial sampling error in STE and SFE
      bse_STE <- BSE ("STE", newdataset)
      bse_SFE <- BSE ("SFE", newdataset)
        
      ##estimate standardised variance
      stdvar_fit.comp1 <- STDVAR (fit.comp1, newdataset)
      stdvar_fit.comp2 <- STDVAR (fit.comp2, newdataset)
        
      if (fit.comp1=="STE") { stdvar_fit.comp1 = stdvar_fit.comp1 - bse_STE }
      if (fit.comp2=="STE") { stdvar_fit.comp2 = stdvar_fit.comp2 - bse_STE }
      if (fit.comp1=="SFE") { stdvar_fit.comp1 = stdvar_fit.comp1 - bse_SFE }
      if (fit.comp2=="SFE") { stdvar_fit.comp2 = stdvar_fit.comp2 - bse_SFE }
      stor.data[i,1] = stdvar_fit.comp1
      stor.data[i,2] = stdvar_fit.comp2
    }
  }
  return (min((2*sum(stor.data[,1]>stor.data[,2]))/iteration, (2*sum(stor.data[,1]<stor.data[,2]))/iteration))
}


## Functions to transform the fitness data to limit data skewness
TRANS_mRS <- function(x) { sqrt(x) }
TRANS_F   <- function(x) { sqrt(x+0.5) }
TRANS_MS  <- function(x) { x }
TRANS_STE <- function(x) { sqrt(x) }
TRANS_SFE <- function(x) { log10(x+1) }


## Function to compute repeatability of individual success across the three mating groups
RPT <- function (fit.comp, TRANS, dataset, iteration){
  dataset$trans.rel.fit.comp <- REL(TRANS(dataset[,fit.comp]))
  rpt_fit.comp <- rpt( trans.rel.fit.comp ~ (1|focal), grname="focal", CI=0.95, datatype="Gaussian", nboot=iteration, npermut=iteration, data=dataset)
  return (rpt_fit.comp)
}

## Function to compute the repeatable variance in a fitness component
RPTVAR <- function(fit.comp, TRANS, dataset_focal, dataset_group) {

  #compute standard variance
  dataset_focal <- dataset_focal[complete.cases(dataset_focal[,fit.comp]),]
  v_fit.comp <- STDVAR (fit.comp, dataset_focal)
  
  #compute repeatability
  rows <- sapply(dataset_focal$focal, function(x) which(x==dataset_group$focal))
  dataset_group = dataset_group[as.vector(rows),]
  rpt_fit.comp <- RPT (fit.comp, TRANS, dataset_group, 0)
  
  #compute repeatable variance
  rptvar_fit.comp <- v_fit.comp * rpt_fit.comp$R

  return(rptvar_fit.comp)
}


## Function to bootstrapp the repeatable variance to assess their 95 CI
## Outcomes are here expressed as percentages of total variance
BOOT_RPTVAR <- function(dataset_focal, dataset_group, iteration){
  stor.data <- data.frame(stdrptvar_mRS=numeric(iteration), stdrptvar_F=numeric(iteration), stdrptvar_MS=numeric(iteration), stdrptvar_STE=numeric(iteration), stdrptvar_SFE=numeric(iteration))
  
  for(i in 1:iteration){
    # estimate variance of mRS, F, MS and STE in boostrapped data (n=150)
    N <- NROW(dataset_focal) 
    rand = sample(1:N,N,replace=TRUE)
    dataset_focal$nb <- 1:N # assign new ID to each focal
    newdataset_focal = dataset_focal[match(rand, dataset_focal$nb), ]
    newdataset_focal$focal = rep(1:N)
    dataset_group$nb <- rep(1:N, each=3)
    rows <- sapply(rand, function(x) which(x==dataset_group$nb))
    newdataset_group = dataset_group[as.vector(rows),]
    newdataset_group$focal = rep(1:N,each=3)
    
    rptvar_mRS <- RPTVAR ("mRS", TRANS_mRS, newdataset_focal, newdataset_group)
    rptvar_F   <- RPTVAR ("F",   TRANS_F,   newdataset_focal, newdataset_group)
    rptvar_MS  <- RPTVAR ("MS",  TRANS_MS,  newdataset_focal, newdataset_group)
    rptvar_STE <- RPTVAR ("STE", TRANS_STE, newdataset_focal, newdataset_group)
    totvar150 <- TOTVAR (newdataset_focal)

    # estimate variance of SFE other boostrapped data (n=147)
    newdataset_focal <- dataset_focal[complete.cases(dataset_focal),]
    N <- NROW(newdataset_focal) 
    rand = sample(1:N, N, replace=TRUE)
    newdataset_focal$new.ID <- 1:N
    newdataset_focal = newdataset_focal[match(rand, newdataset_focal$new.ID), ]
    rows <- sapply(newdataset_focal$focal, function(x) which(x==dataset_group$focal))
    newdataset_group = dataset_group[as.vector(rows),]
    newdataset_focal$focal = rep(1:N)
    newdataset_group$focal = rep(1:N,each=3)

    rptvar_SFE <- RPTVAR ("SFE", TRANS_SFE, newdataset_focal, newdataset_group)
    totvar147  <- TOTVAR (newdataset_focal)
    
    stor.data[i,1] = rptvar_mRS
    stor.data[i,2] = rptvar_F / totvar150
    stor.data[i,3] = rptvar_MS / totvar150
    stor.data[i,4] = rptvar_STE / totvar150
    stor.data[i,5] = rptvar_SFE / totvar147
  }
  output <- data.frame(source=numeric(5), mean=numeric(5), "Q.025"=numeric(5), "Q.500"=numeric(5), "Q.975"=numeric(5))
  output[,1]  <- names(stor.data)
  
  for (i in 1:5){
    output[i,2] <- mean(stor.data[,i])
    output[i,3] <- quantile(stor.data[,i], .025)
    output[i,4] <- quantile(stor.data[,i], .500)
    output[i,5] <- quantile(stor.data[,i], .975)
  }
  return(output)
}


PWCOMP_RPTVAR <- function (fit.comp1, fit.comp2, TRANS1, TRANS2, dataset_focal, dataset_group, iteration){
  
  ##use the n=150 dataset if SFE is not involved in the pairwise comparison
  if (fit.comp1!="SFE" && fit.comp2!="SFE"){  
    
    N <- NROW(dataset_focal)
    stor.data <- data.frame( rptvar_fit.comp1=numeric(iteration), rptvar_fit.comp2=numeric(iteration) )
    
    for(i in 1:iteration){
      rand = sample(1:N, N, replace=TRUE)
      dataset_focal$new.ID <- 1:N
      newdataset_focal = dataset_focal[match(rand, dataset_focal$new.ID), ]
      rows <- sapply(newdataset_focal$focal, function(x) which(x==dataset_group$focal))
      newdataset_group = dataset_group[as.vector(rows),]
      newdataset_focal$focal = rep(1:N)
      newdataset_group$focal = rep(1:N,each=3)
      
      rptvar_fit.comp1 <- RPTVAR (fit.comp1, TRANS1, newdataset_focal, newdataset_group)
      rptvar_fit.comp2 <- RPTVAR (fit.comp2, TRANS2, newdataset_focal, newdataset_group)

      stor.data[i,1] = rptvar_fit.comp1
      stor.data[i,2] = rptvar_fit.comp2
    }
    
  ##use the n=147 dataset if SFE is involved in the pairwise comparison
  } else if(fit.comp1=="SFE" | fit.comp2=="SFE"){
    dataset_focal <- dataset_focal[complete.cases(dataset_focal),]
    N <- NROW(dataset_focal)
    stor.data <- data.frame( rptvar_fit.comp1=numeric(iteration), rptvar_fit.comp2=numeric(iteration) )
    
    for(i in 1:iteration){
      rand = sample(1:N, N, replace=TRUE)
      dataset_focal$new.ID <- 1:N
      newdataset_focal = dataset_focal[match(rand, dataset_focal$new.ID), ]
      rows <- sapply(newdataset_focal$focal, function(x) which(x==dataset_group$focal))
      newdataset_group = dataset_group[as.vector(rows),]
      newdataset_focal$focal = rep(1:N)
      newdataset_group$focal = rep(1:N,each=3)
      
      rptvar_fit.comp1 <- RPTVAR (fit.comp1, TRANS1, newdataset_focal, newdataset_group)
      rptvar_fit.comp2 <- RPTVAR (fit.comp2, TRANS2, newdataset_focal, newdataset_group)

      stor.data[i,1] = rptvar_fit.comp1
      stor.data[i,2] = rptvar_fit.comp2
    }
  }
  list <- list("P value (repeatable variance)"=min((2*sum(stor.data[,1]>stor.data[,2]))/iteration, (2*sum(stor.data[,1]<stor.data[,2]))/iteration))
  return(list)
}