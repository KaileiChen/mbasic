#' @name MBASIC.MADBayes_zeroInflated
#' @title MAD-Bayes method to fit the MBASIC model with zero and one inflation for non-enriched state (s=0).
#' @param Y An N by I matrix containing the original data (before transformation, if there is) from N experiments across I observation units (loci).
#' @param Gamma An N by I matrix for the prior estimated mean for the background state, for N experiments across the I observation units (loci).
#' @param fac A vector of length N denoting the experimental condition for each replicate.
#' @param lambdaw,lambda Tuning parameters.
#' @param family The distribution of family to be used. Either 'lognormal' or 'negbin'. See details for more information.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param S The number of different states.
#' @param verbose Boolean variable for whether the model fitting messages are printed.
#' @param para A list of true paramters.
#' @param zeroInflated whether zeros and ones are inflated
#' @details
#' TODO.
#' @useDynLib MBASIC madbayes_zeroInflated
#' @useDynLib MBASIC madbayes_theta_zeroInflated
#' @return A list object.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @export
MBASIC.MADBayes_zeroInflated <- function(Y, Gamma, fac, lambdaw = 0.2, lambda = 200, maxitr = 100, S = 2, tol = 1e-6, verbose = TRUE, para = NULL, initialize = "kmeans", family = "lognormal", zeroInflated = TRUE) {
	if (zeroInflated) {
  		inflation = getZeroOneTwo(Y)
  	} else {
		inflation = NULL
	}
	
  ## Transform and normalize Data

  if (family == "lognormal"){
	Y <- log(Y+1)
	NormalizeData()
  } else if (family == "poisson"){
	Y <- sqrt(Y + 0.25)
	NormalizeData()
  } else if (family == "binom" & is.null(Gamma) == FALSE) {
	Y <- sqrt(Gamma + 0.5)*asin(sqrt((Y + 0.375) / (Gamma + 0.75) ))
	Gamma <- sqrt(Gamma + 0.5)
	Gamma <- matrix(rep(t(Gamma), S), ncol = ncol(Y), byrow =TRUE)
	zeroInflated = FALSE
  }
 
  
  	    
#	print("internal")
  fit <- MBASIC.MADBayes.internal_zeroInflated(Y = Y, Gamma = Gamma, fac = fac, lambdaw = lambdaw, lambda = lambda, maxitr = maxitr, S = S, tol = tol, verbose = verbose, para = para, initialize = initialize, Theta.init = NULL, Mu.init = NULL, Sigma.init = NULL, clusterLabels.init = NULL, scaleFactor = scaleFactor, family = family, zeroInflated = zeroInflated, inflation = inflation, inflated.prop = NULL)
  return(fit)
}




#' @importFrom cluster silhouette
MBASIC.MADBayes.internal_zeroInflated <- function(Y, Gamma, fac, lambdaw = NULL, lambda, maxitr = 20, S, tol = 1e-8, verbose = TRUE, para = NULL, initialize = "kmeans", Theta.init = NULL, Mu.init = NULL, Sigma.init = NULL, clusterLabels.init = NULL, scaleFactor = NULL, J = NULL, family = "lognormal", zeroInflated = TRUE, inflation = NULL, inflated.prop = NULL) {

#  print("interna")
  if (is.null(inflation) & zeroInflated){
	print("Inflation counts are missing!")
	zeroInflated = FALSE
  }

  GetModelStructure()

  facNames <- as.character(unique(fac))
  facMap <- seq(K)
  names(facMap) <- facNames
  fac <- as.character(fac)
  fac <- facMap[as.character(fac)]
  
  ## Initialize Mu, Sigma, Theta

  if(is.null(Mu.init) | is.null(Sigma.init) | is.null(Theta.init)) {
    InitializeTheta.MADBayes(zeroInflated)
  } else {
    Theta <- Theta.init
    Mu <- Mu.init
    Sigma <- Sigma.init
  }

#  print(inflated.prop)

  if(is.null(lambdaw)) {
    lambdaw <- 2
  }
  
  lambdaw <- lambdaw * mean(na.omit(Sigma))
  lambda <- lambda * lambdaw
  
  ## Initialize cluster
  if(is.null(clusterLabels.init)) {
    InitializeClusters.MADBayes()
  } else {
    clusterLabels <- clusterLabels.init
  }
  
  b <- rep(0, I)
  
#  print(dim(Mu))
  ##  zeta <- mean(b)
  if (zeroInflated) {
    storage.mode(inflation) <- "integer"
	p0 <- inflated.prop[,1]
	p1 <- inflated.prop[,2]
	ret <- .Call("madbayes_zeroInflated", clusterLabels, Theta, Mu, D, Gamma, Y, lambdaw, lambda, maxitr, tol, inflation, p0, p1, package = "MBASIC")
#	ret <- .Call("madbayes", clusterLabels, Theta, Mu, D, Gamma, Y, lambdaw, lambda, maxitr, tol, inflation, p0, p1, package = "MBASIC")
	inflated = cbind(p0, p1)
  } else {
	ret <- .Call("madbayes", clusterLabels, Theta, Mu, D, Gamma, Y, lambdaw, lambda, maxitr, tol, package = "MBASIC")
  }


#  browser()
  Mu = ret$Mu
  Theta = ret$Theta + 1
  if (zeroInflated) {
	inflated.prop = cbind(ret$p0, ret$p1)
  } else {
	inflated.prop = NULL
  }
#  save(ret, file = "Cppresult.RData")


  t0 <- Sys.time()
  if(verbose)
    message("Finished iterations.")

  if(maxitr <= ret$Iter) {
    warning("MADBAYES procedure not converged.")
  }

#  print("begin")
  J <- max(ret$clusterLabels) + 1
  Z <- matrix(0, nrow = I, ncol = J)
  Z[cbind(seq(I), ret$clusterLabels + 1)] <- 1


  W <- ret$W[, seq(max(ret$clusterLabels) + 1), drop = FALSE]

#  print("JZW")  
#  browser()
  Theta.err <- W.err <- Mu.err <- ari <- mcr <- errPreTransformation <- numeric(0)
  if(!is.null(para)) {
    Theta.err <- sqrt(2 * sum(para$Theta != (ret$Theta + 1)) / I / K / S)
    W.f <- matrix(0, nrow = K * S, ncol = J)
    for(s in seq_len(S))
      W.f[ s + S * seq(0, K - 1), ] <- W[ seq_len(K) + K * (s - 1), ]

#	browser()
    mc <- matchCluster(W.f, para$W, Z, para$Z, ret$b, para$non.id)
    W.err <- mc$W.err
    mcr <- mc$mcr
    ## recompute ARI
    trueLabels <- apply(para$Z, 1, which.max)
    trueLabels[para$non.id] <- max(trueLabels) + seq(length(para$non.id))
    ari <- adjustedRandIndex(ret$clusterLabels, trueLabels)
  }

#	print("ari finish")
  ## Compute the loss of each term
  Theta.aug <- W.aug <- matrix(0, nrow = K * S, ncol = I)
  for(s in seq(S)) {
    Theta.aug[seq(K) + K * (s-1), ] <- as.integer(ret$Theta == s - 1)
  }
  W.aug <- tcrossprod(W, Z)
  loss.w <- mean((Theta.aug - W.aug) ^ 2)
  

  ## compute the loss of data fitting
  Theta.Y <- designMat %*% ret$Theta
  Mu.Y <- Sigma.Y <- Y - Y
  for(s in seq(S)) {
    id <- which(Theta.Y == s - 1)
    Mu.Y[id] <- (Gamma[seq(N) + N * (s - 1), ] * rep(Mu[, s], I))[id]
    Sigma.Y[id] <- rep(Sigma[, s], I)[id]
  }
  if (zeroInflated) {
	n0 = apply(ret$Theta==0, 1, sum)
	n0 = round(n0[fac]  * inflated.prop[,1])
	n1 = apply(ret$Theta==1, 1, sum)
	n1 = round(n1[fac] * inflated.prop[,2])
	tempY <- unlist(lapply(1:nrow(Y), function(rowI) {
									ans = Y[rowI,]
									index0 <- which(inflation[rowI, ] == 1)
									index1 <- which(inflation[rowI, ] == 2)
									if (n0[rowI] < length(index0))
										index0 = index0[1:n0[rowI]]
									if (n1[rowI] < length(index1))
										index1 = index1[1:n1[rowI]]	
									ans[index0] <- Mu.Y[rowI, index0]
									ans[index1] <- Mu.Y[rowI, index1]
									return(ans)
								})
			)
	tempY = matrix(tempY, nrow = nrow(Y), byrow = TRUE)
	loss.y <- mean((tempY - Mu.Y) ^ 2)
  } else {
	loss.y <- mean((Y - Mu.Y) ^ 2)
  }
#print("loss finished")
  ## Add normalized data
  ## Y: I * N, each column has mean 1
  Theta.norm <- matrix(0, ncol = I, nrow = K * S)
  PDF <- matrix(0, ncol = I, nrow = N * S)
  if (zeroInflated)
	  pPoisson = 1 - ret$p0 - ret$p1
  for(s in seq(S)) {
    if (zeroInflated & s == 1){
		denY <- lapply(1:N, function(set) {
								index0 <- which(inflation[set, ] == 1)
								index1 <- which(inflation[set, ] == 2)
								ans = dnorm(Y[set,], mean = Mu[set, s]*Gamma[set,], sd = sqrt(Sigma[set, s]), log = FALSE)
								ans[index0] <- pPoisson[set] * ans[index0] + ret$p0[set]
								ans[index1] <- pPoisson[set] * ans[index1] + ret$p1[set]
								ans = log(ans) 
								return(ans)
							}
						)
		denY <- matrix(unlist(denY), ncol = I, byrow = TRUE)
	} else {
		denY <- dnorm(Y, mean = Mu[, s]*Gamma[1:N,], sd = sqrt(Sigma[, s]), log = TRUE)
#		denY <- dbinom(para$Y, para$X, para$Mu[,s], log = TRUE)
	}
    
    PDF[seq(N) + N * (s - 1), ] <- denY
    Theta.norm[seq(K) + K * (s - 1), ] <- crossprod(designMat, denY) / apply(designMat, 2, sum)
  }

  ## for any i k, max_s Theta.norm[i, k, s] = 0
  Theta.norm <- Theta.norm - t(matrix(rep(apply(matrix(t(Theta.norm), nrow = I * K), 1, max), S), nrow = I))
  ## avoid getting Inf exponenant
  Theta.norm[Theta.norm > 5] <- 5
  Theta.norm <- exp(Theta.norm)
  Theta.total <- matrix(0, ncol = I, nrow = K)
  for(s in seq(S)) {
    Theta.total <- Theta.total + Theta.norm[(s - 1) * K + seq(K), ]
  }
  Theta.norm <- Theta.norm / t(matrix(rep(t(Theta.total), S), nrow = I))

  if(J == I | var(ret$clusterLabels) == 0) {
    sil.norm <- 0
  } else {
    dist.norm <- dist(t(Theta.norm), method = "manhattan")
    sil.norm <- mean(cluster::silhouette(ret$clusterLabels, dist.norm)[, 3])
#	browser()
#	sil.norm <- mean(adjusted.silhouette(ret$clusterLabels+1, dist.norm, ret$Theta+1, W))
  }
  print("sil finished")

  P <- matrix(1 / S, ncol = S, nrow = I)
  V <- matrix(1, nrow = N, ncol = S)
  probz <- apply(Z, 2, mean)
  loglik <- .Call("loglik", W, P, V, 1e-10, probz, PDF, fac - 1, seq(S) - 1, package = "MBASIC")
  npars <- ncol(Z) - 1 + prod(dim(W)) * (S - 1) / S + N * S * 2

  Theta <- ret$Theta + 1
  rownames(Theta) <- facNames
  rownames(W) <- rep(facNames, S)

  if (family == "lognormal"){
	Mu <- ret$Mu * scaleFactor
  	Sigma <- ret$Sigma * scaleFactor * scaleFactor
  }
  rownames(Mu) <- rownames(Sigma) <- rownames(V) <- facNames[fac]
#  browser()
  if (!is.null(para)){
	  if (family == "lognormal"){
		Mu.err <- sqrt(mean((Mu - para$Mu) ^ 2))
		errPreTransformation = NA
  	} else if (family == "poisson"){
		Mu.err <- sqrt(mean((Mu - para$Mu) ^ 2))
		poisson_mean <- Mu * Mu - .25
    	errPreTransformation <- sqrt(mean((para$Lambda - poisson_mean) ^ 2))
	  } else if (family == "binom"){
		Mu.err = numeric(0)
		temp = Gamma[1:N,]
		prob = Y / temp
		temp  = temp* temp - 0.5
		prob = sin(prob) * sqrt(temp + 0.75)
		prob = prob * prob / temp
		binom_mean <- lapply(1:N, function(set) {
								sapply(1:S, function(s){
												return(mean(prob[set, Theta[fac[set], ] == s]))
								})
						})
		binom_mean <- matrix(unlist(binom_mean), ncol = S, byrow = TRUE)
		errPreTransformation <- sqrt(mean((para$Mu - binom_mean) ^ 2))
  	}
  }
  
  ans <- new("MBASICFit",
      Theta = Theta,
      W = W,
      clustProb = cbind(0, Z),
      alllik = ret$loss,
      Mu = Mu,
      Sigma = Sigma,
      converged = (ret$Iter <= maxitr),
      Z = Z,
      Iter = ret$Iter,
	  Mu.err = Mu.err,
      Theta.err = Theta.err,
      W.err = W.err,
      ARI = ari,
      MisClassRate = mcr,
      Loss = list(
		inflated.prop = inflated.prop,
		errPreTransformation = errPreTransformation,
        lambdaw = lambdaw,
        lambda = lambda,
        Y = loss.y,
        W = loss.w,
        Silhouette = sil.norm,
        loglik = loglik,
        bic = -2 * loglik + log(N * I) * npars)
    )
#	browser()

	return(ans)
}

#' @name MBASIC.MADBayes.full_zeroInflated
#' @title MAD-Bayes method to fit the MBASIC model.
#' @param Y An N by I matrix containing the data from N experiments across I observation units (loci).
#' @param Gamma An N by I matrix for the prior estimated mean for the background state, for N experiments across the I observation units (loci).
#' @param fac A vector of length N denoting the experimental condition for each replicate.
#' @param lambdap,lambdaw,lambda Tuning parameters.
#' @param family The distribution of family to be used. Either "lognormal" or "negbin". See details for more information.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param S The number of different states.
#' @param ncores The number of CPUs to be used for parallelization.
#' @param nfits The number of random restarts of the model.
#' @details
#' TODO.
#' @useDynLib MBASIC
#' @return A list object including the following fields:
#' \tabular{ll}{
#' allFits \tab A list of \linkS4class{MBASICFit} objects for the best model fit with each lambda.\cr
#' lambda \tab A vector of all lambdas corresponding to \code{allFits}.\cr
#' Loss \tab A vector for the loss corresponding to \code{allFits}.\cr
#' BestFit \tab The \linkS4class{MBASICFit} object with largest Silhouette score.\cr
#' Iter \tab Number of iterations for \code{BestFit}.\cr
#' Time \tab Time in seconds used to fit the model.\cr
#' }
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @import foreach
#' @export
MBASIC.MADBayes.full_zeroInflated <- function(Y, Gamma = NULL, fac, lambdaw = NULL, lambda = NULL, maxitr = 30, S = 2, tol = 1e-10, ncores = 15, nfits = 1, nlambdas = 30, para = NULL, initialize = "kmeans", family = "lognormal", zeroInflated = TRUE, init = NULL) {
  t0 <- Sys.time()
  if(!is.null(lambda)) {
    ncores <- min(c(ncores, length(lambda) * nfits))
  }
  
  startParallel(ncores)
  
  if(!is.null(lambda)) {
    message("given lambdas")
    lambdas <- sort(unique(lambda))
    alllambdas <- rep(lambdas, each = nfits)
    results <- foreach(i = seq_along(alllambdas)) %dopar% {
      set.seed(i + Sys.time())
      fit <- MBASIC.MADBayes_zeroInflated(Y, Gamma, fac, lambdaw = lambdaw, lambda = alllambdas[i], maxitr = maxitr, S = S, tol = tol, verbose = TRUE, para = para, initialize = initialize, family = family, zeroInflated = zeroInflated)
      list(fit = fit, lambda = alllambdas[i])
    }
    initLosses <- NULL
  } else {

    if(!is.numeric(nlambdas)) {
      stop("Error: 'nlambdas' must take a numeric value.")
    }

	if (zeroInflated) {
    	inflation = getZeroOneTwo(Y)
    } else {
		inflation = NULL
    }
	
  ## Transform and normalize Data

    if (family == "lognormal"){
		Y <- log(Y+1)
		NormalizeData()
    } else if (family == "poisson"){
		Y <- sqrt(Y + 0.25)
		NormalizeData()
    } else if (family == "binom" & is.null(Gamma) == FALSE) {
		Y <- sqrt(Gamma + 0.5)*asin(sqrt((Y + 0.375) / (Gamma + 0.75) ))
		Gamma <- sqrt(Gamma + 0.5)
		Gamma <- matrix(rep(t(Gamma), S), ncol = ncol(Y), byrow =TRUE)
		zeroInflated = FALSE
    }

    
    GetModelStructure()
    ## Initialize Theta, Mu, Sigma
    InitializeTheta.MADBayes(zeroInflated)

#    print(inflated.prop)

    ## Initialize clusters and pick a range of lambda values
    if(initialize == "madbayes") {
      ret <- foreach(i = seq(as.integer(sqrt(I)) + 1)) %dopar% {
        .Call("madbayes_init", Theta, 0, S, i, package = "MBASIC")
      }
    } else if(initialize == "kmeans++") {
      ret <- foreach(i = seq(as.integer(sqrt(I)) + 1)) %dopar% {
        .Call("madbayes_init_kmeanspp", Theta, S, i, package = "MBASIC")
      }
    } else if(initialize == "kmeans") {
      Theta.aug <- matrix(0, nrow = K * S, ncol = I)
      for(s in seq(S)) {
        Theta.aug[seq(K) + (s - 1) * K, ] <- (Theta == s)
      }
      ret <- foreach(i = seq(as.integer(sqrt(I)) + 1)) %dopar% {
        fit.kmeans <- kmeans(t(Theta.aug), i)
        return(list(loss = fit.kmeans$tot.withinss * 2, clusterLabels = fit.kmeans$cluster - 1))
      }
    } else {
      stop("Error: a vector for 'lambda' values must be provided.")
    }

    endParallel()
    
    message("Initialized clusters")
    
    initLosses <- rep(0, length(ret))
    allClusterLabels <- list()
    for(i in seq_along(ret)) {
      initLosses[i] = ret[[i]]$loss
      allClusterLabels[[i]] <- ret[[i]]$clusterLabels
    }

    slopes <- abs(diff(sort(initLosses, decreasing = TRUE)))
    slopes <- slopes[-1]
    allLambdas <- (slopes[-1] + slopes[-length(slopes)]) / 2
    allLambdas <- sort(allLambdas, decreasing = TRUE)
    minLambda <- min(allLambdas)
    maxLambda <- max(allLambdas)
    ## lambdas <- seq(minLambda, maxLambda, length = nlambdas)
    lambdas <- unique(quantile(allLambdas, seq(nlambdas) / (nlambdas + 2)))
    lambdas <- lambdas[-c(1, length(lambdas))]
    alllambdas <- rep(lambdas, each = nfits)

    initClusterLabels <- list()
    usedids <- numeric(0)
    freeids <- seq_along(allLambdas)
    allJs <- seq_along(alllambdas)
    for(i in seq_along(alllambdas)) {
      j <- freeids[which.min(abs(allLambdas[freeids] - alllambdas[i]))[1]]
      initClusterLabels[[i]] <- allClusterLabels[[j + 2]]
      usedids <- c(usedids, j)
      freeids <- setdiff(freeids, usedids)
      j <- which.min(abs(allLambdas - alllambdas[i]))[1]
      allJs[i] <- max(allClusterLabels[[j + 2]]) + 1
    }
#	print(allJs)
    message("Selected lambda values: ", paste(allLambdas, collapse = ", "))
    
    results <- foreach(i = seq_along(alllambdas)) %dopar% {
      set.seed(i + Sys.time())
      fit <- MBASIC.MADBayes.internal_zeroInflated(Y, Gamma, fac, lambdaw = lambdaw, lambda = alllambdas[i], maxitr = maxitr, S = S, tol = tol, verbose = TRUE, para = para, initialize = initialize, Theta.init = Theta, Mu.init = Mu, Sigma.init = Sigma, scaleFactor = scaleFactor, J = allJs[i], zeroInflated = zeroInflated, inflation = inflation, inflated.prop = inflated.prop, family = family)
      list(fit = fit, lambda = alllambdas[i])
    }
  }
  

  if(.Platform$OS.type != "unix") {
    stopCluster(cl)
  }
  
  message("Finished individual models")
  ## Within the same lambda, choose the model that minimizes the loss
  bestLosses <- rep(Inf, length(lambdas))
  bestFits <- list()
  bestIters <- rep(0, length(lambdas))
  for(i in seq_along(results)) {
    lambdaid <- which(lambdas == alllambdas[i])[1]
    if(bestLosses[lambdaid] > tail(results[[i]]$fit@alllik, 1)) {
      bestLosses[lambdaid] <- tail(results[[i]]$fit@alllik, 1)
      bestFits[[lambdaid]] <- results[[i]]$fit
      bestIters[lambdaid] <- results[[i]]$fit@Iter
    }
  }
  
  ## Between different lambdas, choose the model with the largest Silhouette score
  bestSil <- -Inf
  bestFit <- NULL
  for(fit in bestFits) {
    if(fit@Loss$Silhouette > bestSil) {
      bestFit <- fit
      bestSil <- fit@Loss$Silhouette
    }
  }

  bestbic <- Inf
  bestFit.bic <- NULL
  for(fit in bestFits) {
    if(fit@Loss$bic < bestbic) {
      bestFit.bic <- fit
      bestbic <- fit@Loss$bic
    }
  }
   
  return(list(allFits = bestFits,
              BestFit = bestFit,
              BestFit.bic = bestFit.bic,
              Iter = bestIters,
              Loss = bestLosses,
              lambda = lambdas,
              Time = as.numeric(Sys.time() - t0, units = "secs"),
              InitLoss = initLosses)
         )
}

InitializeTheta.MADBayes <- function(zeroInflated = FALSE) {
  Inherit(c("N", "S", "K", "I", "Y", "Gamma", "designMat", "D", "maxitr", "tol", "inflation", "family"))
  Mu <- Sigma <- matrix(0, nrow = N, ncol = S)
  Theta <- matrix(0, nrow = K, ncol = I)
  storage.mode(Theta) <- "integer"
  foldChange <- Y / Gamma[seq(N), ]
  if (family != "binom") 
	foldChange[Gamma[seq(N), ] == 0] <- 1
  if (family == "binom"){
	foldChange[Gamma[seq(N), ] < sqrt(10+0.5)] <- asin(sqrt(0.5))
	Y[Gamma[seq(N), ] <= sqrt(5+0.5)] = Gamma[seq(N),][Gamma[seq(N), ] <= sqrt(5+0.5)]*asin(sqrt(0.5))
  }
	

  avgFoldChange <- crossprod(foldChange, designMat) / rep(apply(designMat, 2, sum), each = I)
    

#  for(k in seq(K)) {
#    for(s in seq(S, 1)) {
#      Theta[k, rank(avgFoldChange[, k]) <= s / S * I] <- s - 1
#    }
#  }

  
  for(k in seq(K)) {
    up <- max(avgFoldChange[,k])
    down <- min(avgFoldChange[,k])
    for(s in seq(S, 1)) {
      cutoff = (up-down) * s / S + down
      Theta[k, avgFoldChange[, k] <= cutoff] <- s - 1
    }
  }

  ZeroOne = 0
  if (zeroInflated) {
  	ZeroOne <- crossprod(designMat, (inflation == 1 | inflation == 2) * 1)

	# set loci with ones and zeros apart 
  	Theta[ZeroOne > 0] <- -1
  }

  ## DTheta: N by I
  DTheta <- designMat %*% Theta
  for(n in seq(N)) {
    for(s in seq(S)) {
      Mu[n, s] <- mean(foldChange[n, DTheta[n, ] == s - 1])
      Sigma[n, s] <- var(foldChange[n, DTheta[n, ] == s - 1])
    }
  }
  Sigma[Sigma <= 0] <- min(Sigma[Sigma > 0])

  p0 = p1 = 0
  if (zeroInflated) {
	p0 <- (apply(inflation == 1, 1, sum) + 1) / (apply(DTheta == 0, 1, sum) + 8)
	p1 <- (apply(inflation == 2, 1, sum) + 1) / (apply(DTheta == 0, 1,  sum) + 8)
	p0 <- p0 / (1 + p0) * 0.95
	p1 <- p1 / (1 + p1) * 0.95
  }
  
# put the zero inflated back
  Theta[ZeroOne > 0] <- 0

  storage.mode(Theta) <- "integer"
  
  ## initialize Theta
  
  if (zeroInflated) {
	ret <- .Call("madbayes_theta_zeroInflated", Theta, Mu, D, Gamma, Y, maxitr, tol, inflation, p0, p1, package = "MBASIC")
#	ret <- .Call("madbayes_theta", Theta, Mu, D, Gamma, Y, maxitr, tol, Zero, One, two, p0, p1, package = "MBASIC")
  } else {
	ret <- .Call("madbayes_theta", Theta, Mu, D, Gamma, Y, maxitr, tol, package = "MBASIC")
  }

  Theta <- ret$Theta
  Mu <- ret$Mu
  Sigma <- ret$Sigma
  Sigma[Sigma <= 0] <- min(Sigma[Sigma > 0])

#  Mu[,1] <- max(Mu[,1], 1.41)

  if (zeroInflated) {
    p0 = ret$p0
    p1 = ret$p1
    inflated.prop = matrix(c(p0, p1), ncol = 2)
  } else {
    inflated.prop = NULL
  }
  Return(c("Theta", "Mu", "Sigma", "inflated.prop"))
}

InitializeClusters.MADBayes <- function() {
  Inherit(c("verbose", "I", "initialize", "K", "S", "Theta", "lambda", "lambdaw", "J"))
  if(verbose)
    message("Initialize clusters...")
  
  ## Sample J from an exponential distribution with the median sqrt(I)/4
  if(is.null(J)) {
    J <- sample(seq(I)[-1], 1, prob = exp(-seq(I)/sqrt(I)*4*log(2))[-1])
  } else  {
    J <- as.integer(J)
  }

  if(initialize == "kmeans") {
    ## use K-means to initialize clusters
    Theta.aug <- matrix(0, nrow = K * S, ncol = I)
    for(s in seq(S)) {
      Theta.aug[seq(K) + (s - 1) * K, ] <- (Theta == s)
    }
    ##  clusterLabels <- sample(seq(J), I, replace = TRUE) - 1
    clusterLabels <- kmeans(t(Theta.aug), centers = J)$cluster - 1
  } else if(initialize == "kmeans++" ) {
    ret <- .Call("madbayes_init_kmeanspp", Theta, S, J, package = "MBASIC")
    clusterLabels <- ret$clusterLabels
  } else if(initialize == "madbayes") {
    ret <- .Call("madbayes_init", Theta, lambda / lambdaw, S, I, package = "MBASIC")
    clusterLabels <- ret$clusterLabels
  } else {
    clusterLabels <- sample(seq(J), I, replace = TRUE)
  }
  Return("clusterLabels")
}

GetModelStructure <- function() {

  Inherit(c("fac", "Y", "Gamma"))
  ## prespecified
  K <- length(unique(fac))
  I <- ncol(Y)
  N <- nrow(Y)
  if(length(fac) != N)
    stop("Error: total number of replicates do not match with the number of rows in Y")

  ## design matrix is N by K
  designMat <- matrix(0, ncol = K, nrow = N)
  for(k in 1:K) {
    designMat[fac == unique(fac)[k], k] <- 1
  }
  D <- apply(designMat, 1, function(x) which(x == 1)) - 1
  storage.mode(D) <- "integer"

  Return(c("K", "I", "N", "designMat", "D"))
}

NormalizeData <- function() {
  Inherit(c("Y", "Gamma", "S"))
  ## Scale the data from different replicates
  scaleFactor <- apply(Y, 1, mean)
  N <- length(scaleFactor)
  I <- ncol(Y)
  
  if(is.null(Gamma)) {
    Gamma <- matrix(1, nrow = nrow(Y), ncol = ncol(Y))
  }
  
  if(prod(dim(Y) == dim(Gamma)) != 1) {
    stop("Error: dimensions for Y and Gamma must be the same.")
  }
  
  ## normalize the Gamma
  Gamma <- Gamma / rep(apply(Gamma, 1, mean), ncol(Gamma))
  Y <- Y / scaleFactor
  Gamma <- rbind(Gamma, matrix(0, ncol = ncol(Gamma), nrow = nrow(Gamma) * (S - 1)))
  for(s in seq(S)[-1]) {
    Gamma[(s - 1) * N + seq(N), ] <- rep(apply(Gamma[seq(N), ], 1, mean), I)
  }

  Return(c("Y", "Gamma", "scaleFactor"))

}


getZeroOneTwo <- function(Y) {
  Zero <- (Y == 0) * 1

  One <- (Y == 1) * 2

  Two <- (Y == 2) * 4

  ans  = Zero + One + Two

  return(ans)
}


Inherit <- function(vnames = NULL) {
  if(is.null(vnames)) {
    vnames <- ls(envir = parent.frame(2))
  }
  for(v in vnames) {
    ## do not use get() since it will cause error for non-defined variables
    assign(v, parent.frame(2)[[v]], envir = parent.frame())
  }
}

Return <- function(vnames = NULL) {
  if(is.null(vnames)) {
    vnames <- ls(envir = parent.frame())
  }
  for(v in vnames) {
    ## do not use get() since it will cause error for non-defined variables
    assign(v, parent.frame()[[v]], envir = parent.frame(2))
  }
}

	


adjusted.silhouette <- function(labels, dist.norm, Theta, W){
	n <- attr(dist.norm, "Size")
	K <- nrow(Theta)
	m <- as.matrix(dist.norm)

#	browser()


	ai <- sapply(1:n, function(index){
							distVec <- m[index, labels == labels[index]]
#							distVec <- distVec[distVec > 0]
							proportion <- prod(W[(K*(Theta[,index]-1) + 1:K),labels[index]])
							if (proportion < 5 / length(distVec)){
									ans = median(distVec)
							} else {
									ans = quantile(distVec, proportion / 2) 
							}
							return(ans)
				})

	bi <- lapply(1:n, function(index){
							lowest <- sapply(1:max(labels), function(cl){
													if (cl == labels[index]){
														ans = Inf
													} else { 
														distVec <- m[index, labels == cl]
														proportion <- prod(W[(K*(Theta[,index]-1) + 1:K),cl])
														if (proportion < 5 / sum(labels == cl)){
															ans = mean(distVec)
														} else {
															ans = quantile(distVec, proportion / 2) 
														}
													}
													return(ans)
											})
							ans = min(lowest)
							return(ans)
				})

	si <- rep(0, n)
	den = pmax(as.numeric(ai), as.numeric(bi))
	si[den > 1.0e-6]  = (as.numeric(bi)-as.numeric(ai))[den > 1.0e-6] / den[den > 1.0e-6]
	return(si)
}
#summary(adjusted.silhouette(labels, dist.norm, Theta, W))


binarizeTheta <- function(vec){ return(sum((vec-1)*S^(1:K-1)))}

# ans = apply(Theta, 2, binarizeTheta) 

	
