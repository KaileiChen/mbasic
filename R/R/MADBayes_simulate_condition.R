#' @name madbayes.sim
#' @title Simulate data for the general MBASIC model.
#' @param xi Parameter for the magnitude of each observations. See details for more information.
#' @param family A parameter for the family of distribution to be used, must be "lognormal", "negbin" or "binom". Default: "lognormal".
#' @param f A numerical value that determine the difference of the means between different states. See details for more information. Default: 5.
#' @param I An integer for the total number of units.
#' @param fac A vector of length N denoting the experimental condition for each replicate.
#' @param J An integer for the number of clusters.
#' @param S An integer for the number of states. Default: 2.
#' @param struct An K by J integer matrix. The j-th column denotes the levels for the cluster level parameter. See details for more information. Default: NULL.
#' @param delta A vector of length S, or NULL. This is the dirichlet prior parameter to simulate the probability across the S states for each CLUSTERED unit and each experiment. If NULL, rep(0.1,S) is used.
#' @param delta.non  A vector of length S, or NULL. This is the dirichlet prior parameter to simulate the probability across the S states for each UNCLUSTERED unit and each experiment. If NULL, rep(0.1,S) is used.
#' @param zeta The probability that each unit does not belong to any cluster. Default: 0.1.
#' @details
#' MBASIC.sim allows three types of distributions:\cr
#' For the "lognormal" family, entries in the matrix Y follows distribution: log(Y[n,i] + 1) | Theta[n,i]=s ~ N(Mu[n,s], stdev[s]).\cr
#' For the "negbin" family, entries in the matrix Y follows distribution: Y[n,i] | Theta[n,i]=s ~ NB(Mu[n,s], stdev[s]).\cr
#' For the "binom" family, entries in the matrices X and Y follows distribution: X[n,i] ~ Pois(xi), Y[n,i]| Theta[n,i]=s,X[n,i] ~ Binom(X[n,i],Mu[n,s]). In this package, NB(mu,size) denotes a Negative binomial distribution with mean mu and variance mu(1+mu/size).\cr
#' For "lognormal" or "negbin" families, Mu[n,s]~N(prior.mean[s],prior.sd[s]). Hyper paramters prior.mean and prior.sd are set differently under the two distributional families. For the "lognormal" family, where prior.mean[s] = xi+log((s-1)(f-1)+1), and prior.sd=log(f)/30. For the "negbin" family, prior.mean[s]=xi*((s-1)(f-1)+1), and prior.sd=(f-1)*xi/6. In general, xi is the mean for the state S=1, and f is roughly the ratio between the means from state S=2 and S=1.\cr
#' For the "binom" family, Mu[n,s] ～ Beta(s * f, (S + 1 - s) * f).
#' @return A list containing:
#' \tabular{ll}{
#'  Y \tab An N by I matrix. The (n,i)-th entry is the observed value at the i-th unit for the n-th experiment. \cr
#'  X \tab An N by I matrix for the "binom" family, or NULL otherwise. The (n,i)-th entry is the observed size parameter for the i-th unit in the n-th experiment. \cr
#'  fac \tab Same as the input argument \code{fac}.\cr
#'  Theta \tab A K by I matrix. The (k,i)-th element is the indicator of whether the i-th unit is binding for condition k.\cr
#'  W \tab A K by (J*S) matrix. The (k,J*(s-1)+j)-th entry is the probability of units in the j-th cluster have state s under the k-th experimental condition. \cr
#'  Z \tab An I by J matrix. Each column is the indicator for an individual loci set.\cr
#'  delta \tab Same as the input argument \code{delta}.\cr
#'  zeta \tab Same as the input argumnent \code{zeta}.\cr
#'  non.id \tab A vector of length I indicating the number of each unit not to belong to any cluster.\cr
#'  prior.mean \tab A vector of length S, the hyper means for the S states for each experiment.\cr
#'  prior.sd \tab A vector of length S, the hyper parameters for the dispersions for the S states for each experiment.\cr
#'  Mu \tab A N by S matrix. The (n,s)-th entry is the mean values of the response for the s-th state in the n-th experiment.\cr
#'  stdev \tab A vector of length S. The dispersion parameter for the S states for the observed data.\cr
#' }
#' @seealso \code{\link{MBASIC.sim.state}}
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' dat.sim <- MBASIC.sim(xi = 2, I = 100, fac = rep(1:5, each = 2), J = 3)
#' @importFrom gtools rdirichlet
#' @importFrom msm rtnorm
#' @export

madbayes.sim<- function(xi, family = "lognormal", struct = NULL, I, fac, J, S = 2, f = 5, delta = NULL, delta.non = NULL, zeta = 0.1, condition = TRUE, statemap = NULL) {
  
  if(!family %in% c("lognormal", "negbin", "binom")) {
    stop("Error: 'family' must be one of 'lognormal', 'negbin' or 'binom'.")
  }

  if(is.null(statemap)) {
    statemap <- seq(S)
  }
 
  if(prod(sort(unique(statemap)) == seq(S)) != 1) {
    stop("Error: statemap must be consisted of 1, 2, ... S.")
  }

  M <- length(statemap)
  K <- length(unique(fac))
  N <- length(fac)

  para.theta <- MBASIC.sim.state(I=I, K=K, J=J, S=S, delta=delta, delta.non=delta.non, zeta=zeta, struct = struct, statemap = statemap, condition = condition)
   V <- matrix(0, nrow = N, ncol = M)
  for(s in seq(S)) {
    ids <- which(statemap == s)
    if(length(ids) == 1)
      V[, ids] <- 1
    else
      V[, ids] <- t(matrix(rdirichlet(N, delta[ids]), nrow = length(ids)))
  }
  
  Delta <- matrix(0, nrow = N, ncol = I)
  for(n in seq(N)) {
    for(s in seq(S)) {
      ids <- which(para.theta$Theta[fac[n], ] == s)
      if(length(ids) == 0)
        next
      if(sum(statemap == s) == 1) {
        Delta[n, ids] <- which(statemap == s)
      } else {
        Delta[n, ids] <- sample(which(statemap == s), length(ids), replace = TRUE, prob = V[n, statemap == s])
      }
    }
  }
  
  if(family == "lognormal") {
    prior_mean <- xi + log ((seq(M) - 1) * (f - 1) + 1)
#    prior_mean <- xi*seq(M)
    prior_sd <- 0.05
    ## sdev <- diff(exp(prior_mean))[1]
    ## stdev <- sqrt((log(sdev ^ 2 + exp(2 * prior_mean)) - 2 * prior_mean) / 2)
    stdev <- 0.5
  } else {
    prior_mean <- xi * ((seq(M) - 1) * (f - 1) + 1)
    prior_sd <- 0.5
    sdev <- prior_mean + sqrt(prior_mean)
    stdev <- prior_mean / (sdev / prior_mean - 1)
    stdev[stdev < 0] <- 5
    stdev[stdev > 5] <- 5
  }

  if(family != "binom") {  
    prior.Mu <- t(matrix(rtnorm(M * N, mean = prior_mean, sd = prior_sd, lower = 0), nrow = M))
  } else {
    prior.Mu <- t(matrix(rbeta(M * N, f * seq(M), f * rev(seq(M))), nrow = M))
  }

  ## order the means within each replicate
  for(n in seq(N)) {
    od <- order(prior.Mu[n, ])
    prior.Mu[n, ] <- prior.Mu[n, od]
  }
  
  Mu <- matrix(0, nrow = N, ncol = I)
  for(n in 1:N)
    for(m in seq(M))
      Mu[ n, Delta[n,] == m ] <- prior.Mu[ n, m ]

  X <- NULL
  if(family == "lognormal") {
    Y <- matrix(as.integer(exp(rtnorm(N * I, mean = Mu, sd = stdev, upper = max(Mu) + 1))), nrow = N)
  } else if(family == "negbin") {
    Y <- matrix(rnbinom(N * I, mu = Mu, size = stdev), nrow = N)
  } else {
    X <- matrix(rpois(N * I, lambda = xi), nrow = N)
    Y <- matrix(rbinom(N * I, size = X, prob = Mu), nrow = N)
  }
  
  bkng <- mean(Y[ Delta == 1 ])
  snr <- mean(Y[ Delta == 2 ]) / mean(Y[ Delta == 1 ])

  return(list(Theta = para.theta$Theta, Y = Y, X = X, fac = fac, W = para.theta$W, Z = para.theta$Z, V = V, delta = para.theta$delta, zeta = para.theta$zeta, prior.mean = prior_mean, prior.sd = prior_sd, stdev = stdev, Mu = prior.Mu, bkng = bkng, snr = snr, non.id = para.theta$non.id))
  
}
