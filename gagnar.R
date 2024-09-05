# this code is almost entirely from the supplementary material of:
#  Ren, Y., Zhu, X., Lu, X., and Hu, G. (2024). Graphical assistant grouped network 
# autoregression model: A Bayesian nonparametric recourse.
# Journal of Business and Economic Statistics, 42:49–63.

# # packages needed: "mvtnorm","MASS"
# ipak <- function(pkg){
#   new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#   if (length(new.pkg))
#     install.packages(new.pkg, dependencies = TRUE)
#   try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
# }
# packages <- c("MCMCpack", "distr")
# ipak(packages)

###########################################################################
## GAGNAR MCMC Collapsed Gibbs sampler  
###########################################################################

GAGNAR <- function(Y_tilde, X_tilde, dmat, alpha=1, h=1, 
                  a0=0.01, b0=0.01, tau0=c(0, 0, 0), sigma0=100*diag(3), initNClusters=3, 
                  verbose=TRUE, burn_in=2000, niterations=30000){
  ## Model: z|(alpha, h) \sim gaCRP(alpha, h)
  ##        sigma2_k \sim IG(a0, b0)
  ##        theta_k|sigma2_k \sim MVN(tau0, sigma2_k*sigma0)
  ##        \tilde{Y_i}|z,theta,sigma2 \sim MVN(\tilde{X_i}*theta_{z_i}, sigma2_{z_i}*I)
  #===========================================================#
  ## Input: 
  ##        Y_tilde, T*N matrix ##
  ##        X_tilde, T*N*(p+3) array ##
  ##        dmat, N*N graph distance matrix ##
  ##        alpha = the parameter in Dirichlet distribution that 
  ##                controls the relative size of clusters ##
  ##        h = decay coefficient for graph distance ##
  ##        niterations = the total number of iterations ##
  ##        a0 = hyperparameters for the prior on sigma2 in inverse gamma ##
  ##        b0 = hyperparameters for the prior on sigma2 in inverse gamma ##
  ##        tau0 = hyperparameters for the prior on 
  ##               theta|sigma2 in multivariate normal ##
  ##        sigma0 = hyperparameters for the prior on 
  ##                 theta|sigma2 in multivariate normal ##
  ##        initNClusters = the initial number of clusters ##
  
  ## Output: 
  ##         zout = clustering configuration, n by 1 vector ##
  ##         sigma2out = estimate sigma2, K*1 matrix ##
  ##         thetaout = estimate theta, K*d matrix ##
  #===========================================================#

  t <- dim(Y_tilde)[1]
  N <- dim(Y_tilde)[2]
  d <- dim(X_tilde)[3]  # d = p + 3
  
  dmat[dmat == Inf] <- 100
  W <- exp(-dmat*h)
  diag(W) = 0
  W[dmat==1] <- 1
  W <- W * N * (N-1) / sum(W)  # weight matrix
  
  sigma0_inv <- ginv(sigma0)
  E <- diag(1, nrow = t, ncol = t)  # t by t identity matrix 
  cluster_history <- matrix(0, N, niterations-burn_in)
  params_history <- vector("list", niterations-burn_in)
  sigma_history <- vector("list", niterations-burn_in)
  
  #======== first initialize clustering and beta ========#
  # assign initial cluster for N observations
  clusterAssign <- c(sample(1:initNClusters, size=initNClusters, replace=FALSE),
                     sample(1:initNClusters, size=N-initNClusters, replace=TRUE)) # dim = n*1
  # generate theta and sigma from base distribution
  sigma2 <- rinvgamma(initNClusters, shape = a0, scale = b0)  # dim(sigma2) = [initNClusters, 1]
  theta <- tryCatch({
    sapply(sigma2, function(x) {rmvnorm(1, mean = tau0, sigma = (x * sigma0), checkSymmetry = F)}) %>% t()
  }, error = function(err) {
    print(paste0("ERROR:  ",err))
    set.seed(1111)
    sigma2 <- rinvgamma(initNClusters, shape = a0, scale = b0)  # dim(sigma2) = [initNClusters, 1]
    sapply(sigma2, function(x) {rmvnorm(1, mean = tau0, sigma = (x * sigma0), checkSymmetry = F)}) %>% t() # dim(theta) = [initNClusters, (p+3)]
  })

  #======== start Gibb's sampling ========#
  for (iter in 1:niterations)
  {
    if(verbose) {
      cat(iter, "\r")
    }

    ##======== update z ========##
    clusterSizes = table(as.factor(clusterAssign))  # named vector: {name} cluster; {value} cluster size
    nClusters = length(clusterSizes)  # updated number of clusters
    # cat("number of clusters at the begining of the this iteration: ", nClusters, "\r")
    
    for (i in 1:N)
    { # determine whether ith component is a singleton
      clustersize_of_iclu <- clusterSizes[clusterAssign[i]]
      # cat(i, "th sample's cluster is not a singleton: ", clustersize_of_iclu > 1, "\n")
      if (clustersize_of_iclu > 1){
        # if not a singleton, we have (nClusters + 1) choices
        clusterSizes[clusterAssign[i]] = clusterSizes[clusterAssign[i]] - 1
        # the probs for choosing exsiting clusters
        
        clusterProbs_log = sapply(1:nClusters, function(x){
          sum(W[i, clusterAssign == clusterAssign[x]]) *
            prod(sapply(2:t, function(y){
              exp(-(Y_tilde[y, i] - X_tilde[y,i,] %*% as.matrix(theta[x, ]))^2/(2*sigma2[x])) /
                ((2*pi*sigma2[x]) ^ (1/2))}))
        })

        # the prob for choosing a new cluster
        ## sigma*
        sigma_star <- ginv(sigma0_inv + t(X_tilde[2:t,i,]) %*% X_tilde[2:t,i,])
        ## tau*
        tau_star <- ginv(sigma0_inv + t(X_tilde[2:t,i,]) %*% X_tilde[2:t,i,]) %*%
          (sigma0_inv %*% tau0 + t(X_tilde[2:t,i,]) %*% Y_tilde[2:t, i])
        ## (sigma*)^-1
        sigmastar_inv <- sigma0_inv + t(X_tilde[2:t,i,]) %*% X_tilde[2:t,i,]
        # full conditional distribution (density) for a new cluster
        clusterProbs_log[nClusters+1] <- alpha * 
          (b0^a0 * gamma(a0 + (t-1)/2) * sqrt(det(sigma_star)) * 
          (b0 + (t(tau0) %*% sigma0_inv %*% tau0 + 
                   t(Y_tilde[2:t, i]) %*% Y_tilde[2:t, i] - 
                   t(tau_star) %*% sigmastar_inv %*% tau_star)/2)^(-(a0 + (t-1)/2)) /
          ((2*pi)^((t-1)/2) * gamma(a0) * sqrt(det(sigma0))))
        clusterProbs <- clusterProbs_log

        # choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1, prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i > nClusters) {
          # if the update generate a new cluster label
          new_sigma2 <- rinvgamma(1, shape = a0, scale = b0)
          sigma2 <- c(sigma2, new_sigma2)  # [nCluster + 1, 1]
          new_theta <- tryCatch({
            rmvnorm(1, mean = tau0, sigma = (new_sigma2 * sigma0))
            }, error = function(err) {
            print(paste0("ERROR:  ",err))
            # generate a new theta from MVN(tau0, sigma0)
            rmvnorm(1, mean = tau0, sigma = 100000 * sigma0)
          })
          theta <- rbind(theta, new_theta)  # [nCluster + 1, p + 3]
        }
        
      } 
      if (clustersize_of_iclu <= 1) {
        # if ith component is singleton, we have nClusters choices
        clusterAssign[clusterAssign > clusterAssign[i]] <- 
          clusterAssign[clusterAssign > clusterAssign[i]] - 1

        if (dim(theta)[1] > 2){  # nCluster>2
          sigma2 <- sigma2[-clusterAssign[i]]
          theta <- theta[-clusterAssign[i], ]
        } else {  # process as a matrix with only 1 row
          sigma2 <- sigma2[-clusterAssign[i]]
          theta <- matrix(theta[-clusterAssign[i], ], nrow = nClusters-1)
        }
        
        # the probs for choosing existing clusters
        clusterProbs_log = sapply(1:(nClusters-1), function(x){
          sum(W[i, clusterAssign == clusterAssign[x]]) *
            prod(sapply(2:t, function(y){
              exp(- (Y_tilde[y, i] - X_tilde[y,i,] %*% as.matrix(theta[x, ]))^2/(2*sigma2[x])) /
                ((2*pi*sigma2[x])^(1/2))}))
        })
        
        # the prob for choosing a new cluster
        ## sigma*
        sigma_star <- ginv(sigma0_inv + t(X_tilde[2:t,i,]) %*% X_tilde[2:t,i,])
        ## tau*
        tau_star <- sigma_star %*%
          (sigma0_inv %*% tau0 + t(X_tilde[2:t,i,]) %*% Y_tilde[2:t, i])
        ## (sigma*)^-1
        sigmastar_inv <- sigma0_inv + t(X_tilde[2:t,i,]) %*% X_tilde[2:t,i,]
        
        clusterProbs_log[nClusters] <- alpha * 
          (b0^a0 * gamma(a0 + (t-1)/2) * sqrt(det(sigma_star)) * 
             (b0 + (t(tau0) %*% sigma0_inv %*% tau0 + 
                      t(Y_tilde[2:t, i]) %*% Y_tilde[2:t, i] - 
                      t(tau_star) %*% sigmastar_inv %*% tau_star)/2)^(-(a0 + (t-1)/2)) /
             ((2*pi)^((t-1)/2) * gamma(a0) * sqrt(det(sigma0))))
        
        clusterProbs <- clusterProbs_log

                # choose the cluster number for ith observation
        cluster.i <- sample(1:nClusters, size = 1, prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i == nClusters) {
          new_sigma2 <- rinvgamma(1, shape = a0, scale = b0)
          while (new_sigma2 == Inf) {
            new_sigma2 <- rinvgamma(1, shape = a0, scale = b0)
          }
          sigma2 <- c(sigma2, new_sigma2)  # [nCluster + 1, 1]
          new_theta <- tryCatch({
            rmvnorm(1, mean = tau0, sigma = (new_sigma2 * sigma0))
          }, error = function(err) {
            print(paste0("ERROR:  ",err))
            # generate a new theta from MVN(tau0, sigma0)
            rmvnorm(1, mean = tau0, sigma = 100000 * sigma0)
          })
          theta <- rbind(theta, new_theta)  # [nCluster + 1, p + 3]
        }
      }
      
      # update cluster sizes corresponding to each cluster name
      clusterSizes <- table(as.factor(clusterAssign))
      nClusters <- length(clusterSizes)
    }

    ## update theta and sigma2 ##
    for (r in 1:nClusters){
      ## sigma*
      X_group_r <- X_tilde[ 2:t, which(clusterAssign == r), ]

      if (sum(clusterAssign == r) == 1) {  # if there is single node, then X_group_r appears to be a matrix
        X_tilde_t_sum <- t(X_group_r) %*% X_group_r
      } else if (sum(clusterAssign == r) > 1) { 
        X_list_bygroup <- lapply(1:dim(X_group_r)[2], function(rr) {X_group_r[,rr,]})
        X_tilde_t <- lapply(X_list_bygroup, function(x){t(x) %*% x})
        X_tilde_t_sum <- Reduce('+', X_tilde_t)
      }

      post_sigma_star <- ginv(sigma0_inv + X_tilde_t_sum)
      
      ## tau*
      Y_group_r <- Y_tilde[ 2:t, clusterAssign == r]  # dim(Y_group_r) = [T, Nr]
      
      if (sum(clusterAssign == r) == 1) {  # if there is single node, then Y_group_r appears to be a vector
        X_tilde_t_Y_sum <- t(X_group_r) %*% Y_group_r  # dim = [d, 1]
      } else if (sum(clusterAssign == r) > 1) {
        X_tilde_t_Y <- lapply(1:dim(Y_group_r)[2], function(x) {t(X_list_bygroup[[x]]) %*% Y_group_r[,x]})
        X_tilde_t_Y_sum <- Reduce('+', X_tilde_t_Y)
      }
      
      post_tau_star <- post_sigma_star %*% (sigma0_inv %*% tau0 + X_tilde_t_Y_sum)
      
      ## a*
      post_a_star <- a0 + (t-1) * sum(clusterAssign == r)/2
      
      ## b*
      Y_tilde_sum <- sum(apply(as.matrix(Y_tilde[2:t, clusterAssign == r]), 2, function(x){t(x) %*% x}))
      post_b_star <- (b0 + (t(tau0) %*% (sigma0_inv) %*% tau0 + Y_tilde_sum - t(post_tau_star) %*% (sigma0_inv + X_tilde_t_sum) %*% post_tau_star)/2) %>% as.numeric()
      
      # update sigma2 given Y_tilde. Note that sigma2|Y_tilde \sim IG(a*, b*)
      sigma2[r] <- rinvgamma(1, shape = post_a_star, scale = post_b_star)
      # update theta given (sigma2, Y_tilde). Note that theta|sigma2,Y_tilde \sim N(tau*, sigma2 sigma*)
      theta[r,] <- rmvnorm(1, mean = post_tau_star, sigma = sigma2[r] * post_sigma_star, checkSymmetry = F)

    }
    if (iter > burn_in){
      cluster_history[, iter-burn_in] <- clusterAssign
      params_history[[iter-burn_in]] <- t(theta)
      sigma_history[[iter-burn_in]] <- sigma2
    }
    # History[[iter]] <- list(zout = clusterAssign, sigma2out = sigma2, thetaout = theta)
  }
  return(list("ch"=cluster_history, 
              "p"=params_history,
              "sigma"=sigma_history))
}

