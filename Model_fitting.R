DDP <- function(y, X, aff_mat, alpha=1, max_iters=10000, burn_in=1000, 
                  shape=1, rate=1, delta=c(0.1, 0.05, 0.05), log_lin=F){
    # compulsory inputs:
        # aff_mat: affinity matrix of h_ij as per (Dahl, 2008)
        # y: a matrix of dimension T*N
        # X: a matrix of dimension T*N*p
    # optional arguments:
        # max_iters: length of time for the sampler to run
        # burn_in: number of samples to discard as burn-in
        # shape, rate: parameters for the Gamma distribution
        # delta: length 3 vector for step size in random-walk MH step
    # outputs: a list containing 2 objects:
        # ch: a N*(max_iters-burn_in) matrix where columns are cluster labels
        # p: a list of vectors, each element is a 3*max(ch[, t]) matrix
    
    # setup storage and parameter initialisation
    T_max <- dim(X)[1]; N <- dim(X)[2]; num_params <- dim(X)[3]
    cluster_history <- matrix(0, N, max_iters-burn_in)
    # indexing: params_history[[iter]][params_index, cluster_index]
    params_history <- vector("list", max_iters-burn_in)
    
    cluster_labels <- sample.int(3, size=N, replace=T)
    
    K <- length(table(cluster_labels))
    cluster_params <- matrix(c(rgamma(3*K, shape=shape, rate=rate)), nrow=3, byrow=T)
    for (iter in 1:max_iters){
        # Gibbs sampling for cluster labels
        for (i in 1:N){
            cluster_sizes <- table(cluster_labels)
            c_del <- 0
            K <- length(cluster_sizes)
            old_label <- cluster_labels[i]
            
            # Computing marginal likelihoods node i came from cluster k
            # sort out underflow error
            if(K > 1){
                # compute fitted lambdas for observations 2:T
                fitted <- X[, i, ] %*% cluster_params
                if (log_lin){
                    fitted <- exp(fitted)
                }
                l_cond_lik <- apply(dpois(y[, i], fitted, log=T), 2, sum)
            } else{
                l_cond_lik <- c(sum(dpois(y[, i], fitted, log=T)))
            }
            
            l_proposal_probs <- sapply(1:K, function(k){ l_cond_lik[k] +
                    log(sum(aff_mat[i, cluster_labels == k]))})
            # prevent overflow
            max_prop_prob <- max(l_proposal_probs)
            l_proposal_probs <- l_proposal_probs - max_prop_prob
            
            # if node i is the only element in its cluster
            if (cluster_sizes[old_label] == 1){
                new_label <- sample.int(K, size=1, prob=exp(l_proposal_probs))
                # when the singleton cluster is deleted
                if (old_label != new_label){
                    # relabel so that labels are contiguous
                    cluster_labels[i] <- new_label
                    cluster_labels <- sapply(cluster_labels, function(x) x - (x > old_label))
                    cluster_params <- cluster_params[, -old_label]
                }
            } else {
                # more than one element, so propose a new cluster
                if (log_lin){
                    new_params <- rnorm(3)
                    new_l_marg <- sum(dpois(y[, i], exp((X[, i, ] %*% new_params)), log=TRUE))
                } else{
                    new_params <- rgamma(3, shape=shape, rate=rate)
                    new_l_marg <- sum(dpois(y[, i], (X[, i, ] %*% new_params), log=TRUE))
                }
                l_proposal_probs[K+1] <- log(alpha) + new_l_marg - max_prop_prob
                l_proposal_probs <- l_proposal_probs - max(l_proposal_probs)
                new_label <- sample.int(K+1, size=1, prob=exp(l_proposal_probs))
                if(new_label==(K+1)){
                    cluster_params <- cbind(cluster_params, new_params)
                }
                cluster_labels[i] <- new_label
            }
        }
        K <- length(table(cluster_labels))
        
        # fix array dim for consistency when single cluster
        if(K == 1){
            cluster_params <- array(cluster_params, dim=c(3, 1))
        }
        # Metropolis-Hastings step for cluster parameters
        for (k in 1:K){
            k_nodes <- which(cluster_labels == k)
            y_k <- y[, k_nodes]; X_k <- X[, k_nodes, ]
            cluster_params[, k] <- pr_coef_sampler(y_k, X_k, cluster_params[, k], shape=shape, 
                                                   rate=rate, delta=delta, log_lin=log_lin)
        }
        
        # storing samples
        if (iter > burn_in){
            cluster_history[, iter-burn_in] <- cluster_labels
            params_history[[iter-burn_in]] <- cluster_params
        }
    }
    return(list("ch"=cluster_history, "p"=params_history))
}


pr_coef_sampler <- function(y, X, coefs, delta=c(0.2, 0.05, 0.05), log_lin=FALSE, 
                            shape=1, rate=1){
    # compulsory inputs:
        # y is the (T-1)*N matrix of observations for t = 2, ..., T
        # X is the (T-1)*N*3 array of predictors
        # coefs is a length 3 vector of previously sampled coefficients 
    # outputs: 
        # coefficients from a random walk Metropolis-Hastings step
    
    # deal with dimensional edge cases
    if (length(dim(X)) == 2){
        # R convention drops singleton to matrix, so fix dimension consistency 
        T_max <- length(y)
        num_params <- dim(X)[2]
        y <- array(y, dim=c(T_max, 1))
        X <- array(X, dim=c(T_max, 1, num_params))
    } 
    
    # when the cluster is empty, assuming 3 parameters
    if (dim(y)[2] == 0){
        if (log_lin==TRUE){
            return(rnorm(3))
        }
        return(rgamma(3, shape=shape, rate=rate))
    }
    
    num_params <- dim(X)[3]
    fitted <- t(apply(X, 1, function(x) x %*% coefs))
    if (log_lin==TRUE){
        fitted <- exp(fitted)
    }
    old_llik <- sum(dpois(y, fitted, log=TRUE))
    
    for (p in 1:num_params){
        new_coefs <- coefs 
        new_coefs[p] <- runif(1, min=coefs[p]-delta[p], max=coefs[p]+delta[p])
        if (new_coefs[p] > 0 | log_lin==TRUE){
            new_fitted <- t(apply(X, 1, function(x) x %*% new_coefs))
            if(log_lin==TRUE){
                new_fitted <- exp(new_fitted)
                l_prior_ratio <- dnorm(new_coefs[p], log=TRUE) - dnorm(coefs[p], log=TRUE)
            } else{
                l_prior_ratio <- dgamma(new_coefs[p], shape=shape, rate=rate, log=TRUE) - 
                    dgamma(coefs[p], shape=shape, rate=rate, log=TRUE)
            }
            new_llik <- sum(dpois(y, new_fitted, log=TRUE))
            l_lik_ratio = new_llik - old_llik
            
            l_acceptance_prob <- min(0, l_prior_ratio + l_lik_ratio)
            u <- runif(1)
            if (log(u) < l_acceptance_prob){
                old_llik <- new_llik
                coefs <- new_coefs 
            }
        } 
    }
    return(coefs)
}


FMM <- function(y, X, alpha=c(1, 1, 1, 1), max_iters=10000, burn_in=1000,
                log_lin=FALSE, shape=1, rate=1, delta=c(0.2, 0.1, 0.1)){
    
    # compulsory inputs:
        # y: a matrix of dimension T*N
        # X: a tensor of dimension T*N*num_params
        # alpha: Dirichlet distribution concentration parameter
    # outputs: a list containing 2 objects:
        # ch: a N*(max_iters-burn_in) matrix where columns are cluster label samples
        # p: a (max_iters-burn_in)*num_params*length(alpha) array of cluster coefficients
        # prob: a length(alpha)*(max_iters-burn_in) matrix of mixture proportion samples
    
    # setup storage history and param initialisation
    T_max <- dim(X)[1]; N <- dim(X)[2]; num_params <- dim(X)[3]
    if (alpha == c(1)){
        params_history <- matrix(0, max_iters - burn_in, 3)
        params <- rgamma(3, shape=shape, rate=rate)
        for (iter in 1:max_iters){
            params <- pr_coef_sampler(y, X, params, log_lin=log_lin, 
                                      shape=shape, rate=rate, delta=delta)
            if (iter>burn_in){
                params_history[iter-burn_in, ] <- params
            }
        }
        return(params_history)
        
    } else{
        K <- length(alpha)
        
        cluster_history <- matrix(0, N, max_iters-burn_in)
        params_history <- array(0, dim=c(max_iters-burn_in, num_params, K))
        probs_history <- matrix(0, K, max_iters-burn_in)
        cluster_labels <- c(1:K, sample.int(K, size=N-K, replace=T))
        cluster_params <- matrix(c(rgamma(num_params*K, shape=shape, rate=rate)), 
                                 nrow=3, byrow=T)
        
        
        for (iter in 1:max_iters){
            # Gibbs sampling for cluster allocation probs
            label_counts <- sapply(1:K, function(k) length(which(cluster_labels==k)))
            probs <- rdirichlet(n=1, alpha + label_counts)
            # Gibbs sampling for cluster labels
            for (i in 1:N){
                fitted <- X[, i, ] %*% cluster_params
                if (log_lin){
                    fitted <- exp(fitted)
                }
                l_cond_lik <- apply(dpois(y[, i], fitted, log=T), 2, sum)
                l_proposal_probs <- sapply(1:K, function(k) {l_cond_lik[k] + log(probs[k])})
                max_prop_prob <- max(l_proposal_probs)
                l_proposal_probs <- l_proposal_probs - max_prop_prob
                cluster_labels[i] <- sample.int(K, size=1, prob=exp(l_proposal_probs))
            }
            
            # Metropolis-Hastings step for cluster parameters
            for (k in 1:K){
                k_nodes <- which(cluster_labels == k)
                y_k <- y[, k_nodes]; X_k <- X[, k_nodes, ]
                cluster_params[, k] <- pr_coef_sampler(y_k, X_k, cluster_params[, k],
                                                       log_lin=log_lin, shape=shape, 
                                                       rate=rate, delta=delta)
            }
        
            # storing samples
            if (iter > burn_in){
                cluster_history[, iter-burn_in] <- cluster_labels
                params_history[iter-burn_in, ,] <- cluster_params
                probs_history[, iter-burn_in] <- probs
            }
        }
    }
    return(list("ch"=cluster_history, "p"=params_history, "prob"=probs_history))
}
                                   

PNAR_X <- function(y, X, log_linear=FALSE){
    # Inputs:
        # y: a T_max * N matrix, where columns correspond to time series
        # X: a T_max * N * num_params tensor where X[t, i, ] contains the predictors
        #    for Y[t, i]
        # log_linear: set to TRUE for the canonical log-linear link
    # Outputs: 
        # coefficients of the fitted model
    
    T_max <- dim(y)[1]; N <- dim(y)[2]; num_params <- dim(X)[3]
    pnar_mat <- matrix(0, nrow=T_max*N, ncol=num_params+1)
    index <- 1
    for (t in 1:T_max){
        for (i in 1:N){
            pnar_mat[index, ] <- c(y[t, i], X[t, i, ])
            index <- index + 1
        }
    }
    colnames(pnar_mat) <- c("y", "i", "x", "y_prev")
    if (log_linear){
        pnar_mat <- as.data.frame(pnar_mat)
        start_coefs <- lm(y ~ 0 +., data=pnar_mat)$coefficients
        model <- glm(y ~ 0 +., family=poisson,
                     data=as.data.frame(pnar_mat), start=start_coefs)
    } else{
        obj <- function(b) -sum(dpois(pnar_mat[, 1], pnar_mat[, 2:(1+num_params)] %*% exp(b), 
                                      log = TRUE))
        coefs <- exp(optim(c(0, 0, 0), obj)$par)
    }
    return(coefs)
}


cond_clusters <- function(y, X, cluster_labels, burn_in=1000, max_iters=10000){
    K <- length(unique(cluster_labels))
    param_hist <- array(0, dim=c(K, 3, max_iters))
    mean_params <- array(0, dim=c(3, K))
    for (k in unique(cluster_labels)){
        k_nodes <- which(cluster_labels == k)
        y_k <- y[, k_nodes]; X_k <- X[, k_nodes, ]
        param_hist[k, , 1] <- c(0.5, 0.1, 0.1)
        for (iter in 1:(max_iters-1)){
            param_hist[k, , iter+1] <- pr_coef_sampler(y_k, X_k, param_hist[k, , iter])
        }
        mean_params[, k] <- apply(param_hist[k, , (burn_in+1):max_iters], 1, mean)
    }
    return(list("p"=param_hist[, ,(burn_in+1):max_iters], "mean_p"=mean_params))
}
