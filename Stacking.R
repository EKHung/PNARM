stacking_weights <- function(y, X_pred, mcmc_list, lambda=1.001){
    C <- length(mcmc_list); N <- dim(mcmc_list[[1]]$ch)[1]
    pred_dens <- matrix(0, dim(y)[1], C) 
    S_eff <- rep(0, C)
    for (c in 1:C){
        post_preds <- posterior_predictive(y, X_pred, mcmc_list[[c]], thin=10)$pmf
        pred_dens[, c] <- rowSums(log(post_preds))
        ess <- matrix(0, N, dim(mcmc_list[[1]]$p)[2])
        for (i in 1:N){
            samples <- single_node_samples(mcmc_list[[c]], node=i, thin=10)
            ess[i, ] <- coda::effectiveSize(as.mcmc(t(samples)))
        }
        S_eff[c] <- median(ess)
    }
    pred_dens <- exp(pred_dens - apply(pred_dens, 1, max))
    print(pred_dens)
    print(S_eff)
    alpha <- lambda * S_eff / sum(S_eff)
    
    obj <- function(w){
        w_all <- c(1-sum(w), w)
        return(-sum(log(pred_dens %*% w_all)) -log(ddirichlet(w_all, alpha)))
    }
    gradient <- function(w){
        w_all <- c(1-sum(w), w)
        grad <- rep(0, C-1)
        for (c in 2:C){
            for (t in 1:dim(y)[1]){
                grad[c-1] <- grad[c-1] + (pred_dens[t, c] - pred_dens[t, 1]) / 
                    sum(w_all * pred_dens[t, ])
            }
            grad[c-1] <- grad[c-1] + (alpha[c]-1)/w_all[c]
        }
        return(-grad)
    }
    print(obj(rep(1/C, C-1)))
    ui <- rbind(rep(-1, C - 1), diag(C - 1))  # K-1 simplex constraint matrix
    ci <- c(-1, rep(0, C - 1))
    w <- constrOptim(rep(1/C, C-1), obj, gradient, ui, ci)
    opt_weights <- c(1-sum(w$par), w$par)
    return(opt_weights)
}


FMM_stack <- function(y, X, alpha=c(1, 1, 1, 1), max_iters=6000, burn_in=1000,
                      log_lin=FALSE, shape=1, rate=1, delta=c(0.2, 0.1, 0.1),
                      C=10, test_set=0.2, lambda=1.001){
    # inputs:
    # y: a matrix of dimension T*N
    # X: a tensor of dimension T*N*p
    # alpha: Dirichlet distribution shape parameter
    # outputs:
    #
    train_points <- dim(y)[1] - floor(dim(y)[1]*test_set)
    y_train <- y[1:train_points, ]; X_train <- X[1:train_points, , ]
    y_test <- y[(train_points+1): dim(y)[1], ]
    X_test <- X[(train_points+1): dim(y)[1], , ]
    
    FMM_list <- list()
    for (c in 1:C){
        FMM_list [[c]] <- FMM(y_train, X_train, alpha=alpha, max_iters=max_iters, 
                              burn_in=burn_in, log_lin=log_lin, shape=shape,
                              rate=rate, delta=delta)
    }
    weights <- stacking_weights(y_test, X_test, FMM_list, lambda=lambda)
    return(list("chains"=FMM_list, "weights"=weights))
}


stack_posterior_pmf <- function(y, X_pred, mcmc_list, weights, thin=1){
    C <- length(mcmc_list)
    post_pred <- posterior_predictive(y, X_pred, mcmc_list[[1]], thin=thin)
    lower <- post_pred$lower * weights[1]
    pmf <- post_pred$pmf * weights[1]
    for (c in 2:C){
        post_pred <- posterior_predictive(y, X_pred, mcmc_list[[c]], thin=thin)
        lower <- lower + post_pred$lower * weights[c]
        pmf <- pmf + post_pred$pmf * weights[c]
    }
    return(list("lower"=lower, "pmf"=pmf))
}


stack_score <- function(y, X_pred, stack, thin=1){
    post_pred <- stack_posterior_pmf(y, X_pred, stack$chains, stack$weights, thin=thin)
    score <- sum(-log(post_pred$pmf))
    PIT <- post_pred$lower + runif(length(post_pred))*post_pred$pmf
    return(list("score"=score/length(y), "PIT"=PIT))
}


stack_mean_coefs <- function(mcmc_list, weights, thin=1){
    N <- dim(mcmc_list[[1]]$ch)[1]; C <- length(mcmc_list)
    coefs <- matrix(0, N, 3)
    for (c in 1:C){
        for (i in 1:N){
            samples <- single_node_samples(mcmc_list[[c]], node=i, thin=thin)
            samples <- apply(samples, 1, mean)
            coefs[i, ] <- coefs[i, ] + samples*weights[c]
        }
    }
    return(coefs)
}


stack_co_occurence <- function(mcmc_list, weights, thin=10){
    C <- length(mcmc_list)
    cc <- co_occurence(mcmc_list[[1]]$ch, thin=thin) * weights[1]
    for (c in 2:C){
        cc <- cc + co_occurence(mcmc_list[[c]]$ch, thin=thin)*weights[c]
    }
    return(cc)
}
