co_occurence <- function(clusts, thin=10){
    cc <- matrix(0, nrow(clusts), nrow(clusts))
    for (i in 1:nrow(clusts)) {
        for (j in 1:nrow(clusts)) {
            for (m in 1:(ncol(clusts) %/% thin)) {
                cc[i,j]=cc[i,j]+(clusts[i, m*thin]==clusts[j, m*thin])
            }
        }
    }
    return(cc/(ncol(clusts) %/% thin))
}

cluster_MAP <- function(cluster_mcmc, thin=1, crp=FALSE, cc=NULL){
    N <- nrow(cluster_mcmc)
    optimal_partition <- cluster_mcmc[, 1]
    opt <- 1
    same_optimal <- c()
    if (is.null(cc)){
        cc <- co_occurence(cluster_mcmc)
    }
    min_loss = N^2/2
    for (its in 1:(ncol(cluster_mcmc) %/% thin)){
        loss <- 0
        iter <- its*thin
        for(i in 1:(N-1)){
            for(j in (i+1):N){
                loss = loss + ((cluster_mcmc[i, iter] == cluster_mcmc[j, iter]) - cc[i, j])^2
            }
        }
        if (loss == min_loss){
            if (crp){
                if (cluster_mcmc[, same_optimal[1]] == cluster_mcmc[, iter]){
                    same_optimal <- c(same_optimal, iter)
                }
            } else{
                same_optimal <- c(same_optimal, iter)
            }
        }
        if (loss < min_loss){
            min_loss <- loss; optimal_partition <- cluster_mcmc[, iter]
            same_optimal <- c(iter)
        }
    }
    return(optimal_partition)
}


one_step_forecast <- function(X_pred, mcmc_out, crp=FALSE,
                              pop_adjust=NULL, last_step=NULL, thin=10){
    # Inputs:
    # X_pred: a N * 3 matrix of the predictor variables
    # pop_adjust: length N vector of population adjustments
    # 
    N <- dim(mcmc_out$ch)[1]
    y_hat <- rep(0, N)
    HDP_mean <- matrix(0, nrow=N, ncol=2)
    ESS <- rep(0, N)
    for (i in 1:N){
        sampled_params <- single_node_samples(mcmc_out, node=i, crp=crp, thin=thin)
        pred <- t(sampled_params) %*% X_pred[i, ]
        if (!is.null(last_step)){
            # for reversing lag1-differencing and population adjustment
            pred <- last_step[i] + pop_adjust[i]/100000 * pred
        }
        ESS[i] <- coda::effectiveSize(pred)
        HDP_mean[i, ] <- coda::HPDinterval(as.mcmc(pred))
        y_hat[i] <- mean(pred)
    }
    HDP <- mapply(function(x, y) paste0("[", round(x, digits=2), ", ", 
                                        round(y, digits=2), "]"), 
                  HDP_mean[, 1], HDP_mean[, 2], SIMPLIFY = TRUE)
    return(list("y_hat"=y_hat, "HPD" = HDP, "ESS"=ESS))
}


scaled_error <- function(y, y_hat){
    T_max <- dim(y)[1]
    se <- rep(0, length(y_hat))
    for (i in 1:length(y_hat)){
        se[i] <- abs(y[T_max, i] - y_hat[i]) / 
            mean(abs(y[2:T_max, i] - y[1:(T_max-1), i]))
    }
    return(se)
}


model_fit <- function(true_vals, mcmc_out=NULL, J=10, monte_carlo=NULL, crp=FALSE,
                      X_pred=NULL, poisson_means=NULL, thin=1, gauss=FALSE, 
                      ldensity=FALSE){
    # Inputs:
        # true_vals: a T*N matrix of true values, unless monte_carlo  
        # choice of 3 ways of representing the 'null distribution':
            #1 monte_carlo: M x B vector of Monte Carlo distribution draws
            #2 mcmc_out: MCMC outputs from DDP() or GAGNAR()
            #  X_pred: tensor where [t, i, ] indexes relevant covariates 
            #3 poisson_means: a T*N matrix of Poisson means
    score <- 0
    nonrand_PIT <- rep(0, J)
    if (!is.null(monte_carlo)){
        # assumes distribution is approximated using Monte Carlo simulations
        N <- length(true_vals); B <- dim(monte_carlo)[2]
        emp_dists <- apply(monte_carlo, 1, table)
        for (i in 1:N){
            prob_mass <- emp_dists[[i]][as.character(true_vals[i])]
            if (is.na(prob_mass)){
                # approximation if observed value not simulated
                prob_mass <- dpois(true_vals[i], mean(monte_carlo[i, ]))
            }
            score <- score - log(prob_mass/ B)
        }
        return(unname(score))
    } else {
        # increase dimension by 1 with T_max=1
        if (is.null(dim(true_vals))){
            true_vals <- matrix(true_vals, nrow=1)
            if (!is.null(X_pred)){
                X_pred <- array(X_pred, dim=c(1, dim(X_pred)))
            } else{
                poisson_means <- matrix(poisson_means, nrow=1)
            }
        }
        T_max <- dim(true_vals)[1]; N <- dim(true_vals)[2]
        rPIT <- matrix(0, nrow=T_max, ncol=N)
        cdf_lower <- matrix(0, nrow=T_max, ncol=N)
        cdf_upper <- matrix(0, nrow=T_max, ncol=N)
        if (!is.null(mcmc_out)){
            B <- dim(mcmc_out$ch)[2] %/% thin
            for (i in 1:N){
                sampled_params <- single_node_samples(mcmc_out, node=i, crp=crp, thin=thin)
                sampled_params <- t(sampled_params)
                if (gauss){
                    # get the variances of the Gaussian mixture components
                    sampled_sigmas <- rep(0, B)
                    for (iter in thin*(1:B)){
                        cluster_index <- mcmc_out$ch[i, iter]
                        sampled_sigmas[iter/thin] <- mcmc_out$sigma[[iter]][cluster_index]
                    }
                    for (t in 1:T_max){
                        fitted_vals <- c(sampled_params %*% X_pred[t, i, ])
                        cdf_lower[t, i] <- sum(pnorm(true_vals[t, i]-0.5, mean=fitted_vals,
                                                     sd=sampled_sigmas)) / B
                        cdf_upper[t, i] <- sum(pnorm(true_vals[t, i]+0.5, mean=fitted_vals,
                                                     sampled_sigmas)) / B
                        rPIT[t, i] <- sum(pnorm(true_vals[t, i], mean=fitted_vals,
                                                sampled_sigmas)) / B
                    }
                } else{
                    for (t in 1:T_max){
                        fitted_vals <- sampled_params %*% X_pred[t, i, ]
                        cdf_lower[t, i] <- sum(sapply(fitted_vals, 
                                                      function(x) ppois(true_vals[t, i]-1, x))) / B
                        cdf_upper[t, i] <- sum(sapply(fitted_vals, 
                                                      function(x) ppois(true_vals[t, i], x))) / B
                    }
                }
            }
            if (ldensity==TRUE){
                pmf <- cdf_upper - cdf_lower
                pmf[which(pmf<10^(-12))] <- 10^(-12)
                print(pmf)
                return(sum(log(pmf)))
            }
        } else{  # assume poisson distributed with parameter poisson_means
            for (i in 1:N){
                for (t in 1:T_max){
                    cdf_lower[t, i] <- ppois(true_vals[t, i] - 1, poisson_means[t, i])
                    cdf_upper[t, i] <- ppois(true_vals[t, i], poisson_means[t, i])
                }
            }
        }
        if (gauss){
            score <- -sum(log(cdf_upper - cdf_lower + 10^(-12)))
            return(list("score"=score/(N*T_max), "PIT"=rPIT))
        }
        for (i in 1:N){
            for (t in 1:T_max){
                # prevent underflow error
                score <- score - log(cdf_upper[t, i] - cdf_lower[t, i] + 10^(-12))
                # Brockwell's (2007) randomised prob. integral transform
                if (gauss==FALSE){
                    rPIT[t, i] <- cdf_lower[t, i] + 
                        runif(1)*(cdf_upper[t, i] - cdf_lower[t, i])
                    # Czado et al. (2009) non-randomised PIT
                    for (j in 1:J){
                        if (cdf_lower[t, i] <= j/J){
                            if (cdf_upper[t, i] >= j/J){
                                nonrand_PIT[j] <- nonrand_PIT[j] + 
                                    (j/J-cdf_lower[t, i])/(cdf_upper[t, i] - cdf_lower[t, i] +10^(-12))
                            } else{
                                nonrand_PIT[j] <- nonrand_PIT[j] + 1
                            }
                        }
                    }
                }
            }
        }
    }
    PIT_hist <- nonrand_PIT - c(0, nonrand_PIT[1:(J-1)])
    return(list("score"=score/(N*T_max), "nrPIT"=PIT_hist/(N*T_max), "PIT"=rPIT))
}


mcmc_diagnostics <- function(mcmc_out){
    N <- dim(mcmc_out$ch)[1]
    geweke <- matrix(0, N, 3); ESS <- matrix(0, N, 3)
    for (i in 1:N){
        samples <- as.mcmc(t(single_node_samples(mcmc_out, crp=T, node=i)))
        geweke[i, ] <- geweke.diag(samples)$z
        ESS[i, ] <- effectiveSize(samples)
    }
    return(list("geweke"=geweke, "ESS"=ESS))
}


cluster_number <- function(title, mcmc_out){
    clusts <- apply(mcmc_out$ch, 2, max)
    png(paste("clusts_", title, ".png", sep=""))
    hist(clusts, xlab="Number of clusters",
         ylab="Posterior distribution")
    dev.off()
    return(list("prop"<-table(clusts), 
                "geweke"<- geweke.diag(clusts)$z,
                "ESS" <- effectiveSize(clusts)))
}

