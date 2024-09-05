affinity <- function(dist_mat, fn = function(x) x, standardise=TRUE){
    hmat <- fn(dist_mat)
    diag(hmat) <- 0
    # normalise rows to sum to N-1
    if (standardise){
        hmat <- hmat / rowSums(hmat)
        hmat <- hmat * (nrow(dist_mat)-1)
    }
    return (hmat)
}


pois_mix_sim <- function(X_pred, mcmc_out, nsims=10000, crp=FALSE){
    # simulates per iteration in case of dependencies induced by partition
    N <- dim(mcmc_out$ch)[1]; B <- dim(mcmc_out$ch)[2]
    fitted_vals <- matrix(0, N, B)
    for (i in 1:N){
        sampled_params <- t(single_node_samples(mcmc_out, node=i, crp=crp))
        fitted_vals[i, ] <- sampled_params %*% X_pred[i, ]
    } 
    rand_means <- fitted_vals[, sample.int(B, nsims, replace=TRUE)]
    return(matrix(rpois(N*nsims, rand_means), nrow=N))
}


posterior_predictive <- function(true_vals, X_pred, mcmc_out, crp=FALSE, thin=1){
    if (is.null(dim(true_vals))){ # dimension consistency
        true_vals <- matrix(true_vals, nrow=1)
        X_pred <- array(X_pred, dim=c(1, dim(X_pred)))
    } 
    T_max <- dim(true_vals)[1]; N <- dim(true_vals)[2]
    cdf_lower <- matrix(0, nrow=T_max, ncol=N)
    cdf_upper <- matrix(0, nrow=T_max, ncol=N)
    B <- dim(mcmc_out$ch)[2] %/% thin
    for (i in 1:N){
        sampled_params <- single_node_samples(mcmc_out, node=i, crp=crp, thin=thin)
        sampled_params <- t(sampled_params)
        for (t in 1:T_max){
            fitted_vals <- sampled_params %*% X_pred[t, i, ]
            cdf_lower[t, i] <- sum(sapply(fitted_vals, 
                                          function(x) ppois(true_vals[t, i]-1, x))) / B
            cdf_upper[t, i] <- sum(sapply(fitted_vals, 
                                          function(x) ppois(true_vals[t, i], x))) / B
        }
    }
    pmf <- cdf_upper - cdf_lower
    pmf[which(pmf<10^(-12))] <- 10^(-12)
    return(list("lower"=cdf_lower, "pmf"=pmf))
}


X_from_ts <- function(ts, adj_mat, lags=1, pop_adjust=NULL,
                      log_lin=FALSE, intercept_scale=100000,
                      pop_adjust_intercept=FALSE){
    # inputs: 
    # ts: T*N time series matrix 
    # adj_mat: (unweighted) N*N adjacency matrix 
    # pop_adjust: (optional) a vector of weight adjustments
    # ouputs:
    # X: design tensor 
    max_T <- nrow(ts)
    N <- ncol(ts)
    W <- adj_mat / colSums(adj_mat)
    intercept <- matrix(1, max_T-lags, N)
    if (!is.null(pop_adjust)){
        W <- diag(1/pop_adjust) %*% W %*% diag(pop_adjust)
        if (pop_adjust_intercept){
            intercept <- intercept %*% diag(pop_adjust/intercept_scale)
        }
    }
    lagged <- ts[1:(max_T-lags), ]
    if (log_lin){
        lagged <- log(1 + lagged)
    }
    X <- array(c(intercept, lagged %*% W, lagged), dim=c(max_T-lags, N, 3))
    return(X)
}


single_node_samples <- function(mcmc_out, node=1, crp=FALSE, thin=1){
    # return all the samples corresponding to a single node
    num_samples <- dim(mcmc_out$ch)[2]; T_max <- num_samples %/% thin
    samples <- matrix(0, 3, T_max)
    for (t in thin*(1:T_max)){
        cluster_index <- mcmc_out$ch[node, t]
        if (crp){
            samples[, t/thin] <- mcmc_out$p[[t]][, cluster_index]
        } else{
            samples[, t/thin] <- mcmc_out$p[t, , cluster_index]
        }
    }
    return(samples)
}

all_node_samples <- function(mcmc_out, node_coef=1, crp=FALSE, thin=1){
    # return all samples of coef values, grouped by 'type' of coefficient
    num_samples <- dim(mcmc_out$ch)[2] %/% thin
    N <- dim(mcmc_out$ch)[1]
    samples <- c()
    if (crp){
        for (i in 1:N){
            for (t in 1:num_samples){
                s <- mcmc_out$p[[t*thin]][node_coef, mcmc_out$ch[i, t*thin]]
                samples <- c(samples, s)
            }
        }
    } else{
        for (i in 1:N){
            cluster_inds <- mcmc_out$ch[i, thin*(1:num_samples)]
            samples <- c(samples, mcmc_out$p[cbind(thin*(1:num_samples), 
                                                   rep(node_coef, num_samples), cluster_inds)])
        }
    }
    return(samples)
}


FMM_concat_random_starts <- function(mcmc_list, thin=1){
    sample_inds <- thin*(1:(dim(mcmc_list[[1]]$ch)[2] %/% thin))
    cluster_samples <- mcmc_list[[1]]$ch[, sample_inds]
    param_samples <- mcmc_list[[1]]$p[sample_inds, , ]
    for (l in 2:length(mcmc_list)){
        cluster_samples <- cbind(cluster_samples, mcmc_list[[l]]$ch[, sample_inds])
        param_samples <- abind(param_samples, mcmc_list[[l]]$p[sample_inds, ,], along=1)
    }
    return(list('ch'=cluster_samples, 'p'=param_samples))
}


MCMC_compare <- function(mcmc_samples, trace=FALSE, labels=NULL){
    params_df <- melt(mcmc_samples)
    names(params_df) <- c("row", "column", "value")
    if (trace){
        ggplot(params_df, aes(x = column, y=value, color = factor(row))) +
            geom_line(linewidth=0.2) + theme(legend.position="none") + 
            labs(x="Iteration", y="Value") + 
            theme(axis.text = element_text(size = 18),  
                  axis.title = element_text(size = 18))
    } else{
        xmax <- sort(mcmc_samples, decreasing=TRUE)[50]
        plt <- ggplot(params_df, aes(x = value, color = factor(row))) +
            geom_density(alpha=0.6, linewidth=0.7, linetype='dashed') + 
            xlim(0, xmax) +
            labs(x = "Value",
                 y = "Density",
                 color = "Prior no.") + 
            theme(axis.text = element_text(size = 18),  
                axis.title = element_text(size = 18),
                legend.title = element_text(size = 18),
                legend.text = element_text(size = 18))
        if (!is.null(labels)){
            plt <- plt + scale_color_manual(labels=labels, values=hue_pal()(length(labels)))
        } else{
            plt <- plt + theme(legend.position="none")
        }
        plt
    }
}


param_summary <- function(model_name, mcmc_out=NULL, crp=FALSE, eco=TRUE,
                          coef_means=NULL, col_seq=NULL, breaks=NA){
    if (is.null(coef_means)){
        N <- dim(mcmc_out$ch)[1]
        summary_stats <- matrix(0, N, 6)
        colnames(summary_stats) <- c("theta1", "theta2", "theta3",
                                     "sigma1", "sigma2", "sigma3")
        for (i in 1:N){
            node_samples <- single_node_samples(mcmc_out, node=i, crp=crp)
            summary_stats[i, 1:3] <- apply(node_samples, 1, mean)
            summary_stats[i, 4:6] <- apply(node_samples, 1, sd)
        }
    } else{
        summary_stats <- coef_means
    }
    if (eco){
        row.names(summary_stats) <- counties_info$County
    } 
    if (is.null(col_seq)){
        col_seq <- cluster_MAP(mcmc_out$ch)
    }
    clusters <- split(seq_along(col_seq), col_seq)
    cluster_sizes <- lapply(clusters, length)
    c_order <- as.vector(unlist(clusters))
    summary_stats <- summary_stats[c_order, ]
    pheatmap(summary_stats[, 1:3], cluster_rows=F, cluster_cols=F, 
             angle_col=0, cell_width=4, gaps_row=cumsum(cluster_sizes),
             display_numbers=T, fontsize=15, width=3.5, breaks=breaks, 
             filename=paste(model_name, "param_sum", "png", sep="."),
             color=colorRampPalette(brewer.pal(n=4, name = "YlGnBu"))
             (length(breaks)),
             labels_col=expression(theta[1], theta[2], theta[3]))
}


node_coef_compare <- function(mcmc_list, node_coef=1, crp=FALSE, thin=1){
    L <- length(mcmc_list); num_samples <- dim(mcmc_list[[1]]$ch)[2] %/% thin
    N <- dim(mcmc_list[[1]]$ch)[1]
    all_samples <- matrix(0, L, N*num_samples)
    for (l in 1:L){
        all_samples[l, ] <- all_node_samples(mcmc_list[[l]], node_coef=node_coef,
                                             crp=crp, thin=thin)
    } 
    MCMC_compare(all_samples, trace=FALSE)
}


ARI_diff <- function(cluster_samples, MAP_cluster){
    num_samples <- dim(cluster_samples)[2]
    ARI_series <- rep(0, num_samples)
    for (t in 1:num_samples){
        ARI_series[t] <- adjustedRandIndex(cluster_samples[,t], MAP_cluster)
    }
    return(ARI_series)
}


simulation_boxplot <- function(true_membs, true_params, mean_coefs_array, title){
    N <- length(true_membs); num_params <- dim(true_params)[1]
    K <- dim(true_params)[2]; L <- dim(mean_coefs_array)[1]
    colnames(true_params) <- 1:K; rownames(true_params) <- 1:num_params
    density_list <- matrix(NA, ncol=3)
    for (l in 1:L){
        for (i in 1:N){
            z <- true_membs[i]
            data <- matrix(c(1:3, rep(z, 3), mean_coefs_array[l, i, ]), nrow=3)
            density_list <- rbind(density_list, data)
        }
    }
    df <- as.data.frame(density_list[2:nrow(density_list), ])
    colnames(df) <- c("param", "clust", "diff")
    ggplot(data=df, mapping=aes(x=as.factor(param), y=diff, color=as.factor(clust))) + 
        geom_boxplot() + 
        geom_point(data=as.data.frame(as.table(true_params)), 
                   mapping=aes(x=as.factor(Var1), y=Freq, color=as.factor(Var2)), 
                   position=position_dodge(width=0.75),
                   shape=4, size=3, stroke=1.5) +
        labs(x = "Coefficient", color="Cluster", y="Posterior mean node coefficient") +
        scale_x_discrete(labels=expression(theta[1], theta[2], theta[3])) +
        theme(legend.position = c(0.85, 0.7),
              axis.text = element_text(size = 18),  
              axis.title = element_text(size = 18),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 18))
    ggsave(paste(title, "png", sep="."))
}


sim_mean_coefs <- function(true_membs, true_params, mcmc_out, crp=FALSE){
    N <- length(true_membs); nparams <- dim(true_params)[1]
    mean_node_coefs <- matrix(0, N, nparams); mse <- 0; hpd_prop <- 0; 
    ess <- matrix(0, N, nparams)
    for (i in 1:N){
        samples <- single_node_samples(mcmc_out, node=i, crp=crp)
        hpd <- coda::HPDinterval(as.mcmc(t(samples)), prob=0.90)
        mean_node_coefs[i, ] <- apply(samples, 1, mean)
        z <- true_membs[i]
        for (p in 1:nparams){
            if ((hpd[p, 1] < true_params[p, z]) & (true_params[p, z] < hpd[p, 2])){
                hpd_prop <- hpd_prop + 1
            }
        }
        ess[i, ] <- coda::effectiveSize(as.mcmc(t(samples)))
        mse <- mse + sum((mean_node_coefs[i, ] - true_params[, z])^2)
    }
    return(list("mean_coefs"=mean_node_coefs, "mse"=mse/(N*nparams),
                "hpd"= hpd_prop / (N*nparams), "ess"=median(ess)))
}

mean_node_coefs <- function(mcmc_out, crp=T){
    N <- dim(mcmc_out$ch)[1]; nparams <- dim(mcmc_out$p[[1]])[1]
    mean_node_coefs <- matrix(0, N, nparams)
    for (i in 1:N){
        samples <- single_node_samples(mcmc_out, node=i, crp=crp)
        mean_node_coefs[i, ] <- apply(samples, 1, mean)
    }
    return(mean_node_coefs)
}


coef_quantiles <- function(mcmc_out, quantile=0.75, crp=FALSE){
    N <- dim(mcmc_out$ch)[1]; 
    coef_quantiles <- matrix(0, N, 3); 
    for (i in 1:N){
        samples <- single_node_samples(mcmc_out, node=i, crp=crp)
        coef_quantiles[i, ] <- apply(samples, 1, function(x) {
            quantile(x, probs=c(quantile))})
    }
    return(coef_quantiles)
}


PIT_hist <- function(modfit, title, gauss=F){
    if (gauss){
        xlab <- "PIT"
    } else{
        xlab <- "Randomised PIT"
    }
    png(filename=paste("PIT", title, "png", sep="."))
    par(mar=c(5, 5, 0.5, 1.5))
    hist(as.vector(modfit$PIT), ylim=c(0, 400),xlim=c(0, 1),
         main=NULL, ylab="Count", xlab=xlab, cex.lab=2, cex.axis=2)
    dev.off()
}

