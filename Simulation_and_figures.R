prior_sensitivity <- function(y, X, filename, graph, true_membs, true_params, 
                              max_iters=11000, thin=10, prior9comps=3){
    aff_expdist <- affinity(as.matrix(distances(graph)), function(x) exp(-x))
    crp_aff <- matrix(1, nrow=dim(y)[2], ncol=dim(y)[2])
    diag(crp_aff) <- 0
    aff_expdist4 <- affinity(as.matrix(distances(graph)), function(x) exp(-x/4))
    prior1 <- GPNAR(y, X, aff_expdist, max_iters=max_iters)
    prior2 <- GPNAR(y, X, crp_aff, max_iters=max_iters)
    prior3 <- GPNAR(y, X, aff_expdist4, max_iters=max_iters)
    prior4 <- GPNAR(y, X, aff_expdist, alpha=1/4, max_iters=max_iters)
    prior5 <- GPNAR(y, X, aff_expdist, shape=1/4, rate=1/4, max_iters=max_iters)
    prior6 <- GPNAR(y, X, aff_expdist, shape=1/4, rate=1/2, max_iters=max_iters)
    prior7 <- GPNAR(y, X, aff_expdist, shape=4, rate=16, max_iters=max_iters)
    prior8 <- FMM(y, X, alpha=c(1, 1, 1, 1), max_iters=max_iters)
    prior9 <- FMM(y, X, alpha=rep(1, prior9comps), max_iters=max_iters)
    
    mcmc_list <- list(prior1, prior2, prior3, prior4, prior5, prior6,
                      prior7, prior8, prior9)
    crp <- c(rep(T, 7), F, F)
    plot1 <- c(1, 2, 3, 4, 8, 9)
    plot2 <- c(1, 5, 6, 7)
    num_samples <- dim(mcmc_list[[1]]$ch)[2] %/% thin
    N <- dim(mcmc_list[[1]]$ch)[1]
    coef_metrics_mat <- matrix(0, 3, 9)
    # for (p in mcmc_list){
    #     coef_metrics <- sim_mean_coefs()
    # }
    coef_xlab <- c(expression(theta[1]), expression(theta[2]), expression(theta[3]))
    all_samples <- matrix(0, length(plot1), N*num_samples)
    for (coef in 1:3){
        for (l in plot1){
            all_samples[which(plot1==l), ] <- all_node_samples(mcmc_list[[l]], 
                                                               node_coef=coef, crp=crp[l], thin=thin)
        }
        plt <- MCMC_compare(all_samples, trace=FALSE, labels=plot1)
        plt <- plt + labs(y="Posterior density estimate",
                          x=coef_xlab[coef]) 
        ggsave(filename=paste(filename, 1, coef, "png", sep="."), plt)
    }
    all_samples <- matrix(0, length(plot2), N*num_samples)
    for (coef in 1:3){
        for (l in plot2){
            all_samples[which(plot2==l), ] <- all_node_samples(mcmc_list[[l]], 
                                                               node_coef=coef, crp=crp[l], thin=thin)
        }
        plt <- MCMC_compare(all_samples, trace=FALSE, labels=plot2)
        plt <- plt + labs(y="Posterior density estimate",
                          x=coef_xlab[coef]) 
        ggsave(filename=paste(filename, 2, coef, "png", sep="."), plt)
        
    }
}



coocc_mat_plot <- function(cc, clusters, filename, names=NULL){
    c_order <- as.vector(unlist(split(seq_along(clusters), clusters)))
    if (is.null(names)){
        rownames(cc) <- 1:(dim(cc)[2])
    } else {
        rownames(cc) <- names
    }
    co_occ <- cc[c_order, c_order]
    png(paste(filename, "heatmap", "png", sep="."))
    pheatmap::pheatmap(co_occ, cluster_rows=F, cluster_cols=F, 
                       breaks=seq(0, 1, length.out=101), fontsize=18)
    dev.off()
}



sim_GPNAR <- function(adj_mat, cluster_memberships, params, T_max=25,
                      population=NULL, burn_in=100, log_lin=FALSE){
    N <- nrow(adj_mat)
    if (is.null(population)){
        population <- rep(1, N)
    } else{
        population <- population / 100000
    }
    intercept <- population
    network_mat <- diag(population) %*% (adj_mat / rowSums(adj_mat)) %*% diag(1/population)
    self_mat <- diag(rep(1, N))
    K <- max(cluster_memberships)
    ts <- matrix(0, T_max + burn_in, N)
    for (k in 1:K){
        ki <- which(cluster_memberships == k)
        intercept[ki] <- intercept[ki]*params[1, k]
        network_mat[ki, ] <- network_mat[ki, ]*params[2, k]
        self_mat[ki, ki] <- self_mat[ki, ki]*params[3, k]
        ts[1, ki] <- rpois(length(ki), params[1, k])
    }
    for (t in 2:(T_max+burn_in)){
        if (log_lin){
            expected <- exp(intercept + (network_mat + self_mat) %*% log(ts[t-1, ]+1))
        } else{
            expected <- intercept + (network_mat + self_mat) %*% ts[t-1, ]
        }
        ts[t, ] <- rpois(N, expected)
    }
    return(ts[(burn_in+1):(T_max+burn_in), ])
}

sim_results <- function(title, memberships, params, graph, ts_burn_in=100, 
                        nsims=3, T_max=25, population=NULL, max_iters=11000, 
                        C=5, stack=FALSE, alpha=c(1, 1, 1, 1), sim_adj=NULL){
    png(paste(title, "graph.png", sep="_"))
    par(mar=c(0,0,0,0)+.1)
    plot(graph, vertex.color=hue_pal()(dim(params)[2])[memberships])
    dev.off()
    N <- length(memberships); nparams <- dim(params)[1]
    DDP_means <- array(0, dim=c(nsims, N, nparams))
    FMM_means <- array(0, dim=c(nsims, N, nparams))
    metrics <- matrix(0, nrow=10, ncol=nsims)
    rownames(metrics) <- c("DDP ARI", "DDP mse", "DDP hpd", 
                           "DDP ess", "DDP mixing",
                           "FMM ARI", "FMM mse", "FMM hpd", 
                           "FMM ess", "FMM mixing")
    if (stack){
        stack_iters <- max_iters
        max_iters <- (max_iters-1000)*C + 1000
        stack_metrics <- matrix(0, nrow=2, ncol=nsims)
        rownames(stack_metrics) <- c("Stack ARI", "Stack mse")
        metrics <- rbind(metrics, stack_metrics)
        stack_means <- array(0, dim=c(nsims, N, nparams))
    }
    if (is.null(sim_adj)){
        adj_mat <- igraph::as_adjacency_matrix(graph, sparse=F)
    } else{
        adj_mat <- sim_adj
    }
    aff_mat <- affinity(as.matrix(distances(graph)), function(x) exp(-x))
    for (s in 1:nsims){
        sim_ts <- sim_GPNAR(adj_mat, memberships, params, T_max=T_max, 
                            burn_in=ts_burn_in, population=population)
        if (is.null(population)){
            X <- X_from_ts(sim_ts, adj_mat)
        } else{
            X <- X_from_ts(sim_ts, adj_mat, pop_adjust=population, pop_adjust_intercept=T)
        }

        DDP_fit <- DDP(sim_ts[2:T_max, ], X, aff_mat, max_iters=max_iters)
        MAP_cluster <- cluster_MAP(DDP_fit$ch)
        metrics["DDP ARI", s] <- adjustedRandIndex(memberships, MAP_cluster)
        metrics["DDP mixing", s] <- mean(ARI_diff(DDP_fit$ch, MAP_cluster))
        coef_metrics <- sim_mean_coefs(memberships, params, DDP_fit, crp=T)
        metrics["DDP mse", s] <- coef_metrics$mse
        metrics["DDP hpd", s] <- coef_metrics$hpd
        metrics["DDP ess", s] <- coef_metrics$ess
        DDP_means[s, , ] <- coef_metrics$mean_coefs
        
        FMM_fit <- FMM(sim_ts[2:T_max, ], X, alpha=alpha, max_iters=max_iters)
        MAP_cluster <- cluster_MAP(FMM_fit$ch)
        metrics["FMM ARI", s] <- adjustedRandIndex(memberships, MAP_cluster)
        metrics["FMM mixing", s] <- mean(ARI_diff(FMM_fit$ch, MAP_cluster))
        coef_metrics <- sim_mean_coefs(memberships, params, FMM_fit)
        metrics["FMM mse", s] <- coef_metrics$mse
        metrics["FMM hpd", s] <- coef_metrics$hpd
        metrics["FMM ess", s] <- coef_metrics$ess
        FMM_means[s, , ] <- coef_metrics$mean_coefs
        
        if (stack){
            stack_fit <- FMM_stack(sim_ts[2:T_max, ], X, alpha=alpha,
                                   max_iters=stack_iters, C=C)
            cc <- stack_co_occurence(stack_fit$chains, stack_fit$weights)
            opt_clusts <- cluster_MAP(FMM_concat_random_starts(stack_fit$chains)$ch, cc=cc)
            metrics["Stack ARI", s] <- adjustedRandIndex(memberships, opt_clusts)
            stack_means[s, , ] <- stack_mean_coefs(stack_fit$chains, stack_fit$weights)
            metrics["Stack mse", s] <- mean((stack_means[s, , ] - 
                                                 t(params)[memberships, ])^2)
        }
    }
    print(apply(metrics, 1, mean))
    simulation_boxplot(memberships, params, DDP_means, paste(title, "DDP"))
    simulation_boxplot(memberships, params, FMM_means, paste(title, "FFM"))
    if (stack){
        simulation_boxplot(memberships, params, stack_means, paste(title, "stack"))
        return(list("metrics"=metrics, "DDP_means"=DDP_means, 
                    "FMM_means"=FMM_means, "stack_means"=stack_means))
    } else{
        return(list("metrics"=metrics, "DDP_means"=DDP_means, 
                    "FMM_means"=FMM_means))
    }
}

eco_graph <- function(error, memberships, filename, maxval=NULL){
    if (is.null(maxval)){
        maxval <- max(error)
    }
    normalised_error <- (error/maxval)
    g <- eco_igraph
    num_clusters <- max(memberships)
    cols <- hue_pal()(num_clusters)[memberships]
    V(g)$color <- mapply(function(col, alpha) {
        adjustcolor(col, alpha.f=alpha)}, cols, normalised_error)
    layout=cbind(counties_info$Long, counties_info$Lat)
    par(mar=c(0,0,0,0)+0.01)
    png(paste("errorgraph", filename, "png", sep="."), height=700)
    plot(eco_igraph, vertex.color=V(g)$color, vertex.label.cex=1.8, 
         layout=layout, vertex.frame.color=cols,
         vertex.frame.width=2, asp=1.5, vertex.label.color="black", 
         size=3, margin=rep(0, 4))
    dev.off()
}
