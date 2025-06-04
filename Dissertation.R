
# load required libraries
library('dplyr'); library('igraph'); library('tidyverse'); library('matrixStats')
library('MASS'); library('mvtnorm'); library('MCMCpack'); library('distr')
library('pheatmap'); library('gridExtra')
library('reshape'); library('mclust'); library('coda');
library('xtable'); library('Metrics'); library('scales'); library('abind');
library('RColorBrewer')

# load functions from other R scripts
source("Analysis.R")
source("Stacking.R")
source("Simulation_and_figures.R")
source("Model_fitting.R")
source("Helper_functions.R")
source("gagnar.R")


##### DATA PRE-PROCESSING #####

# Extract datasets for the weekly case counts
covid <- read.csv("COVID-19.csv")
covid$TimeStamp <- as.POSIXct(covid$TimeStamp)

# extract information about counties
counties_info <- covid[1:26, c("CountyName", "PopulationCensus16", "Lat", "Long")]
names(counties_info) <- c("County", "Pop", "Lat", "Long")

# keep weekly observations starting from 1 Mar 2020
covid <- covid[, c("CountyName", "TimeStamp", "ConfirmedCovidCases")] %>% 
    filter(weekdays(TimeStamp) == 'Sunday')
covid_ts <- spread(covid, TimeStamp, ConfirmedCovidCases)
# Turn cumulative cases into increase in the number of cases
covid_ts[, 3:ncol(covid_ts)] <- covid_ts[, 3:ncol(covid_ts)] - 
    covid_ts[, 2:(ncol(covid_ts)-1)]
covid_ts <- t(as.matrix(covid_ts[, 2:ncol(covid_ts)]))

# Obtain the economics hub igraph with alphabetically ordered vertices
load('eco_hubs.Rdata')
eco_adj <- igraph::as_adjacency_matrix(covid_net_eco_hubs_igraph)
eco_adj <- as.matrix(eco_adj)
o <- order(county_index_eco_hubs$CountyName)
eco_adj <- eco_adj[o, o]
row.names(eco_adj) <- 1:26
colnames(eco_adj) <- 1:26
eco_igraph <- igraph::graph_from_adjacency_matrix(eco_adj, mode='undirected')
V(eco_igraph)$name <- counties_info$County

X_count <- X_from_ts(covid_ts, eco_adj, pop_adjust=counties_info$Pop,
                     pop_adjust_intercept=TRUE)
plot(eco_igraph, layout=cbind(counties_info$Long, counties_info$Lat))
eco_distances <- as.matrix(distances(eco_igraph))
eco_aff_expdist <- affinity(eco_distances, function(x) exp(-x))

# obtain differenced series from the covid paper
covid_lag1 <- read.csv("covid_lag1.csv")
covid_lag1 <- covid_lag1[, c("CountyName", "yw", "weeklyCases")]
covid_lag1_agg <- spread(covid_lag1, yw, weeklyCases)
covid_lag1_agg <- t(as.matrix(covid_lag1_agg[, 2:ncol(covid_lag1_agg)]))
X_diff <- X_from_ts(covid_lag1_agg, eco_adj)

# obtain the predictors for the PNAR model
X_pnar <- X_from_ts(covid_ts[1:100, ], eco_adj)
X_log <- X_from_ts(covid_ts[1:100, ], eco_adj, log_lin=TRUE)


##### INITIAL PLOTS and MCMC RUN ######
# line plot of of log case counts
plot_cases <- log(X_count[1:80, 1:6, 3]+1)
colnames(plot_cases) <- counties_info$County[1:6]
df <- data.frame(time = seq(from = as.Date("2020-03-01"), by = "week", length.out = 80), 
                 plot_cases)
df_long <- melt(df, id.vars = "time")
ggplot(df_long, aes(x = time, y = value, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "ln(weekly cases +1)", 
         color = "County") +
    geom_vline(xintercept=df$time[25], linetype='dashed') + 
    geom_vline(xintercept=df$time[44], linetype='dashed') + 
    geom_vline(xintercept=df$time[63], linetype='dashed') + 
    theme(text=element_text(size=12)) + 
    scale_x_continuous(breaks=c(df$time[1], df$time[25], df$time[44], df$time[63]))
ggsave("case_plot.png", width=6, height=4)

# plot the economic hubs network
eco_hubs <- c("Cork", "Dublin", "Galway", "Limerick", "Waterford")

png("eco_hubs_network.png", height=700)
hub_indicator <- rep(0, 26)
hub_indicator[which(counties_info$County %in% eco_hubs)] <- 1
plot(eco_igraph, vertex.label.cex=1.8, 
     vertex.color= c("grey", "orange")[hub_indicator+1],
     layout=cbind(counties_info$Long, counties_info$Lat), 
     asp=1.5, vertex.label.color=c("black","blue")[hub_indicator+1], 
     vertex.size=7)
dev.off()


# initial try for fitting FMM4 model to segment 1 data
X <- X_count[1:23, ,]; y <- X_count[2:24, , 3];
y_true <- X_count[25, , 3]; X_pred <- X_count[24, , ];
y_ts <- X_count[1:25, , 3]; 

seg1_FMM4 <- FMM(y, X)
clusts4 <- cluster_MAP(seg1_FMM4$ch)
# generates heatmap of mean node coefficients, grouped by clusters
param_summary('seg1_FMM4', mcmc_out=seg1_FMM4)

# checking MCMC chain mixing
ddp_mixing <- list()
for (i in 1:10){
    print(i)
    s <- DDP(y, X, eco_aff_expdist)
    ddp_mixing[[i]] <- s
}

FMM4_mixing <- list()
for (i in 1:10){
    print(i)
    s <- FMM(y, X)
    FMM4_mixing[[i]] <- s
}

FMM5_mixing <- list()
for (i in 1:10){
    s <- FMM(y, X, alpha=c(1, 1, 1, 1, 1))
    FMM5_mixing[[i]] <- s
}

# trace plots and posterior density plots 
samples <- matrix(0, 5, 9000)
for (i in 1:5){
    samples[i, ] <- single_node_samples(ddp_mixing[[i]], node=3, crp=TRUE)[2, ]
}
ggsave(filename="ddp_mix_trace.png", MCMC_compare(samples, trace=TRUE))
ggsave(filename="ddp_mix_density.png", MCMC_compare(samples))

samples <- matrix(0, 5, 9000)
for (i in 1:5){
    samples[i, ] <- single_node_samples(FMM4_mixing[[i]], node=3)[2, ]
}
ggsave(filename="FMM4_mix_trace.png", MCMC_compare(samples, trace=TRUE))
ggsave(filename="FMM4_mix_density.png", MCMC_compare(samples))

samples <- matrix(0, 5, 9000)
for (i in 1:5){
    samples[i, ] <- single_node_samples(FMM5_mixing[[i]], node=3)[2, ]
}
ggsave(filename="FMM5_mix_trace.png", MCMC_compare(samples, trace=TRUE))
ggsave(filename="FMM5_mix_density.png", MCMC_compare(samples))

# thin to speed up computations
ggsave(filename="ddp_coef2.png", node_coef_compare(ddp_mixing, crp=T, thin=50, node_coef=2))
ggsave(filename="fmm4_coef2.png", node_coef_compare(FMM4_mixing, node_coef=2, thin=50))
ggsave(filename="fmm5_coef2.png", node_coef_compare(FMM5_mixing, node_coef=2, thin=50))

ddp_partitions <- lapply(ddp_mixing, function(x) mean(ARI_diff(x$ch, cluster_MAP(x$ch))))
print(paste("DDP - mean ARI: samples to ls = ", mean(unlist(ddp_partitions))))
fmm5_partitions <- lapply(FMM5_mixing, function(x) mean(ARI_diff(x$ch, cluster_MAP(x$ch))))
print(paste("FMM5 - mean ARI: samples to ls = ", mean(unlist(fmm5_partitions))))
fmm4_partitions <- lapply(FMM4_mixing, function(x) mean(ARI_diff(x$ch, cluster_MAP(x$ch))))
print(paste("FMM4 - mean ARI: samples to ls = ", mean(unlist(fmm4_partitions))))


##### SIMULATIONS + PRIOR SENSITIVITY #####
pref.probs <- matrix(0.04, nrow=4, ncol=4)
pref.probs <- pref.probs + diag(0.36, nrow=4)
sbm_igraph <- igraph::sample_sbm(40, pref.probs, c(10, 10, 10, 10))
sbm_adj <- as_adjacency_matrix(sbm_igraph, sparse=F)
sbm_memberships <- c(rep(1, 10), rep(2, 10), rep(3, 10), rep(4, 10))

ws_igraph <- sample_smallworld(1, size=30, nei=5, p=0.05)
ws_adj <- as_adjacency_matrix(ws_igraph, sparse=F)
ws_memberships <- c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 15))

params1 <- matrix(c(0.5, 0.8, 1.3, 1.7,
                    0.1, 0.2, 0.3, 0.4,
                    0.2, 0.3, 0.4, 0.5), nrow=3, byrow=T)
sim_sbm1 <- sim_results("sim_sbm1", sbm_memberships, params1, sbm_igraph)
sim_ws1 <- sim_results("sim_ws1", ws_memberships, params1, ws_igraph)
sim_eco1 <- sim_results("sim_eco1", clusts4, params1, eco_igraph, 
                        population=counties_info$Pop)

params2 <- matrix(c(1, 1.5, 1, 1.5,
                    0.3, 0.1, 0.4, 0.2,
                    0.1, 0.8, 0.2, 0.5), nrow=3, byrow=T)
sim_sbm2 <- sim_results("sim_sbm2", sbm_memberships, params2, sbm_igraph)
sim_ws2 <- sim_results("sim_ws2", ws_memberships, params2, ws_igraph)
sim_eco2 <- sim_results("sim_eco2", clusts4, params2, eco_igraph, 
                        population=counties_info$Pop, ts_burn_in=0)

# cluster coefs roughly equal to those estimated
eco_params <- matrix(0, nrow=3, ncol=4)
eco_params[, clusts4[5]] <- c(0.9, 0.01, 0.9)
eco_params[, clusts4[6]] <- c(2.9, 0.01, 0.9)
eco_params[, clusts4[1]] <- c(1.3, 0.7, 0.01)
eco_params[, clusts4[3]] <- c(1.8, 1.5, 0.1)
sim_sbm_eco <- sim_results("sim_sbm_eco", sbm_memberships, eco_params, 
                           sbm_igraph, ts_burn_in=0)
sim_ws_eco <- sim_results("sim_ws_eco", ws_memberships, eco_params, ws_igraph, 
                          ts_burn_in=0)
sim_eco_eco <- sim_results("sim_eco_eco", clusts4, eco_params, eco_igraph,
                           ts_burn_in=0, population=counties_info$Pop)

# simulations for chain stacking
sim_eco_eco_stack <- sim_results("sim_eco_eco_stack", clusts4, eco_params, 
                                 eco_igraph, ts_burn_in=0, stack=T, C=5, max_iters=6000,
                                 population=counties_info$Pop)
set.seed(1)   # for nonzero rows
eco_adj_missing <- eco_adj
edge_list <- which(eco_adj == 1, arr.ind=T)
samples <- sample.int(nrow(edge_list), size=nrow(edge_list)*0.2)
eco_adj_missing[edge_list[samples, ]] <- 0
sim_eco_stack_net <- sim_results("sim_eco_stack_net", clusts4, eco_params, 
                                 eco_igraph, ts_burn_in=0, stack=T, C=5, max_iters=6000,
                                 sim_adj=eco_adj_missing, population=counties_info$Pop)

sim_eco_stack3 <- sim_results("sim_eco_stack3", sbm_memberships, eco_params, 
                              sbm_igraph, ts_burn_in=0, stack=T, C=5, max_iters=6000,
                              alpha=c(1, 1, 1))

# Prior sensitivity
sim_ts <- sim_GPNAR(eco_adj, clusts4, eco_params, population=counties_info$Pop, burn_in=0)
prior_sensitivity(sim_ts[2:25, ], X_from_ts(sim_ts, eco_adj), graph=eco_igraph,
                  filename="eco_sens")
prior_sensitivity(y, X, graph=eco_igraph, filename="true_sens", max_iters=11000,
                  prior9comp=5)


##### MODEL FITTING TO COVID19 DATA #####
times <- c()
start <- Sys.time()
ddp <- DDP(y, X, eco_aff_expdist, max_iters=100000, burn_in=5000)
times <- c(times, Sys.time()-start)
start <- Sys.time()
FMM4 <- FMM_stack(y, X, C=5, max_iters=25000, burn_in=5000)
times <- c(times, Sys.time()-start)
start <- Sys.time()
FMM5 <- FMM_stack(y, X, C=5, max_iters=25000, burn_in=5000, 
                  alpha=c(1, 1, 1, 1, 1))
times <- c(times, Sys.time()-start)

DDP_forecast <- one_step_forecast(X_pred, ddp, crp=TRUE)
FMM4_mean_coefs <- stack_mean_coefs(FMM4$chains, FMM4$weights)
FMM4_yhat <- rowSums(X_pred * FMM4_mean_coefs)
FMM5_mean_coefs <- stack_mean_coefs(FMM5$chains, FMM5$weights)
FMM5_yhat <- rowSums(X_pred * FMM5_mean_coefs)

# fitting the PNAR linear model 
PNAR_coefs <- PNAR::lin_estimnarpq(X_count[1:24, , 3], eco_adj, p=1, uncons=T)
PNAR_yhat <- X_pnar[24, , ] %*% PNAR_coefs$coefs$Estimate

# fitting the PNAR log-linear model
PNARl_coefs <- PNAR::log_lin_estimnarpq(X_count[1:24, , 3], eco_adj, p=1, uncons=T)
PNARl_yhat <- exp(X_log[24, , ] %*% PNARl_coefs$coefs$Estimate)

PNARa_coefs <- PNAR_X(y, X)
PNARa_yhat <- X_pred %*% PNARa_coefs

# fitting the GAGNAR model for 1-lag differenced series
start <- Sys.time()
gagnar_diff <- GAGNAR(covid_lag1_agg[2:24, ], X_diff[1:23, ,], eco_distances,
                      niterations=100000, burn_in=5000)
times <- c(times, Sys.time()-start)
gagnar_diff_all <- gagnar_diff; dnum_samples <- dim(gagnar_diff$ch)[2]

# fitting the GAGNAR model for count series
start <- Sys.time()
gagnar <- GAGNAR(y, X, eco_distances, niterations=150000, burn_in=5000)
times <- c(times, Sys.time()-start)
gagnar_all <- gagnar; num_samples <- dim(gagnar$ch)[2]


# discard more samples based on Geweke's test
gagnar_diff <- list("ch"=gagnar_all$ch[, 5000:dnum_samples],
                    "p"=gagnar_all$p[5000:dnum_samples],
                    "sigma"=gagnar_all$sigma[5000:dnum_samples])
gstart <- 60000
gagnar <- list("ch"=gagnar_all$ch[, gstart:num_samples],
               "p"=gagnar_all$p[gstart:num_samples],
               "sigma"=gagnar_all$sigma[gstart:num_samples])
gagnar_diag <- mcmc_diagnostics(gagnar)

# stationarity and ESS
ddp_diag <- mcmc_diagnostics(ddp)
gagnar_diff_diag <- mcmc_diagnostics(gagnar_diff)
gagnar_diag <- mcmc_diagnostics(gagnar)

# histograms of Geweke's statistic
png("Geweke_DDP.png")
par(mar=c(5, 5, 0.5, 3.5))
hist(ddp_diag$geweke, xlab="Geweke's statistic", ylab="Count",
     main=NULL, cex.lab=2, cex.axis=2)
dev.off()
png("Geweke_gagnar_diff.png")
par(mar=c(5, 5, 0.5, 3.5))
hist(gagnar_diff_diag$geweke, xlab="Geweke's statistic", ylab="Count",
     main=NULL, cex.lab=2, cex.axis=2)
dev.off()
png("Geweke_gagnar.png")
par(mar=c(5, 5, 0.5, 3.5))
hist(gagnar_diag$geweke, xlab="Geweke's statistic", ylab="Count",
     main=NULL, cex.lab=2, cex.axis=2)
dev.off()

# create a table of the effective sample sizes 
ess_tab <- matrix(0, nrow=26, ncol=9)
ess_tab[, 1:3] <- ddp_diag$ESS
ess_tab[, 4:6] <- gagnar_diff_diag$ESS
ess_tab[, 7:9] <- gagnar_diag$ESS
ess_df <- data.frame("County"=counties_info$County)
ess_df[, 2:10] <- ess_tab
print(xtable(ess_df, digits=0))

gagnar_diff_forecast <- one_step_forecast(X_diff[24, ,], gagnar_diff, crp=TRUE, 
                                          pop_adjust=counties_info$Pop,
                                          last_step=X_count[24, , 3])
gagnar_forecast <- one_step_forecast(X_pred, gagnar, crp=TRUE)

# Table for model comparisons
point_forecast_df <- data.frame("County"=counties_info$County, 
                                "Y_25"= X_count[25, , 3],
                                "DDP"=DDP_forecast$y_hat, 
                                "FMM 5"=FMM5_yhat, 
                                "FMM 4"=FMM4_yhat,
                                "GAGNAR diff"=gagnar_diff_forecast$y_hat,
                                "GAGNAR" = gagnar_forecast$y_hat,
                                "PNAR"=PNAR_yhat, 
                                "PNAR ll"=PNARl_yhat,
                                "Adj PNAR"=PNARa_yhat)

se_indiv <- apply(point_forecast_df[, 3:ncol(point_forecast_df)], 2, 
                  function(y) scaled_error(y_ts, y))

mase <- apply(point_forecast_df[, 3:ncol(point_forecast_df)], 2, 
              function(y) mean(scaled_error(y_ts, y)))
point_forecast_df[27, ] <- c(NA, NA, mase)

# 'sharpness' diagnostic -log P(Y=y) and randomised PITs
thin <- 10
DDP_modfit <- model_fit(y, mcmc_out=ddp, crp=TRUE, X_pred=X, thin=thin)
FMM5_modfit <- stack_score(y, X, FMM5, thin=thin)
FMM4_modfit <- stack_score(y, X, FMM4, thin=thin)
GAGNAR_diff_modfit <- model_fit(y, mcmc_out=gagnar_diff, X_pred=X_diff[1:23, , ], 
                                thin=thin, crp=T, gauss=T)
GAGNAR_modfit <- model_fit(y, mcmc_out=gagnar, X_pred=X, thin=thin, crp=T, gauss=T)
PNAR_means <- apply(sweep(X_pnar[1:23, ,], 3, PNAR_coefs$coefs$Estimate, '*'), 
                    c(1, 2), sum)
PNAR_modfit <- model_fit(y, poisson_means=PNAR_means)
PNARl_means <- exp(apply(sweep(X_log[1:23, ,], 3, PNARl_coefs$coefs$Estimate, '*'), 
                         c(1, 2), sum))
PNARl_modfit <- model_fit(y, poisson_means=PNARl_means)
PNARa_modfit <- model_fit(y, poisson_means=t(apply(X, 1, function(x) x %*% PNARa_coefs)))
point_forecast_df[28, ] <- c(NA, NA, DDP_modfit$score, FMM5_modfit$score, FMM4_modfit$score,
                             GAGNAR_diff_modfit$score, GAGNAR_modfit$score, PNAR_modfit$score,
                             PNARl_modfit$score, PNARa_modfit$score)

# prediction scores for time-point 25
DDP_predfit <- model_fit(y_true, mcmc_out=ddp, crp=TRUE, X_pred=X_pred, thin=thin)
FMM5_predfit <- stack_score(y_true, X_pred, FMM5, thin=thin)
FMM4_predfit <- stack_score(y_true, X_pred, FMM4, thin=thin)
GAGNAR_diff_predfit <- model_fit(y_true, mcmc_out=gagnar_diff, crp=TRUE, 
                                 X_pred=X_diff[24, , ], thin=thin, gauss=T)
GAGNAR_predfit <- model_fit(y_true, mcmc_out=gagnar, crp=TRUE, 
                            X_pred=X_pred, thin=thin, gauss=T)
PNAR_predfit <- model_fit(y_true, poisson_means=PNAR_yhat)
PNARl_predfit <- model_fit(y_true, poisson_means=PNARl_yhat)
PNARa_predfit <- model_fit(y_true, poisson_means=PNARa_yhat)
point_forecast_df[29, ] <- c(NA, NA, DDP_predfit$score, FMM5_predfit$score, 
                             FMM4_predfit$score, GAGNAR_diff_predfit$score, GAGNAR_predfit$score,
                             PNAR_predfit$score, PNARl_predfit$score, PNARa_predfit$score)
point_forecast_df[27:29, 1] <- c("MASE", "Training score", "Test score")
point_forecast_df[1:26, 3:ncol(point_forecast_df)] <- se_indiv

xtable(point_forecast_df, digits=c(0, 0, 0, rep(2, 8)))

# histograms of (randomised) probability integral transforms
modfit_list <- list(DDP_modfit, FMM5_modfit, FMM4_modfit, PNAR_modfit, PNARl_modfit, 
                    PNARa_modfit, GAGNAR_modfit, GAGNAR_diff_modfit)
PIT_titles <- c("DDP", "FMM 5 clusters", "FMM 4 clusters", "PNAR", "PNAR log-linear",
                "Adjusted PNAR", "GAGNAR", "GAGNAR_diff")
gauss_tf <- c(rep(F, 6), T, T)
mapply(PIT_hist, modfit_list, PIT_titles, gauss_tf)


# obtain mean node coefficients and co-occurence matrices 
DDP_cc <- co_occurence(ddp$ch)
DDP_cMAP <- cluster_MAP(ddp$ch, cc=DDP_cc)
DDP_mean_coefs <- mean_node_coefs(ddp)
FMM5_cc <- stack_co_occurence(FMM5$chains, FMM5$weights)
FMM5_cMAP <- cluster_MAP(FMM_concat_random_starts(FMM5$chains)$ch, cc=FMM5_cc)
FMM5_mean_coefs <- stack_mean_coefs(FMM5$chains, FMM5$weights)
FMM4_cc <- stack_co_occurence(FMM4$chains, FMM4$weights)
FMM4_cMAP <- cluster_MAP(FMM_concat_random_starts(FMM4$chains)$ch, cc=FMM4_cc)
FMM4_mean_coefs <- stack_mean_coefs(FMM4$chains, FMM4$weights)
GAGNAR_diff_cc <- co_occurence(gagnar_diff$ch, thin=10)
GAGNAR_diff_cMAP <- cluster_MAP(gagnar_diff$ch, cc=GAGNAR_diff_cc)
GAGNAR_diff_mean_coefs <- mean_node_coefs(gagnar_diff)
GAGNAR_cc <- co_occurence(gagnar$ch)
GAGNAR_cMAP <- cluster_MAP(gagnar$ch, cc=GAGNAR_cc)
GAGNAR_mean_coefs <- mean_node_coefs(gagnar)

# create visualisations
max_col <- max(GAGNAR_mean_coefs); min_col <- min(GAGNAR_mean_coefs)
colour_breaks <- c(seq(min_col, -0.05, length.out=10),
                   seq(0, 2.5, length.out=50),
                   seq(2.55, max_col, length.out=20))
param_summary('ddp', coef_means=DDP_mean_coefs, col_seq=DDP_cMAP, breaks=colour_breaks)
param_summary('fmm5', coef_means=FMM5_mean_coefs, col_seq=FMM5_cMAP,
              breaks=colour_breaks)
param_summary('fmm4', coef_means=FMM4_mean_coefs, col_seq=FMM4_cMAP,
              breaks=colour_breaks)
param_summary('gagnar_diff', coef_means=GAGNAR_diff_mean_coefs, col_seq=GAGNAR_diff_cMAP,
              breaks=colour_breaks)
param_summary('gagnar', coef_means=GAGNAR_mean_coefs, col_seq=GAGNAR_cMAP,
              breaks=colour_breaks)

eco_graph(point_forecast_df[1:26, 3], DDP_cMAP, "DDP_predfit",
          maxval=max(point_forecast_df[1:26, 10]))
eco_graph(point_forecast_df[1:26, 4], FMM5_cMAP, "FMM5_predfit",
          maxval=max(point_forecast_df[1:26, 10]))
eco_graph(point_forecast_df[1:26, 7], GAGNAR_cMAP, "GAGNAR_predfit",
          maxval=max(point_forecast_df[1:26, 10]))
eco_graph(point_forecast_df[1:26, 7], GAGNAR_diff_cMAP, "GAGNAR_diff_predfit",
          maxval=max(point_forecast_df[1:26, 10]))
eco_graph(point_forecast_df[1:26, 10], rep(1, 26), "AdjPNAR_predfit",
          maxval=max(point_forecast_df[1:26, 10]))

coocc_mat_plot(DDP_cc, DDP_cMAP, "coocc_DDP", names=counties_info$County)
coocc_mat_plot(FMM4_cc, FMM4_cMAP, "coocc_FMM4", names=counties_info$County)
coocc_mat_plot(FMM5_cc, FMM5_cMAP, "coocc_FMM5", names=counties_info$County)
coocc_mat_plot(GAGNAR_diff_cc, GAGNAR_diff_cMAP, "coocc_GAGNAR_diff",
               names=counties_info$County)
coocc_mat_plot(GAGNAR_cc, GAGNAR_cMAP, "coocc_GAGNAR", names=counties_info$County)


# compare least-squares partition obtained by different models
modnames <- c("DDP", "FMM5", "FMM4", "GAGNAR_diff", "GAGNAR")
cluster_list <- list(DDP_cMAP, FMM5_cMAP, FMM4_cMAP, GAGNAR_diff_cMAP, GAGNAR_cMAP)
models_ARI <- matrix(NA, 5, 5)
for (i in 1:5){
    for (j in i:5){
        models_ARI[i, j] <- adjustedRandIndex(cluster_list[[i]], cluster_list[[j]])
    }
}
rownames(models_ARI) <- modnames
colnames(models_ARI) <- modnames
xtable(models_ARI)

                                                   
#### REVIEWER COMMENTS ####
# checking network sensitivity 
railway_adj <- igraph::as_adjacency_matrix(covid_net_train_igraph)
railway_adj <- as.matrix(railway_adj)
or <- order(county_index_train$CountyName)
railway_adj <- railway_adj[or, or]
row.names(railway_adj) <- 1:26
colnames(railway_adj) <- 1:26
railway_igraph <- igraph::graph_from_adjacency_matrix(railway_adj, 
                                                      mode='undirected')
V(railway_igraph)$name <- counties_info$County
railway_distances <- as.matrix(distances(railway_igraph))
railway_aff_expdist <- affinity(railway_distances, function(x) exp(-x))
png("railway_network.png", height=700)
plot(railway_igraph, vertex.label.cex=1.8, 
     layout=cbind(counties_info$Long, counties_info$Lat), 
     asp=1.5, vertex.size=7)

ts_railway <- X_from_ts(covid_ts, railway_adj, pop_adjust=counties_info$Pop,
                       pop_adjust_intercept=TRUE)
X_railway <- ts_railway[1:23, , ]

complete_adj <- matrix(1, nrow=26, ncol=26) - diag(26)
complete_aff_expdist <- affinity(complete_adj, function(x) exp(-x))
ts_complete <- X_from_ts(covid_ts, complete_adj, pop_adjust=counties_info$Pop,
                       pop_adjust_intercept=TRUE)
X_complete <- ts_complete[1:23, , ]

ddp_railway <- DDP(y, X_railway, railway_aff_expdist, max_iters=100000, 
                   burn_in=5000)
ddp_complete <- DDP(y, X_complete, complete_aff_expdist, max_iters=100000, 
                    burn_in=5000)

railway_modfit <- model_fit(y, mcmc_out=ddp_railway, crp=TRUE, 
                            X_pred=X_railway, thin=10)
railway_yhat <- one_step_forecast(ts_railway[24, , ], ddp_railway, crp=TRUE)
railway_MASE <- mean(scaled_error(y_ts, railway_yhat$y_hat))
complete_modfit <- model_fit(y, mcmc_out=ddp_complete, crp=TRUE, 
                                X_pred=X_complete, thin=10)
complete_yhat <- one_step_forecast(ts_complete[24, , ], ddp_complete, crp=TRUE)
complete_MASE <- mean(scaled_error(y_ts, complete_yhat$y_hat))
railway_predfit <- model_fit(y_true, mcmc_out=ddp_railway, crp=TRUE, 
                            X_pred=ts_railway[24, , ], thin=10)
complete_predfit <- model_fit(y_true, mcmc_out=ddp_complete, crp=TRUE, 
                             X_pred=ts_complete[24, , ], thin=10)

# Bayesian non-mixture model.                                 
FMM1 <- FMM(y, X, max_iters=100000, burn_in=5000, alpha=c(1))
FMM1_modfit <- single_cluster_eval(y, X, sampled_params=FMM1)
FMM1_predfit <- single_cluster_eval(y_true, X_pred, sampled_params=FMM1)
FMM1_MASE <- mean(scaled_error(y_ts, apply(FMM1 %*% t(X_pred), 2, mean)))
