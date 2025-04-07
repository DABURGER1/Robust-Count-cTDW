library(COUNT)
library(runjags)
library(coda)
library(crch)
library(loo)
library(ggplot2)
library(dplyr)
library(scales)
library(DHARMa)
library(xtable)

NCores <- 4

set.seed(12345)

setwd("C:/Users/a239866/OneDrive - Syneos Health/Divan @ Syneos Health - Linked Files/Research/Robust Count")

########################
### Data preparation ###
########################

data(azpro)

azpro <- azpro %>%
  mutate(ProcedureLabel = factor(procedure, levels = c(0, 1),
                                 labels = c("PTCA", "CABG")),
         SexLabel = factor(sex, levels = c(0, 1),
                           labels = c("Female", "Male")),
         AdmitLabel = factor(admit, levels = c(0, 1),
                             labels = c("Elective", "Urgent/Emerg"))
  )

azpro <- azpro %>%
  mutate(Facet = interaction(ProcedureLabel, AdmitLabel, SexLabel, sep = ", "))

LowerBoundary <- floor(min(azpro$los))
UpperBoundary <- ceiling(max(azpro$los)) + 1

facet_levels <- levels(azpro$Facet)

for(facet_name in facet_levels) {
  plot_data <- filter(azpro, Facet == facet_name)
  
  p <- ggplot(plot_data, aes(x = los)) +
    geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))),
                   binwidth = 1,
                   boundary = LowerBoundary,
                   closed = "left",
                   fill = "steelblue1",
                   color = "black") +
    labs(x = "Length of Hospital Stay", y = "Percentage") +
    scale_x_continuous(breaks = seq(LowerBoundary, UpperBoundary, by = 10)) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic(base_size = 15) +
    theme(axis.line = element_line(linewidth = 1, color = "black"),
          panel.grid = element_blank(),
          strip.background = element_blank())
  
  print(p)
  
  safe_facet <- gsub("[^A-Za-z0-9]", "", facet_name)
  
  ggsave(filename = paste0("Manuscript/Output/LOS_Histogram_", safe_facet, 
                           ".pdf"),
         plot = p,
         device = "pdf",
         width = 8,
         height = 6)
}

X_mat <- model.matrix(~ factor(procedure)*factor(admit)*factor(sex), 
                      data = azpro)
N <- nrow(azpro)
Y <- azpro$los
ones_vec <- rep(1, N)

JAGSData <- list(
  N = N,
  Y = Y,
  X = X_mat,
  P = ncol(X_mat),
  ones = ones_vec,
  C = 1e6
)

#######################################
### Contaminated truncated DW model ###
#######################################

cTDWCode <- "
  model {
    for (p in 1:P) {
      beta[p] ~ dnorm(0, 0.001)
    }
    alpha ~ dgamma(0.001, 0.001)
    eta ~ dgamma(0.001, 0.001)T(1, )
    delta ~ dbeta(1, 9)
    alpha2 <- alpha*eta
    for (i in 1:N) {
      log_mStarMinus1[i] <- inprod(beta[], X[i, ])
      mStar[i] <- 1 + exp(log_mStarMinus1[i])
      q1[i] <- exp(log(0.5)/(mStar[i]^(1/alpha) - 1))
      q2[i] <- exp(log(0.5)/(mStar[i]^(1/alpha2) - 1))
      DW1_trunc[i] <- (q1[i]^(Y[i]^(1/alpha)) - 
                      q1[i]^((Y[i] + 1)^(1/alpha)))/q1[i]
      DW2_trunc[i] <- (q2[i]^(Y[i]^(1/alpha2)) - 
                      q2[i]^((Y[i] + 1)^(1/alpha2)))/q2[i]
      pmf_mix[i] <- delta*DW1_trunc[i] + (1 - delta)*DW2_trunc[i]
      p[i] <- pmf_mix[i]/C
      ones[i] ~ dbern(p[i])
    }
  }
"

cTDWInits <- replicate(
  NCores,
  list(
    alpha = 2,
    eta = 2,
    delta = 0.2,
    beta = rep(0, ncol(X_mat)),
    .RNG.name = "base::Mersenne-Twister",
    .RNG.seed = sample.int(n = 1e5, size = 1)
  ),
  simplify = FALSE
)

cTDWParams <- c("beta", "alpha", "eta", "delta")

cTDWRun <- FALSE
cTDWFile <- "Manuscript/Output/LOS_cTDWModelFit.rds"

if (file.exists(cTDWFile) && !cTDWRun) {
  cat("Loading existing cTDW model from:\n", cTDWFile, "\n")
  cTDWFit <- readRDS(cTDWFile)
} else {
  cat("Running cTDW model...\n")
  
  cTDWFit <- run.jags(
    model = cTDWCode,
    data = JAGSData,
    inits = cTDWInits,
    monitor = cTDWParams,
    n.chains = NCores,
    adapt = 2000,
    burnin = 4000,
    sample = 5000,
    thin = 5,
    method = "parallel",
    modules = "glm"
  )
  
  saveRDS(cTDWFit, cTDWFile)
  cat("Model run complete. Saved to file:\n", cTDWFile, "\n")
}

cTDWSummary <- summary(cTDWFit, vars = cTDWParams)
cTDWSummary

cTDWMCMCList <- as.mcmc(cTDWFit)
summary(cTDWMCMCList)

cTDWPost <- as.matrix(cTDWMCMCList, chains = TRUE, combine = TRUE)
dim(cTDWPost)
head(cTDWPost)

##############################
### LOO and K-L divergence ###
##############################

cTDWDiagnostic <- function(post, data_list, cores = 1) {
  N = data_list$N
  Y = data_list$Y
  X = data_list$X
  M = nrow(post)
  
  beta_mat = post[, grep("^beta\\[", colnames(post)), drop = FALSE]
  alpha_vec = post[, "alpha"]
  eta_vec = post[, "eta"]
  delta_vec = post[, "delta"]
  
  fixed_part = beta_mat%*%t(X)
  mStar_mat = 1 + exp(fixed_part)
  
  alpha_mat = matrix(alpha_vec, nrow = M, ncol = N)
  eta_mat = matrix(eta_vec, nrow = M, ncol = N)
  delta_mat = matrix(delta_vec, nrow = M, ncol = N)
  
  alpha2_mat = alpha_mat*eta_mat
  Y_mat = matrix(Y, nrow = M, ncol = N, byrow = TRUE)
  
  q1_mat = exp(log(0.5)/(mStar_mat^(1/alpha_mat) - 1))
  q2_mat = exp(log(0.5)/(mStar_mat^(1/alpha2_mat) - 1))
  
  exponent1 = Y_mat^(1/alpha_mat)
  exponent1p = (Y_mat + 1)^(1/alpha_mat)
  pmf1_mat = (q1_mat^exponent1 - q1_mat^exponent1p)/q1_mat
  
  exponent2 = Y_mat^(1/alpha2_mat)
  exponent2p = (Y_mat + 1)^(1/alpha2_mat)
  pmf2_mat = (q2_mat^exponent2 - q2_mat^exponent2p)/q2_mat
  
  pmf_mix = delta_mat*pmf1_mat + (1 - delta_mat)*pmf2_mat
  log_lik_matrix = log(pmf_mix)
  
  loo_result = loo(log_lik_matrix, cores = cores)
  
  like_matrix = exp(log_lik_matrix)
  like_matrix_t = t(like_matrix)
  
  inv_like_sum = rowSums(1/like_matrix_t)
  log_like_sum = rowSums(t(log_lik_matrix))
  kl_divergences = log(inv_like_sum/M) + (log_like_sum/M)
  
  cpo_inv = rowMeans(1/like_matrix_t)
  lpml = -sum(log(cpo_inv))
  
  return(list(
    log_lik_matrix = log_lik_matrix,
    LOO = loo_result,
    KL_divergences = kl_divergences,
    LPML = lpml
  ))
}

cTDWChecks <- cTDWDiagnostic(
  post = cTDWPost,
  data_list = JAGSData,
  cores = 4
)

cTDWChecks$LOO
cTDWChecks$LPML

KLThreshold <- function(kl) {
  0.5*(1 + sqrt(1 - exp(-2*kl))) >= 0.8
}

cTDWKLVals <- cTDWChecks$KL_divergences

cTDWKLPlotDat <- data.frame(
  Observation = seq_along(cTDWKLVals),
  KL = cTDWKLVals,
  KLCapped = pmin(cTDWKLVals, 1),
  Influential = ifelse(KLThreshold(cTDWKLVals), "Influential", "Not Influential")
)

cTDWKLPlot <- ggplot(cTDWKLPlotDat, aes(x = Observation, y = KLCapped, color = Influential)) +
  geom_segment(aes(xend = Observation, yend = 0), linewidth = 0.5) +
  scale_color_manual(values = c("Influential" = "red", "Not Influential" = "blue")) +
  labs(x = "Observation Index", y = "K-L Divergence") +
  scale_x_continuous(breaks = seq(0, max(cTDWKLPlotDat$Observation), by = 350)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 17),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

cTDWKLPlot

ggsave("Manuscript/Output/LOS_cTDW_KL_Divergence.pdf",
       plot = cTDWKLPlot, device = "pdf", width = 8, height = 6)

###############
### Medians ###
###############

cTDWMedianOne <- function(mStar, alpha, eta, delta, maxY = 1000) {
  alpha2 = alpha*eta
  
  q1 = exp(log(0.5)/(mStar^(1/alpha) - 1))
  q2 = exp(log(0.5)/(mStar^(1/alpha2) - 1))
  
  cdfVal = 0
  for (y in seq_len(maxY)) {
    pmf1 = (q1^(y^(1/alpha)) - q1^((y + 1)^(1/alpha)))/q1
    pmf2 = (q2^(y^(1/alpha2)) - q2^((y + 1)^(1/alpha2)))/q2
    pmf_mix = delta*pmf1 + (1 - delta)*pmf2
    
    cdfVal = cdfVal + pmf_mix
    if (cdfVal >= 0.5) return(y)
  }
  
  return(NA)
}

GetcTDWMedians <- function(post, X_mat, maxY = 1000) {
  beta_cols = grep("^beta\\[", colnames(post))
  alpha_col = "alpha"
  eta_col = "eta"
  delta_col = "delta"
  
  beta_mat = post[, beta_cols, drop = FALSE]
  alpha_vec = post[, alpha_col]
  eta_vec = post[, eta_col]
  delta_vec = post[, delta_col]
  
  combos = unique(X_mat)
  K = nrow(combos)
  M = nrow(beta_mat)
  
  out_list = vector("list", K)
  
  for (k in seq_len(K)) {
    x_row = combos[k, ]
    linpred_vec = beta_mat%*%x_row
    mStar_vec = 1 + exp(linpred_vec)
    
    median_draws = sapply(seq_len(M), function(m) {
      cTDWMedianOne(
        mStar = mStar_vec[m],
        alpha = alpha_vec[m],
        eta = eta_vec[m],
        delta = delta_vec[m],
        maxY = maxY
      )
    })
    
    median_est = median(median_draws, na.rm = TRUE)
    ci_bounds = HPDinterval(as.mcmc(median_draws), prob = 0.95)
    
    fitted_mStar_median = median(mStar_vec, na.rm = TRUE)
    fitted_mStar_hpd = HPDinterval(as.mcmc(mStar_vec), prob = 0.95)
    
    out_list[[k]] = data.frame(
      Model = "cTDW",
      combo_index = k,
      median_post_median = median_est,
      median_CI_lower = ci_bounds[1],
      median_CI_upper = ci_bounds[2],
      fitted_mStar_median = fitted_mStar_median,
      fitted_mStar_lower = fitted_mStar_hpd[1],
      fitted_mStar_upper = fitted_mStar_hpd[2]
    )
  }
  
  results_data = do.call(rbind, out_list)
  return(results_data)
}

cTDWMediansResults = GetcTDWMedians(cTDWPost, X_mat, maxY = 2000)
print(cTDWMediansResults)

#######################
### Residual checks ###
#######################

QValTDW <- function(mStar_vec, alpha) {
  exp(log(0.5)/(mStar_vec^(1/alpha) - 1))
}

QTDW <- function(u, mStar_vec, alpha) {
  qVal = QValTDW(mStar_vec, alpha)
  inside = 1 + log(1 - u)/log(qVal)
  val = inside^alpha - 1
  y = ceiling(val)
  y[y < 1] = 1L
  return(y)
}

rTDW <- function(mStar_vec, alpha) {
  N = length(mStar_vec)
  u = runif(N)
  QTDW(u, mStar_vec, alpha)
}

rCTDW <- function(mStar_vec, alpha, alpha2, delta) {
  N = length(mStar_vec)
  pick1 = (runif(N) < delta)
  out = integer(N)
  
  if (any(pick1)) out[pick1] = rTDW(mStar_vec[pick1], alpha)
  if (any(!pick1)) out[!pick1] = rTDW(mStar_vec[!pick1], alpha2)
  
  return(out)
}

SampleTDW <- function(post, X, nsim = 100) {
  beta_cols = grep("^beta\\[", colnames(post))
  alpha_col = "alpha"
  
  M = nrow(post)
  N = nrow(X)
  
  sim_matrix = matrix(NA_integer_, nrow = N, ncol = nsim)
  
  for (rep_i in seq_len(nsim)) {
    draw_idx = sample.int(M, size = 1)
    
    alpha_draw = post[draw_idx, alpha_col]
    beta_draw = post[draw_idx, beta_cols, drop = FALSE]
    
    linpred = as.numeric(beta_draw%*%t(X))
    mStar_vec = 1 + exp(linpred)
    
    sim_counts = rTDW(mStar_vec, alpha_draw)
    
    sim_matrix[, rep_i] = sim_counts
  }
  
  return(sim_matrix)
}

SamplecTDW <- function(post, X, nsim = 100) {
  beta_cols = grep("^beta\\[", colnames(post))
  alpha_col = "alpha"
  eta_col = "eta"
  delta_col = "delta"
  
  M = nrow(post)
  N = nrow(X)
  sim_matrix = matrix(NA_integer_, nrow = N, ncol = nsim)
  
  for (rep_i in seq_len(nsim)) {
    draw_idx = sample.int(M, size = 1)
    
    alpha_draw = post[draw_idx, alpha_col]
    eta_draw = post[draw_idx, eta_col]
    delta_draw = post[draw_idx, delta_col]
    beta_draw = post[draw_idx, beta_cols, drop = FALSE]
    
    alpha2_draw = alpha_draw*eta_draw
    linpred = as.numeric(beta_draw%*%t(X))
    mStar_vec = 1 + exp(linpred)
    
    sim_counts = rCTDW(mStar_vec, alpha_draw, alpha2_draw, delta_draw)
    sim_matrix[, rep_i] = sim_counts
  }
  
  return(sim_matrix)
}

simcTDW = SamplecTDW(cTDWPost, X_mat, nsim = 500)

cTDWDHARMa = createDHARMa(
  simulatedResponse = simcTDW,
  observedResponse = Y,
  integerResponse = TRUE
)

plot(cTDWDHARMa)
testDispersion(cTDWDHARMa)
testOutliers(cTDWDHARMa, type = "bootstrap")

cTDWError <- residuals(cTDWDHARMa)

cTDWQQData <- data.frame(
  Expected = seq(0, 1, length.out = length(cTDWError)),
  Observed = sort(cTDWError)
)

pdf(file = "Manuscript/Output/LOS_cTDW_UnifQQ.pdf", width = 6, height = 5)

par(mfrow = c(1,1), oma = c(0.5, 0.5, 0.5, 0.5), mar = c(4, 4, 0.5, 0.5))  

qqplot(cTDWQQData$Expected, cTDWQQData$Observed, 
       main = NULL, xlab = "", ylab = "", pch = 16, col = "darkblue", cex = 0.7, bty = "n")

mtext("Expected Residual", side = 1, line = 3, cex = 1.5)
mtext("Observed Residual", side = 2, line = 3, cex = 1.5)
abline(0, 1, col = "red", lwd = 2)

dev.off()

##########################
### Truncated DW model ###
##########################

TDWCode <- "
  model {
    for (p in 1:P) {
      beta[p] ~ dnorm(0, 0.001)
    }
    alpha ~ dgamma(0.001, 0.001)
    for (i in 1:N) {
      log_mStarMinus1[i] <- inprod(beta[], X[i, ])
      mStar[i] <- 1 + exp(log_mStarMinus1[i])
      q[i] <- exp(log(0.5)/(mStar[i]^(1/alpha) - 1))
      DW_trunc[i] <- (q[i]^(Y[i]^(1/alpha)) - q[i]^((Y[i] + 1)^(1/alpha)))/q[i]
      p[i] <- DW_trunc[i]/C
      ones[i] ~ dbern(p[i])
    }
  }
"

TDWInits <- replicate(
  NCores,
  list(
    alpha = 2,
    beta = rep(0, ncol(X_mat)),
    .RNG.name = "base::Mersenne-Twister",
    .RNG.seed = sample.int(n = 1e5, 1)
  ),
  simplify = FALSE
)

TDWParams <- c("beta", "alpha")

TDWRun <- FALSE
TDWFile <- "Manuscript/Output/LOS_TDWModelFit.rds"

if (file.exists(TDWFile) && !TDWRun) {
  cat("Loading existing TDW model from:\n", TDWFile, "\n")
  TDWFit <- readRDS(TDWFile)
} else {
  cat("Running TDW model...\n")
  
  TDWFit <- run.jags(
    model = TDWCode,
    data = JAGSData,
    inits = TDWInits,
    monitor = TDWParams,
    n.chains = NCores,
    adapt = 2000,
    burnin = 4000,
    sample = 5000,
    thin = 5,
    method = "parallel"
  )
  
  saveRDS(TDWFit, TDWFile)
  cat("Model run complete. Saved to file:\n", TDWFile, "\n")
}

TDWSummary <- summary(TDWFit, vars = TDWParams)
TDWSummary

TDWMCMCList <- as.mcmc(TDWFit)
summary(TDWMCMCList)

TDWPost <- as.matrix(TDWMCMCList, chains = TRUE, combine = TRUE)
dim(TDWPost)
head(TDWPost)

##############################
### LOO and K-L divergence ###
##############################

TDWDiagnostic <- function(post, data_list, cores = 1) {
  N = data_list$N
  Y = data_list$Y
  X = data_list$X
  M = nrow(post)
  
  beta_mat = post[, grep("^beta\\[", colnames(post)), drop = FALSE]
  alpha_vec = post[, "alpha"]
  
  fixed_part = beta_mat%*%t(X)
  mStar_mat = 1 + exp(fixed_part)
  
  alpha_mat = matrix(alpha_vec, nrow = M, ncol = N)
  Y_mat = matrix(Y, nrow = M, ncol = N, byrow = TRUE)
  
  q_mat = exp(log(0.5)/(mStar_mat^(1/alpha_mat) - 1))
  
  exponent = Y_mat^(1/alpha_mat)
  exponent_plus = (Y_mat + 1)^(1/alpha_mat)
  pmf_mat = (q_mat^exponent - q_mat^exponent_plus)/q_mat
  
  log_lik_matrix = log(pmf_mat)
  
  loo_result = loo(log_lik_matrix, cores = cores)
  
  like_matrix = exp(log_lik_matrix)
  like_matrix_t = t(like_matrix)
  
  inv_like_sum = rowSums(1/like_matrix_t)
  log_like_sum = rowSums(t(log_lik_matrix))
  kl_divergences = log(inv_like_sum/M) + (log_like_sum/M)
  
  cpo_inv = rowMeans(1/like_matrix_t)
  lpml = -sum(log(cpo_inv))
  
  return(list(
    log_lik_matrix = log_lik_matrix,
    LOO = loo_result,
    KL_divergences = kl_divergences,
    LPML = lpml
  ))
}

TDWChecks <- TDWDiagnostic(
  post = TDWPost,
  data_list = JAGSData,
  cores = 4
)

TDWChecks$LOO
TDWChecks$LPML

TDWKLVals <- TDWChecks$KL_divergences

TDWKLPlotDat <- data.frame(
  Observation = seq_along(TDWKLVals),
  KL = TDWKLVals,
  KLCapped = pmin(TDWKLVals, 1),
  Influential = ifelse(KLThreshold(TDWKLVals), "Influential", "Not Influential")
)

TDWKLPlot <- ggplot(TDWKLPlotDat, aes(x = Observation, y = KLCapped, color = Influential)) +
  geom_segment(aes(xend = Observation, yend = 0), linewidth = 0.5) +
  scale_color_manual(values = c("Influential" = "red", "Not Influential" = "blue")) +
  labs(x = "Observation Index", y = "K-L Divergence") +
  scale_x_continuous(breaks = seq(0, max(TDWKLPlotDat$Observation), by = 350)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 17),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

TDWKLPlot

ggsave("Manuscript/Output/LOS_TDW_KL_Divergence.pdf",
       plot = TDWKLPlot, device = "pdf",
       width = 8, height = 6)

###############
### Medians ###
###############

TDWMedianOne <- function(mStar, alpha, maxY = 1000) {
  q = exp(log(0.5)/(mStar^(1/alpha) - 1))
  
  cdfVal = 0
  for (y in seq_len(maxY)) {
    pmf_y = (q^(y^(1/alpha)) - q^((y + 1)^(1/alpha)))/q
    cdfVal = cdfVal + pmf_y
    if (cdfVal >= 0.5) return(y)
  }
  
  return(NA)
}

GetTDWMedians <- function(post, X_mat, maxY = 2000) {
  beta_cols = grep("^beta\\[", colnames(post))
  alpha_col = "alpha"
  
  beta_mat = post[, beta_cols, drop = FALSE]
  alpha_vec = post[, alpha_col]
  
  combos = unique(X_mat)
  K = nrow(combos)
  M = nrow(beta_mat)
  
  out_list = vector("list", K)
  
  for (k in seq_len(K)) {
    x_row = combos[k, ]
    linpred_vec = beta_mat%*%x_row
    mStar_vec = 1 + exp(linpred_vec)
    
    median_draws = sapply(seq_len(M), function(m) {
      TDWMedianOne(mStar = mStar_vec[m], alpha = alpha_vec[m], maxY = maxY)
    })
    
    median_est = median(median_draws, na.rm = TRUE)
    ci_bounds = HPDinterval(as.mcmc(median_draws), prob = 0.95)
    
    fitted_mStar_median = median(mStar_vec, na.rm = TRUE)
    fitted_mStar_hpd = HPDinterval(as.mcmc(mStar_vec), prob = 0.95)
    
    out_list[[k]] = data.frame(
      Model = "TDW",
      combo_index = k,
      median_post_median = median_est,
      median_CI_lower = ci_bounds[1],
      median_CI_upper = ci_bounds[2],
      fitted_mStar_median = fitted_mStar_median,
      fitted_mStar_lower = fitted_mStar_hpd[1],
      fitted_mStar_upper = fitted_mStar_hpd[2]
    )
  }
  
  results_data = do.call(rbind, out_list)
  return(results_data)
}

TDWMediansResults = GetTDWMedians(TDWPost, X_mat, maxY = 2000)
print(TDWMediansResults)

#######################
### Residual checks ###
#######################

simTDW = SampleTDW(TDWPost, X_mat, nsim = 500)

TDWDHARMa = createDHARMa(
  simulatedResponse = simTDW,
  observedResponse = Y,
  integerResponse = TRUE
)

plot(TDWDHARMa)
testDispersion(TDWDHARMa)
testOutliers(TDWDHARMa, type = "bootstrap")

TDWError <- residuals(TDWDHARMa)

TDWQQData <- data.frame(
  Expected = seq(0, 1, length.out = length(TDWError)),
  Observed = sort(TDWError)
)

pdf(file = "Manuscript/Output/LOS_TDW_UnifQQ.pdf", width = 6, height = 5)

par(mfrow = c(1,1), oma = c(0.5, 0.5, 0.5, 0.5), mar = c(4, 4, 0.5, 0.5))  

qqplot(TDWQQData$Expected, TDWQQData$Observed, 
       main = NULL, xlab = "", ylab = "", pch = 16, col = "darkblue", cex = 0.7, bty = "n")

mtext("Expected Residual", side = 1, line = 3, cex = 1.5)
mtext("Observed Residual", side = 2, line = 3, cex = 1.5)
abline(0, 1, col = "red", lwd = 2)

dev.off()

#################################
### Median estimates combined ###
#################################

FacetLabels <- azpro %>%
  select(procedure, admit, sex) %>%
  distinct() %>%
  mutate(
    procedure_label = ifelse(procedure == 0, "PTCA", "CABG"),
    admit_label = ifelse(admit == 0, "Elective", "Urgent/Emerg"),
    sex_label = ifelse(sex == 0, "Female", "Male"),
    FacetLabel = paste0(procedure_label, 
                        ", ", admit_label, 
                        ", ", sex_label)
  ) %>%
  select(procedure, admit, sex, FacetLabel)

X_mat_df <- as.data.frame(X_mat) %>%
  mutate(
    procedure = ifelse(`factor(procedure)1` == 1, 1, 0),
    admit = ifelse(`factor(admit)1` == 1, 1, 0),
    sex = ifelse(`factor(sex)1` == 1, 1, 0)
  ) %>%
  select(procedure, admit, sex) %>%
  distinct()

FacetMapping <- X_mat_df %>%
  left_join(FacetLabels, by = c("procedure", "admit", "sex")) %>%
  mutate(combo_index = row_number()) %>%
  select(combo_index, FacetLabel)

AllResults <- bind_rows(TDWMediansResults, cTDWMediansResults) %>%
  left_join(FacetMapping, by = "combo_index") %>%
  mutate(FacetLabel = factor(FacetLabel))

AllResults <- AllResults %>%
  mutate(Model = factor(Model, levels = c("TDW", "cTDW")))

pdf("Manuscript/Output/LOS_Median_Bands.pdf", width = 11, height = 6)

ggplot(AllResults, aes(x = Model, y = fitted_mStar_median, color = as.factor(Model))) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(
    aes(ymin = fitted_mStar_lower, ymax = fitted_mStar_upper),
    width = 0.2,
    position = position_dodge(width = 0.3)
  ) +
  facet_wrap(~ FacetLabel, nrow = 2, scales = "free_y") +
  labs(
    x = "Model",
    y = "Median LOS",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11, angle = 30, hjust = 1),
    legend.position = "bottom",
    legend.key.width = unit(1.2, "cm"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  scale_y_continuous(breaks = seq(0, max(ceiling(AllResults$fitted_mStar_median), na.rm = TRUE), 0.5)) +
  scale_color_manual(
    values = c("#1f77b4", "#ff7f0e"),
    guide = guide_legend(override.aes = list(linetype = "solid"))
  )

dev.off()

TDW_stats <- apply(TDWPost, 2, function(x) {
  c(
    Median = median(x),
    LCI_2.5 = quantile(x, probs = 0.025),
    UCI_97.5 = quantile(x, probs = 0.975)
  )
})

rownames(TDW_stats) <- c("Median", "LCI 2.5%", "UCI 97.5%")
print(TDW_stats)

###########################
### Parameter estimates ###
###########################

cTDW_stats <- apply(cTDWPost, 2, function(x) {
  c(
    Median = median(x),
    LCI_2.5 = quantile(x, probs = 0.025),
    UCI_97.5 = quantile(x, probs = 0.975)
  )
})

rownames(cTDW_stats) <- c("Median", "LCI 2.5%", "UCI 97.5%")
print(cTDW_stats)

TDW_df <- data.frame(
  Parameter = colnames(TDW_stats),
  TDW_median = TDW_stats["Median", ],
  TDW_lci = TDW_stats["LCI 2.5%", ],
  TDW_uci = TDW_stats["UCI 97.5%", ],
  stringsAsFactors = FALSE
)

cTDW_df <- data.frame(
  Parameter = colnames(cTDW_stats),
  cTDW_median = cTDW_stats["Median", ],
  cTDW_lci = cTDW_stats["LCI 2.5%", ],
  cTDW_uci = cTDW_stats["UCI 97.5%", ],
  stringsAsFactors = FALSE
)

merged_df <- merge(TDW_df, cTDW_df, by = "Parameter", all = TRUE)

merged_df$sort_order <- sapply(merged_df$Parameter, function(x) {
  if (grepl("^beta\\[\\d+\\]$", x)) {
    as.numeric(sub("beta\\[(\\d+)\\]", "\\1", x))
  } else if (x == "alpha") {
    1000
  } else if (x == "eta") {
    1001
  } else if (x == "delta") {
    1002
  } else {
    9999
  }
})

merged_df <- merged_df[order(merged_df$sort_order), ]
merged_df$sort_order <- NULL

param2latex <- function(x) {
  if (x == "alpha") return("$\\alpha$")
  if (x == "delta") return("$\\delta$")
  if (x == "eta") return("$\\eta$")
  if (grepl("^beta\\[\\d+\\]$", x)) {
    idx <- as.numeric(sub("beta\\[(\\d+)\\]", "\\1", x)) - 1
    return(paste0("$\\beta_{", idx, "}$"))
  }
  return(x)
}

final_df <- data.frame(
  Parameter = sapply(merged_df$Parameter, param2latex),
  TDW_est = sprintf("%.3f", merged_df$TDW_median),
  TDW_left = sprintf("[%.3f;", merged_df$TDW_lci),
  TDW_right = sprintf("%.3f]", merged_df$TDW_uci),
  blank = "",
  cTDW_est = sprintf("%.3f", merged_df$cTDW_median),
  cTDW_left = sprintf("[%.3f;", merged_df$cTDW_lci),
  cTDW_right = sprintf("%.3f]", merged_df$cTDW_uci),
  stringsAsFactors = FALSE
)

xtab <- xtable(
  final_df,
  align = rep("l", ncol(final_df) + 1)
)

sink("Manuscript/Output/LOS_Posterior_Summaries.tex")

print(
  xtab,
  include.rownames = FALSE,
  sanitize.text.function = identity,
  caption.placement = "top",
  hline.after = c(-1, 0, nrow(final_df))
)

sink()