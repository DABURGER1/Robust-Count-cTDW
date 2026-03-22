# =========================
# Libraries and settings
# =========================

library(COUNT)
library(rjags)
library(runjags)
library(coda)
library(crch)
library(loo)
library(ggplot2)
library(dplyr)
library(scales)
library(DHARMa)
library(xtable)

setwd("C:/Users/divan.burger/OneDrive - Cytel/Research/Robust Count")

NCores <- 4
OutputDir <- "Manuscript/Output"
RunLabel <- "LOS"
ConstantC <- 1e6
if (!dir.exists(OutputDir)) dir.create(OutputDir, recursive = TRUE)

set.seed(12345)

# =========================
# Common data
# =========================

data(azpro)

azpro <- azpro %>%
  mutate(
    ProcedureLabel = factor(procedure,
                            levels = c(0, 1),
                            labels = c("PTCA", "CABG")
    ),
    SexLabel = factor(sex,
                      levels = c(0, 1),
                      labels = c("Female", "Male")
    ),
    AdmitLabel = factor(admit,
                        levels = c(0, 1),
                        labels = c("Elective", "Urgent/Emerg")
    )
  )

azpro <- azpro %>%
  mutate(Facet = interaction(ProcedureLabel, AdmitLabel, SexLabel, sep = ", "))

LowerBoundary <- floor(min(azpro$los))
UpperBoundary <- ceiling(max(azpro$los)) + 1

facet_levels <- levels(azpro$Facet)

for (facet_name in facet_levels) {
  plot_data <- filter(azpro, Facet == facet_name)
  
  p <- ggplot(plot_data, aes(x = los)) +
    geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))),
                   binwidth = 1,
                   boundary = LowerBoundary,
                   closed = "left",
                   fill = "steelblue1",
                   color = "black"
    ) +
    labs(x = "Length of Hospital Stay", y = "Percentage") +
    scale_x_continuous(breaks = seq(LowerBoundary, UpperBoundary, by = 10)) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic(base_size = 15) +
    theme(
      axis.line = element_line(linewidth = 1, color = "black"),
      panel.grid = element_blank(),
      strip.background = element_blank()
    )
  
  print(p)
  
  safe_facet <- gsub("[^A-Za-z0-9]", "", facet_name)
  
  ggsave(
    filename = paste0(
      "Manuscript/Output/LOS_Histogram_", safe_facet,
      ".pdf"
    ),
    plot = p,
    device = "pdf",
    width = 8,
    height = 6
  )
}

FormulaObj <- ~ factor(procedure)*factor(admit)*factor(sex)
ResponseVar <- "los"
DataFrame <- azpro

XMat <- model.matrix(FormulaObj, data = DataFrame)
YVec <- as.numeric(DataFrame[[ResponseVar]])
NObs <- length(YVec)

JAGSData <- list(
  N = NObs,
  Y = YVec,
  X = XMat,
  P = ncol(XMat),
  ones = rep(1, NObs),
  C = ConstantC
)

X_mat <- XMat
Y <- YVec

# =========================
# Shared utilities
# =========================

KLThreshold <- function(kl) {
  0.5*(1 + sqrt(1 - exp(-2*kl))) >= 0.8
}

# TDW quantile sampler functions
QValTDW <- function(MStarVec, Alpha) {
  exp(log(0.5)/(MStarVec^(1/Alpha) - 1))
}
QTDW <- function(U, MStarVec, Alpha) {
  QVal <- QValTDW(MStarVec, Alpha)
  Inside <- 1 + log(1 - U)/log(QVal)
  Val <- Inside^Alpha - 1
  YSim <- ceiling(Val)
  YSim[YSim < 1] <- 1L
  YSim
}
RTDW <- function(MStarVec, Alpha) {
  U <- runif(length(MStarVec))
  QTDW(U, MStarVec, Alpha)
}
RCTDW <- function(MStarVec, Alpha, Alpha2, Delta) {
  Pick1 <- (runif(length(MStarVec)) < Delta)
  Out <- integer(length(MStarVec))
  if (any(Pick1)) Out[Pick1] <- RTDW(MStarVec[Pick1], Alpha)
  if (any(!Pick1)) Out[!Pick1] <- RTDW(MStarVec[!Pick1], Alpha2)
  Out
}

SampleTDW <- function(Post, X, nsim = 100) {
  BetaCols <- grep("^beta\\[", colnames(Post))
  M <- nrow(Post)
  N <- nrow(X)
  Sim <- matrix(NA_integer_, nrow = N, ncol = nsim)
  for (r in seq_len(nsim)) {
    i <- sample.int(M, 1)
    Alpha <- Post[i, "alpha"]
    Beta <- Post[i, BetaCols, drop = FALSE]
    LinPred <- as.numeric(Beta%*%t(X))
    MStar <- 1 + exp(LinPred)
    Sim[, r] <- RTDW(MStar, Alpha)
  }
  Sim
}

SampleCTDW <- function(Post, X, nsim = 100) {
  BetaCols <- grep("^beta\\[", colnames(Post))
  M <- nrow(Post)
  N <- nrow(X)
  Sim <- matrix(NA_integer_, nrow = N, ncol = nsim)
  for (r in seq_len(nsim)) {
    i <- sample.int(M, 1)
    Alpha <- Post[i, "alpha"]
    Eta <- Post[i, "eta"]
    Delta <- Post[i, "delta"]
    Alpha2 <- Alpha*Eta
    Beta <- Post[i, BetaCols, drop = FALSE]
    LinPred <- as.numeric(Beta%*%t(X))
    MStar <- 1 + exp(LinPred)
    Sim[, r] <- RCTDW(MStar, Alpha, Alpha2, Delta)
  }
  Sim
}

rTNB <- function(MuVec, Alpha) {
  Prob <- Alpha/(Alpha + MuVec)
  F0 <- dnbinom(0, size = Alpha, prob = Prob)
  U <- runif(length(MuVec))
  Target <- U*(1 - F0) + F0
  as.integer(qnbinom(Target, size = Alpha, prob = Prob))
}
SampleTNB <- function(Post, X, nsim = 100) {
  BetaCols <- grep("^beta\\[", colnames(Post))
  M <- nrow(Post)
  N <- nrow(X)
  Sim <- matrix(NA_integer_, nrow = N, ncol = nsim)
  for (r in seq_len(nsim)) {
    i <- sample.int(M, 1)
    Alpha <- Post[i, "alpha"]
    Beta <- Post[i, BetaCols, drop = FALSE]
    LinPred <- as.numeric(Beta%*%t(X))
    Mu <- exp(LinPred)
    Sim[, r] <- rTNB(Mu, Alpha)
  }
  Sim
}

# Medians
TDWMedianOne <- function(MStar, Alpha, MaxY = 2000) {
  q <- exp(log(0.5)/(MStar^(1/Alpha) - 1))
  CDF <- 0
  for (y in seq_len(MaxY)) {
    pmf <- (q^(y^(1/Alpha)) - q^((y + 1)^(1/Alpha)))/q
    CDF <- CDF + pmf
    if (CDF >= 0.5) {
      return(y)
    }
  }
  NA_integer_
}
GetTDWMedians <- function(Post, X, MaxY = 2000) {
  BetaCols <- grep("^beta\\[", colnames(Post))
  BetaMat <- Post[, BetaCols, drop = FALSE]
  AlphaVec <- Post[, "alpha"]
  Combos <- unique(X)
  K <- nrow(Combos)
  M <- nrow(BetaMat)
  Out <- vector("list", K)
  for (k in seq_len(K)) {
    XRow <- Combos[k, ]
    Lin <- BetaMat%*%XRow
    MStar <- 1 + exp(Lin)
    Draws <- sapply(seq_len(M), function(m) TDWMedianOne(MStar[m], AlphaVec[m], MaxY))
    Med <- median(Draws, na.rm = TRUE)
    CI <- HPDinterval(as.mcmc(Draws), prob = 0.95)
    MMed <- median(MStar, na.rm = TRUE)
    MHPD <- HPDinterval(as.mcmc(MStar), prob = 0.95)
    Out[[k]] <- data.frame(
      Model = "TDW",
      ComboIndex = k,
      MedianPostMedian = Med,
      MedianCILower = CI[1],
      MedianCIUpper = CI[2],
      FittedMStarMedian = MMed,
      FittedMStarLower = MHPD[1],
      FittedMStarUpper = MHPD[2]
    )
  }
  do.call(rbind, Out)
}

CTDWMedianOne <- function(MStar, Alpha, Eta, Delta, MaxY = 2000) {
  Alpha2 <- Alpha*Eta
  q1 <- exp(log(0.5)/(MStar^(1/Alpha) - 1))
  q2 <- exp(log(0.5)/(MStar^(1/Alpha2) - 1))
  CDF <- 0
  for (y in seq_len(MaxY)) {
    pmf1 <- (q1^(y^(1/Alpha)) - q1^((y + 1)^(1/Alpha)))/q1
    pmf2 <- (q2^(y^(1/Alpha2)) - q2^((y + 1)^(1/Alpha2)))/q2
    CDF <- CDF + Delta*pmf1 + (1 - Delta)*pmf2
    if (CDF >= 0.5) {
      return(y)
    }
  }
  NA_integer_
}
GetCTDWMedians <- function(Post, X, MaxY = 2000) {
  BetaCols <- grep("^beta\\[", colnames(Post))
  BetaMat <- Post[, BetaCols, drop = FALSE]
  AlphaVec <- Post[, "alpha"]
  EtaVec <- Post[, "eta"]
  DeltaVec <- Post[, "delta"]
  Combos <- unique(X)
  K <- nrow(Combos)
  M <- nrow(BetaMat)
  Out <- vector("list", K)
  for (k in seq_len(K)) {
    XRow <- Combos[k, ]
    Lin <- BetaMat%*%XRow
    MStar <- 1 + exp(Lin)
    Draws <- sapply(seq_len(M), function(m) {
      CTDWMedianOne(MStar[m], AlphaVec[m], EtaVec[m], DeltaVec[m], MaxY)
    })
    Med <- median(Draws, na.rm = TRUE)
    CI <- HPDinterval(as.mcmc(Draws), prob = 0.95)
    MMed <- median(MStar, na.rm = TRUE)
    MHPD <- HPDinterval(as.mcmc(MStar), prob = 0.95)
    Out[[k]] <- data.frame(
      Model = "cTDW",
      ComboIndex = k,
      MedianPostMedian = Med,
      MedianCILower = CI[1],
      MedianCIUpper = CI[2],
      FittedMStarMedian = MMed,
      FittedMStarLower = MHPD[1],
      FittedMStarUpper = MHPD[2]
    )
  }
  do.call(rbind, Out)
}

# =========================
# Model runners
# =========================

RunCTDWPrepared <- function(
    DeltaLower, DeltaUpper,
    EtaPrior = c("TruncatedGamma", "Uniform"),
    EtaGammaShape = 0.001,
    EtaGammaRate = 0.001,
    EtaUniformLower = 1,
    EtaUniformUpper = 10,
    Adapt = 2000, Burnin = 4000, Sample = 5000, Thin = 5,
    NSimDHARMa = 500,
    PerDeltaSubdir = FALSE, Overwrite = TRUE, Seed = 12345) {
  set.seed(Seed)
  
  EtaPrior <- match.arg(EtaPrior)
  
  if (DeltaLower < 0 || DeltaUpper > 1 || DeltaLower >= DeltaUpper) {
    stop("Require 0 <= DeltaLower < DeltaUpper <= 1.")
  }
  
  if (EtaPrior == "Uniform") {
    if (EtaUniformLower < 1 || EtaUniformLower >= EtaUniformUpper) {
      stop("For EtaPrior = 'Uniform', require 1 <= EtaUniformLower < EtaUniformUpper.")
    }
  }
  
  # JAGS data
  JagsDataLocal <- JAGSData
  JagsDataLocal$DeltaLower <- DeltaLower
  JagsDataLocal$DeltaUpper <- DeltaUpper
  JagsDataLocal$EtaGammaShape <- EtaGammaShape
  JagsDataLocal$EtaGammaRate <- EtaGammaRate
  JagsDataLocal$EtaUniformLower <- EtaUniformLower
  JagsDataLocal$EtaUniformUpper <- EtaUniformUpper
  
  # Names
  LTag <- gsub("\\.", "p", formatC(DeltaLower, format = "f", digits = 3))
  UTag <- gsub("\\.", "p", formatC(DeltaUpper, format = "f", digits = 3))
  IntervalTag <- paste0("L", LTag, "_U", UTag)
  
  EtaPriorTag <- if (EtaPrior == "TruncatedGamma") {
    "EtaTG"
  } else {
    paste0(
      "EtaU",
      gsub("\\.", "p", formatC(EtaUniformLower, format = "f", digits = 3)),
      "_",
      gsub("\\.", "p", formatC(EtaUniformUpper, format = "f", digits = 3))
    )
  }
  
  DirOut <- if (PerDeltaSubdir) {
    file.path(OutputDir, paste0(RunLabel, "_", IntervalTag, "_", EtaPriorTag))
  } else {
    OutputDir
  }
  if (!dir.exists(DirOut)) dir.create(DirOut, recursive = TRUE)
  
  BaseTag <- paste0(RunLabel, "_cTDW_", IntervalTag, "_", EtaPriorTag)
  
  FitFile <- file.path(DirOut, paste0(BaseTag, "_ModelFit.rds"))
  LooFile <- file.path(DirOut, paste0(BaseTag, "_LOO.txt"))
  KLPlotFile <- file.path(DirOut, paste0(BaseTag, "_KL_Divergence.pdf"))
  KLDataFile <- file.path(DirOut, paste0(BaseTag, "_KL_Data.csv"))
  DHARMaPlotFile <- file.path(DirOut, paste0(BaseTag, "_DHARMa.pdf"))
  DHARMaTxtFile <- file.path(DirOut, paste0(BaseTag, "_DHARMa_Tests.txt"))
  QQPlotFile <- file.path(DirOut, paste0(BaseTag, "_UnifQQ.pdf"))
  MediansFile <- file.path(DirOut, paste0(BaseTag, "_Medians.csv"))
  
  EtaPriorLine <- if (EtaPrior == "TruncatedGamma") {
    "eta ~ dgamma(EtaGammaShape, EtaGammaRate)T(1, )"
  } else {
    "eta ~ dunif(EtaUniformLower, EtaUniformUpper)"
  }
  
  # -------- fit or load --------
  if (file.exists(FitFile) && !Overwrite) {
    message("Loading existing fit: ", FitFile)
    Fit <- readRDS(FitFile)
  } else {
    message(
      "Running cTDW with delta in (", DeltaLower, ", ", DeltaUpper,
      ") and eta prior = ", EtaPrior, " ..."
    )
    
    ModelCode <- paste0("
            model {
                for (p in 1:P) { beta[p] ~ dnorm(0, 0.001) }
                alpha ~ dgamma(0.001, 0.001)
                ", EtaPriorLine, "
                delta ~ dunif(DeltaLower, DeltaUpper)

                alpha2 <- alpha*eta

                for (i in 1:N) {
                    log_mStarMinus1[i] <- inprod(beta[], X[i, ])
                    mStar[i] <- 1 + exp(log_mStarMinus1[i])

                    q1[i] <- exp(log(0.5)/(mStar[i]^(1/alpha) - 1))
                    q2[i] <- exp(log(0.5)/(mStar[i]^(1/alpha2) - 1))

                    DW1_trunc[i] <- (q1[i]^(Y[i]^(1/alpha)) - q1[i]^((Y[i] + 1)^(1/alpha)))/q1[i]
                    DW2_trunc[i] <- (q2[i]^(Y[i]^(1/alpha2)) - q2[i]^((Y[i] + 1)^(1/alpha2)))/q2[i]

                    pmf_mix[i] <- delta*DW1_trunc[i] + (1 - delta)*DW2_trunc[i]
                    p[i] <- pmf_mix[i]/C
                    ones[i] ~ dbern(p[i])
                }
            }
        ")
    
    PCols <- ncol(XMat)
    DeltaInit <- 0.5*(DeltaLower + DeltaUpper)
    
    EtaInit <- if (EtaPrior == "Uniform") {
      0.5*(EtaUniformLower + EtaUniformUpper)
    } else {
      2
    }
    
    Inits <- replicate(
      NCores,
      list(
        alpha = 2,
        eta = EtaInit,
        delta = DeltaInit,
        beta = rep(0, PCols),
        .RNG.name = "base::Mersenne-Twister",
        .RNG.seed = sample.int(1e6, 1)
      ),
      simplify = FALSE
    )
    
    Params <- c("beta", "alpha", "eta", "delta")
    
    Fit <- run.jags(
      model = ModelCode,
      data = JagsDataLocal,
      inits = Inits,
      monitor = Params,
      n.chains = NCores,
      adapt = Adapt,
      burnin = Burnin,
      sample = Sample,
      thin = Thin,
      method = "parallel",
      modules = "glm"
    )
    
    saveRDS(Fit, FitFile)
    message("Saved fit to: ", FitFile)
  }
  
  # Posterior
  MCMC <- as.mcmc(Fit)
  Post <- as.matrix(MCMC)
  
  # Diagnostics
  Diagnostic <- function(PostMat, DataList, Cores = 1) {
    N <- DataList$N
    X <- DataList$X
    Y <- DataList$Y
    M <- nrow(PostMat)
    BetaCols <- grep("^beta\\[", colnames(PostMat))
    Beta <- PostMat[, BetaCols, drop = FALSE]
    Alpha <- PostMat[, "alpha"]
    Eta <- PostMat[, "eta"]
    Delta <- PostMat[, "delta"]
    Fixed <- Beta%*%t(X)
    MStar <- 1 + exp(Fixed)
    A <- matrix(Alpha, nrow = M, ncol = N)
    E <- matrix(Eta, nrow = M, ncol = N)
    D <- matrix(Delta, nrow = M, ncol = N)
    A2 <- A*E
    YMat <- matrix(Y, nrow = M, ncol = N, byrow = TRUE)
    Q1 <- exp(log(0.5)/(MStar^(1/A) - 1))
    Q2 <- exp(log(0.5)/(MStar^(1/A2) - 1))
    E1 <- YMat^(1/A)
    E1P <- (YMat + 1)^(1/A)
    E2 <- YMat^(1/A2)
    E2P <- (YMat + 1)^(1/A2)
    PMF1 <- (Q1^E1 - Q1^E1P)/Q1
    PMF2 <- (Q2^E2 - Q2^E2P)/Q2
    PMF <- D*PMF1 + (1 - D)*PMF2
    LogLik <- log(PMF)
    LooRes <- loo(LogLik, cores = Cores)
    Like <- exp(LogLik)
    LikeT <- t(Like)
    InvLikeSum <- rowSums(1/LikeT)
    LogLikeSum <- rowSums(t(LogLik))
    KLD <- log(InvLikeSum/M) + (LogLikeSum/M)
    CPOInv <- rowMeans(1/LikeT)
    LPML <- -sum(log(CPOInv))
    list(LogLik = LogLik, LOO = LooRes, KLDivergences = KLD, LPML = LPML)
  }
  
  Checks <- Diagnostic(Post, JagsDataLocal, Cores = NCores)
  
  capture.output(
    {
      print(Checks$LOO)
      cat("\nLPML: ", Checks$LPML, "\n", sep = "")
    },
    file = LooFile
  )
  
  # KL plot
  KLVals <- Checks$KLDivergences
  KLData <- data.frame(
    Observation = seq_along(KLVals),
    KL = KLVals,
    KLCapped = pmin(KLVals, 1),
    Influential = ifelse(KLThreshold(KLVals), "Influential", "Not Influential")
  )
  write.csv(KLData, KLDataFile, row.names = FALSE)
  
  KLPlot <- ggplot(KLData, aes(x = Observation, y = KLCapped, color = Influential)) +
    geom_segment(aes(xend = Observation, yend = 0), linewidth = 0.5) +
    scale_color_manual(values = c("Influential" = "red", "Not Influential" = "blue")) +
    labs(x = "Observation Index", y = "K-L Divergence") +
    scale_x_continuous(breaks = seq(0, max(KLData$Observation), by = 350)) +
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
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  ggsave(filename = KLPlotFile, plot = KLPlot, device = "pdf", width = 8, height = 6)
  
  # DHARMa residuals
  Sim <- SampleCTDW(Post, XMat, nsim = NSimDHARMa)
  DHObj <- createDHARMa(simulatedResponse = Sim, observedResponse = YVec, integerResponse = TRUE)
  pdf(file = DHARMaPlotFile, width = 7, height = 6)
  plot(DHObj)
  dev.off()
  capture.output(
    {
      print(testDispersion(DHObj))
      print(testOutliers(DHObj, type = "bootstrap"))
    },
    file = DHARMaTxtFile
  )
  
  # QQ plot
  DHError <- residuals(DHObj)
  QQData <- data.frame(Expected = seq(0, 1, length.out = length(DHError)), Observed = sort(DHError))
  pdf(file = QQPlotFile, width = 6, height = 5)
  par(mfrow = c(1, 1), oma = c(0.5, 0.5, 0.5, 0.5), mar = c(4, 4, 0.5, 0.5))
  qqplot(
    QQData$Expected,
    QQData$Observed,
    main = NULL,
    xlab = "",
    ylab = "",
    pch = 16,
    col = "darkblue",
    cex = 0.7,
    bty = "n"
  )
  mtext("Expected Residual", side = 1, line = 3, cex = 1.5)
  mtext("Observed Residual", side = 2, line = 3, cex = 1.5)
  abline(0, 1, col = "red", lwd = 2)
  dev.off()
  
  # Medians
  MediansRes <- GetCTDWMedians(Post, XMat, MaxY = 2000)
  write.csv(MediansRes, MediansFile, row.names = FALSE)
  
  list(
    Fit = Fit,
    Posterior = Post,
    Diagnostics = list(
      LOO = Checks$LOO,
      LPML = Checks$LPML,
      KLDivergences = Checks$KLDivergences
    ),
    Medians = MediansRes,
    Files = list(
      FitRDS = FitFile,
      LOOText = LooFile,
      KLPlotPDF = KLPlotFile,
      KLDataCSV = KLDataFile,
      DHARMaPDF = DHARMaPlotFile,
      DHARMaText = DHARMaTxtFile,
      QQPlotPDF = QQPlotFile,
      MediansCSV = MediansFile
    ),
    Settings = list(
      DeltaLower = DeltaLower,
      DeltaUpper = DeltaUpper,
      EtaPrior = EtaPrior,
      EtaGammaShape = EtaGammaShape,
      EtaGammaRate = EtaGammaRate,
      EtaUniformLower = EtaUniformLower,
      EtaUniformUpper = EtaUniformUpper
    )
  )
}

RunTDWPrepared <- function(
    Adapt = 2000, Burnin = 4000, Sample = 5000, Thin = 5,
    NSimDHARMa = 500, Overwrite = TRUE, Seed = 12345) {
  set.seed(Seed)
  
  BaseTag <- paste0(RunLabel, "_TDW")
  FitFile <- file.path(OutputDir, paste0(BaseTag, "_ModelFit.rds"))
  LooFile <- file.path(OutputDir, paste0(BaseTag, "_LOO.txt"))
  KLPlotFile <- file.path(OutputDir, paste0(BaseTag, "_KL_Divergence.pdf"))
  KLDataFile <- file.path(OutputDir, paste0(BaseTag, "_KL_Data.csv"))
  DHARMaPlotFile <- file.path(OutputDir, paste0(BaseTag, "_DHARMa.pdf"))
  DHARMaTxtFile <- file.path(OutputDir, paste0(BaseTag, "_DHARMa_Tests.txt"))
  QQPlotFile <- file.path(OutputDir, paste0(BaseTag, "_UnifQQ.pdf"))
  MediansFile <- file.path(OutputDir, paste0(BaseTag, "_Medians.csv"))
  
  # -------- fit or load --------
  if (file.exists(FitFile) && !Overwrite) {
    message("Loading existing fit: ", FitFile)
    Fit <- readRDS(FitFile)
  } else {
    message("Running TDW ...")
    ModelCode <- "
      model {
        for (p in 1:P) { beta[p] ~ dnorm(0, 0.001) }
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
    
    Inits <- replicate(
      NCores,
      list(
        alpha = 2,
        beta = rep(0, ncol(XMat)),
        .RNG.name = "base::Mersenne-Twister",
        .RNG.seed = sample.int(1e6, 1)
      ),
      simplify = FALSE
    )
    Params <- c("beta", "alpha")
    
    Fit <- run.jags(
      model = ModelCode,
      data = JAGSData,
      inits = Inits,
      monitor = Params,
      n.chains = NCores,
      adapt = Adapt,
      burnin = Burnin,
      sample = Sample,
      thin = Thin,
      method = "parallel",
      modules = "glm"
    )
    saveRDS(Fit, FitFile)
    message("Saved fit to: ", FitFile)
  }
  
  MCMC <- as.mcmc(Fit)
  Post <- as.matrix(MCMC)
  
  # Diagnostics
  Diagnostic <- function(Post, DataList, Cores = 1) {
    N <- DataList$N
    X <- DataList$X
    Y <- DataList$Y
    M <- nrow(Post)
    Beta <- Post[, grep("^beta\\[", colnames(Post)), drop = FALSE]
    Alpha <- Post[, "alpha"]
    Fixed <- Beta%*%t(X)
    MStar <- 1 + exp(Fixed)
    A <- matrix(Alpha, nrow = M, ncol = N)
    YMat <- matrix(Y, nrow = M, ncol = N, byrow = TRUE)
    q <- exp(log(0.5)/(MStar^(1/A) - 1))
    e <- YMat^(1/A)
    ep <- (YMat + 1)^(1/A)
    PMF <- (q^e - q^ep)/q
    LogLik <- log(PMF)
    LooRes <- loo(LogLik, cores = Cores)
    Like <- exp(LogLik)
    LikeT <- t(Like)
    InvLikeSum <- rowSums(1/LikeT)
    LogLikeSum <- rowSums(t(LogLik))
    KLD <- log(InvLikeSum/M) + (LogLikeSum/M)
    CPOInv <- rowMeans(1/LikeT)
    LPML <- -sum(log(CPOInv))
    list(LOO = LooRes, KLDivergences = KLD, LPML = LPML, LogLik = LogLik)
  }
  
  Checks <- Diagnostic(Post, JAGSData, Cores = NCores)
  
  capture.output(
    {
      print(Checks$LOO)
      cat("\nLPML: ", Checks$LPML, "\n", sep = "")
    },
    file = LooFile
  )
  
  # KL
  KLVals <- Checks$KLDivergences
  KLData <- data.frame(
    Observation = seq_along(KLVals),
    KL = KLVals,
    KLCapped = pmin(KLVals, 1),
    Influential = ifelse(KLThreshold(KLVals), "Influential", "Not Influential")
  )
  write.csv(KLData, KLDataFile, row.names = FALSE)
  KLPlot <- ggplot(KLData, aes(x = Observation, y = KLCapped, color = Influential)) +
    geom_segment(aes(xend = Observation, yend = 0), linewidth = 0.5) +
    scale_color_manual(values = c("Influential" = "red", "Not Influential" = "blue")) +
    labs(x = "Observation Index", y = "K-L Divergence") +
    scale_x_continuous(breaks = seq(0, max(KLData$Observation), by = 350)) +
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
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  ggsave(filename = KLPlotFile, plot = KLPlot, device = "pdf", width = 8, height = 6)
  
  # DHARMa
  Sim <- SampleTDW(Post, XMat, nsim = NSimDHARMa)
  DHObj <- createDHARMa(simulatedResponse = Sim, observedResponse = YVec, integerResponse = TRUE)
  pdf(file = DHARMaPlotFile, width = 7, height = 6)
  plot(DHObj)
  dev.off()
  capture.output(
    {
      print(testDispersion(DHObj))
      print(testOutliers(DHObj, type = "bootstrap"))
    },
    file = DHARMaTxtFile
  )
  
  # QQ
  DHError <- residuals(DHObj)
  QQData <- data.frame(Expected = seq(0, 1, length.out = length(DHError)), Observed = sort(DHError))
  pdf(file = QQPlotFile, width = 6, height = 5)
  par(mfrow = c(1, 1), oma = c(0.5, 0.5, 0.5, 0.5), mar = c(4, 4, 0.5, 0.5))
  qqplot(QQData$Expected, QQData$Observed,
         main = NULL, xlab = "", ylab = "",
         pch = 16, col = "darkblue", cex = 0.7, bty = "n"
  )
  mtext("Expected Residual", side = 1, line = 3, cex = 1.5)
  mtext("Observed Residual", side = 2, line = 3, cex = 1.5)
  abline(0, 1, col = "red", lwd = 2)
  dev.off()
  
  # Medians
  MediansRes <- GetTDWMedians(Post, XMat, MaxY = 2000)
  write.csv(MediansRes, MediansFile, row.names = FALSE)
  
  list(
    Fit = Fit, Posterior = Post, Medians = MediansRes,
    Diagnostics = list(LOO = Checks$LOO, LPML = Checks$LPML, KLDivergences = Checks$KLDivergences)
  )
}

RunTNBPrepared <- function(
    Adapt = 2000, Burnin = 4000, Sample = 5000, Thin = 5,
    NSimDHARMa = 500, Overwrite = TRUE, Seed = 12345) {
  set.seed(Seed)
  
  BaseTag <- paste0(RunLabel, "_TNB")
  FitFile <- file.path(OutputDir, paste0(BaseTag, "_ModelFit.rds"))
  LooFile <- file.path(OutputDir, paste0(BaseTag, "_LOO.txt"))
  KLPlotFile <- file.path(OutputDir, paste0(BaseTag, "_KL_Divergence.pdf"))
  KLDataFile <- file.path(OutputDir, paste0(BaseTag, "_KL_Data.csv"))
  DHARMaPlotFile <- file.path(OutputDir, paste0(BaseTag, "_DHARMa.pdf"))
  DHARMaTxtFile <- file.path(OutputDir, paste0(BaseTag, "_DHARMa_Tests.txt"))
  QQPlotFile <- file.path(OutputDir, paste0(BaseTag, "_UnifQQ.pdf"))
  
  # -------- fit or load --------
  if (file.exists(FitFile) && !Overwrite) {
    message("Loading existing fit: ", FitFile)
    Fit <- readRDS(FitFile)
  } else {
    message("Running TNB ...")
    ModelCode <- "
      model {
        for (p in 1:P) { beta[p] ~ dnorm(0, 0.001) }
        alpha ~ dgamma(0.001, 0.001)
        for (i in 1:N) {
          LogMu[i] <- inprod(beta[], X[i, ])
          Mu[i] <- exp(LogMu[i])
          Prob[i] <- alpha/(alpha + Mu[i])
          OneMinusProb[i] <- 1 - Prob[i]
          LogPMF_NB[i] <- loggam(Y[i] + alpha) - loggam(alpha) - loggam(Y[i] + 1) +
                          alpha*log(Prob[i]) + Y[i]*log(OneMinusProb[i])
          LogDenom[i] <- log(1 - pow(Prob[i], alpha))
          LogLik[i] <- LogPMF_NB[i] - LogDenom[i]
          Lik[i] <- exp(LogLik[i])
          p[i] <- Lik[i]/C
          ones[i] ~ dbern(p[i])
        }
      }
    "
    
    Inits <- replicate(
      NCores,
      list(
        alpha = 1.5,
        beta = rep(0, ncol(XMat)),
        .RNG.name = "base::Mersenne-Twister",
        .RNG.seed = sample.int(1e6, 1)
      ),
      simplify = FALSE
    )
    Params <- c("beta", "alpha")
    
    Fit <- run.jags(
      model = ModelCode,
      data = JAGSData,
      inits = Inits,
      monitor = Params,
      n.chains = NCores,
      adapt = Adapt,
      burnin = Burnin,
      sample = Sample,
      thin = Thin,
      method = "parallel",
      modules = "glm"
    )
    saveRDS(Fit, FitFile)
    message("Saved fit to: ", FitFile)
  }
  
  MCMC <- as.mcmc(Fit)
  Post <- as.matrix(MCMC)
  
  Diagnostic <- function(Post, DataList, Cores = 1) {
    N <- DataList$N
    X <- DataList$X
    Y <- DataList$Y
    M <- nrow(Post)
    Beta <- Post[, grep("^beta\\[", colnames(Post)), drop = FALSE]
    Alpha <- Post[, "alpha"]
    Fixed <- Beta%*%t(X)
    Mu <- exp(Fixed)
    A <- matrix(Alpha, nrow = M, ncol = N)
    Prob <- A/(A + Mu)
    YMat <- matrix(Y, nrow = M, ncol = N, byrow = TRUE)
    LogPMF <- lgamma(YMat + A) - lgamma(A) - lgamma(YMat + 1) + A*log(Prob) + YMat*log1p(-Prob)
    LogDen <- log(-expm1(A*log(Prob)))
    LogLik <- LogPMF - LogDen
    LooRes <- loo(LogLik, cores = Cores)
    Like <- exp(LogLik)
    LikeT <- t(Like)
    InvLikeSum <- rowSums(1/LikeT)
    LogLikeSum <- rowSums(t(LogLik))
    KLD <- log(InvLikeSum/M) + (LogLikeSum/M)
    CPOInv <- rowMeans(1/LikeT)
    LPML <- -sum(log(CPOInv))
    list(LOO = LooRes, LPML = LPML, KLDivergences = KLD, LogLik = LogLik)
  }
  
  Checks <- Diagnostic(Post, JAGSData, Cores = NCores)
  
  capture.output(
    {
      print(Checks$LOO)
      cat("\nLPML: ", Checks$LPML, "\n", sep = "")
    },
    file = LooFile
  )
  
  # KL
  KLVals <- Checks$KLDivergences
  KLData <- data.frame(
    Observation = seq_along(KLVals),
    KL = KLVals,
    KLCapped = pmin(KLVals, 1),
    Influential = ifelse(KLThreshold(KLVals), "Influential", "Not Influential")
  )
  write.csv(KLData, KLDataFile, row.names = FALSE)
  KLPlot <- ggplot(KLData, aes(x = Observation, y = KLCapped, color = Influential)) +
    geom_segment(aes(xend = Observation, yend = 0), linewidth = 0.5) +
    scale_color_manual(values = c("Influential" = "red", "Not Influential" = "blue")) +
    labs(x = "Observation Index", y = "K-L Divergence") +
    scale_x_continuous(breaks = seq(0, max(KLData$Observation), by = 350)) +
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
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  ggsave(filename = KLPlotFile, plot = KLPlot, device = "pdf", width = 8, height = 6)
  
  # DHARMa
  Sim <- SampleTNB(Post, XMat, nsim = NSimDHARMa)
  DHObj <- createDHARMa(simulatedResponse = Sim, observedResponse = YVec, integerResponse = TRUE)
  pdf(file = DHARMaPlotFile, width = 7, height = 6)
  plot(DHObj)
  dev.off()
  capture.output(
    {
      print(testDispersion(DHObj))
      print(testOutliers(DHObj, type = "bootstrap"))
    },
    file = DHARMaTxtFile
  )
  
  # QQ
  DHError <- residuals(DHObj)
  QQData <- data.frame(Expected = seq(0, 1, length.out = length(DHError)), Observed = sort(DHError))
  pdf(file = QQPlotFile, width = 6, height = 5)
  par(mfrow = c(1, 1), oma = c(0.5, 0.5, 0.5, 0.5), mar = c(4, 4, 0.5, 0.5))
  qqplot(QQData$Expected, QQData$Observed,
         main = NULL, xlab = "", ylab = "",
         pch = 16, col = "darkblue", cex = 0.7, bty = "n"
  )
  mtext("Expected Residual", side = 1, line = 3, cex = 1.5)
  mtext("Observed Residual", side = 2, line = 3, cex = 1.5)
  abline(0, 1, col = "red", lwd = 2)
  dev.off()
  
  list(Fit = Fit, Posterior = Post, Diagnostics = Checks)
}

# =========================
# Run two cTDW fits
# =========================

CTDWResMixed <- list()

CTDWResMixed$L0p500_U1p000_EtaTG <- RunCTDWPrepared(
  DeltaLower = 0.50,
  DeltaUpper = 1.00,
  EtaPrior = "TruncatedGamma",
  Adapt = 2000,
  Burnin = 4000,
  Sample = 5000,
  Thin = 5,
  NSimDHARMa = 500,
  PerDeltaSubdir = FALSE,
  Overwrite = FALSE,
  Seed = 12345
)

CTDWResMixed$L0p000_U1p000_EtaU1_10 <- RunCTDWPrepared(
  DeltaLower = 0.00,
  DeltaUpper = 1.00,
  EtaPrior = "Uniform",
  EtaUniformLower = 1,
  EtaUniformUpper = 10,
  Adapt = 2000,
  Burnin = 4000,
  Sample = 5000,
  Thin = 5,
  NSimDHARMa = 500,
  PerDeltaSubdir = FALSE,
  Overwrite = FALSE,
  Seed = 12346
)

TDWRes <- RunTDWPrepared(
  Adapt = 2000, Burnin = 4000, Sample = 5000, Thin = 5,
  NSimDHARMa = 500,
  Overwrite = FALSE,
  Seed = 12345
)

TNBRes <- RunTNBPrepared(
  Adapt = 2000, Burnin = 4000, Sample = 5000, Thin = 5,
  NSimDHARMa = 500,
  Overwrite = FALSE,
  Seed = 12345
)

# =========================
# TDW vs cTDW median plots
# =========================

FacetLabels <- DataFrame %>%
  select(procedure, admit, sex) %>%
  distinct() %>%
  mutate(
    procedure_label = ifelse(procedure == 0, "PTCA", "CABG"),
    admit_label = ifelse(admit == 0, "Elective", "Urgent/Emerg"),
    sex_label = ifelse(sex == 0, "Female", "Male"),
    FacetLabel = paste0(procedure_label, ", ", admit_label, ", ", sex_label)
  ) %>%
  select(procedure, admit, sex, FacetLabel)

XMatDF <- as.data.frame(XMat) %>%
  mutate(
    procedure = ifelse(`factor(procedure)1` == 1, 1, 0),
    admit = ifelse(`factor(admit)1` == 1, 1, 0),
    sex = ifelse(`factor(sex)1` == 1, 1, 0)
  ) %>%
  select(procedure, admit, sex) %>%
  distinct()

FacetMapping <- XMatDF %>%
  left_join(FacetLabels, by = c("procedure", "admit", "sex")) %>%
  mutate(ComboIndex = dplyr::row_number()) %>%
  select(ComboIndex, FacetLabel)

MakeMedianPlot <- function(cTDWMediansResults, TDWMediansResults, FacetMapping, OutputFile) {
  AllResults <- bind_rows(TDWMediansResults, cTDWMediansResults) %>%
    left_join(FacetMapping, by = c("ComboIndex")) %>%
    mutate(
      FacetLabel = factor(FacetLabel),
      Model = factor(Model, levels = c("TDW", "cTDW"))
    )
  
  pdf(OutputFile, width = 11, height = 6)
  print(
    ggplot(AllResults, aes(x = Model, y = FittedMStarMedian, color = as.factor(Model))) +
      geom_point(position = position_dodge(width = 0.3), size = 3) +
      geom_errorbar(aes(ymin = FittedMStarLower, ymax = FittedMStarUpper),
                    width = 0.2, position = position_dodge(width = 0.3)
      ) +
      facet_wrap(~FacetLabel, nrow = 2, scales = "free_y") +
      labs(x = "Model", y = "Median LOS", color = "Model") +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 30, hjust = 1),
        legend.position = "bottom",
        legend.key.width = grid::unit(1.2, "cm"),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")
      ) +
      scale_y_continuous(breaks = seq(0, max(ceiling(AllResults$FittedMStarUpper), na.rm = TRUE), 0.4)) +
      scale_color_manual(
        values = c("#1f77b4", "#ff7f0e"),
        guide = guide_legend(override.aes = list(linetype = "solid"))
      )
  )
  dev.off()
}

MakeMedianPlot(
  cTDWMediansResults = CTDWResMixed$L0p000_U1p000_EtaU1_10$Medians,
  TDWMediansResults = TDWRes$Medians,
  FacetMapping = FacetMapping,
  OutputFile = file.path(OutputDir, paste0(RunLabel, "_L0p000_U1p000_EtaU1_10_Median_Bands.pdf"))
)

MakeMedianPlot(
  cTDWMediansResults = CTDWResMixed$L0p500_U1p000_EtaTG$Medians,
  TDWMediansResults = TDWRes$Medians,
  FacetMapping = FacetMapping,
  OutputFile = file.path(OutputDir, paste0(RunLabel, "_L0p500_U1p000_EtaTG_Median_Bands.pdf"))
)

cTDWCompare <- bind_rows(
  CTDWResMixed$L0p500_U1p000_EtaTG$Medians %>%
    mutate(Source = "DeltaU05_1_EtaTG"),
  CTDWResMixed$L0p000_U1p000_EtaU1_10$Medians %>%
    mutate(Source = "DeltaU0_1_EtaU1_10")
) %>%
  left_join(FacetMapping, by = c("ComboIndex")) %>%
  mutate(
    FacetLabel = factor(FacetLabel),
    Source = factor(Source, levels = c("DeltaU05_1_EtaTG", "DeltaU0_1_EtaU1_10"))
  )

LevelOrder <- c("DeltaU05_1_EtaTG", "DeltaU0_1_EtaU1_10")

LabelExpr <- expression(
  delta %~% "Uniform" * "(" * 0.5 * ", " * 1 * ")" * ", " ~~
    eta %~% "Gamma" * "(" * 0.001 * ", " * 0.001 * ")" * plain(I)["(1, " * infinity * ")"],
  delta %~% "Uniform" * "(" * 0 * ", " * 1 * ")" * ", " ~~
    eta %~% "Uniform" * "(" * 1 * ", " * 10 * ")"
)

ColorMap <- c(
  "DeltaU05_1_EtaTG" = "#ff7f0e",
  "DeltaU0_1_EtaU1_10" = "black"
)

grDevices::cairo_pdf(
  file.path(OutputDir, paste0(RunLabel, "_cTDW_Median_Comparison_MixedPriors.pdf")),
  width = 11,
  height = 5
)

print(
  ggplot(cTDWCompare, aes(x = Source, y = FittedMStarMedian, color = Source)) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = FittedMStarLower, ymax = FittedMStarUpper),
      width = 0.2
    ) +
    facet_wrap(~FacetLabel, nrow = 2, scales = "free_y") +
    labs(x = NULL, y = "Median LOS", color = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.key.width = grid::unit(1.2, "cm"),
      axis.line = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black")
    ) +
    scale_y_continuous(
      breaks = seq(0, max(ceiling(cTDWCompare$FittedMStarUpper), na.rm = TRUE), 0.25)
    ) +
    scale_x_discrete(
      limits = LevelOrder,
      drop = FALSE
    ) +
    scale_color_manual(
      values = ColorMap,
      breaks = LevelOrder,
      limits = LevelOrder,
      labels = LabelExpr,
      drop = FALSE,
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(linetype = 1)
      )
    )
)

dev.off()

# =========================
# Posterior summary tables
# =========================

BuildPosteriorTable <- function(TDWPost, cTDWPost, OutputFile) {
  SummarizePosterior <- function(Post) {
    Stats <- apply(Post, 2, function(x) {
      c(
        Median = median(x),
        LCI_2.5 = quantile(x, 0.025),
        UCI_97.5 = quantile(x, 0.975)
      )
    })
    rownames(Stats) <- c("Median", "LCI 2.5%", "UCI 97.5%")
    return(Stats)
  }
  
  TDWStats <- SummarizePosterior(TDWPost)
  cTDWStats <- SummarizePosterior(cTDWPost)
  
  TDW_DF <- data.frame(
    Parameter = colnames(TDWStats),
    TDW_median = TDWStats["Median", ],
    TDW_lci = TDWStats["LCI 2.5%", ],
    TDW_uci = TDWStats["UCI 97.5%", ],
    stringsAsFactors = FALSE
  )
  cTDW_DF <- data.frame(
    Parameter = colnames(cTDWStats),
    cTDW_median = cTDWStats["Median", ],
    cTDW_lci = cTDWStats["LCI 2.5%", ],
    cTDW_uci = cTDWStats["UCI 97.5%", ],
    stringsAsFactors = FALSE
  )
  
  MergedDF <- merge(TDW_DF, cTDW_DF, by = "Parameter", all = TRUE)
  
  MergedDF$SortOrder <- sapply(MergedDF$Parameter, function(x) {
    if (grepl("^beta\\[\\d+\\]$", x)) {
      as.numeric(sub("beta\\[(\\d+)\\]", "\\1", x))
    } else if (x == "alpha") 1000 else if (x == "eta") 1001 else if (x == "delta") 1002 else 9999
  })
  MergedDF <- MergedDF[order(MergedDF$SortOrder), ]
  MergedDF$SortOrder <- NULL
  
  ParamToLatex <- function(x) {
    if (x == "alpha") {
      return("$\\alpha$")
    }
    if (x == "delta") {
      return("$\\delta$")
    }
    if (x == "eta") {
      return("$\\eta$")
    }
    if (grepl("^beta\\[\\d+\\]$", x)) {
      idx <- as.numeric(sub("beta\\[(\\d+)\\]", "\\1", x)) - 1
      return(paste0("$\\beta_{", idx, "}$"))
    }
    x
  }
  
  FinalDF <- data.frame(
    Parameter = sapply(MergedDF$Parameter, ParamToLatex),
    TDW_est = sprintf("%.3f", MergedDF$TDW_median),
    TDW_left = sprintf("[%.3f;", MergedDF$TDW_lci),
    TDW_right = sprintf("%.3f]", MergedDF$TDW_uci),
    blank = "",
    cTDW_est = sprintf("%.3f", MergedDF$cTDW_median),
    cTDW_left = sprintf("[%.3f;", MergedDF$cTDW_lci),
    cTDW_right = sprintf("%.3f]", MergedDF$cTDW_uci),
    stringsAsFactors = FALSE
  )
  
  XTab <- xtable(FinalDF, align = rep("l", ncol(FinalDF) + 1))
  sink(OutputFile)
  print(
    XTab,
    include.rownames = FALSE,
    sanitize.text.function = identity,
    caption.placement = "top",
    hline.after = c(-1, 0, nrow(FinalDF))
  )
  sink()
}

BuildPosteriorTable(
  TDWPost = TDWRes$Posterior,
  cTDWPost = CTDWResMixed$L0p500_U1p000_EtaTG$Posterior,
  OutputFile = file.path(OutputDir, paste0(RunLabel, "_Posterior_Summaries_L0p500_U1p000_EtaTG.tex"))
)

BuildPosteriorTable(
  TDWPost = TDWRes$Posterior,
  cTDWPost = CTDWResMixed$L0p000_U1p000_EtaU1_10$Posterior,
  OutputFile = file.path(OutputDir, paste0(RunLabel, "_Posterior_Summaries_L0p000_U1p000_EtaU1_10.tex"))
)

SummarizePosterior <- function(Post) {
  Stats <- apply(Post, 2, function(x) {
    c(
      Median = median(x),
      LCI_2.5 = quantile(x, 0.025),
      UCI_97.5 = quantile(x, 0.975)
    )
  })
  rownames(Stats) <- c("Median", "LCI 2.5%", "UCI 97.5%")
  return(Stats)
}

cTDWStatsDelta05EtaTG <- SummarizePosterior(CTDWResMixed$L0p500_U1p000_EtaTG$Posterior)
cTDWStatsDelta01EtaU <- SummarizePosterior(CTDWResMixed$L0p000_U1p000_EtaU1_10$Posterior)

cTDW_05_DF <- data.frame(
  Parameter = colnames(cTDWStatsDelta05EtaTG),
  cTDW_05_median = cTDWStatsDelta05EtaTG["Median", ],
  cTDW_05_lci = cTDWStatsDelta05EtaTG["LCI 2.5%", ],
  cTDW_05_uci = cTDWStatsDelta05EtaTG["UCI 97.5%", ],
  stringsAsFactors = FALSE
)

cTDW_1_DF <- data.frame(
  Parameter = colnames(cTDWStatsDelta01EtaU),
  cTDW_1_median = cTDWStatsDelta01EtaU["Median", ],
  cTDW_1_lci = cTDWStatsDelta01EtaU["LCI 2.5%", ],
  cTDW_1_uci = cTDWStatsDelta01EtaU["UCI 97.5%", ],
  stringsAsFactors = FALSE
)

MergedDF <- merge(cTDW_05_DF, cTDW_1_DF, by = "Parameter", all = TRUE)

MergedDF$SortOrder <- sapply(MergedDF$Parameter, function(x) {
  if (grepl("^beta\\[\\d+\\]$", x)) {
    as.numeric(sub("beta\\[(\\d+)\\]", "\\1", x))
  } else if (x == "alpha") 1000 else if (x == "eta") 1001 else if (x == "delta") 1002 else 9999
})
MergedDF <- MergedDF[order(MergedDF$SortOrder), ]
MergedDF$SortOrder <- NULL

ParamToLatex <- function(x) {
  if (x == "alpha") {
    return("$\\alpha$")
  }
  if (x == "delta") {
    return("$\\delta$")
  }
  if (x == "eta") {
    return("$\\eta$")
  }
  if (grepl("^beta\\[\\d+\\]$", x)) {
    idx <- as.numeric(sub("beta\\[(\\d+)\\]", "\\1", x)) - 1
    return(paste0("$\\beta_{", idx, "}$"))
  }
  x
}

FinalDF <- data.frame(
  Parameter = sapply(MergedDF$Parameter, ParamToLatex),
  cTDW_05_est = sprintf("%.3f", MergedDF$cTDW_05_median),
  cTDW_05_left = sprintf("[%.3f;", MergedDF$cTDW_05_lci),
  cTDW_05_right = sprintf("%.3f]", MergedDF$cTDW_05_uci),
  blank = "",
  cTDW_1_est = sprintf("%.3f", MergedDF$cTDW_1_median),
  cTDW_1_left = sprintf("[%.3f;", MergedDF$cTDW_1_lci),
  cTDW_1_right = sprintf("%.3f]", MergedDF$cTDW_1_uci),
  stringsAsFactors = FALSE
)

XTab <- xtable(FinalDF, align = rep("l", ncol(FinalDF) + 1))
sink(file.path(OutputDir, paste0(RunLabel, "_Posterior_Summaries_cTDW_Comparison_MixedPriors.tex")))
print(
  XTab,
  include.rownames = FALSE,
  sanitize.text.function = identity,
  caption.placement = "top",
  hline.after = c(-1, 0, nrow(FinalDF))
)
sink()
