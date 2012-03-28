endorse <- function(Y,
                    data,
                    data.village = NA,
                    village = NA,
                    treat = NA,
                    na.strings = 99,
                    identical.lambda = TRUE,
                    covariates = FALSE,
                    formula.indiv = NA,
                    hierarchical = FALSE,
                    formula.village = NA,
                    x.start = 0,
                    s.start = 0,
                    beta.start = 1,
                    tau.start = NA,
                    lambda.start = 0,
                    omega2.start = .1,
                    theta.start = 0,
                    phi2.start = .1,
                    kappa.start = 0,
                    psi2.start = 1,
                    delta.start = 0,
                    zeta.start = 0,
                    rho2.start = 1,
                    mu.beta = 0,
                    mu.x = 0,
                    mu.theta = 0,
                    mu.kappa = 0,
                    mu.delta = 0,
                    mu.zeta = 0,
                    precision.beta = 0.04,
                    precision.x = 1,
                    precision.theta = 0.04,
                    precision.kappa = 0.04,
                    precision.delta = 0.04,
                    precision.zeta = 0.04,
                    s0.omega2 = .2,
                    nu0.omega2 = 30,
                    s0.phi2 = .2,
                    nu0.phi2 = 30,
                    s0.psi2 = .2,
                    nu0.psi2 = 30,
                    s0.sig2 = 1,
                    nu0.sig2 = 400,
                    s0.rho2 = .2,
                    nu0.rho2 = 30,
                    MCMC = 20000,
                    burn = 1000,
                    thin = 1,
                    mh = TRUE,
                    prop = 0.001,
                    x.sd = TRUE,
                    tau.out = FALSE,
                    s.out = FALSE,
                    omega2.out = TRUE,
                    phi2.out = TRUE,
                    psi2.out = TRUE, 
                    verbose = TRUE,
                    seed.store = FALSE,
                    update = FALSE,
                    update.start = NULL) {

  mda = FALSE

  if (!identical.lambda & hierarchical)
    stop("Options are not consistent. If 'identical.lambda = TRUE', 'hierarchical' must be set at
         'FALSE.'")

  if (covariates == FALSE) formula.indiv <- ~ 1

  cov.mat <- model.matrix(formula.indiv, data)
  M <- ncol(cov.mat)

  var.names.indiv <- colnames(cov.mat)

  #############################################
  ## NEED TO MODIFY
  #############################################
  data <- data[complete.cases(cov.mat),] ## do NOT work. NOT DROP NA obs
  cov.mat <- cov.mat[complete.cases(cov.mat), ]

  N <- nrow(data)
  #############################################
  J <- length(Y)

  if (hierarchical) {
    village <- as.integer(as.factor(eval(parse(text = paste("data$", village,
                                                 sep = ""))))) - 1
    G <- length(unique(village))
    if (!is.na(data.village[1, 1])) {
      cov.village.mat <- model.matrix(formula.village, data.village)
      P <- ncol(cov.village.mat)
      var.names.village <- colnames(cov.village.mat)
    }
  } else {
    G <- 1
    P <- 1
    cov.village.mat <- double(1)
    village <- rep(0, times = N)
  }
  
  response <- matrix(NA, nrow = N, ncol = J)
  temp.Qs <- paste("Y$Q", 1:J, sep ="")

  if (is.na(treat[1])) {
    endorse <- matrix(0, nrow = N, ncol = J)
  } else {
    endorse <- treat
  }

  K <- rep(NA, times = J)
  
  if (is.na(na.strings[1])) {
    for (i in 1:J) {
      temp.group <- eval(parse(text = paste("length(", temp.Qs[i], ")", sep ="")))
      K[i] <- temp.group
      
      for (j in 1:temp.group) {
        varname <- eval(parse(text = paste(temp.Qs[i], "[j]", sep = "")))
        response[, i] <- eval(parse(text = paste("ifelse(!is.na(data$", varname,
                                      "), data$", varname, ", response[, i])",
                                      sep = "")))

        if (is.na(treat[1])) {
          endorse[, i] <- eval(parse(text = paste("ifelse(!is.na(data$", varname,
                                       "), j - 1, endorse[, i])",
                                       sep = "")))
        }
      }
    }
  } else {
    for (i in 1:J) {
      temp.group <- eval(parse(text = paste("length(", temp.Qs[i], ")", sep ="")))
      K[i] <- temp.group

      for (j in 1:temp.group) {
        varname <- eval(parse(text = paste(temp.Qs[i], "[j]", sep = "")))
        response[, i] <- eval(parse(text = paste("ifelse(!is.na(data$", varname,
                                      ") & !(data$", varname,
                                      " %in% na.strings), data$", varname,
                                      ", response[, i])", sep = "")))

        if (is.na(treat[1])) {
          endorse[, i] <- eval(parse(text = paste("ifelse(!is.na(data$", varname,
                                       "), j - 1, endorse[, i])",
                                       sep = "")))
        }
      }
    }
  }

  for (i in 1:J){
    response[, i] <- as.integer(as.ordered(response[, i]))
  }

  L <- apply(response, 2, max, na.rm = TRUE)
  max.L <- max(L)
  
  response <- response - 1

  for (i in 1:J) {
    response[, i] <- ifelse(is.na(response[, i]), -1, response[, i])
  }
  
  if (is.na(treat[1])) {
    K <- max(K) - 1
  } else {
    K <- max(treat)
  }

  if (update) {
    beta.start <- update.start$beta.start
    
    tau.start <- matrix(-99, nrow = J, ncol = max.L)
    for (j in 1:J) {
      tau.start[j, 1:(max.L - 1)] <- update.start$tau.start[((max.L - 1) *
                                                             (j - 1) + 1):((max.L -
                                                                            1) * j)]
      tau.start[j, L[j]] <- max(tau.start[j, ]) + 1000
    }

    x.start <- update.start$x.start
    s.start <- update.start$s.start
    lambda.start <- update.start$lambda.start
    omega2.start <- update.start$omega2.start
    theta.start <- update.start$theta.start
    phi2.start <- update.start$phi2.start
    if (covariates) delta.start <- update.start$delta.start
    if (hierarchical) {
      kappa.start <- update.start$kappa.start
      psi2.start <- update.start$psi2.start
      zeta.start <- update.start$zeta.start
      rho2.start <- update.start$rho2.start
    }

    .Random.seed <- update.start$seed
  } else {
    if (is.na(tau.start[1])) {
      tau.start <- matrix(-99, nrow = J, ncol = max.L)

      for (j in 1:J){
        temp.tau <- seq(from = 0, to = .5 * (L[j] - 2), by = .5)
        for (i in 1:(L[j] - 1)) {
          tau.start[j, i] <- temp.tau[i]
        }
        tau.start[j, L[j]] <- max(temp.tau) + 1000
      }
    } else if (class(tau.start) != "matrix" | dim(tau.start)[1] != J | dim(tau.start)[2] != max.L) {
      stop(paste("Incorrect input for tau.start. It should be a ",
                 J, "-by-", max.L, " matrix.", sep = ""))
    }
    
    if (length(x.start) == 1) {
      x.start <- rep(x.start, times = N)
    } else if (length(x.start) != N) {
      stop(paste("The length of x.start should be", N))
    }

    if (length(s.start) == 1) {
      s.start <- as.numeric(matrix(s.start, nrow = N, ncol = J))
      s.start[as.integer(endorse) == 0] <- 0
    } else if (class(s.start) != "matrix" | dim(s.start)[1] != N | dim(s.start)[2] != J) {
      stop(paste("Incorrect input for s.start. It should be a ",
                 N, "-by-", J, " matrix.", sep = ""))
    } else {
      s.start <- as.double(s.start)
      s.start[as.integer(endorse) == 0] <- 0
    }
      
    if (length(beta.start) == 1) {
      beta.start <- matrix(beta.start, nrow = J, ncol = 2)
    } else if (class(beta.start) != "matrix" | dim(beta.start)[1] != J | dim(beta.start)[2] != 2) {
      stop(paste("Incorrect input for beta.start. It should be a ",
                 J, "-by-", 2, " matrix.", sep = ""))
    }

    if (length(lambda.start) == 1) {
      lambda.start <- matrix(lambda.start, nrow = G + M, ncol = J * K)
    } else if (covariates & !identical.lambda & (class(lambda.start) != "matrix" |
                                                 dim(lambda.start)[1] != M |
                                                 dim(lambda.start)[2]  != J * K)) {
      stop(paste("Incorrect input for lambda.start. It should be a ",
                 M, "-by-", J * K, " matrix.", sep = ""))
    } else if (!covariates & !identical.lambda & length(lambda.start) != J * K) {
      stop(paste("Incorrect input for lambda.start. Its length should be ", J * K,
                 sep = ""))
    } else if (covariates & identical.lambda & hierarchical &
               (class(lambda.start) != "matrix" | dim(lambda.start)[1] != G + M - 1 |
                dim(lambda.start)[2] != K)) {
      stop(paste("Incorrect input for lambda.start. It should be a ",
                 G + M - 1, "-by-", K, " matrix.", sep = ""))
    } else if (!covariates & identical.lambda & hierarchical & length(lambda.start) != G * K) {
      stop(paste("Incorrect input for lambda.start. Its length should be ", G * K,
                 sep = ""))
    }

    if (length(omega2.start) == 1) {
      omega2.start <- rep(omega2.start, times = J * K)
    } else if (identical.lambda & length(omega2.start) != K) {
      stop(paste("Incorrect input for omega2.start. Its length should be ", K,
                 sep = ""))
    } else if (!identical.lambda & length(omega2.start) != J * K) {
      stop(paste("Incorrect input for omega2.start. Its length should be ", J * K,
                 sep = ""))
    }

    if (length(theta.start) == 1) {
      theta.start <- matrix(theta.start, nrow = M, ncol = K)
    }

    if (length(phi2.start) == 1) {
      phi2.start <- rep(phi2.start, times = K * M)
    }

  if (hierarchical) {
    if (length(kappa.start) == 1) {
      kappa.start <- rep(kappa.start, times = K * P)
    } else if (class(kappa.start) != "matrix" | dim(kappa.start)[1] != P |
               dim(kappa.start)[2] != K) {
      stop(paste("Incorrect input for kappa.start. It should be a ", P, "-by", K, " matrix.",
                 sep = ""))
    }

    if (length(psi2.start == 1)) {
      psi2.start <- rep(psi2.start, times = K)
    }

    if (length(zeta.start) == 1) {
      zeta.start <- rep(zeta.start, times = P)
    }

  }

    if (length(delta.start) == 1) {
      if (hierarchical) {
        delta.start <- rep(delta.start, times = G + M - 1)
      } else {
        delta.start <- rep(delta.start, times = M)
      }
    } else {
      if (hierarchical & length(delta.start) != (G + M - 1)) {
      	 stop(paste("Incorrect input for delta.start. Its length should be ", G + M - 1,
                 sep = ""))
      } else if (!hierarchical & length(delta.start) != M) {
      	 stop(paste("Incorrect input for delta.start. Its length should be ", M,
                 sep = ""))
      }
    }
  }
  
  if (length(mu.beta) == 1) {
    mu.beta <- matrix(mu.beta, nrow = J, ncol = 2)
  }

  if (length(mu.theta) == 1) {
    mu.theta <- rep(mu.theta, times = M * K)
  } else if (identical.lambda & hierarchical & (class(mu.theta) != "matrix" |
                                                dim(mu.theta)[1] != (M - 1) |
                                                dim(mu.theta)[2] != K)) {
    stop(paste("Incorrect input for mu.theta. It should be a ", M - 1, "-by-", K, " matrix.",
               sep = ""))
  }

  if (covariates) {
    if (length(mu.delta) == 1) {
      if (hierarchical) {
        mu.delta <- rep(mu.delta, times = M - 1)
      } else {
        mu.delta <- rep(mu.delta, times = M)
      }
    }
  } else {
    if (length(mu.x) == 1) {
      mu.delta <- rep(mu.x, times = N)
    } else {
      mu.delta <- mu.x
    }
  }

  if (hierarchical) {
    if (length(mu.kappa) == 1) {
      mu.kappa <- rep(mu.kappa, times = K * P)
    } else if (class(mu.kappa) != "matrix" | dim(mu.kappa)[1] != P | dim(mu.kappa)[2] != K){
      stop(paste("Incorrect input for mu.kappa. It should be a ", P, "-by-", K, " matrix.",
                 sep = ""))
    }

    if (length(mu.zeta) == 1) {
      mu.zeta <- rep(mu.zeta, times = P)
    } else if (length(mu.zeta) != P){
      stop(paste("Incorrect input for mu.zeta. Its length should be ", P,
                 sep = ""))
    }
  }

  precision.beta <- diag(precision.beta, nrow = 2, ncol = 2)

  if (hierarchical)
    precision.delta <- diag(precision.delta, nrow = M - 1, ncol = M - 1)
  else
    precision.delta <- diag(precision.delta, nrow = M, ncol = M)
  

  if (length(precision.theta) == 1) {
    precision.theta <- rep(precision.theta, times = M)
  }

  if (hierarchical) {
    precision.kappa <- diag(precision.kappa, nrow = P, ncol = P)
    precision.zeta <- diag(precision.zeta, nrow = P, ncol = P)
  }

  if (length(prop) != J) {
    prop <- rep(prop, times = J)
  }

  printout <- floor( (MCMC - burn) / thin )

  temp <- .C("R2endorse",
             as.integer(response),
             as.integer(endorse),
             as.double(cov.mat),
             as.double(cov.village.mat),
             as.integer(village),
             as.integer(c(N, J, M, K, max.L)),
             as.integer(L),
             as.integer(c(G, P)),
             as.double(x.start),
             as.double(s.start),
             as.double(beta.start),
             as.double(tau.start),
             as.double(lambda.start),
             as.double(omega2.start),
             as.double(theta.start),
             as.double(phi2.start),
             as.double(kappa.start),
             as.double(psi2.start),
             as.double(delta.start),
             as.double(zeta.start),
             as.double(rho2.start),
             as.double(mu.beta),
             as.double(precision.beta),
             as.double(precision.x),
             as.double(mu.theta),
             as.double(precision.theta),
             as.double(mu.kappa),
             as.double(precision.kappa),
             as.double(mu.delta),
             as.double(precision.delta),
             as.double(mu.zeta),
             as.double(precision.zeta),
             as.double(c(s0.omega2, nu0.omega2, s0.phi2, nu0.phi2, s0.psi2, nu0.psi2, s0.sig2,
                       nu0.sig2, s0.rho2, nu0.rho2)),
             as.integer(c(mda, mh, x.sd, tau.out, s.out, omega2.out, phi2.out,
                          psi2.out, verbose, seed.store, covariates,
                          identical.lambda, hierarchical, MCMC, burn, thin)),
             as.double(prop),
             betaStore = double(printout * 2 * J),
             tauStore = if (tau.out) double(printout * (max.L - 1) * J) else double(1),
             xStore = if (x.sd) double(printout) else double(printout * N),
             sStore = if (s.out) double(printout * N * J) else double(1),
             lambdaStore = if (!identical.lambda) double(printout * J * K * M) else if (hierarchical) double(printout * K * (G + M - 1)) else double(printout * K * M),
             thetaStore = if (!identical.lambda) double(printout * K * M) else double(1),
             kappaStore = if (hierarchical) double(printout * K * P) else double(1),
             deltaStore = if (hierarchical) double(printout * (G + M - 1)) else if (covariates) double(printout * M) else double(1),
             zetaStore = if (hierarchical) double(printout * P) else double(1),
             omega2Store = if (omega2.out & !identical.lambda) double(printout * J * K) else if (omega2.out) double(printout * K) else double(1),
             phi2Store = if (phi2.out & !identical.lambda) double(printout * K * M) else double(1),
             psi2Store = if (psi2.out & hierarchical) double(printout * K) else double(1),
	     sig2Store = if (covariates | hierarchical) double(printout) else double(1),
             rho2Store = if (hierarchical) double(printout) else double(1), 
             betaLast = if (seed.store) double(2 * J) else double(1),
             tauLast = if (seed.store) double((max.L - 1) * J) else double(1),
             xLast = if (seed.store) double(N) else double(1),
             sLast = if (seed.store) double(N * J) else double(1),
             lambdaLast = if (seed.store & !identical.lambda) double(J * K * M) else if (seed.store) double(K * (M + G - 1)) else double(1),
             thetaLast = if (seed.store & !identical.lambda) double(K * M) else double(1),
             kappaLast = if (seed.store & hierarchical) double(K * P) else double(1),
             deltaLast = if (seed.store & hierarchical) double(G + M - 1) else if (seed.store & !hierarchical) double (M) else double(1),
             zetaLast = if (seed.store & hierarchical) double(P) else double(1),
             omega2Last = if (seed.store & !identical.lambda) double(J * K) else if (seed.store) double(K) else double(1),
             phi2Last = if (seed.store) double(K * M) else double(1),
             psi2Last = double(K),
	     sig2Last = double(1),
             rho2Last = double(1),
             accept.ratio = double(J),
             package = "endorse")

  seedStore <- .Random.seed

  res <- list(beta = matrix(as.double(temp$betaStore), byrow = TRUE, ncol = 2 * J,
                nrow = printout),
              tau = if (tau.out) matrix(as.double(temp$tauStore), byrow = TRUE,
                ncol = (max.L - 1) * J, nrow = printout) else NULL,
              x = if (x.sd) matrix(as.double(temp$xStore)[1:printout], ncol = 1,
                nrow = printout) else matrix(as.double(temp$xStore), byrow = TRUE,
                  ncol = N, nrow = printout),
              s = if (s.out) matrix(as.double(temp$sStore), byrow = TRUE,
                ncol = N * J, nrow = printout) else NULL,
              lambda = if (identical.lambda & hierarchical) matrix(as.double(temp$lambdaStore), byrow = TRUE, ncol = (G + M - 1) * K, nrow = printout) else if (identical.lambda) matrix(as.double(temp$lambdaStore), byrow = TRUE, ncol = M * K, nrow = printout) else matrix(as.double(temp$lambdaStore), byrow = TRUE, ncol = J * M * K, nrow = printout),
              theta = if (identical.lambda) NULL else matrix(as.double(temp$thetaStore), byrow = TRUE, ncol = K * M, nrow = printout),
              kappa = if (hierarchical) matrix(as.double(temp$kappaStore), byrow = TRUE, ncol = P * K, nrow = printout) else NULL,
              delta = if (hierarchical) matrix(as.double(temp$deltaStore), byrow = TRUE, ncol = G + M - 1, nrow = printout) else if (covariates) matrix(as.double(temp$deltaStore), byrow = TRUE, ncol = M, nrow = printout) else NULL,
              zeta = if (hierarchical) matrix(as.double(temp$zetaStore), byrow = TRUE, ncol = P, nrow = printout) else NULL,
              omega2 = if (identical.lambda & omega2.out) matrix(as.double(temp$omega2Store), byrow = TRUE, ncol = K, nrow = printout) else if (omega2.out) matrix(as.double(temp$omega2Store), byrow = TRUE, ncol = J * K, nrow = printout) else NULL,
              phi2 = if (!identical.lambda & phi2.out) matrix(as.double(temp$phi2Store), byrow = TRUE, ncol = K * M, nrow = printout) else NULL,
              psi2 = if (hierarchical & psi2.out) matrix(as.double(temp$psi2Store), byrow = TRUE, ncol = K, nrow = printout) else NULL,
	      sig2 = if (hierarchical | covariates) matrix(as.double(temp$sig2Store), ncol = 1, nrow = printout) else NULL,
              rho2 = if (hierarchical) matrix(as.double(temp$rho2Store), byrow = TRUE, ncol = 1, nrow = printout) else NULL,
              accept.ratio = if (mh) as.double(temp$accept.ratio) else NULL,
              seed = if (seed.store) seedStore else NULL,
              beta.start = if (seed.store) matrix(as.double(temp$betaLast), nrow = J, ncol = 2, byrow = TRUE) else NULL,
              tau.start = if (seed.store) as.double(temp$tauLast) else NULL,
              x.start = if (seed.store) as.double(temp$xLast) else NULL,
              s.start = if (seed.store) matrix(as.double(temp$sLast), nrow = N, ncol = J, byrow = TRUE) else NULL,
              lambda.start = if (seed.store & covariates) as.double(temp$lambdaLast) else NULL,
              theta.start = if (seed.store & !identical.lambda & covariates) as.double(temp$thetaLast) else NULL,
              kappa.start = if (seed.store & hierarchical) as.double(temp$kappaLast) else NULL,
              delta.start = if (seed.store & (covariates | hierarchical)) as.double(temp$deltaLast) else NULL,
              zeta.start = if (seed.store & hierarchical) as.double(temp$zetaLast) else NULL,
              omega2.start = if (seed.store) as.double(temp$omega2Last) else NULL,
              phi2.start = if (seed.store & !identical.lambda) as.double(temp$phi2Last) else NULL,
              psi2.start = if (seed.store & hierarchical) as.double(temp$psi2Last) else NULL,
	      sig2.start = if (seed.store & (covariates | hierarchical)) as.double(temp$sig2.start) else NULL,
              rho2.start = if (seed.store & hierarchical) as.double(temp$rho2Last) else NULL,
              village.indicator = village + 1,
              model.matrix.indiv = cov.mat,
              formula.indiv = formula.indiv,
              model.matrix.village = if (hierarchical) cov.village.mat else NULL,
              formula.village = if (hierarchical) formula.village else NULL,
              hierarchical = hierarchical,
              identical.lambda = identical.lambda,
              num.act = K,
              num.policy = J,
              treat = endorse)

  colnames(res$beta) <- paste(rep(c("alpha", "beta"), times = J),
                              rep(1:J, each = 2), sep = ".")
  res$beta <- mcmc(res$beta, start = burn + 1, end = MCMC, thin = thin)

  if (tau.out) {
    temp.names <- paste("tau", 1:J, sep = "")
    colnames(res$tau) <- paste(rep(temp.names, each = (max.L - 1)),
                               rep(1:(max.L - 1), times = J), sep = ".")
    res$tau <- mcmc(res$tau, start = burn + 1, end = MCMC, thin = thin)
  }

  if (x.sd) {
    colnames(res$x) <- "sd.x"
    res$x <- mcmc(res$x, start = burn + 1, end = MCMC, thin = thin)
  } else {
    colnames(res$x) <- paste("x", 1:N, sep = ".")
    res$x <- mcmc(res$x, start = burn + 1, end = MCMC, thin = thin)
  }

  if (s.out) {
    colnames(res$s) <- paste("s", rep(1:nrow(data), each = J),
                             rep(1:J, times = nrow(data)), sep = ".")
  }

  
  if (identical.lambda) {
    if (hierarchical) {
      if (covariates) {
        colnames(res$lambda) <- paste(rep(c(paste("village", 1:G, sep = "."), var.names.indiv[2:M]), times = K),
                                      rep(1:K, each = G + M - 1), sep = ".")
      } else {
        colnames(res$lambda) <- paste(rep(paste("village", 1:G, sep = "."), times = K),
                                      rep(1:K, each = G), sep = ".")
      }
    } else {
      colnames(res$lambda) <- paste(rep(var.names.indiv, times = K),
                                    rep(1:K, each = M), sep =".")
    }
  }
  res$lambda <- mcmc(res$lambda, start = burn + 1, end = MCMC, thin = thin)


  if (identical.lambda) {
    colnames(res$omega2) <- paste("omega2.", 1:K, sep = "")
  } else {
    colnames(res$omega2) <- paste("omega2.", rep(1:J, each = K), rep(1:K,
                                                        times = J), sep = "")
  }
  res$omega2 <- mcmc(res$omega2, start = burn + 1, end = MCMC, thin = thin)  

  if (identical.lambda == FALSE) { 
    temp.names <- paste("theta", 1:K, sep = "")
    if (covariates) {
      colnames(res$theta) <- paste(rep(temp.names, each = M), rep(1:M, times = K),
                                   sep = ".")
    } else {
      colnames(res$theta) <- temp.names
    }
    res$theta <- mcmc(res$theta, start = burn + 1, end = MCMC, thin = thin)

    temp.names <- paste("phi2.", 1:K, sep = "")
    if (covariates) {
      colnames(res$phi2) <- paste(rep(temp.names, each = M), rep(1:M, times = K),
                                  sep = ".")
    } else {
      colnames(res$phi2) <- temp.names
    }
    res$phi2 <- mcmc(res$phi2, start = burn + 1, end = MCMC, thin = thin)
  }

  if (hierarchical) {
    colnames(res$kappa) <- paste(rep(var.names.village, times = K),
                                 rep(1:K, each = P), sep = ".")
    res$kappa <- mcmc(res$kappa, start = burn + 1, end = MCMC, thin = thin)
    if (psi2.out) {
      colnames(res$psi2) <- paste("psi2", 1:K, sep = ".")
      res$psi2 <- mcmc(res$psi2, start = burn + 1, end = MCMC, thin = thin)
    }
  }
  

  if (hierarchical) {
    if (covariates) {
      colnames(res$delta) <- c(paste("village", 1:G, sep = "."), var.names.indiv[2:M])
    } else {
      colnames(res$delta) <- paste("village", 1:G, sep = ".")
    }
    res$delta <- mcmc(res$delta, start = burn + 1, end = MCMC, thin = thin)  
  } else if (covariates) {
      colnames(res$delta) <- var.names.indiv
      res$delta <- mcmc(res$delta, start = burn + 1, end = MCMC, thin = thin)  
  }

  if (hierarchical) {
    colnames(res$sig2) <- "sig2"
    res$sig2 <- mcmc(res$sig2, start = burn + 1, end = MCMC, thin = thin)
  }

  if (hierarchical) {
    colnames(res$zeta) <- var.names.village
    res$zeta <- mcmc(res$zeta, start = burn + 1, end = MCMC, thin = thin)
    colnames(res$rho2) <- "rho2"
    res$rho2 <- mcmc(res$rho2, start = burn + 1, end = MCMC, thin = thin)
  }

  names(res$accept.ratio) <- paste("Q", 1:J, sep = "")

  return(res)
}
