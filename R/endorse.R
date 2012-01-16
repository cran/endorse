endorse <- function(Y, # = list(Q1 = c(VARIABLE NAMES), ...)
                    data,
                    treat = NA, # Recommended: N * J matrix of treatment status
                    na.strings = NA,
                    covariates = FALSE,
                    formula = NA,
                    x.start = NA,  # scalar or vector of length N
                    s.start = 0,  # scalar or N * J matrix
                    beta.start = 1, # scalar or J * 2 matrix
                    tau.start = NA,   # J * max(L_{j}) matrix
                    lambda.start = 0, # scalar or (J * M) * K matrix
                    omega2.start = 1, #scalar or J * K matrix
                    theta.start = 0,   #scalar or vector of length K * M
                    phi2.start = 1, # scalar or vector of length K * M
                    delta.start = 0, # scalar or vector of length M
                    mu.beta = 0, #scalar or J * 2 matrix
                    mu.x = 0, # scalar or vector of length N
                    mu.theta = 0, # scalar or vector of length M
                    mu.delta = 0, # scalar or vector of length M
                    precision.beta = 0.1, # scalar or 2 * 2 diagonal matrix
                    precision.x = 1,    # scalar
                    precision.theta = 0.1,# scalar or vector of length M
                    precision.delta = 0.1, # scalar or M * M matrix
                    s0.omega2= 1,
                    nu0.omega2 = 1,
                    s0.phi2 = 1,
                    nu0.phi2 = 1,
                    MCMC = 20000,
                    burn = 1000,
                    thin = 1,
                    mda = TRUE,
                    mh = TRUE,
                    prop = 0.001,
                    x.sd = TRUE,
                    tau.out = FALSE
                    ) {

  if (covariates) {
    cov.mat <- model.matrix(formula, data)
    M <- ncol(cov.mat)

    for (i in 1:M) {
      data <- data[!is.na(cov.mat[, i]),]
      cov.mat <- cov.mat[!is.na(cov.mat[, i]), ]
    }

  } else {
    M <- 1
  }

  N <- nrow(data)
  J <- length(Y)
  
  response <- matrix(-1, nrow = N, ncol = J)
  temp.Qs <- paste("Y$Q", 1:J, sep ="")

  if (is.na(treat[1])) {
    endorse <- matrix(0, nrow = N, ncol = J)
  } else {
    endorse <- treat
  }

  K <- rep(NA, times = J)  
  
  if (is.na(na.strings)) {
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
        response[, i] <- eval(parse(text = paste("ifelse(data$", varname,
                                      " != na.strings, data$", varname,
                                      ", response[, i])", sep = "")))

        if (is.na(treat[1])) {
          endorse[, i] <- eval(parse(text = paste("ifelse(data$", varname,
                                        " != na.strings, j - 1, response[, i])",
                                        sep = "")))
        }
      }
    }
  }

  for (i in 1:J){
    response[, i] <- as.integer(as.ordered(response[, i]))
  }

  L <- apply(response, 2, max) - 1
  max.L <- max(L)
  
  response <- response - 2
  
  if (is.na(treat[1])) {
    K <- max(K) - 1
  } else {
    K <- max(treat)
  }

  if (is.na(tau.start[1])) {
    tau.start <- matrix(-99, nrow = J, ncol = max.L)

    for (j in 1:J){
      temp.tau <- seq(from = 0, to = .5 * (L[j] - 2), by = .5)
      for (i in 1:(L[j] - 1)) {
        tau.start[j, i] <- temp.tau[i]
      }
      tau.start[j, L[j]] <- max(temp.tau) + 1000
    }
  }
  
  if (is.na(x.start)) {
    x.start <- rnorm(N)
  } else if (length(x.start) == 1) {
    x.start <- rep(x.start, times = N)
  }

  if (length(s.start) == 1) {
    s.start <- matrix(s.start, nrow = N, ncol = J)
  }

  if (length(beta.start) == 1) {
    beta.start <- matrix(beta.start, nrow = J, ncol = 2)
  }

  if (length(lambda.start) == 1) {
    lambda.start <- matrix(lambda.start, nrow = J * M, ncol = K)
  }

  if (length(omega2.start) == 1) {
    omega2.start <- matrix(omega2.start, nrow = J, ncol = K)
  }

  if (length(theta.start) == 1) {
    theta.start <- rep(theta.start, times = K * M)
  }

  if (length(phi2.start) == 1) {
    phi2.start <- rep(phi2.start, times = K * M)
  }

  if (length(delta.start) == 1) {
    delta.start <- rep(delta.start, times = M)
  }
  
  if (length(mu.beta) == 1) {
    mu.beta <- matrix(mu.beta, nrow = J, ncol = 2)
  }

  if (length(mu.x) == 1) {
    mu.x <- rep(mu.x, times = N)
  }

  if (length(mu.theta) == 1) {
    mu.theta <- rep(mu.theta, times = M)
  }

  if (length(mu.delta) == 1) {
    mu.delta <- rep(mu.delta, times = M)
  }

  precision.beta <- diag(precision.beta, nrow = 2, ncol = 2)

  precision.delta <- diag(precision.delta, nrow = M, ncol = M)

  if (length(precision.theta) == 1) {
    precision.theta <- rep(precision.theta, times = M)
  }

  accept <- 0
  
  printout <- floor( (MCMC - burn) / thin )

  if (covariates) {
    if (x.sd) {
      if (tau.out) {
        temp <- .C("R2endorse",
                   as.integer(response),
                   as.integer(endorse),
                   as.double(cov.mat),
                   as.integer(N),
                   as.integer(J),
                   as.integer(M),
                   as.integer(K),
                   as.integer(L),
                   as.integer(max.L),
                   as.double(x.start),  # vector of length N
                   as.double(s.start),  # N * J matrix
                   as.double(beta.start),
                   as.double(tau.start),
                   as.double(lambda.start), # lambda
                   as.double(omega2.start),  # omega
                   as.double(theta.start),  # theta
                   as.double(phi2.start),   #phi
                   as.double(delta.start),
                   as.double(mu.beta),
                   as.double(precision.beta),
                   as.double(precision.x),
                   as.double(mu.theta),
                   as.double(precision.theta),
                   as.double(mu.delta),
                   as.double(precision.delta),
                   as.double(s0.omega2),
                   as.integer(nu0.omega2),
                   as.double(s0.phi2),
                   as.integer(nu0.phi2),
                   as.integer(MCMC),
                   as.integer(burn),
                   as.integer(thin),
                   as.integer(mda),
                   as.integer(mh),
                   as.double(prop),
                   as.integer(accept),
                   as.integer(x.sd),
                   as.integer(tau.out),
                   betaStore = double(printout * 2 * J),
                   tauStore = double(printout * (max(L)-1) * J),
                   xStore = double(printout),
                   lambdaStore = double(printout * J * M * K),
                   thetaStore = double(printout * K * M),
                   deltaStore = double(printout * M),
                   package = "endorse")
      } else {
        temp <- .C("R2endorse",
                   as.integer(response),
                   as.integer(endorse),
                   as.double(cov.mat),
                   as.integer(N),
                   as.integer(J),
                   as.integer(M),
                   as.integer(K),
                   as.integer(L),
                   as.integer(max.L),
                   as.double(x.start),  # vector of length N
                   as.double(s.start),  # N * J matrix
                   as.double(beta.start),
                   as.double(tau.start),
                   as.double(lambda.start), # lambda
                   as.double(omega2.start),  # omega
                   as.double(theta.start),  # theta
                   as.double(phi2.start),   #phi
                   as.double(delta.start),
                   as.double(mu.beta),
                   as.double(precision.beta),
                   as.double(precision.x),
                   as.double(mu.theta),
                   as.double(precision.theta),
                   as.double(mu.delta),
                   as.double(precision.delta),
                   as.double(s0.omega2),
                   as.integer(nu0.omega2),
                   as.double(s0.phi2),
                   as.integer(nu0.phi2),
                   as.integer(MCMC),
                   as.integer(burn),
                   as.integer(thin),
                   as.integer(mda),
                   as.integer(mh),
                   as.double(prop),
                   as.integer(accept),
                   as.integer(x.sd),
                   as.integer(tau.out),
                   betaStore = double(printout * 2 * J),
                   tauStore = double(1),
                   xStore = double(printout),
                   lambdaStore = double(printout * J * M * K),
                   thetaStore = double(printout * K * M),
                   deltaStore = double(printout * M),
                   package = "endorse")
      }
    } else {
      if (tau.out) {
        temp <- .C("R2endorse",
                   as.integer(response),
                   as.integer(endorse),
                   as.double(cov.mat),
                   as.integer(N),
                   as.integer(J),
                   as.integer(M),
                   as.integer(K),
                   as.integer(L),
                   as.integer(max.L),
                   as.double(x.start),  # vector of length N
                   as.double(s.start),  # N * J matrix
                   as.double(beta.start),
                   as.double(tau.start),
                   as.double(lambda.start), # lambda
                   as.double(omega2.start),  # omega
                   as.double(theta.start),  # theta
                   as.double(phi2.start),   #phi
                   as.double(delta.start),
                   as.double(mu.beta),
                   as.double(precision.beta),
                   as.double(precision.x),
                   as.double(mu.theta),
                   as.double(precision.theta),
                   as.double(mu.delta),
                   as.double(precision.delta),
                   as.double(s0.omega2),
                   as.integer(nu0.omega2),
                   as.double(s0.phi2),
                   as.integer(nu0.phi2),
                   as.integer(MCMC),
                   as.integer(burn),
                   as.integer(thin),
                   as.integer(mda),
                   as.integer(mh),
                   as.double(prop),
                   as.integer(accept),
                   as.integer(x.sd),
                   as.integer(tau.out),
                   betaStore = double(printout * 2 * J),
                   tauStore = double(printout * (max(L)-1) * J),
                   xStore = double(printout * N),
                   lambdaStore = double(printout * J * M * K),
                   thetaStore = double(printout * K * M),
                   deltaStore = double(printout * M),
                   package = "endorse")
      } else {
        temp <- .C("R2endorse",
                   as.integer(response),
                   as.integer(endorse),
                   as.double(cov.mat),
                   as.integer(N),
                   as.integer(J),
                   as.integer(M),
                   as.integer(K),
                   as.integer(L),
                   as.integer(max.L),
                   as.double(x.start),  # vector of length N
                   as.double(s.start),  # N * J matrix
                   as.double(beta.start),
                   as.double(tau.start),
                   as.double(lambda.start), # lambda
                   as.double(omega2.start),  # omega
                   as.double(theta.start),  # theta
                   as.double(phi2.start),   #phi
                   as.double(delta.start),
                   as.double(mu.beta),
                   as.double(precision.beta),
                   as.double(precision.x),
                   as.double(mu.theta),
                   as.double(precision.theta),
                   as.double(mu.delta),
                   as.double(precision.delta),
                   as.double(s0.omega2),
                   as.integer(nu0.omega2),
                   as.double(s0.phi2),
                   as.integer(nu0.phi2),
                   as.integer(MCMC),
                   as.integer(burn),
                   as.integer(thin),
                   as.integer(mda),
                   as.integer(mh),
                   as.double(prop),
                   as.integer(accept),
                   as.integer(x.sd),
                   as.integer(tau.out),
                   betaStore = double(printout * 2 * J),
                   tauStore = double(1),
                   xStore = double(printout * N),
                   lambdaStore = double(printout * J * M * K),
                   thetaStore = double(printout * K * M),
                   deltaStore = double(printout * M),
                   package = "endorse")
      }
    }
  } else {
    if (x.sd) {
      if (tau.out) {
        temp <- .C("R2endorseNoCov",
                   as.integer(response),
                   as.integer(endorse),
                   as.integer(N),
                   as.integer(J),
                   as.integer(K),
                   as.integer(L),
                   as.integer(max.L),
                   as.double(x.start),  # vector of length N
                   as.double(s.start),  # N * J matrix
                   as.double(beta.start),
                   as.double(tau.start),
                   as.double(lambda.start), # lambda
                   as.double(omega2.start),  # omega
                   as.double(theta.start),  # theta
                   as.double(phi2.start),   #phi
                   as.double(mu.beta),
                   as.double(precision.beta),
                   as.double(mu.x),
                   as.double(precision.x),
                   as.double(mu.theta),
                   as.double(precision.theta),
                   as.double(s0.omega2),
                   as.integer(nu0.omega2),
                   as.double(s0.phi2),
                   as.integer(nu0.phi2),
                   as.integer(MCMC),
                   as.integer(burn),
                   as.integer(thin),
                   as.integer(mda),
                   as.integer(mh),
                   as.double(prop),
                   as.integer(accept),
                   as.integer(x.sd),
                   as.integer(tau.out),
                   betaStore = double(printout * 2 * J),
                   tauStore = double(printout * (max(L)-1) * J),
                   xStore = double(printout),
                   lambdaStore = double(printout * J * M * K),
                   thetaStore = double(printout * K * M),
                   package = "endorse")
      } else {
        temp <- .C("R2endorseNoCov",
                   as.integer(response),
                   as.integer(endorse),
                   as.integer(N),
                   as.integer(J),
                   as.integer(K),
                   as.integer(L),
                   as.integer(max.L),
                   as.double(x.start),  # vector of length N
                   as.double(s.start),  # N * J matrix
                   as.double(beta.start),
                   as.double(tau.start),
                   as.double(lambda.start), # lambda
                   as.double(omega2.start),  # omega
                   as.double(theta.start),  # theta
                   as.double(phi2.start),   #phi
                   as.double(mu.beta),
                   as.double(precision.beta),
                   as.double(mu.x),
                   as.double(precision.x),
                   as.double(mu.theta),
                   as.double(precision.theta),
                   as.double(s0.omega2),
                   as.integer(nu0.omega2),
                   as.double(s0.phi2),
                   as.integer(nu0.phi2),
                   as.integer(MCMC),
                   as.integer(burn),
                   as.integer(thin),
                   as.integer(mda),
                   as.integer(mh),
                   as.double(prop),
                   as.integer(accept),
                   as.integer(x.sd),
                   as.integer(tau.out),
                   betaStore = double(printout * 2 * J),
                   tauStore = double(1),
                   xStore = double(printout),
                   lambdaStore = double(printout * J * M * K),
                   thetaStore = double(printout * K * M),
                   package = "endorse")
      }
    } else {
            if (tau.out) {
        temp <- .C("R2endorseNoCov",
                   as.integer(response),
                   as.integer(endorse),
                   as.integer(N),
                   as.integer(J),
                   as.integer(K),
                   as.integer(L),
                   as.integer(max.L),
                   as.double(x.start),  # vector of length N
                   as.double(s.start),  # N * J matrix
                   as.double(beta.start),
                   as.double(tau.start),
                   as.double(lambda.start), # lambda
                   as.double(omega2.start),  # omega
                   as.double(theta.start),  # theta
                   as.double(phi2.start),   #phi
                   as.double(mu.beta),
                   as.double(precision.beta),
                   as.double(mu.x),
                   as.double(precision.x),
                   as.double(mu.theta),
                   as.double(precision.theta),
                   as.double(s0.omega2),
                   as.integer(nu0.omega2),
                   as.double(s0.phi2),
                   as.integer(nu0.phi2),
                   as.integer(MCMC),
                   as.integer(burn),
                   as.integer(thin),
                   as.integer(mda),
                   as.integer(mh),
                   as.double(prop),
                   as.integer(accept),
                   as.integer(x.sd),
                   as.integer(tau.out),
                   betaStore = double(printout * 2 * J),
                   tauStore = double(printout * (max(L)-1) * J),
                   xStore = double(printout * N),
                   lambdaStore = double(printout * J * M * K),
                   thetaStore = double(printout * K * M),
                   package = "endorse")
      } else {
        temp <- .C("R2endorseNoCov",
                   as.integer(response),
                   as.integer(endorse),
                   as.integer(N),
                   as.integer(J),
                   as.integer(K),
                   as.integer(L),
                   as.integer(max.L),
                   as.double(x.start),  # vector of length N
                   as.double(s.start),  # N * J matrix
                   as.double(beta.start),
                   as.double(tau.start),
                   as.double(lambda.start), # lambda
                   as.double(omega2.start),  # omega
                   as.double(theta.start),  # theta
                   as.double(phi2.start),   #phi
                   as.double(mu.beta),
                   as.double(precision.beta),
                   as.double(mu.x),
                   as.double(precision.x),
                   as.double(mu.theta),
                   as.double(precision.theta),
                   as.double(s0.omega2),
                   as.integer(nu0.omega2),
                   as.double(s0.phi2),
                   as.integer(nu0.phi2),
                   as.integer(MCMC),
                   as.integer(burn),
                   as.integer(thin),
                   as.integer(mda),
                   as.integer(mh),
                   as.double(prop),
                   as.integer(accept),
                   as.integer(x.sd),
                   as.integer(tau.out),
                   betaStore = double(printout * 2 * J),
                   tauStore = double(1),
                   xStore = double(printout * N),
                   lambdaStore = double(printout * J * M * K),
                   thetaStore = double(printout * K * M),
                   package = "endorse")
      }
    }
  }

  if (covariates) {
    if (tau.out) {
      if (x.sd) {
        res <- list(beta = matrix(as.double(temp$betaStore), byrow = TRUE,
                       ncol = 2*J, nrow = printout),
                    tau = matrix(as.double(temp$tauStore), byrow = TRUE,
                      ncol = (max.L-1)*J, nrow = printout),
                    x = matrix(as.double(temp$xStore), byrow = TRUE, ncol = 1,
                      nrow = printout),
                    lambda = matrix(as.double(temp$lambdaStore), byrow = TRUE,
                      ncol = J * M * K, nrow = printout),
                    theta = matrix(as.double(temp$thetaStore), byrow = TRUE,
                      ncol = K * M, nrow = printout),
                    delta = matrix(as.double(temp$deltaStore), byrow = TRUE,
                      ncol = M, nrow = printout))
      } else {
        res <- list(beta = matrix(as.double(temp$betaStore), byrow = TRUE,
                       ncol = 2*J, nrow = printout),
                    tau = matrix(as.double(temp$tauStore), byrow = TRUE,
                      ncol = (max.L-1)*J, nrow = printout),
                    x = matrix(as.double(temp$xStore), byrow = TRUE, ncol = N,
                      nrow = printout),
                    lambda = matrix(as.double(temp$lambdaStore), byrow = TRUE,
                      ncol = J * M * K, nrow = printout),
                    theta = matrix(as.double(temp$thetaStore), byrow = TRUE,
                      ncol = K * M, nrow = printout),
                    delta = matrix(as.double(temp$deltaStore), byrow = TRUE,
                      ncol = M, nrow = printout))
      }
    } else {
      if (x.sd) {
        res <- list(beta = matrix(as.double(temp$betaStore), byrow = TRUE,
                       ncol = 2*J, nrow = printout ),
                    x = matrix(as.double(temp$xStore), byrow = TRUE, ncol = 1,
                      nrow = printout),
                    lambda = matrix(as.double(temp$lambdaStore), byrow = TRUE,
                      ncol = J * M * K, nrow = printout),
                    theta = matrix(as.double(temp$thetaStore), byrow = TRUE,
                      ncol = K * M, nrow = printout),
                    delta = matrix(as.double(temp$deltaStore), byrow = TRUE,
                      ncol = M, nrow = printout))
      } else {
        res <- list(beta = matrix(as.double(temp$betaStore), byrow = TRUE,
                       ncol = 2*J, nrow = printout ),
                    x = matrix(as.double(temp$xStore), byrow = TRUE, ncol = N,
                      nrow = printout),
                    lambda = matrix(as.double(temp$lambdaStore), byrow = TRUE,
                      ncol = J * M * K, nrow = printout),
                    theta = matrix(as.double(temp$thetaStore), byrow = TRUE,
                      ncol = K * M, nrow = printout),
                    delta = matrix(as.double(temp$deltaStore), byrow = TRUE,
                      ncol = M, nrow = printout))
      }
    }
  } else {
    if (tau.out) {
      if (x.sd) {
        res <- list(beta = matrix(as.double(temp$betaStore), byrow = TRUE,
                       ncol = 2*J, nrow = printout ),
                    tau = matrix(as.double(temp$tauStore), byrow = TRUE,
                      ncol = (max.L-1)*J, nrow = printout),
                    x = matrix(as.double(temp$xStore), byrow = TRUE, ncol = 1,
                      nrow = printout),
                    lambda = matrix(as.double(temp$lambdaStore), byrow = TRUE,
                      ncol = J * M * K, nrow = printout),
                    theta = matrix(as.double(temp$thetaStore), byrow = TRUE,
                      ncol = K * M, nrow = printout))
      } else {
        res <- list(beta = matrix(as.double(temp$betaStore), byrow = TRUE,
                       ncol = 2*J, nrow = printout ),
                    tau = matrix(as.double(temp$tauStore), byrow = TRUE,
                      ncol = (max.L-1)*J, nrow = printout),
                    x = matrix(as.double(temp$xStore), byrow = TRUE, ncol = N,
                      nrow = printout),
                    lambda = matrix(as.double(temp$lambdaStore), byrow = TRUE,
                      ncol = J * M * K, nrow = printout),
                    theta = matrix(as.double(temp$thetaStore), byrow = TRUE,
                      ncol = K * M, nrow = printout))
      }
    } else {
      if (x.sd) {
        res <- list(beta = matrix(as.double(temp$betaStore), byrow = TRUE,
                       ncol = 2*J, nrow = printout ),
                    x = matrix(as.double(temp$xStore), byrow = TRUE, ncol = 1,
                      nrow = printout),
                    lambda = matrix(as.double(temp$lambdaStore), byrow = TRUE,
                      ncol = J * M * K, nrow = printout),
                    theta = matrix(as.double(temp$thetaStore), byrow = TRUE,
                      ncol = K * M, nrow = printout))
      } else {
        res <- list(beta = matrix(as.double(temp$betaStore), byrow = TRUE,
                       ncol = 2*J, nrow = printout ),
                    x = matrix(as.double(temp$xStore), byrow = TRUE, ncol = N,
                      nrow = printout),
                    lambda = matrix(as.double(temp$lambdaStore), byrow = TRUE,
                      ncol = J * M * K, nrow = printout),
                    theta = matrix(as.double(temp$thetaStore), byrow = TRUE,
                      ncol = K * M, nrow = printout))
      }
    }
  }

  colnames(res$beta) <- paste(rep(c("alpha", "beta"), times = J),
                              rep(1:J, each = 2), sep = ".")
  res$beta <- as.mcmc(res$beta)

  if (tau.out) {
    temp.names <- paste("tau", 1:J, sep = "")
    colnames(res$tau) <- paste(rep(temp.names, each = (max.L - 1)),
                               rep(1:(max.L - 1), times = J), sep = ".")
    res$tau <- as.mcmc(res$tau)
  }

  if (x.sd) {
    colnames(res$x) <- "sd.x"
  } else {
    colnames(res$x) <- paste("x", 1:N, sep = ".")
    res$x <- as.mcmc(res$x)
  }

  temp.names <- paste("lambda", rep(1:J, each = K), rep(1:K, times = J), sep = "")
  colnames(res$lambda) <- paste(rep(temp.names, each = M), rep(1:M, times = (J * K)),
                                sep = ".")
  res$lambda <- as.mcmc(res$lambda)
  
  temp.names <- paste("theta", 1:K, sep = "")
  colnames(res$theta) <- paste(rep(temp.names, each = M), rep(1:M, times = K), sep = ".")
  res$theta <- as.mcmc(res$theta)

  if (covariates) {
    colnames(res$delta) <- paste("delta", 1:M, sep = "")
    res$delta <- as.mcmc(res$delta)
  }

  return(res)
}
