library("party")
library("GeneNet")
library("mvtnorm")
library("MASS")
library("corpcor")
library("networktree")
library("strucchange")
library("Formula")
library("tidyverse")
library("IsingSampler")
library("profvis")

# -----------------------------------------------------------------------------------------
source("~/Documents/Studium/Master_UniversityAmsterdam/Masterthesis/03_Programming/Ising_Fitting_MOB.R")
source("~/Documents/Studium/Master_UniversityAmsterdam/Masterthesis/03_Programming/Sources_Ising/B Ising - EMVS - Functions.R")

#### STEP 0: Simulate Data

IsingDataSimulation <- function(n.nds, n.obs, n.zvar, 
                                delta.cor, delta.all = TRUE, 
                                mean.dif = FALSE, mean.change = .6, 
                                part.type = "numeric"){
  
  # Number of nodes to include in the network
  p <- n.nds    
  # Number of splitting variables to generate
  n.zvar <- n.zvar
  # Number of observations to generate
  n <- n.obs  
  
  # Step A: 
  # Randomly generate an underlying interaction matrix
  sigma <- matrix(0, p, p)                      #matrix of pairwise interactions
  for(s in 1:(p-1)) {
    for(t in (s+1):p) {
      sigma[s, t] <- sigma[t, s] <- runif(1, -1, 1) / 6
    }
  }
  
  # Step B:
  # Modify the interaction matrix
  delta.cor <- delta.cor
  delta.cor.matrix <- matrix(sample(c(delta.cor, - delta.cor), p*p, replace = TRUE), nrow = p)
  
  if (delta.all == TRUE){
    sigma2 <- sigma + delta.cor.matrix
  } else {
    sigma2 <- sigma
    sigma2[1, 2] <- sigma2[1, 2] + delta.cor
    sigma2[2, 1] <- sigma2[2, 1] + delta.cor
  }
  
  sigma2[sigma2 < -.9999999] <- -.9999999
  sigma2[sigma2 > .9999999] <- .9999999
  
  # Step C:
  # Specify the means
  diag(sigma) <- runif(p, -1, 1)                  #main effects on diagonal of sigma
  
  # Modify the means if wanted
  if (mean.dif == TRUE){
    #  modify the mean vector to simulate mean difference
    diag(sigma2)[1] <- diag(sigma)[1] + mean.change
  }
  
  # Step D: 
  # Creating (several) partitioning variables
  
  # In case no change should occur
  
  if (delta.cor == 0){
    if (part.type == "numeric") {
      z1 <- rnorm(n.obs, 0, 1)
    } else if (part.type == "discrete") {
      z1 <- sample(c(0,1, 2, 3), n.obs, replace = TRUE)
    } else {
      z1 <- sample(c(0,1), n.obs, replace = TRUE)
    }
  } else {
    # In case change should occur, sort the part variable
    if (part.type == "numeric") {
      z1 <- sort(rnorm(n.obs, 0, 1))
    } else if (part.type == "discrete") {
      z1 <- rep(c(0,1, 2, 3), each = n.obs/4)
    } else {
      z1 <- rep(c(0,1), each = n.obs/2)
    }
  }
  
  # Create any additional partitioning variables
  if (part.type == "numeric") {
    dsplitvars <- matrix(rnorm(n.obs*(n.zvar-1), 0, 1), nrow = n.obs, ncol = (n.zvar-1))
  } else if (part.type == "discrete") {
    dsplitvars <- matrix(sample(c(0,1, 2, 3), n.obs*(n.zvar-1), replace = TRUE), nrow = n.obs, ncol = (n.zvar-1))
  } else {
    dsplitvars <- matrix(sample(c(0,1), n.obs*(n.zvar-1), replace = TRUE), nrow = n.obs, ncol = (n.zvar-1))
  }

  dsplitvars <- cbind(as.data.frame(z1) , as.data.frame(dsplitvars))
  
  
  # Step E: 
  # Generate data from the determined main effects- interaction - matrix
  # First with the unaltered matrix for n/2 obs
  if(delta.cor == 0) {
    x1 <- matrix(0, n, p)
    for(it in 1:n) {
      for(s in 1:p) 
        x1[,s] <- 1*(rlogis(n) <= sigma[s, s] + x1[,-s] %*% sigma[-s, s])
    }
    dnodevars <- as.data.frame(x1)
  } else {
    x1 <- matrix(0, n/2, p)
    for(it in 1:n/2) {
      for(s in 1:p) 
        x1[,s] <- 1*(rlogis(n/2) <= sigma[s, s] + x1[,-s] %*% sigma[-s, s])
      # if(it / 50 == round(it / 50)) {
      #   plot.ecdf(rowSums(x1))
      #   print(paste("Data simulation; iteration", it))
      # }
    }
    # Second with the altered matrix for n/2 obs
    x2 <- matrix(0, n/2, p)
    for(it in 1:n/2) {
      for(s in 1:p) 
        x2[,s] <- 1*(rlogis(n/2) <= sigma2[s, s] + x2[,-s] %*% sigma2[-s, s])
      # if(it / 50 == round(it / 50)) {
      #   plot.ecdf(rowSums(x2))
      #   print(paste("Data simulation; iteration", it))
      # }
    }
    # rowbind both matrixes 
    dnodevars <- as.data.frame(rbind(x1, x2)) 
  }

  
  
  # Step F:
  # Generate output
  # Combine both node and split variables into one dataframe
  d <- cbind(dsplitvars, dnodevars)
  
  # Specify column names
  colnames(d)[1:n.zvar] <- paste0("z", 1:n.zvar)
  colnames(d)[(n.zvar+1):(n.zvar+p)] <- paste0("y", 1:p)
  
  return(d)
}

# data <- IsingDataSimulation(n.nds = 15, n.obs = 500, n.zvar = 1,
#                             delta.cor = 0.5, delta.all = TRUE, mean.dif = FALSE,
#                             mean.change = .6)

# -----------------------------------------------------------------------------------------

#### STEP 1: Set up Ising fit function 
####### - needs to have objects 
####### - (1) coef (coefficients, estimated parameters), 
####### - (2) logLik (maximized log-likelihood),
####### - (3) emprical estimating function (score function)


IsingFit <- function(y, x = NULL, start = NULL, weights = NULL,
                     offset = NULL, model = c("correlation", "main"), ...,
                     estfun = FALSE, object = FALSE) {
  
  nodevars <- y 
  
  p <- ncol(nodevars) # number of nodes
  n <- nrow(nodevars) # number of observations
  
  nodevars <- as.matrix(nodevars)
  
  ### put dots in a list
  dotlist <- list(...)
  
  ### check if correlation matrix is identified
  if(n <= p*(p-1)/2) {
    stop("Isingfit: n < k*(k-1)/2, correlation matrix is not identified.")
  }
  
  
  # --- STEP 1: 
  # Fitting the Ising Model
  # Using Ising Sampler package
  # fit <- EstimateIsing(nodevars, responses = c(0, 1L), method = c("pl"))
  # R <- fit$graph
  # mean <- fit$thresholds
  
  # -------------------------
  # Maartens Function 
  fit <- maximum.pseudoposterior(x = nodevars, prior.variance = Inf)
  
  # main effect of variables on diagonal
  main <- diag(fit$sigma)
  
  # interactions on off-diagonal
  R <- fit$sigma
  inter <- as.vector(fit$sigma)[as.vector(lower.tri(fit$sigma))]
  diag(R) <- 0
  
  # ----------------
  
  # --- STEP 2: Calculate score functions 
  nodevars <- t(nodevars)
  # Step 2a: 
  # Compute matrixes to make actual computation easier
  M <- R%*%nodevars # compute the sums 
  A <- exp(main + M)
  D <- 1 + A
  E <- A/ D
  
  
  # Step 2b: 
  ## Score Functions
  
  # Score Function for main effect
  score_main <- nodevars - E # 0 up to the 3rd decimal
  
  
  # Score Function for interaction 
  score_rho <- combn(p, 2,
                     function(x) (2*(nodevars[x[1] ,] * nodevars[x[2], ]) - (E[x[1], ] * nodevars[x[2], ]) - (E[x[2], ] * nodevars[x[1], ])) )
  
  
  # --- STEP 3: Store objects depending which measures to consider
  # Set-Up
  scores <- c()
  coef <- c()
  objnames <- c()
  ynam <- if(is.null(colnames(nodevars))) 1L:p else colnames(nodevars)
  
  # Main 
  if ("main" %in% model) {
    scores <- cbind(scores, t(score_main))
    coef <- c(coef, main)
    objnames <- c(objnames, paste0("main_", ynam))
  }
  # Correlation
  if ("correlation" %in% model) {
    scores <- cbind(scores, score_rho)
    coef <- c(coef, R[upper.tri(R)])
    objnames <- c(objnames, combn(p, 2, function(x) paste0("cor_", x[1], "_", x[2])))
  }
  
  #  colnames(scores) <- objnames
  names(coef) <- objnames
  colnames(scores) <- objnames
  
  # --- STEP 4: Compute Log Likelihood of Ising function
  S <- main + M
  loglik <- sum(nodevars*S) - sum(log(D))
  
  if(estfun == FALSE) {
    scores <- NULL
  }
  
  # --- STEP 5: Compute inverse of Hessian
  # index <- rbind(1: (p * (p - 1)/2), combn(p, 2))
  # inverse.hessian <- invert.hessian(sigma.map = fit$sigma, index = t(index), x = t(nodevars), prior.variance = Inf)
  # vc <- - inverse.hessian
  vc <- NULL
  
  res <- list(coefficients = coef,
              objfun = -loglik,
              estfun = scores,
              object = vc)
  
  return(res)
}

#fit <- IsingFit(y = data[, 2:16], estfun = TRUE, model = c("main", "correlation"))

#-----------------------------------------------------

#### STEP 2: Feed function to mob 
####### - fit is the self-defined fit function 
####### - control is pre-defined control function 
####### - 

IsingNetworkTree <- function(nodevars, splitvars, #Data 
                             type = c("cor", "pcor", "glasso"),
                             model = c("main", "correlation"), # Define which parameters to split by
                             alpha = .05, bonferroni = TRUE, minsplit = "simple" ) { # Control for MOB  
  
  # Rename column names of the node and split variables
  nodevars <- nodevars
  splitvars <- splitvars
  
  if(is.null(colnames(nodevars))){
    if(ncol(nodevars) == 1){
      colnames(nodevars) <- "nodevars1"
    } else {
      colnames(nodevars) <- paste('nodevars',1:ncol(nodevars),sep="")
    }
  }
  
  if(is.null(colnames(splitvars))){
    if(ncol(splitvars) == 1){
      colnames(splitvars) <- "splitvars1"
    } else {
      colnames(splitvars) <- paste('splitvars',1:ncol(splitvars),sep="")
    }
  }
  
  # Create dataframe used by the MOB function 
  d <- data.frame(nodevars,splitvars)
  
  # Write the formula; everything behind the ~ is used as a partitioning variable
  form <- paste(paste(colnames(nodevars), collapse=" + "), "~",paste(colnames(splitvars), collapse=" + "))
  form <- as.formula(form)
  
  # Define control variables for the MOB
  if(minsplit == "simple") {
    p <- ncol(nodevars)
    minsplit <- (p*(p-1)) # !!!! Further define what the minimum sample size needed is for a reliable estimate
  }
  control <- partykit::mob_control(alpha = alpha, bonferroni = bonferroni, minsplit = minsplit)
  control$ytype <- "matrix"
  
  # Run MOB function
  res <- partykit::mob(formula = form, data = d, fit = IsingFit,  model = model, control = control)
  
  # Change class of 
  class(res) <- c("networktree", "mob_networktree", type[1], class(res))
  
  # Print results
  return(res)
}
# 
# data <- IsingDataSimulation(n.nds = 15, n.obs = 500, n.zvar = 1, delta.cor = 0.3,
#                             delta.all = FALSE, mean.dif = FALSE, mean.change = .6,
#                             part.type = "discrete")
# 
# nodevars <- as.data.frame(data[, 2:16])
# splitvars <- as.data.frame(data[, 1])
# 
# test <- ISINGnetworktree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
# test


