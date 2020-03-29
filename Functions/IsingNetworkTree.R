library("party")
library("mvtnorm")
library("MASS")
library("corpcor")
library("networktree")
library("strucchange")
library("Formula")
library("tidyverse")

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
