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
