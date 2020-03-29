# Functions necessary for a computationally inexpensive Pseudolikelihood-Estimation
# Credits for those functions to Maarten Marsman 

one.dimensional.NR <- function(sigma.map, x, suf.stat, prior.variance) {
  p <- ncol(x)

  ##############################################################################
  # MAIN EFFECTS; SIGMA[S,S]
  ##############################################################################
  for(s in 1:p) {
    phi <- sigma.map[s, s] + x[, -s] %*% sigma.map[s, -s]
    pr <- exp(phi) / (1 + exp(phi))
    
    d <- suf.stat[s, s] - sum(pr) - sigma.map[s, s] / prior.variance
    dd <- -sum(pr * (1 - pr)) - 1 / prior.variance
    
    sigma.map[s, s] <- sigma.map[s, s] - d / dd
  }
  
  ##############################################################################
  # INTERACTION EFFECTS; SIGMA[S, T]
  ##############################################################################
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      phi <- sigma.map[s, s] + x[, -s] %*% sigma.map[-s, s]
      prS <- exp(phi) / (1 + exp(phi))
      
      phi <- sigma.map[t, t] + x[, -t] %*% sigma.map[-t, t]
      prT <- exp(phi) / (1 + exp(phi))
      
      d <- 2 * suf.stat[s, t] -  x[, t] %*% prS - x[, s] %*% prT - sigma.map[s, t] / prior.variance
      dd <- -sum(x[, t] * (prS * (1 - prS))) - sum(x[, s] * (prT * (1 - prT))) - 1 / prior.variance

      sigma.map[t, s] <- sigma.map[s, t] - d / dd
      sigma.map[s, t] <- sigma.map[s, t] - d / dd
    }
  }
  
  return(sigma.map)
}

invert.hessian <- function(sigma.map, index, x, prior.variance) {
  # The Hessian is build up as 
  # H = (A   , B)
  #     (B^T , D)
  # Where A contains the derivatives w.r.t the main effects.
  # Where D contains the derivatives w.r.t the interaction effects.
  
  # Since inverse of A is cheap, we use schur complement to invert H.
  #  schur <- solve(D - t(B) %*% iA %*% B)
  # See below. 
  
  #This code can be improved with better memory allocation. Now, everything is stored twice.
  #This code can be improved with better inversion methods. Cholesky?
  #This code can be improved with more efficient computing of matrix D.
  #This code can be improved with implementation in C for speed and accuracy / precision.
  
  p <- ncol(x)
  n <- nrow(x)
  
  #P(x_{is} = 1 | x[i,-s]) * P(x_{is} = 0 | x[i,-s])
  pq <- matrix(0, n, p)
  for(s in 1:p) {
    phi <- sigma.map[s, s] + x[, -s] %*% sigma.map[-s, s]
    pq[, s] <- exp(phi) / (1 + exp(phi)) ^ 2
  }
  
  #iA is inverse of A
  iA <- matrix(0, nrow = p, ncol = p)
  for(s in 1:p)
    iA[s, s] <- 1 / (-sum(pq[, s]) - 1 / prior.variance)
  
  # The matrix D contains derivatives w.r.t. interactions. 
  # The interactions are indexed with the matrix index.
  
  D <- matrix(0, p * (p - 1) / 2, p * (p-1) / 2)
  for(row in 1:(p * (p - 1) / 2)) {
    s <- index[row, 2]
    t <- index[row, 3]
    
    #d ^ 2 / d sigma.map[s, t] ^ 2               #with (s < t)
    D[row, row] <- -x[, t] %*% pq[, s] - x[, s] %*% pq[, t] - 1 / prior.variance
    
    #d ^ 2 / d sigma.map[s, t] / d sigma.map[s, r]   #with (s < t)
    I <- which(index[, 2] == s & index[, 3] != t)
    for(i in I) {
      r <- index[i, 3]
      col <- index[i, 1]
      D[row, col] <- -(x[, t] * x[, r]) %*% pq[, s]
    }
    #d^2 / d sigma.map[s, t] / d sigma.map[r, s]   #with (s < t)
    I <- which(index[, 3] == s)
    for(i in I) {
      r <- index[i, 2]
      col <- index[i, 1]
      D[row, col] <- -(x[, t] * x[, r]) %*% pq[, s]
    }
    
    #d^2 / d sigma.map[s, t] / d sigma.map[c, t]   #with (s < t)
    I <- which(index[, 3] == t & index[, 2] != s)
    for(i in I) {
      c <- index[i, 2]
      col <- index[i, 1]
      D[row, col] <- -(x[, s] * x[, c]) %*% pq[, t]
    }
    #d^2 / d sigma.map[s, t] / d sigma.map[t, c]   #with (s < t)
    I <- which(index[, 2] == t)
    for(i in I) {
      c <- index[i, 3]
      col <- index[i, 1]
      D[row, col] <- -(x[, s] * x[, c]) %*% pq[, t]
    }
  }
  
  B <- matrix(0, nrow = p, ncol = p * (p - 1) / 2)
  for(s in 1:p) {
    I <- which(index[, 2] == s)
    for(i in I) {
      col <- index[i, 1]
      t <- index[i, 3]
      B[s, col] <- - x[, t] %*% pq[, s]
    }
    
    I <- which(index[, 3] == s)
    for(i in I) {
      col <- index[i, 1]
      t <- index[i, 2]
      B[s, col] <- - x[, t] %*% pq[, s]
    }
  }
  
  #now compute inverse
  inverse.hessian <- matrix(data = 0, nrow = p * (p + 1) / 2, ncol = p * (p + 1) / 2)
  
  M <- 1:p                       #index for main effects
  I <- (1:(p * (p + 1) / 2))[-M] #index for interaction effects
  
  #compute schur complement; schur = solve (D - t(B) %*% iA %*% B)
  inverse.hessian[I, I] <- D - t(B) %*% iA %*% B
  inverse.hessian[I, I] <- solve(inverse.hessian[I, I]) 
  #-iA %*% B %*% schur
  inverse.hessian[M, I] <- -iA %*% B %*% inverse.hessian[I, I]
  #-schur %*% t(B) %*% iA
  inverse.hessian[I, M] <- t(inverse.hessian[M, I])
  #iA + iA %*% B %*% schur %*% t(B) %*% iA
  inverse.hessian[M, M] <- iA - inverse.hessian[M, I] %*% t(B) %*% iA  
  
  #One can also do it numerically:
  #  inverse.hessian[M, M] <- solve(iA)
  #  inverse.hessian[M, I] <- B
  #  inverse.hessian[I, M] <- t(B)
  #  inverse.hessian[I, I] <- D
  #  inverse.hessian <- solve(inverse.hessian)
  
  return(inverse.hessian)
}

multi.dimensional.NR <- function(sigma.map, index, x, suf.stat, prior.variance) {
  n <- nrow(x)
  p <- ncol(x)
  
  delta <- eta <- vector(length = p * (p + 1) / 2)
  prob <- matrix(0, nrow = n, ncol = p)
  
  #Compute D, the vector of first-order derivatives
  for(s in 1:p) {
    phi <- sigma.map[s, s] + x[, -s] %*% sigma.map[-s, s]
    prob[, s] <- exp(phi) / (1 + exp(phi))
    delta[s] <- suf.stat[s, s] - sum(prob[, s]) - sigma.map[s, s] / prior.variance
  }
  for(s in 1:(p - 1)) {
    for(t in (s + 1): p) {
      row <- p + index[index[, 2] == s & index[, 3] == t, 1]
      delta[row] <- 2 * suf.stat[s, t] - x[, t] %*% prob[, s] - x[, s] %*% prob[, t] - sigma.map[s, t] / prior.variance
    }
  }
  
  #Compute inverse Hessian
  inverse.hessian <- invert.hessian(sigma.map = sigma.map, index = index, x = x, prior.variance = prior.variance)
  
  #Assign eta values
  eta[1:p] <- diag(sigma.map)
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      row <- index[index[,2] == s & index[, 3] == t, 1] + p
      eta[row] <- sigma.map[s, t]
    }
  }
  
  #Newton - Raphson step
  eta <- eta - inverse.hessian %*% delta
  
  #Assign sigma values
  diag(sigma.map) <- eta[1:p]
  for(s in 1:(p-1)) {
    for(t in (s+1):p) {
      row <- p + index[index[,2] == s & index[, 3] == t, 1]
      sigma.map[s, t] <- eta[row]
      sigma.map[t, s] <- eta[row]
    }
  }
  return(sigma.map)
}

proportional.log.pseudoposterior <- function(sigma.map, x, prior.variance) {
  p <- ncol(x)
  q <- 0
  for(s in 1:p) {
    phi <- sigma.map[s, s] + x[, -s] %*% sigma.map[-s, s]
    q <- q + x[, s] %*% phi - sum(log(1 + exp(phi)))
    if(is.finite(prior.variance))
      q <- q + dnorm(x = sigma.map[s, s], sd = sqrt(prior.variance), log = TRUE)
  }
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      if(is.finite(prior.variance))
        q <- q + dnorm(x = sigma.map[s, t], sd = sqrt(prior.variance), log = TRUE)
    }
  }
  return(q)
}

invert.hessian.emvs <- function(sigma, 
                                gamma, 
                                index, 
                                x, 
                                prior.variance.intercepts, 
                                nu0, 
                                nu1) {
  # The Hessian is build up as 
  # H = (A   , B)
  #     (B^T , D)
  # Where A contains the derivatives w.r.t the main effects.
  # Where D contains the derivatives w.r.t the interaction effects.
  
  # Since inverse of A is cheap, we use schur complement to invert H.
  #  schur <- solve(D - t(B) %*% iA %*% B)
  # See below. 
  
  #This code can be improved with better memory allocation. Now, everything is stored twice.
  #This code can be improved with better inversion methods. Cholesky?
  #This code can be improved with more efficient computing of matrix D.
  #This code can be improved with implementation in C for speed and accuracy / precision.
  
  p <- ncol(x)
  n <- nrow(x)
  
  #P(x_{is} = 1 | x[i,-s]) * P(x_{is} = 0 | x[i,-s])
  pq <- matrix(data = 0, nrow = n, ncol = p)
  for(s in 1:p) {
    phi <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
    pq[, s] <- exp(phi) / (1 + exp(phi)) ^ 2
  }
  
  #iA is inverse of A: derivatives of intercepts
  iA <- matrix(data = 0, nrow = p, ncol = p)
  for(s in 1:p)
    iA[s, s] <- 1 / (-sum(pq[, s]) - 1 / prior.variance.intercepts)
  
  # The matrix D contains derivatives w.r.t. interactions. 
  # The interactions are indexed with the matrix index.
  
  D <- matrix(0, p * (p - 1) / 2, p * (p-1) / 2)
  for(row in 1:(p * (p - 1) / 2)) {
    s <- index[row, 2]
    t <- index[row, 3]
    
    E <- gamma[s, t] / nu1[s, t] + (1 - gamma[s, t]) / nu0[s, t]
    
    #d^2 / d sigma.emvs[s, t] ^ 2               #with (s < t)
    D[row, row] <- -x[, t] %*% pq[, s] - x[, s] %*% pq[, t] - E
    
    #d^2 / d sigma.map[s, t] / d sigma.map[s, r]   #with (s < t)
    I <- which(index[, 2] == s & index[, 3] != t)
    for(i in I) {
      r <- index[i, 3]
      col <- index[i, 1]
      D[row, col] <- -(x[, t] * x[, r]) %*% pq[, s]
    }
    #d^2 / d sigma.map[s, t] / d sigma.map[r, s]   #with (s < t)
    I <- which(index[, 3] == s)
    for(i in I) {
      r <- index[i, 2]
      col <- index[i, 1]
      D[row, col] <- -(x[, t] * x[, r]) %*% pq[, s]
    }
    
    #d^2 / d sigma.map[s, t] / d sigma.map[c, t]   #with (s < t)
    I <- which(index[, 3] == t & index[, 2] != s)
    for(i in I) {
      c <- index[i, 2]
      col <- index[i, 1]
      D[row, col] <- -(x[, s] * x[, c]) %*% pq[, t]
    }
    #d^2 / d sigma.map[s, t] / d sigma.map[t, c]   #with (s < t)
    I <- which(index[, 2] == t)
    for(i in I) {
      c <- index[i, 3]
      col <- index[i, 1]
      D[row, col] <- -(x[, s] * x[, c]) %*% pq[, t]
    }
  }
  
  B <- matrix(0, nrow = p, ncol = p * (p - 1) / 2)
  for(s in 1:p) {
    I <- which(index[, 2] == s)
    for(i in I) {
      col <- index[i, 1]
      t <- index[i, 3]
      B[s, col] <- - x[, t] %*% pq[, s]
    }
    
    I <- which(index[, 3] == s)
    for(i in I) {
      col <- index[i, 1]
      t <- index[i, 2]
      B[s, col] <- - x[, t] %*% pq[, s]
    }
  }
  
  #now compute inverse
  inverse.hessian <- matrix(data = 0, nrow = p * (p + 1) / 2, ncol = p * (p + 1) / 2)
  
  M <- 1:p                       #index for main effects
  I <- (1:(p * (p + 1) / 2))[-M] #index for interaction effects
  
  #compute schur complement; schur = solve (D - t(B) %*% iA %*% B)
  inverse.hessian[I, I] <- D - t(B) %*% iA %*% B
  inverse.hessian[I, I] <- solve(inverse.hessian[I, I])               
  #-iA %*% B %*% schur
  inverse.hessian[M, I] <- -iA %*% B %*% inverse.hessian[I, I]
  #-schur %*% t(B) %*% iA
  inverse.hessian[I, M] <- t(inverse.hessian[M, I])
  #iA + iA %*% B %*% schur %*% t(B) %*% iA
  inverse.hessian[M, M] <- iA - inverse.hessian[M, I] %*% t(B) %*% iA  
  
  #One can also do it numerically:
  #  inverse.hessian[M, M] <- solve(iA)
  #  inverse.hessian[M, I] <- B
  #  inverse.hessian[I, M] <- t(B)
  #  inverse.hessian[I, I] <- D
  #  inverse.hessian <- solve(inverse.hessian)
  
  return(inverse.hessian)
}



maximum.pseudoposterior <- function(x, prior.variance = 1) {
  p <- ncol(x)
  n <- nrow(x)
  
  sigma.map <- matrix(0, nrow = p, ncol = p)         #ising parameters
  suf.stat <- t(x) %*% x                             #pre-compute sufficient statistics

  #Compute startvalues
  for(it in 1:5)
    sigma.map <- one.dimensional.NR(sigma.map = sigma.map, 
                                    x = x, 
                                    suf.stat = suf.stat, 
                                    prior.variance = prior.variance)
  Q <- proportional.log.pseudoposterior(sigma.map = sigma.map, 
                                        x = x, 
                                        prior.variance = prior.variance)
  
  #The matrix "index" is needed to convert the matrix "sigma" to the vector 
  #"eta" during optimization.
  index <- matrix(0, nrow = p * (p - 1) / 2, ncol = 3)
  tel <- 0
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      tel <- tel + 1
      index[tel, 1] <- tel
      index[tel, 2] <- s
      index[tel, 3] <- t
    }
  } 

  for(it in 1:1e5) {
    #update ising parameters
    sigma.map <- multi.dimensional.NR(sigma.map = sigma.map, 
                                      index = index, 
                                      x = x, 
                                      suf.stat = suf.stat, 
                                      prior.variance = prior.variance)
    #compute proportional log pseudoposterior
    Q[it + 1] <- proportional.log.pseudoposterior(sigma.map = sigma.map, 
                                                  x = x, 
                                                  prior.variance = prior.variance)
    if(is.null(abs(Q[it+1] - Q[it])) | is.nan(abs(Q[it+1] - Q[it])) | abs(Q[it+1] - Q[it]) < sqrt(.Machine$double.eps))
      break
  }
  return(list(sigma = sigma.map, Q = Q))
}

