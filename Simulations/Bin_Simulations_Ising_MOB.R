### BINARY PARTITIONING VARIABLE
## WATCHOUT: Only combinations of n.obs and n.nds where n.obs is large enough
# -----------------------------------------------------------------------------------------

# SIMULATIONS:
### A - 1: Single delta rho change (i.e., one interaction changes)
### A - 2: Multiple delta rho changes (i.e., all interactions change)
### B: Single delta rho change & spurious partitioning variables (i.e., z.var > 1)
### C: False Positive: delta rho = 0 
### D: False Positive: delta rho = 0 & delta mu neq 0

# -----------------------------------------------------------------------------------------

#### STEP 0: Source Functions
#source("~/Documents/Studium/Master_UniversityAmsterdam/Masterthesis/03_Programming/Ising_Fitting_MOB.R")
#source("~/Documents/Studium/Master_UniversityAmsterdam/Masterthesis/03_Programming/Sources_Ising/B Ising - EMVS - Functions.R")

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

#### SIMULATION A - 1: Single delta rho change (i.e., one interaction changes)

# Node - 5

cor <- c(.1, .3, .5)
obs <- c(100, 500, 1500)
id <- 0
sim <- 50
BSimA_s5 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


# Simulate different amount of correlation changes
for(o in 1:length(obs)){
  for(d in 1:length(cor)){
    for (iter in 1:sim){
      # Keep Track of the ID
      id <- id + 1
      
      # Determine Simulation variables
      n.zvar <- 1
      n.nds <- 5
      n.obs <- obs[o]
      delta.all <- FALSE
      delta.cor <- cor[d]
      mean.dif <- FALSE
      mean.change <- .6
      
      #Simualte Data
      data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                  delta.cor = delta.cor, delta.all = delta.all, 
                                  mean.dif = mean.dif, mean.change = mean.change, 
                                  part.type = "binary")
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, 
                               model = c("correlation"))
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      BSimA_s5$n.zvar[id] <- n.zvar
      BSimA_s5$n.nds[id] <- n.nds
      BSimA_s5$n.obs[id] <- n.obs
      BSimA_s5$delta.all[id] <- delta.all
      BSimA_s5$delta.cor[id] <- delta.cor
      BSimA_s5$mean.dif[id] <- mean.dif
      BSimA_s5$mean.change[id] <- mean.change
      BSimA_s5$split[id] <- split
    }
  }
}


write.csv(BSimA_s5, "BSimA_s5.csv")

# --------------------------------------------------------------------------------
# N.Nodes: 15

cor <- c(.1, .3, .5)
obs <- c(500, 1500)
id <- 0
sim <- 50
BSimA_s15 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


# Simulate different amount of correlation changes
for(o in 1:length(obs)){
  for(d in 1:length(cor)){
    for (iter in 1:sim){
      # Keep Track of the ID
      id <- id + 1
      
      # Determine Simulation variables
      n.zvar <- 1
      n.nds <- 15
      n.obs <- obs[o]
      delta.all <- FALSE
      delta.cor <- cor[d]
      mean.dif <- FALSE
      mean.change <- .6
      
      #Simualte Data
      data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                  delta.cor = delta.cor, delta.all = delta.all, 
                                  mean.dif = mean.dif, mean.change = mean.change, 
                                  part.type = "binary")
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars,
                               model = c("correlation"))
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      BSimA_s15$n.zvar[id] <- n.zvar
      BSimA_s15$n.nds[id] <- n.nds
      BSimA_s15$n.obs[id] <- n.obs
      BSimA_s15$delta.all[id] <- delta.all
      BSimA_s15$delta.cor[id] <- delta.cor
      BSimA_s15$mean.dif[id] <- mean.dif
      BSimA_s15$mean.change[id] <- mean.change
      BSimA_s15$split[id] <- split
    }
  }
}


write.csv(BSimA_s15, "BSimA_s15.csv")



# -----------------------------------------------------------------------------------------

# N.Nodes: 25 
cor <- c(.1, .3, .5)
id <- 0
sim <- 50
BSimA_s25 <- data.frame(id = numeric(sim*length(cor)))


# Simulate different amount of correlation changes

for(d in 1:length(cor)){
  for (iter in 1:sim){
    # Keep Track of the ID
    id <- id + 1
    
    # Determine Simulation variables
    n.zvar <- 1
    n.nds <- 25
    n.obs <- 1500
    delta.all <- FALSE
    delta.cor <- cor[d]
    mean.dif <- FALSE
    mean.change <- .6
    
    #Simualte Data
    data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                delta.cor = delta.cor, delta.all = delta.all, 
                                mean.dif = mean.dif, mean.change = mean.change, 
                                part.type = "binary")
    
    # Run Ising MOB split function
    nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
    splitvars <- as.data.frame(data[, 1:n.zvar])
    tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
    
    # Extract result
    res <- unlist(tree)
    if(is.null(res$node.kids.id)){
      split <- "no"
    } else {
      split <- "yes"
    }
    
    # Write Result to Datafile
    BSimA_s25$n.zvar[id] <- n.zvar
    BSimA_s25$n.nds[id] <- n.nds
    BSimA_s25$n.obs[id] <- n.obs
    BSimA_s25$delta.all[id] <- delta.all
    BSimA_s25$delta.cor[id] <- delta.cor
    BSimA_s25$mean.dif[id] <- mean.dif
    BSimA_s25$mean.change[id] <- mean.change
    BSimA_s25$split[id] <- split
    print(id)
  }
}


BSimA_s <- rbind(BSimA_s5, BSimA_s15, BSimA_s25)
write.csv(BSimA_s, "BSimA_s.csv")

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

#### SIMULATION A - 2: Multiple delta rho changes (i.e., all interactions change)

# N.Nodes: 5
cor <- c(.1, .3, .5)
obs <- c(100, 500, 1500)
id <- 0
sim <- 50
BSimA_m5 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


# Simulate different amount of correlation changes
for(o in 1:length(obs)){
  for(d in 1:length(cor)){
    for (iter in 1:sim){
      # Keep Track of the ID
      id <- id + 1
      
      # Determine Simulation variables
      n.zvar <- 1
      n.nds <- 5
      n.obs <- obs[o]
      delta.all <- TRUE
      delta.cor <- cor[d]
      mean.dif <- FALSE
      mean.change <- .6
      
      #Simualte Data
      data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                  delta.cor = delta.cor, delta.all = delta.all, 
                                  mean.dif = mean.dif, mean.change = mean.change, 
                                  part.type = "binary")
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      BSimA_m5$n.zvar[id] <- n.zvar
      BSimA_m5$n.nds[id] <- n.nds
      BSimA_m5$n.obs[id] <- n.obs
      BSimA_m5$delta.all[id] <- delta.all
      BSimA_m5$delta.cor[id] <- delta.cor
      BSimA_m5$mean.dif[id] <- mean.dif
      BSimA_m5$mean.change[id] <- mean.change
      BSimA_m5$split[id] <- split
    }
  }
}


write.csv(BSimA_m5, "BSimA_m5.csv")

#-----------------------------------------------------------------------------

# N.Nodes: 15
cor <- c(.1, .3, .5)
obs <- c(500, 1500)
id <- 0
sim <- 50
BSimA_m15 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


# Simulate different amount of correlation changes
for(o in 1:length(obs)){
  for(d in 1:length(cor)){
    for (iter in 1:sim){
      # Keep Track of the ID
      id <- id + 1
      
      # Determine Simulation variables
      n.zvar <- 1
      n.nds <- 15
      n.obs <- obs[o]
      delta.all <- TRUE
      delta.cor <- cor[d]
      mean.dif <- FALSE
      mean.change <- .6
      
      #Simualte Data
      data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                  delta.cor = delta.cor, delta.all = delta.all, 
                                  mean.dif = mean.dif, mean.change = mean.change, 
                                  part.type = "binary")
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      BSimA_m15$n.zvar[id] <- n.zvar
      BSimA_m15$n.nds[id] <- n.nds
      BSimA_m15$n.obs[id] <- n.obs
      BSimA_m15$delta.all[id] <- delta.all
      BSimA_m15$delta.cor[id] <- delta.cor
      BSimA_m15$mean.dif[id] <- mean.dif
      BSimA_m15$mean.change[id] <- mean.change
      BSimA_m15$split[id] <- split
      print(id)
    }
  }
}

write.csv(BSimA_m15, "BSimA_m15.csv")

#-----------------------------------------------------------------------------

# N.Nodes: 25
cor <- c(.1, .3, .5)
id <- 0
sim <- 50
BSimA_m25 <- data.frame(id = numeric(sim*length(cor)))


# Simulate different amount of correlation changes

for(d in 1:length(cor)){
  for (iter in 1:sim){
    # Keep Track of the ID
    id <- id + 1
    
    # Determine Simulation variables
    n.zvar <- 1
    n.nds <- 25
    n.obs <- 1500
    delta.all <- TRUE
    delta.cor <- cor[d]
    mean.dif <- FALSE
    mean.change <- .6
    
    #Simualte Data
    data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                delta.cor = delta.cor, delta.all = delta.all, 
                                mean.dif = mean.dif, mean.change = mean.change, 
                                part.type = "binary")
    
    # Run Ising MOB split function
    nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
    splitvars <- as.data.frame(data[, 1:n.zvar])
    tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
    
    # Extract result
    res <- unlist(tree)
    if(is.null(res$node.kids.id)){
      split <- "no"
    } else {
      split <- "yes"
    }
    
    # Write Result to Datafile
    BSimA_m25$n.zvar[id] <- n.zvar
    BSimA_m25$n.nds[id] <- n.nds
    BSimA_m25$n.obs[id] <- n.obs
    BSimA_m25$delta.all[id] <- delta.all
    BSimA_m25$delta.cor[id] <- delta.cor
    BSimA_m25$mean.dif[id] <- mean.dif
    BSimA_m25$mean.change[id] <- mean.change
    BSimA_m25$split[id] <- split
    print(id)
  }
}


write.csv(BSimA_m25, "BSimA_m25.csv")

BSimA_m <- rbind(BSimA_m5, BSimA_m15, BSimA_m25)
write.csv(BSimA_m, "BSimA_m.csv")
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------


#### SIMULATION B: Single delta rho change & spurious partitioning variables (i.e., z.var > 1)


# N.Nodes: 5
cor <- c(.1, .3, .5)
obs <- c(100, 500, 1500)
id <- 0
sim <- 50
BSimB_s5 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


# Simulate different amount of correlation changes
for(o in 1:length(obs)){
  for(d in 1:length(cor)){
    for (iter in 1:sim){
      # Keep Track of the ID
      id <- id + 1
      
      # Determine Simulation variables
      n.zvar <- 5
      n.nds <- 5
      n.obs <- obs[o]
      delta.all <- FALSE
      delta.cor <- cor[d]
      mean.dif <- FALSE
      mean.change <- .6
      
      #Simualte Data
      data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                  delta.cor = delta.cor, delta.all = delta.all, 
                                  mean.dif = mean.dif, mean.change = mean.change, 
                                  part.type = "binary")
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      BSimB_s5$n.zvar[id] <- n.zvar
      BSimB_s5$n.nds[id] <- n.nds
      BSimB_s5$n.obs[id] <- n.obs
      BSimB_s5$delta.all[id] <- delta.all
      BSimB_s5$delta.cor[id] <- delta.cor
      BSimB_s5$mean.dif[id] <- mean.dif
      BSimB_s5$mean.change[id] <- mean.change
      BSimB_s5$split[id] <- split
      
    }
  }
}

write.csv(BSimB_s5, "BSimB_s5.csv")

#---------------------------------------------------------------------------------------

# N.Nodes: 15
cor <- c(.1, .3, .5)
obs <- c(500, 1500)
id <- 0
sim <- 50
BSimB_s15 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


# Simulate different amount of correlation changes
for(o in 1:length(obs)){
  for(d in 1:length(cor)){
    for (iter in 1:sim){
      # Keep Track of the ID
      id <- id + 1
      
      # Determine Simulation variables
      n.zvar <- 5
      n.nds <- 15
      n.obs <- obs[o]
      delta.all <- FALSE
      delta.cor <- cor[d]
      mean.dif <- FALSE
      mean.change <- .6
      
      #Simualte Data
      data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                  delta.cor = delta.cor, delta.all = delta.all, 
                                  mean.dif = mean.dif, mean.change = mean.change, 
                                  part.type = "binary")
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      BSimB_s15$n.zvar[id] <- n.zvar
      BSimB_s15$n.nds[id] <- n.nds
      BSimB_s15$n.obs[id] <- n.obs
      BSimB_s15$delta.all[id] <- delta.all
      BSimB_s15$delta.cor[id] <- delta.cor
      BSimB_s15$mean.dif[id] <- mean.dif
      BSimB_s15$mean.change[id] <- mean.change
      BSimB_s15$split[id] <- split
      print(id)
    }
  }
}

write.csv(BSimB_s15, "BSimB_s15.csv")

# -------------------------------------------------------------------------

# N.Nodes: 25
cor <- c(.1, .3, .5)
id <- 0
sim <- 50
BSimB_s25 <- data.frame(id = numeric(sim*length(cor)))


# Simulate different amount of correlation changes
for(d in 1:length(cor)){
  for (iter in 1:sim){
    # Keep Track of the ID
    id <- id + 1
    
    # Determine Simulation variables
    n.zvar <- 5
    n.nds <- 25
    n.obs <- 1500
    delta.all <- FALSE
    delta.cor <- cor[d]
    mean.dif <- FALSE
    mean.change <- .6
    
    #Simualte Data
    data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                delta.cor = delta.cor, delta.all = delta.all, 
                                mean.dif = mean.dif, mean.change = mean.change, 
                                part.type = "binary")
    
    # Run Ising MOB split function
    nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
    splitvars <- as.data.frame(data[, 1:n.zvar])
    tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
    
    # Extract result
    res <- unlist(tree)
    if(is.null(res$node.kids.id)){
      split <- "no"
    } else {
      split <- "yes"
    }
    
    # Write Result to Datafile
    BSimB_s25$n.zvar[id] <- n.zvar
    BSimB_s25$n.nds[id] <- n.nds
    BSimB_s25$n.obs[id] <- n.obs
    BSimB_s25$delta.all[id] <- delta.all
    BSimB_s25$delta.cor[id] <- delta.cor
    BSimB_s25$mean.dif[id] <- mean.dif
    BSimB_s25$mean.change[id] <- mean.change
    BSimB_s25$split[id] <- split
    print(id)
  }
}


write.csv(BSimB_s25, "BSimB_s25.csv")

BSimB_s <- rbind(BSimB_s5, BSimB_s15, BSimB_s25)
write.csv(BSimB_s, "BSimB_s.csv")
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------


#### SIMULATION C: False Positive: delta rho = 0 

# N.Nodes: 5
obs <- c(100, 500, 1500)
id <- 0
sim <- 50
BSimC_5 <- data.frame(id = numeric(sim*length(obs)))


for(c in 1:length(obs)){
  for (it in 1:sim){
    # Keep Track of the ID
    id <- id + 1
    
    # Determine Simulation variables
    n.zvar <- 1
    n.nds <- 5
    n.obs <- obs[c]
    delta.all <- FALSE
    delta.cor <- 0
    mean.dif <- FALSE
    mean.change <- .6
    
    #Simualte Data
    data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                delta.cor = delta.cor, delta.all = delta.all, 
                                mean.dif = mean.dif, mean.change = mean.change, 
                                part.type = "binary")
    
    # Run Ising MOB split function
    nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
    splitvars <- as.data.frame(data[, 1:n.zvar])
    tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
    
    # Extract result
    res <- unlist(tree)
    if(is.null(res$node.kids.id)){
      split <- "no"
    } else {
      split <- "yes"
    }
    
    # Write Result to Datafile
    #    SimC_5$id[id] <- id 
    BSimC_5$n.zvar[id] <- n.zvar
    BSimC_5$n.nds[id] <- n.nds
    BSimC_5$n.obs[id] <- n.obs
    BSimC_5$delta.all[id] <- delta.all
    BSimC_5$delta.cor[id] <- delta.cor
    BSimC_5$mean.dif[id] <- mean.dif
    BSimC_5$mean.change[id] <- mean.change
    BSimC_5$split[id] <- split
  }
}

write.csv(BSimC_5, "BSimC_5.csv")

# --------------------------------------------------------------------------------

# Simulation C: 
# N.Nodes: 15
obs <- c(500, 1500)
id <- 0
sim <- 50
BSimC_15 <- data.frame(id = numeric(sim*length(obs)))


for(c in 1:length(obs)){
  for (it in 1:sim){
    # Keep Track of the ID
    id <- id + 1
    
    # Determine Simulation variables
    n.zvar <- 1
    n.nds <- 15
    n.obs <- obs[c]
    delta.all <- FALSE
    delta.cor <- 0
    mean.dif <- FALSE
    mean.change <- .6
    
    #Simualte Data
    data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                delta.cor = delta.cor, delta.all = delta.all, 
                                mean.dif = mean.dif, mean.change = mean.change, 
                                part.type = "binary")
    
    # Run Ising MOB split function
    nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
    splitvars <- as.data.frame(data[, 1:n.zvar])
    tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
    
    # Extract result
    res <- unlist(tree)
    if(is.null(res$node.kids.id)){
      split <- "no"
    } else {
      split <- "yes"
    }
    
    # Write Result to Datafile
    #    SimC_5$id[id] <- id 
    BSimC_15$n.zvar[id] <- n.zvar
    BSimC_15$n.nds[id] <- n.nds
    BSimC_15$n.obs[id] <- n.obs
    BSimC_15$delta.all[id] <- delta.all
    BSimC_15$delta.cor[id] <- delta.cor
    BSimC_15$mean.dif[id] <- mean.dif
    BSimC_15$mean.change[id] <- mean.change
    BSimC_15$split[id] <- split
  }
}


write.csv(BSimC_15, "BSimC_15.csv")

# --------------------------------------------------------------------------------

# Simulation C:
# N.Nodes: 25
id <- 0
sim <- 50
BSimC_25 <- data.frame(id = numeric(sim))


for (it in 1:sim){
  # Keep Track of the ID
  id <- id + 1
  
  # Determine Simulation variables
  n.zvar <- 1
  n.nds <- 25
  n.obs <- 1500
  delta.all <- FALSE
  delta.cor <- 0
  mean.dif <- FALSE
  mean.change <- .6
  
  #Simualte Data
  data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                              delta.cor = delta.cor, delta.all = delta.all, 
                              mean.dif = mean.dif, mean.change = mean.change, 
                              part.type = "binary")
  
  # Run Ising MOB split function
  nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
  splitvars <- as.data.frame(data[, 1:n.zvar])
  tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
  
  # Extract result
  res <- unlist(tree)
  if(is.null(res$node.kids.id)){
    split <- "no"
  } else {
    split <- "yes"
  }
  
  # Write Result to Datafile
  #    SimC_5$id[id] <- id 
  BSimC_25$n.zvar[id] <- n.zvar
  BSimC_25$n.nds[id] <- n.nds
  BSimC_25$n.obs[id] <- n.obs
  BSimC_25$delta.all[id] <- delta.all
  BSimC_25$delta.cor[id] <- delta.cor
  BSimC_25$mean.dif[id] <- mean.dif
  BSimC_25$mean.change[id] <- mean.change
  BSimC_25$split[id] <- split
}


write.csv(BSimC_25,"BSimC_25.csv", row.names = FALSE)

BSimC <- rbind(BSimC_5, BSimC_15, BSimC_25)

write.csv(BSimC,"BSimC.csv", row.names = FALSE)

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

