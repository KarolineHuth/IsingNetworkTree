### DISCRETE PARTITIONING VARIABLE
## WATCHOUT: Only combinations of n.obs and n.nds where n.obs is large enough
# -----------------------------------------------------------------------------------------

# SIMULATIONS:
### A - 1: Single delta rho change (i.e., one interaction changes)
### A - 2: Multiple delta rho changes (i.e., all interactions change)
### B: Single delta rho change & spurious partitioning variables (i.e., z.var > 1)
### C: False Positive: delta rho = 0 

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

#### SIMULATION A - 1: Single delta rho change (i.e., one interaction changes)

# Node - 5

cor <- c(.1, .3, .5)
obs <- c(100, 500, 1500)
id <- 0
sim <- 50
DSimA_s5 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
                                  part.type = "discrete")
      
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
      DSimA_s5$n.zvar[id] <- n.zvar
      DSimA_s5$n.nds[id] <- n.nds
      DSimA_s5$n.obs[id] <- n.obs
      DSimA_s5$delta.all[id] <- delta.all
      DSimA_s5$delta.cor[id] <- delta.cor
      DSimA_s5$mean.dif[id] <- mean.dif
      DSimA_s5$mean.change[id] <- mean.change
      DSimA_s5$split[id] <- split
    }
  }
}


write.csv(DSimA_s5, "DSimA_s5.csv")

# --------------------------------------------------------------------------------
# N.Nodes: 15

cor <- c(.1, .3, .5)
obs <- c(500, 1500)
id <- 0
sim <- 50
DSimA_s15 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
                                  part.type = "discrete")
      
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
      DSimA_s15$n.zvar[id] <- n.zvar
      DSimA_s15$n.nds[id] <- n.nds
      DSimA_s15$n.obs[id] <- n.obs
      DSimA_s15$delta.all[id] <- delta.all
      DSimA_s15$delta.cor[id] <- delta.cor
      DSimA_s15$mean.dif[id] <- mean.dif
      DSimA_s15$mean.change[id] <- mean.change
      DSimA_s15$split[id] <- split
      print(id)
    }
  }
}


write.csv(DSimA_s15, "DSimA_s15.csv")



# -----------------------------------------------------------------------------------------

# N.Nodes: 25 
cor <- c(.1, .3, .5)
id <- 0
sim <- 50
DSimA_s25 <- data.frame(id = numeric(sim*length(cor)))


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
                                part.type = "discrete")
    
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
    DSimA_s25$n.zvar[id] <- n.zvar
    DSimA_s25$n.nds[id] <- n.nds
    DSimA_s25$n.obs[id] <- n.obs
    DSimA_s25$delta.all[id] <- delta.all
    DSimA_s25$delta.cor[id] <- delta.cor
    DSimA_s25$mean.dif[id] <- mean.dif
    DSimA_s25$mean.change[id] <- mean.change
    DSimA_s25$split[id] <- split
    print(id)
  }
}

write.csv(DSimA_s25, "DSimA_s25.csv")
DSimA_s <- rbind(DSimA_s5, DSimA_s15, DSimA_s25)
write.csv(DSimA_s, "DSimA_s.csv")

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

#### SIMULATION A - 2: Multiple delta rho changes (i.e., all interactions change)

# N.Nodes: 5
cor <- c(.1, .3, .5)
obs <- c(100, 500, 1500)
id <- 0
sim <- 50
DSimA_m5 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
                                  part.type = "discrete")
      
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
      DSimA_m5$n.zvar[id] <- n.zvar
      DSimA_m5$n.nds[id] <- n.nds
      DSimA_m5$n.obs[id] <- n.obs
      DSimA_m5$delta.all[id] <- delta.all
      DSimA_m5$delta.cor[id] <- delta.cor
      DSimA_m5$mean.dif[id] <- mean.dif
      DSimA_m5$mean.change[id] <- mean.change
      DSimA_m5$split[id] <- split
      print(id)
    }
  }
}


write.csv(DSimA_m5, "DSimA_m5.csv")

#-----------------------------------------------------------------------------

# N.Nodes: 15
cor <- c(.1, .3, .5)
obs <- c(500, 1500)
id <- 0
sim <- 50
DSimA_m15 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
                                  part.type = "discrete")
      
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
      DSimA_m15$n.zvar[id] <- n.zvar
      DSimA_m15$n.nds[id] <- n.nds
      DSimA_m15$n.obs[id] <- n.obs
      DSimA_m15$delta.all[id] <- delta.all
      DSimA_m15$delta.cor[id] <- delta.cor
      DSimA_m15$mean.dif[id] <- mean.dif
      DSimA_m15$mean.change[id] <- mean.change
      DSimA_m15$split[id] <- split
      print(id)
    }
  }
}

write.csv(DSimA_m15, "DSimA_m15.csv")

#-----------------------------------------------------------------------------

# N.Nodes: 25
cor <- c(.1, .3, .5)
id <- 0
sim <- 50
DSimA_m25 <- data.frame(id = numeric(sim*length(cor)))


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
                                part.type = "discrete")
    
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
    DSimA_m25$n.zvar[id] <- n.zvar
    DSimA_m25$n.nds[id] <- n.nds
    DSimA_m25$n.obs[id] <- n.obs
    DSimA_m25$delta.all[id] <- delta.all
    DSimA_m25$delta.cor[id] <- delta.cor
    DSimA_m25$mean.dif[id] <- mean.dif
    DSimA_m25$mean.change[id] <- mean.change
    DSimA_m25$split[id] <- split
    print(id)
  }
}


write.csv(DSimA_m25, "DSimA_m25.csv")

DSimA_m <- rbind(DSimA_m5, DSimA_m15, DSimA_m25)
write.csv(DSimA_m, "DSimA_m.csv")
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------


#### SIMULATION B: Single delta rho change & spurious partitioning variables (i.e., z.var > 1)


# N.Nodes: 5
cor <- c(.1, .3, .5)
obs <- c(100, 500, 1500)
id <- 0
sim <- 50
DSimB_s5 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
                                  part.type = "discrete")
      
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
      DSimB_s5$n.zvar[id] <- n.zvar
      DSimB_s5$n.nds[id] <- n.nds
      DSimB_s5$n.obs[id] <- n.obs
      DSimB_s5$delta.all[id] <- delta.all
      DSimB_s5$delta.cor[id] <- delta.cor
      DSimB_s5$mean.dif[id] <- mean.dif
      DSimB_s5$mean.change[id] <- mean.change
      DSimB_s5$split[id] <- split
      print(id)
    }
  }
}

write.csv(DSimB_s5, "DSimB_s5.csv")

#---------------------------------------------------------------------------------------

# N.Nodes: 15
cor <- c(.1, .3, .5)
obs <- c(500, 1500)
id <- 0
sim <- 50
DSimB_s15 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
                                  part.type = "discrete")
      
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
      DSimB_s15$n.zvar[id] <- n.zvar
      DSimB_s15$n.nds[id] <- n.nds
      DSimB_s15$n.obs[id] <- n.obs
      DSimB_s15$delta.all[id] <- delta.all
      DSimB_s15$delta.cor[id] <- delta.cor
      DSimB_s15$mean.dif[id] <- mean.dif
      DSimB_s15$mean.change[id] <- mean.change
      DSimB_s15$split[id] <- split
      print(id)
    }
  }
}

write.csv(DSimB_s15, "DSimB_s15.csv")

# -------------------------------------------------------------------------

# N.Nodes: 25
cor <- c(.1, .3, .5)
id <- 0
sim <- 50
DSimB_s25 <- data.frame(id = numeric(sim*length(cor)))


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
                                part.type = "discrete")
    
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
    DSimB_s25$n.zvar[id] <- n.zvar
    DSimB_s25$n.nds[id] <- n.nds
    DSimB_s25$n.obs[id] <- n.obs
    DSimB_s25$delta.all[id] <- delta.all
    DSimB_s25$delta.cor[id] <- delta.cor
    DSimB_s25$mean.dif[id] <- mean.dif
    DSimB_s25$mean.change[id] <- mean.change
    DSimB_s25$split[id] <- split
    print(id)
  }
}


write.csv(DSimB_s25, "DSimB_s25.csv")

DSimB_s <- rbind(DSimB_s5, DSimB_s15, DSimB_s25)
write.csv(DSimB_s, "DSimB_s.csv")
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------


#### SIMULATION C: False Positive: delta rho = 0 

# N.Nodes: 5
obs <- c(100, 500, 1500)
id <- 0
sim <- 50
DSimC_5 <- data.frame(id = numeric(sim*length(obs)))


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
                                part.type = "discrete")
    
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
    DSimC_5$n.zvar[id] <- n.zvar
    DSimC_5$n.nds[id] <- n.nds
    DSimC_5$n.obs[id] <- n.obs
    DSimC_5$delta.all[id] <- delta.all
    DSimC_5$delta.cor[id] <- delta.cor
    DSimC_5$mean.dif[id] <- mean.dif
    DSimC_5$mean.change[id] <- mean.change
    DSimC_5$split[id] <- split
  }
}

write.csv(DSimC_5, "DSimC_5.csv")

# --------------------------------------------------------------------------------

# Simulation C: 
# N.Nodes: 15
obs <- c(500, 1500)
id <- 0
sim <- 50
DSimC_15 <- data.frame(id = numeric(sim*length(obs)))


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
                                part.type = "discrete")
    
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
    DSimC_15$n.zvar[id] <- n.zvar
    DSimC_15$n.nds[id] <- n.nds
    DSimC_15$n.obs[id] <- n.obs
    DSimC_15$delta.all[id] <- delta.all
    DSimC_15$delta.cor[id] <- delta.cor
    DSimC_15$mean.dif[id] <- mean.dif
    DSimC_15$mean.change[id] <- mean.change
    DSimC_15$split[id] <- split
  }
}


write.csv(DSimC_15, "DSimC_15.csv")

# --------------------------------------------------------------------------------

# Simulation C:
# N.Nodes: 25
id <- 0
sim <- 50
DSimC_25 <- data.frame(id = numeric(sim))


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
                              part.type = "discrete")
  
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
  DSimC_25$n.zvar[id] <- n.zvar
  DSimC_25$n.nds[id] <- n.nds
  DSimC_25$n.obs[id] <- n.obs
  DSimC_25$delta.all[id] <- delta.all
  DSimC_25$delta.cor[id] <- delta.cor
  DSimC_25$mean.dif[id] <- mean.dif
  DSimC_25$mean.change[id] <- mean.change
  DSimC_25$split[id] <- split
}


write.csv(DSimC_25,"DSimC_25.csv", row.names = FALSE)

DSimC <- rbind(DSimC_5, DSimC_15, DSimC_25)

write.csv(DSimC,"DSimC.csv", row.names = FALSE)

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
