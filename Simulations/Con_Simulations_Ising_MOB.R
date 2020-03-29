### CONTINUOUS PARTITIONING VARIABLE
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
# -----------------------------------------------------------------------------------------

#### SIMULATION A - 1: Single delta rho change (i.e., one interaction changes)

# Node - 5

cor <- c(.1, .3, .5)
obs <- c(100, 500, 1500)
min <- c(35, 220, 720)
id <- 0
sim <- 50
CSimA_s5 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
                                  mean.dif = mean.dif, mean.change = mean.change)
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, 
                               model = c("correlation"), minsplit = min[o])
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      CSimA_s5$n.zvar[id] <- n.zvar
      CSimA_s5$n.nds[id] <- n.nds
      CSimA_s5$n.obs[id] <- n.obs
      CSimA_s5$delta.all[id] <- delta.all
      CSimA_s5$delta.cor[id] <- delta.cor
      CSimA_s5$mean.dif[id] <- mean.dif
      CSimA_s5$mean.change[id] <- mean.change
      CSimA_s5$split[id] <- split
    }
  }
}


write.csv(CSimA_s5, "CSimA_s5.csv")

# --------------------------------------------------------------------------------
# N.Nodes: 15

cor <- c(.1, .3, .5)
obs <- c(500, 1500)
min <- c(220, 720)
id <- 0
sim <- 50
CSimA_s15 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
                                  mean.dif = mean.dif, mean.change = mean.change)
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars,
                               model = c("correlation"), minsplit = min[o])
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      CSimA_s15$n.zvar[id] <- n.zvar
      CSimA_s15$n.nds[id] <- n.nds
      CSimA_s15$n.obs[id] <- n.obs
      CSimA_s15$delta.all[id] <- delta.all
      CSimA_s15$delta.cor[id] <- delta.cor
      CSimA_s15$mean.dif[id] <- mean.dif
      CSimA_s15$mean.change[id] <- mean.change
      CSimA_s15$split[id] <- split
    }
  }
}


write.csv(CSimA_s15, "CSimA_s15.csv")



# -----------------------------------------------------------------------------------------

# N.Nodes: 25 
cor <- c(.1, .3, .5)
id <- 0
sim <- 5
CSimA_s25 <- data.frame(id = numeric(sim*length(cor)))


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
                                mean.dif = mean.dif, mean.change = mean.change)
    
    # Run Ising MOB split function
    nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
    splitvars <- as.data.frame(data[, 1:n.zvar])
    tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars,
                             model = c("correlation"), minsplit = 720)
    
    # Extract result
    res <- unlist(tree)
    if(is.null(res$node.kids.id)){
      split <- "no"
    } else {
      split <- "yes"
    }
    
    # Write Result to Datafile
    CSimA_s25$n.zvar[id] <- n.zvar
    CSimA_s25$n.nds[id] <- n.nds
    CSimA_s25$n.obs[id] <- n.obs
    CSimA_s25$delta.all[id] <- delta.all
    CSimA_s25$delta.cor[id] <- delta.cor
    CSimA_s25$mean.dif[id] <- mean.dif
    CSimA_s25$mean.change[id] <- mean.change
    CSimA_s25$split[id] <- split
  }
}

write.csv(CSimA_s25, "CSimA_s25.csv")

CSimA_s <- rbind(CSimA_s5, CSimA_s15, CSimA_s25)
write.csv(CSimA_s, "CSimA_s.csv")

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

#### SIMULATION A - 2: Multiple delta rho changes (i.e., all interactions change)

# N.Nodes: 5
cor <- c(.1, .3, .5)
obs <- c(100, 500, 1500)
min <- c(35, 220, 720)
id <- 0
sim <- 50
CSimA_m5 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
                                  mean.dif = mean.dif, mean.change = mean.change)
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, 
                               model = c("correlation"), minsplit = min[o])
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      CSimA_m5$n.zvar[id] <- n.zvar
      CSimA_m5$n.nds[id] <- n.nds
      CSimA_m5$n.obs[id] <- n.obs
      CSimA_m5$delta.all[id] <- delta.all
      CSimA_m5$delta.cor[id] <- delta.cor
      CSimA_m5$mean.dif[id] <- mean.dif
      CSimA_m5$mean.change[id] <- mean.change
      CSimA_m5$split[id] <- split
      print(id)
    }
  }
}


write.csv(CSimA_m5, "CSimA_m5.csv")

#-----------------------------------------------------------------------------

# N.Nodes: 15
cor <- c(.1, .3, .5)
obs <- c(500, 1500)
min <- c(220, 720)
id <- 0
sim <- 50
CSimA_m15 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
      
      #Simulate Data
      data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                  delta.cor = delta.cor, delta.all = delta.all, 
                                  mean.dif = mean.dif, mean.change = mean.change)
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, 
                               model = c("correlation"), minsplit = min[o])
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      CSimA_m15$n.zvar[id] <- n.zvar
      CSimA_m15$n.nds[id] <- n.nds
      CSimA_m15$n.obs[id] <- n.obs
      CSimA_m15$delta.all[id] <- delta.all
      CSimA_m15$delta.cor[id] <- delta.cor
      CSimA_m15$mean.dif[id] <- mean.dif
      CSimA_m15$mean.change[id] <- mean.change
      CSimA_m15$split[id] <- split
      print(id)
    }
  }
}



write.csv(CSimA_m15, "CSimA_m15.csv")

#-----------------------------------------------------------------------------

# N.Nodes: 25
cor <- c(.1, .3, .5)
id <- 0
sim <- 50
CSimA_m25_extra <- data.frame(id = numeric(sim))


# Simulate different amount of correlation changes

#for(d in 1:length(cor)){
  for (iter in 1:sim){
    # Keep Track of the ID
    id <- id + 1
    
    # Determine Simulation variables
    n.zvar <- 1
    n.nds <- 25
    n.obs <- 1500
    delta.all <- TRUE
    delta.cor <- .5
    mean.dif <- FALSE
    mean.change <- .6
    
    #Simualte Data
    data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                delta.cor = delta.cor, delta.all = delta.all, 
                                mean.dif = mean.dif, mean.change = mean.change)
    
    # Run Ising MOB split function
    nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
    splitvars <- as.data.frame(data[, 1:n.zvar])
    tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars,
                             model = c("correlation"), minsplit = 720)
    
    # Extract result
    res <- unlist(tree)
    if(is.null(res$node.kids.id)){
      split <- "no"
    } else {
      split <- "yes"
    }
    
    # Write Result to Datafile
    CSimA_m25_extra$n.zvar[id] <- n.zvar
    CSimA_m25_extra$n.nds[id] <- n.nds
    CSimA_m25_extra$n.obs[id] <- n.obs
    CSimA_m25_extra$delta.all[id] <- delta.all
    CSimA_m25_extra$delta.cor[id] <- delta.cor
    CSimA_m25_extra$mean.dif[id] <- mean.dif
    CSimA_m25_extra$mean.change[id] <- mean.change
    CSimA_m25_extra$split[id] <- split
    print(id)
  }
#}

CSimA_m25[101:150, ] <- CSimA_m25_extra


write.csv(CSimA_m25, "CSimA_m25.csv")

CSimA_m <- rbind(CSimA_m5[, -1], CSimA_m15, CSimA_m25)
write.csv(CSimA_m, "CSimA_m.csv")

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------


#### SIMULATION B: Single delta rho change & spurious partitioning variables (i.e., z.var > 1)

# N.Nodes: 5
cor <- c(.1, .3, .5)
obs <- c(100, 500, 1500)
min <- c(35, 220, 720)
id <- 0
sim <- 50
CSimB_s5 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
                                  mean.dif = mean.dif, mean.change = mean.change)
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars,
                               model = c("correlation"), minsplit = min[o])
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      CSimB_s5$n.zvar[id] <- n.zvar
      CSimB_s5$n.nds[id] <- n.nds
      CSimB_s5$n.obs[id] <- n.obs
      CSimB_s5$delta.all[id] <- delta.all
      CSimB_s5$delta.cor[id] <- delta.cor
      CSimB_s5$mean.dif[id] <- mean.dif
      CSimB_s5$mean.change[id] <- mean.change
      CSimB_s5$split[id] <- split
      print(id)
    }
  }
}

write.csv(CSimB_s5, "CSimB_s5.csv")

#---------------------------------------------------------------------------------------

# N.Nodes: 15
cor <- c(.1, .3, .5)
obs <- c(500, 1500)
min <- c(220, 720)
id <- 0
sim <- 50
CSimB_s15 <- data.frame(id = numeric(sim*length(cor)*length(obs)))


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
                                  mean.dif = mean.dif, mean.change = mean.change)
      
      # Run Ising MOB split function
      nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
      splitvars <- as.data.frame(data[, 1:n.zvar])
      tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars,
                               model = c("correlation"), minsplit = min[o])
      
      # Extract result
      res <- unlist(tree)
      if(is.null(res$node.kids.id)){
        split <- "no"
      } else {
        split <- "yes"
      }
      
      # Write Result to Datafile
      CSimB_s15$n.zvar[id] <- n.zvar
      CSimB_s15$n.nds[id] <- n.nds
      CSimB_s15$n.obs[id] <- n.obs
      CSimB_s15$delta.all[id] <- delta.all
      CSimB_s15$delta.cor[id] <- delta.cor
      CSimB_s15$mean.dif[id] <- mean.dif
      CSimB_s15$mean.change[id] <- mean.change
      CSimB_s15$split[id] <- split
      print(id)
    }
  }
}

write.csv(CSimB_s15, "CSimB_s15.csv")

# -------------------------------------------------------------------------

# N.Nodes: 25
cor <- c(.1, .3, .5)
id <- 0
sim <- 50
CSimB_s25 <- data.frame(id = numeric(sim*length(cor)))


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
                                mean.dif = mean.dif, mean.change = mean.change)
    
    # Run Ising MOB split function
    nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
    splitvars <- as.data.frame(data[, 1:n.zvar])
    tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars,
                             model = c("correlation"), minsplit = 720)
    
    # Extract result
    res <- unlist(tree)
    if(is.null(res$node.kids.id)){
      split <- "no"
    } else {
      split <- "yes"
    }
    
    # Write Result to Datafile
    CSimB_s25$n.zvar[id] <- n.zvar
    CSimB_s25$n.nds[id] <- n.nds
    CSimB_s25$n.obs[id] <- n.obs
    CSimB_s25$delta.all[id] <- delta.all
    CSimB_s25$delta.cor[id] <- delta.cor
    CSimB_s25$mean.dif[id] <- mean.dif
    CSimB_s25$mean.change[id] <- mean.change
    CSimB_s25$split[id] <- split
    print(id)
  }
}


write.csv(CSimB_s25, "CSimB_s25.csv")

CSimB_s <- rbind(CSimB_s5, CSimB_s15, CSimB_s25)
write.csv(CSimB_s, "CSimB_s.csv")
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------


#### SIMULATION C: False Positive: delta rho = 0 

# N.Nodes: 5
obs <- c(100, 500, 1500)
min <- c(35, 220, 720)
id <- 0
sim <- 50
CSimC_5 <- data.frame(id = numeric(sim*length(obs)))


for(o in 1:length(obs)){
  for (it in 1:sim){
    # Keep Track of the ID
    id <- id + 1
    
    # Determine Simulation variables
    n.zvar <- 1
    n.nds <- 5
    n.obs <- obs[o]
    delta.all <- FALSE
    delta.cor <- 0
    mean.dif <- FALSE
    mean.change <- .6
    
    #Simualte Data
    data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                delta.cor = delta.cor, delta.all = delta.all, 
                                mean.dif = mean.dif, mean.change = mean.change)
    
    # Run Ising MOB split function
    nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
    splitvars <- as.data.frame(data[, 1:n.zvar])
    tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, 
                             model = c("correlation"), minsplit = min[o])
    
    # Extract result
    res <- unlist(tree)
    if(is.null(res$node.kids.id)){
      split <- "no"
    } else {
      split <- "yes"
    }
    
    # Write Result to Datafile
    CSimC_5$n.zvar[id] <- n.zvar
    CSimC_5$n.nds[id] <- n.nds
    CSimC_5$n.obs[id] <- n.obs
    CSimC_5$delta.all[id] <- delta.all
    CSimC_5$delta.cor[id] <- delta.cor
    CSimC_5$mean.dif[id] <- mean.dif
    CSimC_5$mean.change[id] <- mean.change
    CSimC_5$split[id] <- split
    print(id)
  }
}

write.csv(CSimC_5, "CSimC_5.csv")

# --------------------------------------------------------------------------------

# Simulation C: 
# N.Nodes: 15
obs <- c(500, 1500)
min <- c(220, 720)
id <- 0
sim <- 50
CSimC_15 <- data.frame(id = numeric(sim*length(obs)))


for(o in 1:length(obs)){
  for (it in 1:sim){
    # Keep Track of the ID
    id <- id + 1
    
    # Determine Simulation variables
    n.zvar <- 1
    n.nds <- 15
    n.obs <- obs[o]
    delta.all <- FALSE
    delta.cor <- 0
    mean.dif <- FALSE
    mean.change <- .6
    
    #Simualte Data
    data <- IsingDataSimulation(n.nds = n.nds, n.obs = n.obs, n.zvar = n.zvar, 
                                delta.cor = delta.cor, delta.all = delta.all, 
                                mean.dif = mean.dif, mean.change = mean.change)
    
    # Run Ising MOB split function
    nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
    splitvars <- as.data.frame(data[, 1:n.zvar])
    tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars,
                             model = c("correlation"), minsplit = min[o])
    
    # Extract result
    res <- unlist(tree)
    if(is.null(res$node.kids.id)){
      split <- "no"
    } else {
      split <- "yes"
    }
    
    # Write Result to Datafile
    CSimC_15$n.zvar[id] <- n.zvar
    CSimC_15$n.nds[id] <- n.nds
    CSimC_15$n.obs[id] <- n.obs
    CSimC_15$delta.all[id] <- delta.all
    CSimC_15$delta.cor[id] <- delta.cor
    CSimC_15$mean.dif[id] <- mean.dif
    CSimC_15$mean.change[id] <- mean.change
    CSimC_15$split[id] <- split
    print(id)
  }
}


write.csv(CSimC_15, "CSimC_15.csv")

# --------------------------------------------------------------------------------

# Simulation C:
# N.Nodes: 25
id <- 0
sim <- 50
CSimC_25 <- data.frame(id = numeric(sim))


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
                              mean.dif = mean.dif, mean.change = mean.change)
  
  # Run Ising MOB split function
  nodevars <- as.data.frame(data[, (n.zvar + 1):(n.zvar + n.nds)])
  splitvars <- as.data.frame(data[, 1:n.zvar])
  tree <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, 
                           model = c("correlation"), minsplit = 720)
  
  # Extract result
  res <- unlist(tree)
  if(is.null(res$node.kids.id)){
    split <- "no"
  } else {
    split <- "yes"
  }
  
  # Write Result to Datafile
  CSimC_25$n.zvar[id] <- n.zvar
  CSimC_25$n.nds[id] <- n.nds
  CSimC_25$n.obs[id] <- n.obs
  CSimC_25$delta.all[id] <- delta.all
  CSimC_25$delta.cor[id] <- delta.cor
  CSimC_25$mean.dif[id] <- mean.dif
  CSimC_25$mean.change[id] <- mean.change
  CSimC_25$split[id] <- split
  print(id)
}


write.csv(CSimC_25,"CSimC_25.csv")

CSimC <- rbind(CSimC_5, CSimC_15, CSimC_25)

write.csv(CSimC,"CSimC.csv")

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------



