# Augments data output from AS_dataset_sim

# New file due to corruption in v3.

# Additionally; changing grid.st.complete to be a vector from 1:R instead of sampling R. 

# # Testdata
# source("AS_dataset_sim_v4.3.R")
# AS.simulator()

# Goal is to 
# 1. Load the requisite files from AS.simulator()
# 2. Grab the APPROPRIATE grid items, augment them
# 3. Create zst's
# 4. Create grid.st
# 5. Create Sst. 
# 6. Pad appropriate Y with M-nrow(Y) 0's 

# Test M
# M = 40; m5.unob = 0

AS.augmentor = function(M = NULL){
  
  require(foreach)
  
  # Get files from wd -------------------------------------------------------------------------
  filetemp = "data_"
  files = grep(filetemp,dir(), perl = T, value = TRUE)
  f = length(files)
  
  for(i in 1:f) load(files[i])
  
  Nmat = data.complete$Nmat
  
  max.N = max(colSums(Nmat))
  
  # Extract old setup params ------------------------------------------------------------------
  
  xlim = data.setup$xlim
  ylim = data.setup$ylim
  A = data.setup$A
  R = data.setup$R
  
  # Set M if unset ----------------------------------------------------------------------------
  
  if(is.null(M)){
    message(paste("NOTE: Max N is", max.N))
    M = as.integer(readline(prompt = "M NOT SET, ENTER VALUE: "))
  }
  
  if(M > 3*max.N) message("Warning: M is greater than 3*max.N, may result in inefficiencies")
  
  # Pull out relevant data ------------------------------------------------
  
  # GRIDS
  grid.complete = data.complete$gridlist
  grid.subset = data.subset$newgrid.scr
  grid.m5 = data.m5$newgrid.m5
  
  grid.all = list("grid.complete" = grid.complete, "grid.subset" = grid.subset, "grid.m5" = grid.m5)
  
  # Y's
  Y.complete = data.complete$Y
  Y.subset = data.subset$Yscr
  Y.m5 = data.m5$Y.m5
  
  Y.all = list("Y.complete" = Y.complete, "Y.subset" = Y.subset, "Y.m5" = Y.m5)
  
  # S's
  Slist = data.complete$Slist
  
  # covmat
  covmat = data.complete$covmat
  
  # Robs - number sites sampled by SCR
  Robs = data.subset$Robs
  
  # where covmat > 0
  
  browser()
  
  cov.no.0 = foreach(col = 1:ncol(covmat)) %do%  {covmat[,col] > 0}
  
  # Augment grid items ----------------------------------------------------
  aug.fn = function(grid.obj, y.obj){
    mapply(FUN = function(x,y){c(x, rep(NA, M-nrow(y)))}, x = grid.obj, y = y.obj, SIMPLIFY = F)
  }
  
  # GRID.ID analog --------------------------------------------------------
  aug.grid.id = mapply(FUN = aug.fn, grid.obj = grid.all, y.obj = Y.all, SIMPLIFY = F)
  
  # Functions for GRID.ST --------------------------------------------------------
  
  # want to replace the values of GRID.ID with NA's, and then an assortment of
  # grid #'s where originally were NA's. 
  
  # So if gridid = [1,2,3,3,NA,NA,NA], change to [NA,NA,NA,1,2,3]
  
  # Starting grid for subset data - note, operates as Andy's code
  
  grid.st.sbfn = function(grid.obj){
    grid.st = grid.obj
    
    for(a in 1:A){
      grid.st[[a]][!is.na(grid.obj[[a]])] = NA
      grid.st[[a]][is.na(grid.obj[[a]])] = 1:Robs[[a]]
    }
    return(grid.st)
  }
  
  # Starting grid for complete data
  
  # Original note: "below I don't include plots as starting values if the
  # covariate = 0, causes bad starting value problem"
  
  # This is the utility of cov.no.0
  
  grid.st.compfun = function(grid.obj){ #Accepts a slice s of aug.grid.id[[s]]
    
    grid.st = grid.obj
    
    for(a in 1:A){
      grid.st[[a]][!is.na(grid.obj[[a]])] = NA
      grid.st[[a]][is.na(grid.obj[[a]])] = sample((1:R)[cov.no.0[[a]]], sum(is.na(grid.obj[[a]])), replace = TRUE)
    }
    return(grid.st)
  }
  
  # Starting grid for M5 data - mimics subset data in that it assigns sequence
  # 1:(#sites sampled) instead of sampling random grids. 
  
  # UNNECESSARY, A MISTAKE. THERE ARE ROBS SITES IN M5, THIS LEAVES SITES EMPTY TO START.
  
  # # Analog for Robs for model 5
  # m5max = lapply(grid.m5, max, na.rm = TRUE)
  
  grid.st.m5fn = function(grid.obj){
    grid.st = grid.obj
    
    for(a in 1:A){
      grid.st[[a]][!is.na(grid.obj[[a]])] = NA
      grid.st[[a]][is.na(grid.obj[[a]])] = 1:Robs[[a]]
    }
    return(grid.st)
  }
  
  # GRID.ST ANALOG FINAL ---------------------------------------------------------------
  
  # OUTPUT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  grid.st.complete = suppressWarnings(grid.st.sbfn(aug.grid.id$grid.complete)) # Changed from grid.id.compfn in v.3.1.
  grid.st.subset = suppressWarnings(grid.st.sbfn(aug.grid.id$grid.subset))
  grid.st.m5 = suppressWarnings(grid.st.m5fn(aug.grid.id$grid.m5))
  # OUTPUT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  # Zst analog --------------------------------------------------------------------------------
  
  
  # Zst subset
  zst.method.rand = function(grid.obj, Y.obj, add1){
    
    # tests to see if we're using subset data.
    # if so, use Robs to determine how many random {0,1} to generate
    # else we're in m5 land; find out how many sites have individuals using the lapply function
    
    Map(f = function(x,y,z){c(rep(1,length(x)), # ones for observed individuals
                              rep(1,y), # additional ones for Robs individuals
                              rbinom(M-nrow(z)-y,1,0.5))}, # random starting {0,1} for unobserved
        x = grid.obj, y = add1, z = Y.obj)
  }
  
  zst.method.alt = function(grid.obj, Y.obj, unob.status = 0){
    
    Map(f = function(x,y){c(rep(1, length(x)), # ones for observed individuals
                            rep(unob.status, M-nrow(y)))}, # either 1's or 0's for unobserved individuals (everyone else)
        x = grid.obj, y = Y.obj)
  }
  
  
  # OUTPUT :::::::::::::::::::::::::::
  zst.subset.rand = zst.method.rand(grid.obj = grid.subset, Y.obj = Y.subset, add1 = Robs)
  zst.subset.1 = zst.method.alt(grid.obj = grid.subset, Y.obj = Y.subset, unob.status = 1)
  zst.subset.0 = zst.method.alt(grid.obj = grid.subset, Y.obj = Y.subset, unob.status = 0)
  # OUTPUT :::::::::::::::::::::::::::
  
  # For zst for complete
  
  
  # OUTPUT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  zst.complete.rand = zst.method.rand(grid.obj = grid.complete, Y.obj = Y.complete, add1 = 0)
  zst.complete.1 = zst.method.alt(grid.obj = grid.complete, Y.obj = Y.complete, unob.status = 1)
  zst.complete.0 = zst.method.alt(grid.obj = grid.complete, Y.obj = Y.complete, unob.status = 0)
  # OUTPUT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  # Zst for m5
  
  # OUTPUT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  zst.m5.rand = zst.method.rand(grid.obj = grid.m5, Y.obj = Y.m5, add1 = Robs)
  zst.m5.1 = zst.method.alt(grid.obj = grid.m5, Y.obj = Y.m5, unob.status = 1)
  zst.m5.0 = zst.method.alt(grid.obj = grid.m5, Y.obj = Y.m5, unob.status = 0)
  # OUTPUT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  
  # Sst analog -----------------------------------------------------------------------------
  
  # Taking the Slist and rbinding M - nrow(S) random numbers to it
  
  Slist.aug.fn = function(S){
    
    lapply(X = S, FUN = function(x){rbind(x, cbind(
      runif(M-nrow(x), xlim[1], xlim[2]), 
      runif(M-nrow(x), ylim[1], ylim[2])))})
    
  }
  
  # OUTPUT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  Slist.aug = Slist.aug.fn(Slist)
  # OUTPUT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  # Y augmentation analog ------------------------------------------------------------------
  
  # Taking each Y and rbinding M-nrow(Y[[a]]) all 0 encounter histories
  
  # Apply to each Y
  
  Y.aug.fn = function(Y){
    
    lapply(Y, function(x){rbind(x, matrix(0, nrow = M-nrow(x), ncol = ncol(x)))})
    
  }
  
  # OUTPUT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  Y.complete.aug = Y.aug.fn(Y.complete)
  Y.subset.aug = Y.aug.fn(Y.subset)
  Y.m5.aug = Y.aug.fn(Y.m5)
  
  
  # OUTPUT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  
  
  
  
  # FINAL DATA EXPORT ----------------------------------------------------------------------
  export = list("grid.id.complete" = aug.grid.id$grid.complete, "grid.id.subset" = aug.grid.id$grid.subset,
                "grid.id.m5" = aug.grid.id$grid.m5, 
                "grid.st.complete" = grid.st.complete, "grid.st.subset" = grid.st.subset, "grid.st.m5" = grid.st.m5, 
                "zst.complete.rand" = zst.complete.rand, "zst.complete.1" = zst.complete.1, "zst.complete.0" = zst.complete.0, 
                "zst.subset.rand" = zst.subset.rand, "zst.subset.1" = zst.subset.1, "zst.subset.0" = zst.subset.0,
                "zst.m5.rand" = zst.m5.rand, "zst.m5.1" = zst.m5.1, "zst.m5.0" = zst.m5.0, 
                "Slist.aug" = Slist.aug, 
                "Y.complete.aug" = Y.complete.aug, "Y.subset.aug" = Y.subset.aug, "Y.m5.aug" = Y.m5.aug,
                "Robs" = Robs,
                "Robs.m5" = Robs,
                "M" = M)
  
  complete.aug = export[grep(pattern = "complete|Slist", x = names(export), value = TRUE)]
  subset.aug = export[grep(pattern = ".subset|Slist|Robs$", x = names(export), value = TRUE)]
  m5.aug = export[grep(pattern = ".m5|Slist|Robs.m5", x = names(export), value = TRUE)]
  
  # WHY CAN'T I GET ROBS INTO SUBSET.AUG?!
  
  # save into separate augmentation packages.
  
  save(complete.aug, file = "aug_complete.Rdata")
  save(subset.aug, file = "aug_subset.Rdata")
  save(m5.aug, file = "aug_m5.Rdata")
  
}


