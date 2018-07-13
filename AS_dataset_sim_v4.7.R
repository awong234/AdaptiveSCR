# Version 4.7

# What's new in 4.7:

# Added retention of all positions, not just those encountered.

#### Adaptive Sampling Dataset Simulation ####

# # Test values
# # Small:
# A = 3;R = 30; K = 3; sigma = 1; lam.N = 2; cov.mult = 3; cutoff = 4
# # Large
# A = 1000;R = 100; K = 3; sigma = 1; lam.N = 2; cov.mult = 3; cutoff = 4


#### NOTICE ####

# A is an OUTER simulation replicate for A runs of each model
# R is an INNER replicate of R state spaces ("plots") within 1 replicate of 1 model


# Simulation function ============================================================================================

##### BEGIN COMPLETE DATA SIM ################################################################################
AS.simulator = function(A = 3,
                        R = 10,
                        K = 3,
                        cutoff = 4,
                        sigma = 1, 
                        lam.N = 2, 
                        cov.mult = 3,
                        traps = NULL){ #Defaults to test. Can change later.
  
  
  # Optional
  library(ggplot2)
  library(doParallel)
  # library(microbenchmark)
  # library(scrbook)
  # library(beepr)
  
  # Required
  require(fields)
  
  if(is.null(traps)){ #If alternative trap setups aren't specified, do the original scheme.
    
    # State space
    xlim<- c(-1, 4)
    ylim<- c(-1, 5)
    
    # Sample units
    # line<-rbind( c(1,1), c(1,3), c(2,3), c(2, 1) )
    # plot(line, xlim=xlim,ylim=ylim)
    # lines(line)
    
    # Chop the lines up into 0.1 length units for SCR traps - doesn't change between sims
    delta<- 0.10
    trap1<- cbind(1,seq(1,3,delta))
    trap2<- cbind(seq(1 + delta,2,delta), 3)
    trap3<- cbind(2, seq(1,3-delta,delta))
    traps<- rbind(trap1,trap2,trap3)
    
    # points(traps<- rbind(trap1,trap2,trap3),pch=20)
  }
  
  data.setup = list("xlim" = xlim, "ylim" = ylim, "traps" = traps, "A" = A, "R" = R,
                    "K" = K, "cutoff" = cutoff, "sigma" = sigma, "lam.N" = lam.N, 
                    "cov.mult" = cov.mult)
  
  save(data.setup, file = "data_setup.Rdata")
  
  # Simulate population size of each state-space -------------------------------------------------
  
  
  #**************************************************************
  Nmat = matrix(data = rpois(R*A, lambda = lam.N), nrow = R, ncol = A)
  
  rownames(Nmat) = paste("grid", seq(1:nrow(Nmat))) # rownames, for information
  #**************************************************************
  
  
  # Simulate covariates. ----------------------------------------------------------------------
  
  # There are as many covariates as there are entries in Nmat - A*R. 
  
  #**************************************************************
  covmat = matrix(data = rpois(A*R, lambda = Nmat*cov.mult), nrow = R, ncol = A)
  
  rownames(covmat) = rownames(Nmat)
  #**************************************************************
  
  # ***LIST MAKER*** ====================================================================================
  
  # Use for all lists to be generated. FUN is to be the function generating the
  # list. FUN = expression(`expression`)
  
  #**************************************************************#**************************************************************#
  #**************************************************************#**************************************************************#
  #**************************************************************#**************************************************************#
  
  # VERY IMPORTANT - implemented a lot!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  list.maker = function(FUN){
    rlist = list1 = list()
    for(a in 1:A){
      for(r in 1:R){
        rlist[[r]] = eval(FUN)
      }
      list1[[a]] = rlist
      rlist = list()
    }
    return(list1)
  }
  
  # RETURNS LISTS OF STRUCTURE LIST[[A]][[R]]. CAN RETRIEVE DATA FOR ANY
  # SIMULATION a, AND ANY REPLICATE r. ACCOMMODATES ANY DATA TYPE (animal positions, dmats, etc.)
  
  # Accepts any function & variables wrapped in "expression()".
  
  #**************************************************************#**************************************************************#
  #**************************************************************#**************************************************************#
  #**************************************************************#**************************************************************#
  
  # Simulate animal positions, for each R ------------------------------------------------------
  
  # But don't do it where N == 0. 
  
  # loops not all that bad with lists Slist will have dim Slist[[A]][[R]]. This 
  # is in contrast to the arrays; R moves fastest, and so it should be on the 
  # end.
  
  s.sim.list = function(n){ #accepts a single Nmat[r,a] value
    if(n == 0) return(NA)
    else{
      cbind(runif(n,xlim[1],xlim[2]),
            runif(n,ylim[1],ylim[2]))
    }
  }
  
  # s.sim.list is much closer to the original intended - mostly no pesky NA's
  
  # However, if N == 0, then positions are cbind(NA,NA); else simulate positions
  # as originally simulated. 
  
  # NULLS are DANGEROUS WITHIN FOR LOOPS OVER LISTS!!!!!!
  
  #**************************************************************
  # Analog to Sr
  Slist.byR = list.maker(FUN = expression(s.sim.list(Nmat[r,a])))
  
  # rowID - this function gives an individual ID for each individual in x, 
  # whether that is the dmat or Y or something else. Used primarily to verify 
  # that the functions eliminating unobserved individuals (using ycap) is doing
  # its job correctly
  
  rowID = function(x){
    for(a in 1:A){
      for(r in 1:R){
        if(NA %in% x[[a]][[r]]){next}
        else rownames(x[[a]][[r]]) = paste("ID", a,"-",r,"-",1:nrow(x[[a]][[r]]), sep = "")
      }
    }
    return(x)  # returns ID{a}-{r}{ind#}
  }
  
  Slist.byR = rowID(Slist.byR)
  
  # Analog to S - later we'll use this along with ycap to trim down. 
  Slist = list()
  for(a in 1:A){
    Slist[[a]] = do.call(what = rbind, args = Slist.byR[[a]])
  }
  
  # Remove NA's
  Slist = lapply(X = Slist, FUN = function(x){x[complete.cases(x),]})
  
  
  # Specify DMATS  ---------------------------------------------------------------------
  
  dmat = list.maker(FUN = expression(if(NA %in% Slist.byR[[a]][[r]]){return(NA)}
                                     else(rdist(Slist.byR[[a]][[r]], traps))))
  
  # If there's an NA, do nothing. Else, compute distance matrix between Sr and traps. 
  
  dmat = rowID(dmat)
  
  # For ID's of individuals. I have verified that a particular individual's
  # distance matrices are correctly calculated from their positions in Slist.byR
  
  
  
  # Simulate Y ----------------------------------------------------------------------------
  
  sim.y = function(dmat){ #accepts a single slice of Dmat - that is Dmat[[a]][[r]]
    N = nrow(dmat)
    if(NA %in% dmat){return(NA)} # for NA distance matrices (because can't compute with NA positions), return NA
    else{ #else simulate captures as originally specified.
      col = ncol(dmat)
      pmat = 0.1*exp(-dmat*dmat/(2*sigma*sigma))
      
      # dim(pmat) == dim(dmat)
      
      Y = NULL
      for(j in 1:N){ # note, for loop faster with small N - most of the time, on average faster. 
        y = rbinom(nrow(traps), K,  pmat[j,])
        Y = rbind(Y, y)
      }
      return(Y)}
  }
  
  # browser()
  

  # Each R is separated. Nested list Y[[a]][[r]]
  Y.byR = list.maker(FUN = expression(sim.y(dmat[[a]][[r]])))
  
  Y.byR = rowID(Y.byR)
  
  # Direct analog to Y in first sim. Y.byR[[a]] - R's are condensed.
  Y = list()
  for(a in 1:A){
    Y[[a]] = do.call(what = rbind, args = Y.byR[[a]])
  }
  
  Y = lapply(X = Y, FUN = function(x){x[complete.cases(x),]}) # remove NA's
  

  # Analog to ycap ------------------------------------------------------------------------------
  # recall, TRUE if captured at least once
  ycap = lapply(X = Y, function(x){apply(x,1,sum)>0})
  
  # Analog to nind
  # nind = sum(ycap)
  
  nind = lapply(ycap, sum)
  
  # keep those names for which ycap == T; for use when trimming Y.final below
  keep.names = lapply(ycap, function(x){names(which(x == T))})
  
  
  # Y.final - Apply ycap ----------------------------------------------------------------------------------------------------------------
  
  #**************************************************************#**************************************************************
  Y.final = mapply(FUN = function(x,y){x[y,]}, x = Y, y = keep.names)
  #**************************************************************#**************************************************************
  
 
  # The names of those individuals for whom ycap returns TRUE are used to slice 
  # from Y. Those individuals for whom ycap returns FALSE are not included.
  
  
  # THE GRID OBJECT ------------------------------------------------------------------------------
  
  # 1 gridlist per A. list of A grid objects, then. for each N[,a], do 
  # rep(1:R,N[,a])
  
  gridlist = list()
  for(a in 1:A){
    gridlist[[a]] = rep(1:R,Nmat[,a])
  }
  
  # now subset grid by ycap as done in original code
  #**************************************************************#************************************************
  gridlist.final = mapply(FUN = function(x,y){x[y]}, x = gridlist, y = ycap, SIMPLIFY = FALSE)
  #**************************************************************#************************************************
  
  # Slist.final  ---------------------------------------------------------------------------------------------
  
  Slist.final = mapply(FUN = function(x,y){x[y,]}, x = Slist, y = ycap)
  
  
  # EXPORT DATA ===================================================================================================  
  data.complete = list("Nmat" = Nmat, "covmat" = covmat, "SlistFull" = Slist, "Slist" = Slist.final,
                       "Y" = Y.final, "nind" = nind, "gridlist" = gridlist.final, 
                       "traps" = traps, "K" = K)
  
  save(data.complete, file = "data_complete.Rdata")
  
  # save(Nmat, covmat, Slist.final, Y.final, nind, gridlist.final, traps, K, 
  #      file = "data_complete.Rdata")
  
  # END COMPLETE DATA SECTION ###########################################################################
  
  
##### BEGIN SUBSET DATA SECTION########################################################################
  
  
  # From andy's code:
  
  # Now I subset the data into two parts:
  # The SCR data only for plots where "covariate > cutoff"
  # The index data for the non-SCR plots
  
  # `ss` analog ---------------------------------------------------------------
  
  # Recall, ss is the vector of gridID's for which covariate > cutoff.
  
  # for each a in A, we eventually need a subset (1:R)[covmat[,a] > cutoff] However, 
  # while covmat[,a] > cutoff is a vector of R elements, the resulting vectors 
  # for ss won't be the same across all a's. See why:
  
  # colSums(covmat > cutoff) # The number of TRUE's are not consistent across a's. 
  
  # Therefore we need a list for the ss's.
  
  sslist = foreach(a = 1:A) %do% (1:R)[covmat[,a] > cutoff]
  
  # Subset Nmat ----------------------------------------------------------------
  
  # Plot #'s where covariate > cutoff. Going to have to be a list of A elements 
  # each containing a vector of N's for which the covmat > cutoff.
  
  # originally it was N[ss]. So, analog is Nmat[,a][ss]. The following specifies
  # that for each simulation, draw the appropriate column from Nmat and subset
  # by the appropriate list element in ss.list. 
  
  Nscr.ss = foreach(a = 1:A) %do% Nmat[,a][sslist[[a]]]
  
  # verification - notice that the grid ID's (rownames) accompanying Nmat match ss.list
  
  # Keep object --------------------------------------------------------------
  
  # Recall, in the keep object we're keeping the gridID's for each individual that
  # are in plots where the plot is part of ss. 
  
  # for all the individuals caught > 1, we're checking that gridID is in ss. Returns a
  # vector of length `nind` for each a in A. 
  
  keeplist = mapply(FUN = function(x,y){x %in% y}, x = gridlist.final, y = sslist)
  
  # Since gridlist[[a]] is just 1:R repeated Nmat[,a] times, we're just keeping the
  # indices matching ss.list[[a]]
  
  # NOTE: I changed the use of gridlist to gridlist.final. Why? Because we only 
  # want those individuals that were caught at least once. Check the original
  # code: grid = grid[ycap], AND THEN the subset grid %in% ss is applied. 
  
  # Subset Y.final ------------------------------------------------------------------
  
  # We use Y.final, and keep.list. `keep.list` is a vector of T/F of size `nind`. The 
  # original code was Y = Y[keep,], stating that for all columns(traps), we are 
  # taking the individuals' encounter histories that fall into plots where cov >
  # cutoff.
  
  Y.ss = mapply(FUN = function(x,y){x[y,]}, x = Y.final, y = keeplist)
  
  # Note, changed Y to Y.final. We don't want those individuals who weren't captured.
  
  # Verify - ss.list[[a]] is the grids we're keeping due to the condition. Since
  # the rownames for Y are ID{A}-{R}-{#}, we can see that only those R #'s that
  # were in ss.list[[a]] exist in Y.ss[[a]].
  
  nscrlist = lapply(X = Y.ss, FUN = nrow)
  
  # subset grid.list -------------------------------------------------------
  
  # Same as subsetting Y - the subset operates over `nind` items. Y.final[[a]] has
  # `nind` rows, and grid.list[[a]] has length `nind`. 
  
  gridlist.ss = mapply(FUN = function(x,y){x[y]}, x = gridlist.final, y = keeplist)
  
  # Changed gridlist to gridlist.final
  
  # Subset covariate -------------------------------------------------------
  
  # do just like Nmat
  
  covlist.ss = foreach(a = 1:A) %do% covmat[,a][sslist[[a]]]
  
  # Newgrid ---------------------------------------------------------------------
  
  # step 1: newg = 1:length(ss)
  # step 2: names(newg) = ss
  # step 3: newgrid = newg[as.character(grid)]
 
  newgridfn = function(grid.ls, ss){
    
    newg = foreach(a = 1:A) %do% {
      as.integer(as.factor(ss[[a]]))
    }
    
    for(a in 1:A){
      names(newg[[a]]) = ss[[a]]
    }
    
    # ***************************************************************************
    new.gridlist = foreach(a = 1:A) %do% newg[[a]][as.character(grid.ls[[a]])]
    # ***************************************************************************
    
  }
  
  new.gridlist = newgridfn(gridlist.ss, sslist)
  
  # Robs ----------------------------------------------------------------------------
  
  Robs = lapply(covlist.ss, length)
  
  # covindex & Rindex ------------------------------------------------------------------------
  
  # covmat[covmat <= cutoff]
  
  covindex = foreach(a = 1:A) %do% covmat[,a][covmat[,a] <= cutoff]
  
  Rindex = lapply(covindex, length)
  
  # EXPORT DATA -----------------------------------------------------------------------------
  data.subset = list("Nscr" = Nscr.ss, "Yscr" = Y.ss, "newgrid.scr" = new.gridlist, 
                     "Robs" = Robs, "covscr" = covlist.ss, "covindex" = covindex, "Rindex" = Rindex)
  
  save(data.subset, file = "data_subset.Rdata")
  
  # END SUBSET DATA SECTION -------------------------------------------------
  
#### BEGIN MODEL 5 DATA SECTION ###########################################################################
  
  # Want a subset of data.complete where for each a in A simulations the sample 
  # size of sites selected is equal to the number of sites selected in
  # data.subset; this is Robs. 
  
  # QUESTION! Should we be using gridlist.final, or gridlist to develop the
  # keep.m5 object? 
  
  # Gridlist is the full set of gridID's for all individuals. Gridlist.final is
  # the set of gridID's for those individuals captured at least once. 
  
  # When we subset, we are picking random plots to survey - yes, we want this. 
  # When we subset N, we are grabbing N at those plots - yes. When we subset 
  # gridlist, we are potentially grabbing individuals who we never encountered -
  # we don't want this. So using gridlist.final is correct. 
  
  # When we subset Y, we take those in grids that we are sampling (ss.m5), but 
  # of course, we don't want those who we didn't capture at least once, so using
  # Y.final is appropriate. 
  
  # Subset m5 --------------------------------------------------------------------------------------
  
  ss.m5 = lapply(X = Robs, FUN = function(x){sort(sample(x = c(1:R), size = x))})
  
  # we're using each Robs[[a]] to set the sample size for the sample function. 
  # This accomplishes the goal of sampling a random sample EQUAL to those 
  # selected by adaptive sampling. So, if `n` sites were sampled during adaptive
  # sampling in simulation `a`, `n` RANDOM sites will ALSO be sampled in simulation `a`
  # for model 5. 
  
  # Subset Nmat ------------------------------------------------------------------------------------
  
  N.m5 = foreach(a = 1:A) %do% Nmat[,a][ss.m5[[a]]]
  
  # Subset grid (keep) -----------------------------------------------------------------------------
  
  keep.m5 = mapply(FUN = function(x,y){x %in% y}, x = gridlist.final, y = ss.m5)
  
  # recall: 
  # gridlist.final = mapply(FUN = function(x,y){x[y]}, x = gridlist, y = ycap, SIMPLIFY = FALSE)
  
  # Of those individuals caught, we are assessing whether their grid is within
  # the subset we're sampling. If not, we ignore their capture histories NEXT.
  
  # Subset Y ---------------------------------------------------------------------------------------
  
  Y.m5 = mapply(FUN = function(x,y){x[y,]}, x = Y.final, y = keep.m5)
  
  # Here, among those caught (Y.final), we're only interested in the capture
  # histories of those within the subset of sites we're sampling.
  
  # n.m5 - number of individuals sampled ------------------------------------------------------------
  
  n.m5 = lapply(Y.m5, nrow)
  
  
  # Subset gridlist.final ---------------------------------------------------------------------------
  
  gridlist.m5 = mapply(FUN = function(x,y){x[y]}, x = gridlist.final, y = keep.m5)
  
  # Subset covariate (just in case) -----------------------------------------------------------------
  
  # I don't know whether model 5 will be using covariate data or not. I think
  # not, but I'll include it just in case it's needed. 
  
  # Once again, we're just grabbing covariates at those sites we visited. 
  
  cov.m5 = foreach(a = 1:A) %do% covmat[,a][ss.m5[[a]]]
  
  unvisited = foreach(a = 1:A) %do% (1:R)[!(1:R) %in% ss.m5[[a]]]
  unvisited.cov.m5 = foreach(a = 1:A) %do% covmat[,a][unvisited[[a]]]
  
  # Newgrid -------------------------------------------------------------
  
  # newg = foreach(a = 1:A) %do% {
  #   as.integer(as.factor(sslist[[a]]))
  # }
  # 
  # for(a in 1:A){
  #   names(newg[[a]]) = sslist[[a]]
  # }
  # 
  new.gridlist.m5 = newgridfn(gridlist.m5, ss.m5)
  
  # Robs analog ---------------------------------------------------------
  
  # NOTE: There is no need for an additional object here. Robs == # sites surveyed in m5. 
  
  # ALSO NOTE: lapply(gridlist.m5, max) is NOT going to be # sites surveyed, 
  # because gridlist doesn't include any individuals not surveyed. Use cov.m5
  # names to get sites surveyed whether or not individuals exist.
  
  # Index stuff ---------------------------------------------------------
  
  # Do we need this? We aren't measuring anything at sites we don't visit. Won't include. 
  
  # EXPORT DATA ---------------------------------------------------------
  data.m5 = list("n.m5" = n.m5, "Y.m5" = Y.m5, "newgrid.m5" = new.gridlist.m5, "cov.m5" = cov.m5, "unv.cov.m5" = unvisited.cov.m5)
  
  save(data.m5, file = "data_m5.Rdata")
  
## Write truth values #########################################################  
  
  # Including: 
  
  # Sim # 
  # # ind captured > 1 
  # # ind captured > 1 in SCR subset of A.S.
  # # ind captured > 1 in SRS subset (M5)
  # Total # ind
  # Total # ind in SCR subsets
  # Total # ind in M5 subsets
  
  headers = c("a", "Robs", "nind", "nscrlist", "n.m5", "sumNmat", "sumNscr", "sumNm5")
  
  tv = matrix(NA, nrow = A, ncol = length(headers))
  
  for(a in 1:A){
    tv[a,] = c(a, Robs[[a]], nind[[a]], nscrlist[[a]], n.m5[[a]], sum(Nmat[,a]), sum(Nscr.ss[[a]]), sum(N.m5[[a]]))
  }
  
  colnames(tv) = headers
  
  write.csv(x = tv, file = "truevals.csv", row.names = F)
                 
} # END AS.SIMULATOR
