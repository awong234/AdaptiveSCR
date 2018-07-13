library(doParallel)
library(jagsUI)

# Extract aug files

files = grep(".Rdata", x = dir(), value = TRUE)

for(i in 1:length(files)){
  load(files[i])
}

load("data_setup.Rdata")

# GET SETUP DATA --------------------------------------------------------------------------------------------------------

# Extraction function

extract = function(what){invisible(Map(f = function(x,y){assign(x = x, value = y, envir = .GlobalEnv)}, x = names(what), y = what))}

# Get augmented datasets out 

extract(subset.aug[c("grid.id.subset", "grid.st.subset", "zst.subset.rand", "Slist.aug", "Y.subset.aug", "Robs")])

# Get setup variables out

extract(data.setup[c("xlim", "ylim", "traps", "A", "R", "K")])

# Need cov.scr analog

extract(data.subset[c("covscr", "covindex", "Rindex")])

# What's M? 

M = length(grid.id.subset[[1]])

rm(data.complete, data.m5, data.setup, data.subset, subset.aug, complete.aug)


# Test model 1 run - most complex -----------------------------------------------------------------------------------------------

# Model succeeds in converging

parameters <- c("psi","p0","sigma","beta0","Ntotal","Ntotal.scr","Ntotal.index","beta.cov")

inits <- function(){list (z=zst.subset.rand[[a]],sigma=runif(1,0.2, 1) ,S=Slist.aug[[a]], p0=0.1 ,grid.id=grid.st.subset[[a]], beta.cov=1 ) }

nburn = 500
niter = 10000

# registerDoParallel(cores = detectCores()-1)
registerDoParallel(cores = 2)
getDoParWorkers()

jagsrun = function(a){
  data <- list(Robs = Robs[[a]], grid.id = grid.id.subset[[a]], X = as.matrix(traps), K = K, Y = Y.subset.aug[[a]],
               bigM = M, ntraps = nrow(traps), ylim = ylim, xlim = xlim, cov.scr = covscr[[a]],
               cov.index = covindex[[a]], Rindex = Rindex[[a]])


  # m1.out = jags(data, inits, parameters, "model1.txt", n.thin=1, n.chains=4,
  #               n.burnin = nburn, n.iter = niter, parallel = TRUE)

  m1.out = jags(data, inits, parameters, "model1.txt", n.thin=1, n.chains=4, n.burnin = nburn,
                    n.iter = niter, parallel=TRUE)

  save(m1.out, file = paste("MODEL1_OUT\\MODEL1_",a, ".Rdata",sep = ""))
  
  write(a, file = "MODEL1_OUT\\Sims_completed.txt", ncolumns = 1, append = TRUE)
}

# Check if output folder exists -- if not, create it
if (!dir.exists("MODEL1_OUT")) {
  dir.create("MODEL1_OUT")
  file.create("MODEL1_OUT\\Sims_completed.txt")
}

# Completed simulation #'s; which simulations already completed?

A.del = as.numeric(readLines("MODEL1_OUT\\Sims_completed.txt"))

# Remaining simulation #'s. This gets rid of any analyses that have already been completed (i.e. written to Sims_completed.txt) so that the analyses can be started and stopped at will. 

A.seq = (1:A)
A.seq = A.seq[!A.seq %in% A.del]

foreach(a = A.seq, .packages = "jagsUI") %dopar% jagsrun(a)
