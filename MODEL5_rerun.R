library(doParallel)
library(jagsUI)


files = grep(".Rdata", x = dir(), value = TRUE)

for(i in 1:length(files)){
  load(files[i])
}

load("data_setup.Rdata")

extract = function(what){invisible(Map(f = function(x,y){assign(x = x, value = y, envir = .GlobalEnv)}, x = names(what), y = what))}

extract(m5.aug[c("grid.id.m5", "grid.st.m5", "zst.m5.rand", "Slist.aug", "Y.m5.aug", "Robs.m5")])

extract(data.setup[c("xlim", "ylim", "traps", "A", "R", "K")])

M = length(grid.id.m5[[1]])

suppressWarnings(rm(data.complete, data.m5, data.setup, data.subset, subset.aug, complete.aug))

# Test Model 4 ------------------------------------------------------------------------------------

inits <- function(){list(z=zst.m5.rand[[a]],sigma=runif(1,0.2, 1), S=Slist.aug[[a]], p0=0.1 ,grid.id=grid.st.m5[[a]])}


parameters <- c("psi","p0","sigma","beta0","Ntotal")

# registerDoParallel(cores = detectCores()-1)
# getDoParWorkers()

nburn = 0
niter = 10000

jagsrun = function(a){
  
  data <- list (R=Robs.m5[[a]], grid.id=grid.id.m5[[a]], X= as.matrix(traps), K=K, Y=Y.m5.aug[[a]], bigM=M, ntraps=nrow(traps), 
                ylim=ylim, xlim=xlim)
  
  m5.out = autojags(data, inits, parameters, "model5.txt", n.thin=1, n.chains=4, n.burnin = nburn,
                    iter.increment = 200, save.all.iter = TRUE, parallel=T)
  
  # m5.out = jags(data, inits, parameters, "model5.txt", n.thin=1, n.chains=4, n.burnin = nburn, n.iter = niter, parallel = TRUE)

  
  save(m5.out, file = paste("MODEL5_OUT\\MODEL5_",a, ".Rdata",sep = ""))
  
  write(a, file = "MODEL5_OUT\\Sims_completed.txt", ncolumns = 1, append = TRUE)
}

# Check if output folder exists -- if not, create it
if (!dir.exists("MODEL4_OUT")) {
  dir.create("MODEL4_OUT")
  file.create("MODEL4_OUT\\Sims_completed.txt")
}

# Completed simulation #'s; which simulations already completed?

A.del = as.numeric(readLines("MODEL4_OUT\\Sims_completed.txt"))

# Remaining simulation #'s. This gets rid of any analyses that have already been completed (i.e. written to Sims_completed.txt) so that the analyses can be started and stopped at will. 

A.seq = (1:A)
A.seq = A.seq[!A.seq %in% A.del]

foreach(a = A.seq, .packages = "jagsUI") %dopar% jagsrun(a)
