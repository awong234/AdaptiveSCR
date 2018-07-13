library(doParallel)
library(jagsUI)


files = grep(".Rdata", x = dir(), value = TRUE)

for(i in 1:length(files)){
  load(files[i])
}

digest::digest(complete.aug)
digest::digest(data.complete)

extract = function(what){invisible(Map(f = function(x,y){assign(x = x, value = y, envir = .GlobalEnv)}, x = names(what), y = what))}

extract(complete.aug[c("grid.id.complete", "grid.st.complete", "zst.complete.1", "zst.complete.rand", "zst.complete.0", "Slist.aug", "Y.complete.aug")])

extract(data.setup[c("xlim", "ylim", "traps", "A", "R", "K")])

covariate = data.complete$covmat

M = length(grid.id.complete[[1]])

suppressWarnings(rm(data.complete, data.m5, data.setup, data.subset, subset.aug, complete.aug))

inits <- function(){list(z = zst.complete.1[[a]], sigma = runif(1,0.2, 1), S = Slist.aug[[a]], p0 = 0.1, grid.id = grid.st.complete[[a]], beta.cov = 0, beta0 = 0)}

# -------------------------------------------

parameters <- c("psi","p0","sigma","beta0","Ntotal","beta.cov") 

# registerDoParallel(cores = detectCores()-1)
# getDoParWorkers()

nburn = 1000
niter = 10000

jagsrun = function(a){
  
  data <- list (R=R, grid.id=grid.id.complete[[a]], X= as.matrix(traps), K=K, Y=Y.complete.aug[[a]], bigM=M, ntraps=nrow(traps), 
                ylim=ylim, xlim=xlim, covariate=unname(covariate[,a]))
  
  m3.out = jags(data, inits, parameters, "model3.txt", n.thin=1, n.chains=4,
                n.burnin = nburn, n.iter = niter, parallel = F)
  
  save(m3.out, file = paste("MODEL3_OUT\\MODEL3_",a, ".Rdata",sep = ""))
  write(a, file = "MODEL3_OUT\\Sims_completed.txt", ncolumns = 1, append = TRUE)
}

# foreach(a = 1) %do% jagsrun(a)


# Check if output folder exists -- if not, create it
if (!dir.exists("MODEL3_OUT")) {
  dir.create("MODEL3_OUT")
  file.create("MODEL3_OUT\\Sims_completed.txt")
}

# Completed simulation #'s
A.del = as.numeric(readLines("MODEL3_OUT\\Sims_completed.txt"))

# Remaining simulation #'s
A.seq = (1:A)

A.seq = A.seq[!A.seq %in% A.del]

foreach(a = A.seq, .packages = "jagsUI") %dopar% jagsrun(a)
