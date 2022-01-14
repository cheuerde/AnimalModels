library(pacman)
p_load(brms, cmdstanr, BGLR)

# only have to be done once 
# install_cmdstan(cores = 6)

# load the BGLR data
data(wheat)

# marker covariates
M <- wheat.X

sids <- colnames(M)

# phenotypes
Y <- wheat.Y

mids <- rownames(Y)

# make a data frame
d <- data.frame(
		id = as.factor(mids),
		yOne = Y[,1],
		yTwo = Y[,2]
)

# make G
mus <- colMeans(M)
Mc <- sweep(M, 2, FUN = "-", mus)
vars <- apply(M, 2, var)

G <- tcrossprod(Mc) / sum(vars)
diag(G) <- diag(G) + 0.01
colnames(G) <- rownames(G) <- mids

##################################################
### create stand model - just for sanity check ###
##################################################

# the model formula
brmf <- brmsformula(
		    mvbind(
			   yOne
			   ) ~ 1 + 
		    (1|gr(id, cov = G)) # this is the random animal effect with G as covariance structure
)

stanCode <- make_stancode(
			  formula = brmf,
			  data = d,
			  data2 = list(G = G),
			  family = gaussian()
)

stanData <- make_standata(
			  formula = brmf,
			  data = d,
			  data2 = list(G = G),
			  family = gaussian()
)

# check the stan code
# this is particularly useful to check the priors that brms imposed on the parameters.
# the default prior choice by brms follows the stan recommendations for the most part and
# are excellent!
stanCode

##########################################
### brms run - univariate animal model ###
##########################################

# stan needs significantly fewer iterations - HMC NUTS sampler explores the
# posterior extremely efficiently
niter = 1000
chains = 5
cores = 5
threads = 1

mod <- brm(
	   formula = brmf,
	   data = d, 
	   data2 = list(G = G),
	   family = gaussian(),
	   chains = chains,
	   cores = cores,
	   threads = threading(threads),
	   backend = "cmdstanr",
	   iter = niter
)

# extract posterior of variance components
vc <- VarCorr(mod, summary = FALSE)

# get h2 - note: we square because those are sds (only the case for univariate models)
vA <- vc[["id"]]$sd[,1]**2
vE <- vc[["residual__"]]$sd[,1]**2

# this is the posterior of h2
h2 <- vA / (vA + vE)

# get h2 and sd
h2_mean <- mean(h2)
h2_sd <- sd(h2)

# get the breeding values
up = ranef(mod, summary = FALSE)
upA <- up$id[,,1]

# the posterior means of breeding balues
u <- colMeans(upA)

# the PEVs
PEV <- apply(upA, 2, var)

# the reliabilities
REL <- 1 - (PEV / mean(vA))
