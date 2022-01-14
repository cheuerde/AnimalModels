# Claas Heuer, June 2021
# 
# Illustration of using brms/stan for running multivariate repeatability/animal
# models.
# NOTE: this example assumes identical model for all traits and also identical response
# distributions. That does not have to be the case - the set of fixed and random effects
# can be trait dependent, even in a multivariate setting as well as mutliple response distributions
# in a single multivaraite model. For some overview, see: https://cran.r-project.org/web/packages/brms/vignettes/brms_multivariate.html

library(pacman)
p_load(brms, cmdstanr, BGLR, pedigreemm, tidyverse)

# only have to be done once 
# install_cmdstan(cores = 6)

# load the BGLR data
data(milk)

# get pedigree and make A
P <- pedCows

mids <- P@label

# make subset
midsIn <- as.character(unique(milk$id))[1:300]

L <- as(t(relfactor(P)), "dgCMatrix")

A <- tcrossprod(L[match(midsIn, mids),])
rownames(A) <- colnames(A) <- midsIn

# make subset of data
D <- milk[milk$id %in% midsIn,]

D$idA <- factor(D$id, levels = midsIn)
D$idP <- as.character(D$id)

##################################################
### create stand model - just for sanity check ###
##################################################

# the model formula
brmf <- brmsformula(
		    mvbind(
			   milk,
			   fat
			   ) ~ 1 + 
		    factor(lact) +
		    herd +
		    (1|g|gr(idA, cov = A)) + # additive genetic effect
		    (1|p|gr(id)) # PE effect
)

stanCode <- make_stancode(
			  formula = brmf,
			  data = D,
			  data2 = list(A = A),
	                  # family = poisson("log"), # for a poisson response with log link function
			  family = gaussian()
)

stanData <- make_standata(
			  formula = brmf,
			  data = D,
			  data2 = list(A = A),
	                  # family = poisson("log"), # for a poisson response with log link function
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

# more iterations for a more complex model
niter = 2000
chains = 5
cores = 5
threads = 1

mod <- brm(
	   formula = brmf,
	   data = D, 
	   data2 = list(A = A),
	   family = gaussian(),
	   # family = poisson("log"), # for a poisson response with log link function
	   chains = chains,
	   cores = cores,
	   threads = threading(threads),
	   backend = "cmdstanr",
	   iter = niter
)

# extract posterior of covariance-matrices
vcp <- VarCorr(mod, summary = FALSE)

# additive genetic covariance matrix
G <- vcp$idA$cov

# pe covariance matrix
PE <- vcp$id$cov

# residucal covariance matrix
R <- vcp$residual$cov

# phenotypic
P <- G + PE + R

# extract the variance components for the traits
vA <- matrix(as.numeric(NA), nrow = nrow(G), ncol = 2)
vP <- matrix(as.numeric(NA), nrow = nrow(G), ncol = 2)

colnames(vA) <- colnames(vP) <- c("milk", "fat")

# fill matrices with diagonal elements of the var-covar matrix
for(i in 1:nrow(G)) vA[i,] <- diag(G[i,,])
for(i in 1:nrow(G)) vP[i,] <- diag(P[i,,])

# h2
h2 <- colMeans(vA / vP)
h2_sd <- apply(vA / vP, 2, sd)

h2Out <- data.frame(
		Trait = c("milk", "fat"),
		h2 = h2,
		h2_sd = h2_sd
		)

vAmeans <- colMeans(vA)
vPmeans <- colMeans(vP)


#######################
### Breeding Values ###
#######################

up = ranef(mod, summary = FALSE)
upA = up$idA

ids <- rownames(upA[1,,])
idsIn <- ids

out <- list()

# make some predictions
for(i in 1:length(idsIn)) {

	tmp <- tibble(
		      id = idsIn[i]
	)

	# now get the blups for this animals
	thisUp <- upA[,match(idsIn[i], upIds),,drop = TRUE]

	# make prediction
	tmp <- tmp %>%
		mutate(
		       milkPred = mean(thisUp[,"milk_Intercept"]),
		       milkPred_PEV = var(thisUp[,"milk_Intercept"]),
		       milkPred_REL = 1 - (milkPred_PEV / vAmeans["milk"]),

		       fatPred = mean(thisUp[,"fat_Intercept"]),
		       fatPred_PEV = var(thisUp[,"fat_Intercept"]),
		       fatPred_REL = 1 - (fatPred_PEV / vAmeans["fat"]),
		       )

	out[[i]] <- tmp

}

datPred <- bind_rows(out)


