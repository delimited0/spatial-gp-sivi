library(spBayes)

cov_func <- function(loc, p, sig) sig * exp(-fields::rdist(loc) / p)
cov_dist_func <- function(dists, p, sig) sig * exp(- dists / p)

set.seed(1)
n <- 200
s <- matrix(runif(2*n, 0, 1), n, 2)
s_dists <- fields::rdist(s)

X <- cbind(1, s)
coeff <- c(1, -5, 10)

tau2 <- 1
sigma2 <- 1
phi <- 5

w <- MASS::mvrnorm(1, rep(0, n), cov_func(s, phi, sigma2))
eps <- MASS::mvrnorm(1, rep(0, n), tau2 * diag(n))
y <- X %*% coeff + w + eps

df <- data.frame(y = y, lon = s[, 1], lat = s[, 2])
model_bayes <- spBayes::spLM(y ~ ., data = df, coords = s, 
                             starting = list("phi" = .5, "sigma.sq" = 1, "tau.sq" = .1),
                             tuning = list("phi" = .1, "sigma.sq" = 0, "tau.sq" = .05),
                             priors = list("phi.Unif" = c(.01, 10), 
                                           "sigma.sq.IG" = c(2, 1),
                                           "tau.sq.IG" = c(2, 1), 
                                           "beta.norm" = list(mean = rep(0, 3),
                                                              cov = 5 * diag(3))),
                             cov.model = "exponential",
                             n.samples = 5000, verbose = TRUE,
                             n.report = 1000)
model_rec <- spBayes::spRecover(model_bayes)
summary(model_rec$p.theta.recover.samples)
summary(model_rec$p.beta.recover.samples)

plot(model_rec$p.theta.recover.samples)
plot(model_rec$p.beta.recover.samples)

plot(tau.sq ~ phi, data = model_rec$p.theta.samples)

