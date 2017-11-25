# script for LCA in diagnostic testing
#
# uses Gibbs steps to estimate unknown prevalence, sensitivity and specificity from
# table data.

# Data
#
# number of populations
P <- 2
# number of tests
T <- 2
# counts
#      1++,1+-,1-+,1--,2++,2+-,2-+,2--
#counts <- c(500,100, 50,400,400,120, 10,100)
#           ++,+-,-+,--
counts <- c(500,100, 50,400,400,120, 10,100)

# Priors, assuming Beta distribution with parameters alpha, beta.
p_alpha <- rep(1, P)
p_beta  <- rep(1, P)
s_alpha <- rep(1, T)
s_beta  <- rep(1, T)
c_alpha <- rep(1, T)
c_beta  <- rep(1, T)

# MCMC control
iters  <- 10000
burnin <- 1000
thin   <- 10

# Now run the MCMC model...

# Simple test for data entered correctly
if (length(counts) != P*(2^T)) {
	stop("counts must be of length (number of populations)*2^(number of tests)")
}

# Setup matrices for y and m
y <- matrix(counts, nrow=P, ncol=2^T, byrow=T)
m <- matrix(0, nrow=T, ncol=2^T)
for (i in 1:T)
{
	for (k in 1:(2^T))
	{
		# does 2^T-k have bit T-i set?
		r <- (2^T-k) %/% (2^(T-i));
		if (r %% 2)
			m[i,k] <- 1
	}
}

# initial values
p <- rep(0.5, P) # prevalence
s <- rep(0.5, T) # test sensitivity
c <- rep(0.5, T) # test specificity

# latent true positives for each y
x <- y

# posterior
post <- matrix(NA, iters/thin, P+2*T)
post_names <- rep("",P+2*T)
for (i in 1:P)
	post_names[i] <- paste("Pop",i,"Prev",sep="_")
for (i in 1:T)
	post_names[P+i] <- paste("Test",i,"Sens",sep="_")
for (i in 1:T)
	post_names[P+T+i] <- paste("Test",i,"Spec",sep="_")
colnames(post) <- post_names

# run the MCMC
j <- 1
for (l in 1:(iters + burnin))
{
 	# Gibbs steps:
	# 1. Sample from x_{iM}
	for (i in 1:P) {
		for (k in 1:(2^T)) {
			p1 <- p[i];
			p2 <- (1-p[i]);
			for (t in 1:T) {
				if (m[t,k]) {
					p1 <- p1*s[t];
					p2 <- p2*(1-c[t]);
				} else {
					p1 <- p1*(1-s[t]);
					p2 <- p2*c[t];
				}
			}
			x[i,k] = rbinom(1, y[i,k], p1/(p1+p2))
		}
	}
	# 2. Sample from p_i
	for (i in 1:P) {
#		sx <- 0;
#		sy <- 0;
#		for (k in 1:(2^T)) {
#			sx <- sx + x[i,k];
#			sy <- sy + y[i,k];
#		}
		# optimisation:
		sx <- sum(x[i,])
		sy <- sum(y[i,])
		p[i] <- rbeta(1, p_alpha[i] + sx, p_beta[i] + sy - sx);
	}
	# 3. Sample from s_t
	for (t in 1:T) {
		alpha <- s_alpha[t];
		beta  <- s_beta[t];
		for (i in 1:P) {
			# optimisation
			alpha <- alpha + sum(m[t,]*x[i,])
			beta  <- beta  + sum((1-m[t,])*x[i,])
#			for (k in 1:(2^T)) {
#				if (m[t,k]) {
#					alpha <- alpha + x[i,k];
#				} else {
#					beta <- beta + x[i,k];
#				}
#			}
		}
		s[t] <- rbeta(1, alpha, beta);
	}
	# 4. Sample from c_t
	for (t in 1:T) {
		alpha <- c_alpha[t];
		beta  <- c_beta[t];
		for (i in 1:P) {
			# optimisation
			alpha <- alpha + sum((1-m[t,])*(y[i,]-x[i,]))
			beta  <- beta  + sum(m[t,]*(y[i,]-x[i,]))
#			for (k in 1:(2^T)) {
#				if (m[t,k]) {
#					beta <- beta + y[i,k] - x[i,k];
#				} else {
#					alpha <- alpha + y[i,k] - x[i,k];
#				}
#			}
		}
		c[t] <- rbeta(1, alpha, beta);
	}
	# store sample
	if (l > burnin && (l-burnin) %% thin == 0) {
		post[j,] <- c(p, s, c)
		j <- j + 1;
	}
}

# Plot posteriors
par(mfrow=c(3,2))
plot(post[,1], type="l", ylim=c(0,1), main="Prevalence traces")
if (P > 1) {
	for (i in 2:P)
		lines(post[,i], col=i)
}
d <- list(); xlim = NULL; ylim = NULL
for (i in 1:P) {
	d[[i]] <- density(post[,i])
	xlim = range(xlim, d[[i]]$x)
	ylim = range(ylim, d[[i]]$y)
}
plot(d[[1]], ylim=ylim, xlim=xlim, main="Prevalence density")
if (P > 1) {
	for (i in 2:P) {
		lines(d[[i]], col=i)
	  # and the prior :)
	  x = seq(xlim[1], xlim[2], length.out=100)
	  y = dbeta(x, p_alpha[i], p_beta[i])
	  lines(x, y, col=i, lty="dotted")
	}
}

plot(post[,P+1], type="l", ylim=c(0,1), main="Sensitivity traces")
if (T > 1) {
	for (t in 2:T)
		lines(post[,P+t], col=t)
}
d <- list(); xlim = NULL; ylim = NULL
for (t in 1:T) {
	d[[t]] <- density(post[,P+t])
	xlim = range(xlim, d[[t]]$x)
	ylim = range(ylim, d[[t]]$y)
}
plot(d[[1]], ylim=ylim, xlim=xlim, main="Sensitivity density")
if (T > 1) {
	for (t in 2:T) {
		lines(d[[t]], col=t)
	  # and the prior :)
	  x = seq(xlim[1], xlim[2], length.out=100)
	  y = dbeta(x, s_alpha[t], s_beta[t])
	  lines(x, y, col=t, lty="dotted")
	}
}

plot(post[,P+T+1], type="l", ylim=c(0,1), main="Specificity traces")
if (T > 1) {
	for (t in 2:T)
		lines(post[,P+T+t], col=t)
}
d <- list(); xlim = NULL; ylim = NULL
for (t in 1:T) {
	d[[t]] <- density(post[,P+T+t])
	xlim = range(xlim, d[[t]]$x)
	ylim = range(ylim, d[[t]]$y)
}
plot(d[[1]], ylim=ylim, xlim=xlim, main="Specificity density")
if (T > 1) {
	for (t in 2:T) {
		lines(d[[t]], col=t)
	  x = seq(xlim[1], xlim[2], length.out=100)
	  y = dbeta(x, s_alpha[t], s_beta[t])
	  lines(x, y, col=t, lty="dotted")
	}
}

# summary
tab <- apply(post, 2, function(x) { c(mean(x), sd(x), quantile(x, c(0.025, 0.975))) })
rownames(tab) <- c("mean", "sd", "2.5%", "97.5%")
print(tab)

# NOTES ON HOW IT ALL WORKS...

# Let M = (01001) be a binary string representing a test result across multiple tests, where
# a 1 specifies a positive test and a 0 specifies a negative.

# Let y_{iM} be the test result from pattern M in population i.
# Let p_{i} be the prevalence on populationi i.
# Let s_{t} be the test sensitivity of test t.
# Let c_{t} be the test specificity of test t.

# Then the likelihood of y given p, s and c may be written
#
# L(y | p, s, c) = \prod_{i=1}^P \prod_M \in Z_T [p_i\prod_{t=1}^T s_t^M_t(1-s_t)^(1-M_t) + (1-p_i)\prod_{t=1}^T (1-c_t)^M_tc_t^(1-M_t)]^y_{iM}.
# given that we don't know the actual latent times $x_i$
#
# where Z_{tM} = 0 if M contains a 1 at position t, 1 if M contains a 1 at position t.
#
# To proceed, let x_{iM} be the number of true positives that report test pattern M in population i.
# Then x_{iM} are unknown, but can be incorporated into the likelihood at latent variables to setup an
# appropriate gibbs sampler. The likelihood becomes
#
# L(y | p, s, c) = \prod_{i=1}^P \prod_M \in Z_T \sum_{x_{iM}=0}^{y_{iM}} [p_i\prod_{t=1}^T s_t^{M_t}(1-s_t)^{1-M_t}]^{x_{iM}}[(1-p_i)\prod_{t=1}^T (1-c_t)^{M_t}c_t^{1-M_t}]^{y_{iM}-x_{iM}}
# conditional on the $x_{iM}$'s then,
#
# L(y | p, s, c, x) = \prod_{i=1}^P \prod_M \in Z_T [p_i\prod_{t=1}^T s_t^{M_t}(1-s_t)^{1-M_t}]^{x_{iM}}[(1-p_i)\prod_{t=1}^T (1-c_t)^{M_t}c_t^{1-M_t}]^{y_{iM}-x_{iM}}
#                   = \prod_{i=1}^P \prod_M \in Z_T p_i^{x_iM} \prod_{t=1}^T s_t^{M_tx_{iM}}(1-s_t)^{(1-M_t)^{x_{iM}} (1-p_i)^{y_iM}-x_{iM}} \prod_{t=1}^T (1-c_t)^{M_t(y_{iM}-x_{iM})}c_t^{(1-M_t)(y_{iM}-x_{iM)}.
#                   = \prod_{i=1}^P p_i^{\sum_M x_iM} (1-p_i)^{\sum_M y_{iM}-x_{iM}} \prod_{t=1}^T s_t^{\sum_M M_tx_{iM}}(1-s_t)^{\sum_M (1-M_t)^{x_{iM}} \prod_{t=1}^T (1-c_t)^{\sum_M M_t(y_{iM}-x_{iM})}c_t^{\sum_M (1-M_t)(y_{iM}-x_{iM)}.
#
# Assuming Beta priors, the conditional posterior of p, s, c are also Beta:
#
# p(p_i | s, c, x, y) \sim Beta(\sum_M x_{iM} + \alpha_{p_i}, \sum_M y_{iM}-x_{iM} + \beta_{p_i})
# 
# p(s_t | p, c, x, y) \sim Beta(\sum_M\sum_i M_t x_{iM} + \alpha_{s_t}, \sum_M\sum_i (1-M_t)x_{iM} + \beta_{s_t})
#
# p(c_t | s, c, x, y) \sim Beta(\sum_M\sum_i (1-M_t)(y_{iM}-x_{iM}) + \alpha_{c_t}, \sum_M\sum_i M_t(y_{iM}-x_{iM}) + \beta_{c_t})
#
# Finally, the conditional posterior of x is Binomial:
#
# p(x_{iM} | p, s, c, y) \sim Binomial(y_{iM}, \frac{a_{iM}}{a_{iM} + b_{iM}})
#
# where
#  a_{iM} = p_i\prod_{t=1}^T s_t^{M_t}(1-s_t)^{1-M_t}
#  b_{iM} = (1-p_i)\prod_{t=1}^T (1-c_t)^{M_t}c_t^{1-M_t}
#

