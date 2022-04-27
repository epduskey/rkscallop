# Thin plate spline Bayesian estimation
# h: depth
# px: coordinates for probability of presence
# x: x-coordinate
# y: y-coordinate
# z: z-coordinate (i.e. surface height e.g. density)
# Kd: number of knots for depth portion
# Kc: number of knots for coordinate portion
# nchains: number of chains in Bayesian model run
# nburn: number of iterations to discard from each chain
# niter: number of iterations in each chain
# offset: area offset
# smush.low: for depth spline, location of first point normalized to entire depth range
# smush.hi: for depth spline, location of last point normalized to entire depth range
# plot.k: do you want to plot location of the knots?
# par: do you want to parallelize?
# jt: do you want to jitter your tps knots?
#
library(cmdstanr)
library(gstat)

cubetps = function(h, px, xc, yc, x, y, z, Kd, Kc, nburn, niter, offset, smush.low, smush.hi, plot.k = T, par = T, jt = F) {
	
	# Calculate distance between points for depth
	Mxh = matrix(h, nrow = length(h), ncol = length(h), byrow = T)
	Myh = matrix(h, nrow = length(h), ncol = length(h), byrow = F)
	Dh = abs(Mxh - Myh)
	
	# Get knot values for depth
	dknots = seq(min(h) + smush.low*diff(range(h)), max(h) - smush.hi*diff(range(h)), length.out = Kd)
	
	# Set up the matrix dZ_K where dZ_K[i,j] = ((h[i]-dknots[j])_+)^3
	dZK = matrix(NA, nrow = length(h), ncol = length(dknots))
	for(i in 1:length(h)) {
		for(j in 1:length(dknots)) {
			dZK[i,j] = ifelse(h[i] - dknots[j] < 0, 0, (h[i]-dknots[j])^3) 
		}
	}
	
	# Calculate penalty matrix
	omegad = pcube(dknots, maxdat = max(h))
	svd.omegad = svd(omegad)
	sqrt.omegad = t(svd.omegad$v %*% (t(svd.omegad$u)*sqrt(svd.omegad$d)))
	dZ = t(solve(sqrt.omegad, t(dZK)))
	dX = matrix(c(rep(1, length(h)), h), nrow = length(h), ncol = 2)
	dG = cbind(dX,dZ)
	
	# Calculate distances between points for the coordinates
	Mxc = matrix(x, nrow = length(x), ncol = length(x))
	Myc = matrix(y, nrow = length(y), ncol = length(y))
	Dc = sqrt((Mxc-t(Mxc))^2 + (Myc-t(Myc))^2)
	
	# Set up the knots (big K) for the coordinates, knots will be regularly spaced in the grid
	Kc.df = data.frame(x = xc, y = yc)
	Kc.clust = hclust(dist(Kc.df))
	Kc.tree = cutree(Kc.clust, k = Kc)
	cknots = matrix(NA, nrow = Kc, ncol = 2)
	for(i in 1:Kc) {
		cknots[i, ] = apply(Kc.df[Kc.tree == i, ], 2, mean)
	} 
	if(jt) {
		cknots = jitter(cknots, factor = 1)
	}

	# Set the matrix cZ_K where cZ_K[i,j] = r^2*log(r); r = sqrt((kx[i]-x[j])^2 + (ky[i]-y[j])^2)
	cZK = matrix(NA,  nrow = length(xc), ncol = nrow(cknots))
	for(i in 1:nrow(cZK)) {
		for(j in 1:ncol(cZK)) {
			r = sqrt((xc[i]-cknots[j,1])^2 + (yc[i]-cknots[j,2])^2)
			cZK[i,j] = (r^2)*log(r)
		}
	}
	if(plot.k) {par(mfrow=c(1, 2)); plot(h, z); points(dknots, rep(mean(z), length(dknots)), pch = 16, col = "red"); 
				plot(xc, yc); points(cknots[,1], cknots[,2], pch = 16, col = "red")}

	# Calculate the penalty matrix
	omegac = pplate(cknots)
	svd.omegac = svd(omegac)
	sqrt.omegac = t(svd.omegac$v %*% (t(svd.omegac$u)*sqrt(svd.omegac$d)))
	cZ = t(solve(sqrt.omegac, t(cZK)))
	cX = matrix(c(xc, yc), nrow = length(xc), ncol = 2)
	cG = cbind(cX,cZ)
	
	# Combine both matrices
	G = cbind(dG,cG)
	
	# Write and compile Stan model
	cat("
	data
	{
		// Count observations and area offset
		int<lower=1> n;	
		array[n] int<lower=0> y;
		vector[n] os;
		
		// Depth observations for hurdle model
		vector[n] px;
		
		// Presence
		array[n] int presence;
		
		// Knots
		int<lower=1> kd;
		int<lower=1> kc;
		
		// Design matrix
		array[n] vector[kd+2+kc+2] G;
		
		// Distance matrix
		matrix<lower=0>[n,n] D;
		
	}
	parameters
	{
		// Spline coefficients and variance
		vector[kd+2+kc+2] b;
		real<lower=0> dtau;
		real<lower=0> ctau;
		
		// Variogram parameters
		real<lower=0> nugget;
		real<lower=0> range;	
		real<lower=0> sill;
		
		// Local residual means and global variance
		vector[n] mu_residual;
		
		// Probability of presence depth coefficients
		real<lower=0,upper=1> a;
		real k;
		real s;
		
		// Overdispersion
		real<lower=0> phi;
	}
	transformed parameters
	{	
		// Probability of presence
		vector<lower=0,upper=1>[n] p;
		for(i in 1:n) {
			p[i] = a*exp(-square(px[i] - k)/s);
		}
		
		// Expected mean response and residual
		vector[n] log_mu;
		vector[n] residual;
		for(i in 1:n) {
			log_mu[i] = dot_product(b,G[i]) + os[i];
			residual[i] = y[i] - exp(log_mu[i]);
		}
	}
	model
	{
		// Priors
		
		// Spline model coefficient and variance priors
		b[1:2] ~ normal(0,1);
		b[3:kd+2] ~ normal(0,dtau);
		b[kd+3:kd+4] ~ normal(0,1);
		b[kd+5:kd+5-1+kc] ~ normal(0,ctau);
		dtau ~ normal(0,1) T[0,];
		ctau ~ normal(0,1) T[0,];
		phi ~ gamma(1,0.5);
		
		// Variogram priors
		nugget ~ normal(0,1) T[0,];
		range ~ normal(0,1) T[0,];
		sill ~ normal(0,1000) T[0,];
		
		// Hurdle model priors
		a ~ uniform(0,1);
		k ~ normal(0,1);
		s ~ gamma(1,5);
		
		// Residual mean and variance priors
		mu_residual ~ normal(0,1);
		
		// Residuals
		residual ~ multi_normal(mu_residual,sill*exp(-D/range)+nugget);
		
		// Likelihood
		
		// Count observations with zero-inflation
		for(i in 1:n) {
			if(y[i] == 0) {
				target += log_sum_exp(bernoulli_lpmf(1|1-p[i]), bernoulli_lpmf(0|1-p[i]) + neg_binomial_2_log_lpmf(y[i]|log_mu[i],phi));
			} else {
				target += bernoulli_lpmf(0|1-p[i]) + neg_binomial_2_log_lpmf(y[i]|log_mu[i],phi);
			}
		}
	}
	
	
	", file = "Code/cubetps.stan")
	cubetpsstan = cmdstan_model("Code/cubetps.stan")

	
	# Organize data for JAGS run
	dat = list(y = round(z),		# Surface height z for spline surface likelihood
				presence = ifelse(round(z) > 0, 1, 0),		# Presence
				n = length(z),		# Number of total points
				kc = Kc,			# Number of coordinate knots
				kd = Kd,			# Number of depth knots
				G = G,		# Restructured design matrix
				D = Dc,		# Matrix of distances between points
				os = offset,			# Area offset
				px = px		# Depth at each observation
				)
	
	# Stan model run
	ostan = cubetpsstan$sample(
		data = dat,
		chains = 6,
		parallel_chains = 6,
		refresh = 100,
		adapt_delta = 0.9,
		max_treedepth = 15,
		iter_warmup = nburn,
		iter_sampling = niter,
		init = function() {list(nugget = 10, range = 0.1, sill = 1000, dtau = 0.5, ctau = 0.5, gtau = 0.5, a = 0.9, k = mean(px), s = 0.1)}
	)
	ostan$draws()
	ostan$sampler_diagnostics()
	ostan$cmdstan_diagnose()
	ostan$summary()
	
	# Return model object
	return(list(ostan = ostan, G = G))
	
}