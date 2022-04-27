library(gstat)

# Simulate 3d Gaussian data
#	a: maximum probability of presence
#	k: location where maximum probability of presence occurs
#	s: spread of probability of presence; higher numbers indicate greater spread
#	pcoord: directional coordinate over which presence varies (string "x" or "y")
# 	zmu: mean of spatial function
# 	dmu: mean of depth function
# 	zs: variance of spatial function
# 	ds: variance of depth function
# 	zmag: magnitude/multiplier of spatial function
# 	dmag: magnitude/multiplier of depth function
#	dcoord: directional coordinate over which depth varies (string "x" or "y")
# 	xy: points at which to evaluate the spatial function
# 	xc: points at which to evaluate autocorrelation (may be scaled xy to reduce magnitude of distance values)
#	xd: points at which to evaluate the depth function
# 	psill: partial sill for auto function
# 	rnge: range for auto function
# 	nugget: nugget for auto function
# 	model: variogram model type for auto function
# 	errvg: variance of white noise for auto function
# 	returns simulated density values for each coordinate value or pair in xy
gauss = function(a, k, s, pcoord = "x", zmu, dmu, zs, ds, zmag, dmag, dcoord = "x", xy, xc, xd, psill, rnge, nugget, model = "Exp", errvg) {
	
	# Zero-inflation
	pind = ifelse(pcoord == "x", 1, 2)
	probs = a * exp((-(xy[,pind] - k)^2)/s)
	presence = rbinom(nrow(xy), 1, probs)
	
	# Spatial density
	mu_space = rep(NA, nrow(xy))
	for(i in 1:length(mu_space)) {
		xi = t(t(xy[i,]))
		mu_space[i] = zmag*exp((-1/2)*(t(xi-zmu)%*%solve(zs)%*%(xi-zmu)))
	}
	
	# Depth density
	cind = ifelse(dcoord == "x", 1, 2)
	mu_depth = dmag*exp(-(xy[,cind] - dmu)^2/(2*ds^2))
	
	# Total density
	mu = mu_space*mu_depth
	
	# Mean plus autocorrelated residuals
	zres = exp(auto(xc[,1], xc[,2], psill, rnge, nugget, err = errvg))
	z = presence * rpois(nrow(xy), mu*zres)
	
	return(z)
}



# Autocorrrelated residuals
# 	x: x-coordinate
# 	y: y-coordinate
# 	psill: partial sill of variogram
# 	rnge: range of variogram
# 	nugget: nugget of variogram
# 	model: variogram model, exponential by default
# 	err: variance of white noise added to residuals
# 	returns autocorrelated residual of length = length(x)
auto = function(x, y, psill, rnge, nugget, model = "Exp", err) {
	
	# Create variogram model
	vmod = vgm(psill = psill, range = rnge, nugget = nugget, model = model)
	
	# Calculate distance between points
	Mx = matrix(x, nrow = length(x), ncol = length(x))
	My = matrix(y, nrow = length(y), ncol = length(y))
	D = sqrt((Mx-t(Mx))^2 + (My-t(My))^2)
	
	# Create variance-covariance matrix from variogram
	vcov = matrix(variogramLine(vmod, maxdist = 2*max(D), dist_vector = c(D), covariance = T)$gamma, nrow = length(x), ncol = length(x))
	LL = t(chol(vcov))
	
	# Add autocorrelation to random normal draws
	res = (LL %*% (rnorm(length(x), 0, err)))
	return(res)
}