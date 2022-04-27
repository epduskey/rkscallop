# Compute values of penalty matrix
library(Ryacas0)
x = Sym("x")
f1 = function(x) {return( x^2 )}
f2 = function(x) {return( x^3 )}
fu = function(x,kappa) {return( (x - kappa)^3 )}
f = c(f1, f2)
# pcube = function(knots, maxdat, type = "K") {
	
	# # Set x as a symbolic variable
	# x = Sym("x")
	
	# # Lower and upper integration range
	# low = c(rep(min(knots), 4), knots)
	# up = maxdat
	
	# # Results matrix
	# pmat = matrix(0, nrow = 4 + length(knots), ncol = 4 + length(knots))
	
	# for(i in 3:nrow(pmat)) {
		# for(j in 3:ncol(pmat)) {
			
			# # Set the first function
			# if(i <= 4) {
				# g1 = deriv(deriv(f[i-2][[1]](x), x), x)
			# } else {
				# g1 = deriv(deriv(fu(x, knots[i-4]), x), x)
			# }
			
			# # Set the second function
			# if(j <= 4) {
				# g2 = deriv(deriv(f[j-2][[1]](x), x), x)
			# } else {
				# g2 = deriv(deriv(fu(x, knots[j-4]), x), x)
			# }
			# symbo = as.expression(Integrate(g1*g2, x))
			# x = max(low[i], low[j])
			# intlow = eval(symbo)
			# x = up
			# intup = eval(symbo)
			# pmat[i,j] = intup - intlow
			# x = Sym("x")
		# }
	# }
	# if(type == "K") {return(pmat[5:dim(pmat)[1], 5:dim(pmat)[2]])}
	# if(type == "all") {return(pmat)}		# To return penalty of parametric portion as well
# }

pcube = function(knots, maxdat, type = "K") {
	
	# Set x as a symbolic variable
	x = Sym("x")
	
	# Lower and upper integration range
	low = c(rep(min(knots), 2), knots)
	up = maxdat
	
	# Results matrix
	pmat = matrix(0, nrow = 2 + length(knots), ncol = 2 + length(knots))
	
	for(i in 3:nrow(pmat)) {
		for(j in 3:ncol(pmat)) {
			
			# Set the first function
			g1 = deriv(deriv(fu(x, knots[i-2]), x), x)
			
			# Set the second function
			g2 = deriv(deriv(fu(x, knots[j-2]), x), x)
			
			symbo = as.expression(Integrate(g1*g2, x))
			x = max(low[i], low[j])
			intlow = eval(symbo)
			x = up
			intup = eval(symbo)
			pmat[i,j] = intup - intlow
			x = Sym("x")
		}
	}
	if(type == "K") {return(pmat[3:dim(pmat)[1], 3:dim(pmat)[2]])}
	if(type == "all") {return(pmat)}		# To return penalty of parametric portion as well
}

pplate = function(knots) {
	
	# Set up the matrix
	pmat = matrix(NA, nrow = nrow(knots), ncol = nrow(knots))
	
	# Calculate pairwise distances between knots
	for(i in 1:nrow(pmat)) {
		for(j in 1:ncol(pmat)) {
			pmat[i,j] = sqrt((knots[i,1]-knots[j,1])^2 + (knots[i,2]-knots[j,2])^2)
		}
	}
	
	# Return the penalty (note this matrix is only proportional to the penalty matrix)
	return(pmat)
}