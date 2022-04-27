# Simulate data

# Created: April 21, 2022
# Last modified: April 26, 2022

# Source scripts
source("Code/simfuncs.R")

# Contents (ctrl-f):
#	I. Choose parameters
#	II. Simulate data
#	III. Visualize


########## I. Choose parameters ##########

# Set seed
set.seed(1509)

# Number of unique x-values
nx = 15

# Number of unique y-values	
ny = 15

# Data scale	
sc = 100

# x-coordinates	
x = rep(seq(1,100,length.out=nx), ny)*sc
xscale = scale(x)

# y-coordinates
y = rep(seq(1,100,length.out=ny), each = nx)*sc
yscale = scale(y)

# Depth varies in x direction
depth = -1*unique(c(xscale))^2 + 2*max(xscale)^2 + rnorm(nx*ny)
xd = scale(depth)

# All points
xy = matrix(c(xscale,yscale), nrow = length(xscale), ncol = 2)
xc = matrix(c(xscale,yscale), nrow = length(x), ncol = 2)

# Choose zero-inflation parameters
a = 0.8
k = 0
s = 0.5

# Choose variogram parameters
psill = 200
rnge = sqrt((max(xscale)-min(xscale))^2+(max(yscale)-min(yscale))^2)/20
nugget = 5

# Choose spatial mean
zmag = 5
zmu = matrix(0, nrow = 2, ncol = 1)

# Choose depth mean
dmag = 5
dmu = 0

# Choose sigma values
zs = 0.4*diag(2)
ds = 0.8


########## II. Simulate data ##########

# Simulate
dat = data.frame(x = x, 
	y = y, 
	xs = xscale, 
	ys = yscale, 
	depth = depth,
	dscale = xd,
	z = gauss(a, k, s, "x", zmu, dmu, zs, ds, zmag, dmag, "x", xy, xc, xd, psill, rnge, nugget, errvg = 0.1))

# Add presence
dat$presence = ifelse(dat$z == 0, 0, 1)

# Save data
write.table(dat, "Data/dat.txt")


########## III. Visualize ##########

# Use wireframe to plot counts over space
# wireframe(z ~ x + y, dat, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)

# Plot counts over depth
# plot(z ~ depth, data = dat)