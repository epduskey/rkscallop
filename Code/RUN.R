# Run model on simulated data

# Created: April 21, 2022
# Last modified: April 26, 2022

# Set working directory
setwd(paste(mypath, "rkscallop-main", sep = ""))

# Source scripts
# 	datsim: simulates data with variability in space and over depth
#	omega: creates penalty matrix for splines
#	cubetps_stan: sets up and runs Stan model
source("Code/datsim.R")
source("Code/omega.R")
source("Code/cubetps_stan.R")

# Load packages
library(lattice)
library(mgcv)
library(fields)
library(sp)

# Contents (ctrl-f):
#	I. Load data
#	II. Run Stan models
#	III. Run GAM+OK with hurdle
#	IV. Plot Stan output
#	V. Plot GAM+OK output
#	VI. Compare totals


########## I. Load data ##########

# Load data
dat = read.table("Data/dat.txt", header = T)

# Visualize the data
wireframe(z ~ x + y, dat, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)


########## II. Run Stan models ##########

# Knots
kd = 5
kc = 15

# Run the model on simulated data
mod.stan = cubetps(h = dat$dscale, 
	px = dat$xs, 
	xc = dat$xs, 
	yc = dat$ys, 
	x = dat$xs, 
	y = dat$ys, 
	z = dat$z, 
	Kd = kd, 
	Kc = kc, 
	nburn = 1000, 
	niter = 1000, 
	offset = rep(0,nrow(dat)), 
	smush.low = 0.01, 
	smush.hi = 0.01, 
	plot.k = T, 
	par = T)

# Save model output
# save(object = mod.stan, file = "Output/modstan.rda")

# Load model output
load("Output/modstan.rda")


########## III. Run GAM+OK with hurdle ##########

# Run 
hurdle.gam = gam(presence ~ s(xs), data = dat, family = binomial())
summary(hurdle.gam)
dat.gamok = dat
dat.gamok$ppres = predict(hurdle.gam, type = "response")

# Run GAM
cutoff = 0.5
dat.gam = gam(z ~ s(depth,k=kd) + s(xs,ys,k=kc), data = dat.gamok[dat.gamok$ppres >= cutoff, ], family = poisson())
summary(dat.gam)
dat.gamok$preds = ifelse(dat.gamok$ppres >= cutoff, 1, 0) * predict(dat.gam, newdata = dat.gamok, type = "response")


########## IV. Plot Stan output ##########

# Initialize plot
par(mar = c(3,3,3,1), mfrow = c(2,4))

# Get summary output
out = as.data.frame(mod.stan$ostan$summary())

# Get estimated median and residuals
dat.stan = dat
dat.stan$ppres = out[startsWith(out$variable,"p["), ]$median
dat.stan$preds = exp(out[grepl("log_mu",out$variable), ]$median)
dat.stan$resid = out[startsWith(out$variable,"residual"), ]$median

# Plot Stan - zero-inflation
plot(ppres ~ xs, data = dat.stan, pch = 16, ylim = c(0,1), main = "stan - zeroes")

# Plot Stan - depth
plot(preds ~ depth, data = dat.stan, pch = 16, ylim = c(0,max(dat$z)), main = "stan - depth")

# Plot GAM - space
stan.sp = xyz2img(dat.stan, xcol = 3, ycol = 4, zcol = 10)
image.plot(stan.sp, breaks = seq(0,max(dat.stan$preds),length.out=1001), col = hcl.colors(1000,"Blues",rev=T), axes = F, main = "stan - space")
box()

# Plot variogram
par(mar = c(3,3,3,1))
coordinates(dat.stan) = c("xs","ys")
stan.vars = out[out$variable %in% c("nugget","range","sill"), ]$median
stan.fit = vgm(nugget = stan.vars[1], range = stan.vars[2], psill = stan.vars[3], model = "Exp")
stan.vario = variogramLine(stan.fit, maxdist = 2)
plot(gamma ~ dist, data = stan.vario, type = 'l', lwd = 2, main = "stan - vario")


########## V. Plot GAM+OK output ##########

# Plot GAM - hurdle
plot(ppres ~ xs, data = dat.gamok, pch = 16, ylim = c(0,1), main = "gam - zeroes")

# Plot GAM - depth
plot(preds ~ depth, data = dat.gamok, pch = 16, main = "gam - depth")

# Plot GAM - space
dat.sp = xyz2img(dat.gamok, xcol = 3, ycol = 4, zcol = 10)
image.plot(dat.sp, breaks = seq(0,max(dat.gamok$preds),length.out=1001), col = hcl.colors(1000,"Blues",rev=T), axes = F, main = "gam - space")
box()

# Plot variogram
par(mar = c(3,3,3,1))
dat.gamok$resid = dat.gamok$z - predict(dat.gam, newdata = dat.gamok, type = "response")
coordinates(dat.gamok) = c("xs","ys")
gamok.bins = variogram(resid ~ 1, dat.gamok)
gamok.fit = fit.variogram(gamok.bins, vgm(nugget = min(gamok.bins$gamma), psill = max(gamok.bins$gamma), range = 0.2, model = "Exp"))
gamok.vario = variogramLine(gamok.fit, maxdist = 2)
plot(gamma ~ dist, gamok.vario, type = 'l', lwd = 2, main = "gam - vario")


########## VI. Compare totals ##########

# Totals for GAM+OK
total.gamok = dat.gamok$preds * ifelse(dat.gamok$ppres >= cutoff, 1, 0) + krige(resid ~ 1, dat.gamok, SpatialPoints(cbind(dat[,c("xs","ys")])), model = gamok.fit)$var1.pred
sum(total.gamok)

# Totals for Stan
total.stan = dat.stan$ppres * dat.stan$preds + krige(resid ~ 1, dat.stan, SpatialPoints(cbind(dat[,c("xs","ys")])), model = stan.fit)$var1.pred
sum(total.stan)
