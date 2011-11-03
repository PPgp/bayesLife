library(bayesLife)

start.test <- function(name) cat('\n<=== Starting test of', name,'====\n')
test.ok <- function(name) cat('\n==== Test of', name, 'OK.===>\n')

test.get.wpp.data <- function(wpp.year=2008) {
	test.name <- 'getting WPP data'
	start.test(test.name)
	data <- bayesLife:::get.wpp.e0.data()
	stopifnot(length(dim(data$e0.matrix))==2)
	stopifnot(ncol(data$e0.matrix)==158)
	stopifnot(nrow(data$e0.matrix)==12)
	stopifnot(rownames(data$e0.matrix[12])=='2008')
	test.ok(test.name)
}

test.estimate.mcmc <- function() {
	sim.dir <- tempfile()
    # run MCMC
    test.name <- 'estimating MCMC'
	start.test(test.name)
    m <- run.e0.mcmc(nr.chains=1, iter=10, thin=1, output.dir=sim.dir)
    stopifnot(m$mcmc.list[[1]]$finished.iter == 10)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 10)
	test.ok(test.name)
	
	# continue MCMC
	test.name <- 'continuing MCMC'
	start.test(test.name)
	m <- continue.e0.mcmc(iter=10, output.dir=sim.dir)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 20)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 20)
	stopifnot(!is.element(900, m$meta$regions$country_code)) # 'World' should not be included
	test.ok(test.name)
	
	# run MCMC for an aggregation
	test.name <- 'estimating MCMC for extra areas'
	start.test(test.name)
	data.dir <- file.path(.find.package("bayesLife"), 'data')
	m <- run.e0.mcmc.extra(sim.dir=sim.dir, 
					my.e0.file=file.path(data.dir, 'my_e0_template.txt'), burnin=0)
	stopifnot(is.element(900, m$meta$regions$country_code)) # 'World' should be included
	test.ok(test.name)
	
	# run prediction
	test.name <- 'running projections'
	start.test(test.name)
	pred <- e0.predict(m, burnin=0, verbose=FALSE)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 20)
	stopifnot(!is.element(903, pred$mcmc.set$regions$country_code))
	npred <- dim(pred$e0.matrix.reconstructed)[2]
	test.ok(test.name)
	
	# run MCMC for another aggregation
	test.name <- 'running projections on extra areas'
	start.test(test.name)
	m <- run.e0.mcmc.extra(sim.dir=sim.dir, countries=903, burnin=0)
	# run prediction only for the area 903
	pred <- e0.predict.extra(sim.dir=sim.dir, verbose=FALSE)
	stopifnot(dim(pred$e0.matrix.reconstructed)[2] == npred+1)
	stopifnot(is.element(903, pred$mcmc.set$meta$regions$country_code))
	stopifnot(!is.null(bayesTFR:::get.trajectories(pred, 903)$trajectories))
	test.ok(test.name)
    
	test.name <- 'shifting the median'
	start.test(test.name)
    country <- 'Netherlands'
    country.idx <- get.country.object(country, m$meta)$index
    projs <- summary(pred, country=country)$projections
    e0.median.shift(sim.dir, country=country, shift=5.3, from=2051, to=2080)
	shifted.pred <- get.e0.prediction(sim.dir)
	shifted.projs <- summary(shifted.pred, country=country)$projections
	stopifnot(all(projs[10:15,c(1,3:dim(projs)[2])]+5.3 == shifted.projs[10:15,c(1,3:dim(projs)[2])]))
	stopifnot(all(projs[c(1:9, 16:19),c(1,3:dim(projs)[2])] == shifted.projs[c(1:9, 16:19),
								c(1,3:dim(projs)[2])]))
	test.ok(test.name)
	
	test.name <- 'resetting the median'
	start.test(test.name)
	shifted.pred <- e0.median.shift(sim.dir, country=country, reset = TRUE)
	shifted.projs <- summary(shifted.pred, country=country)$projections
	stopifnot(all(projs[,c(1,3:dim(projs)[2])] == shifted.projs[,c(1,3:dim(projs)[2])]))
	test.ok(test.name)
	
	test.name <- 'setting the median'
	start.test(test.name)
	expert.values <- c(90.5, 91, 93.8)
    shift <- expert.values - pred$quantiles[country.idx, '0.5',4:6] # Netherlands has index 106
	mod.pred <- e0.median.set(sim.dir, country=country, values=expert.values, years=2024)
	mod.projs <- summary(mod.pred, country=country)$projections
	stopifnot(all(mod.projs[4:6, c(1,3:dim(projs)[2])]==projs[4:6, c(1,3:dim(projs)[2])]+shift))
	stopifnot(all(mod.projs[c(1:3,7:19), c(1,3:dim(projs)[2])]==projs[c(1:3,7:19), c(1,3:dim(projs)[2])]))
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}

test.existing.simulation <- function() {
	test.name <- 'retrieving MCMC results'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
	m <- get.e0.mcmc(sim.dir, low.memory=FALSE, chain.ids=c(1,2))
	stopifnot(length(m$mcmc.list)==2)
	stopifnot(dim(m$mcmc.list[[1]]$traces)[1]==26) # because the chains are thinned by two + init value
	m <- get.e0.mcmc(sim.dir)
	summary(m)
	summary(e0.mcmc(m, 1), par.names.cs=NULL)
	stopifnot(bayesTFR:::get.total.iterations(m$mcmc.list) == 150)
	stopifnot(bayesTFR:::get.stored.mcmc.length(m$mcmc.list, burnin=30) == 30)
	test.ok(test.name)
	
	test.name <- 'retrieving projection results'
	start.test(test.name)
	pred <- get.e0.prediction(sim.dir)
	s <- summary(pred, country='Japan')
	stopifnot(s$nr.traj == 30)
	stopifnot(all(dim(s$projections)==c(19,11)))
	mb <- get.thinned.e0.mcmc(m, thin=2, burnin=30)
	s <- summary(mb, meta.only=TRUE)
	stopifnot(s$iters == 30)
	test.ok(test.name)
}

test.DLcurve <- function() {
	test.name <- 'plotting DL curves'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
	m <- get.e0.mcmc(sim.dir)
	filename <- tempfile()
	png(filename=filename)
	e0.DLcurve.plot(m, 'Slovenia')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
}

test.e0trajectories <- function() {
	test.name <- 'plotting e0 trajectories'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
	pred <- get.e0.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	png(filename=filename)
	e0.trajectories.plot(pred, 'Australia')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	test.name <- 'tabulating e0 trajectories'
	start.test(test.name)
	t <- e0.trajectories.table(pred, "Australia", pi=c(90, 80, 70, 52))
	stopifnot(all(dim(t) == c(30, 9)))
	test.ok(test.name)
}

test.plot.density <- function() {
	test.name <- 'plotting parameter density'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
	pred <- get.e0.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	png(filename=filename)
	e0.pardensity.cs.plot('Ireland', pred)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
}

test.plot.map <- function() {
	test.name <- 'creating e0 maps'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
	pred <- get.e0.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	e0.map(pred, projection.year=2098, device='png', device.args=list(filename=filename))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	test.name <- 'creating e0 maps of observed data'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
	pred <- get.e0.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	e0.map(pred, projection.year=1974, device='png', device.args=list(filename=filename))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	test.name <- 'creating parameter maps'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
	pred <- get.e0.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	e0.map(pred, par.name='z.c', device='png', device.args=list(filename=filename))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
}
	
test.get.parameter.traces <- function() {
	test.name <- 'getting parameter traces'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
	m <- get.e0.mcmc(sim.dir, low.memory=TRUE)
	traces <- get.e0.parameter.traces(m$mcmc.list, burnin=10, 
					thinning.index=c(4, 41, 59))
	stopifnot(nrow(traces)==3)
	m.check <- get.e0.mcmc(sim.dir, low.memory=FALSE, burnin=10, chain.ids=c(1,3))
	stopifnot(traces[1,'omega']==m.check$mcmc.list[[1]]$traces[4,'omega'])
	# indices 41 and 59 in the collapsed traces correspond to indices 1 and 19, respectively, in chain 3
	stopifnot(all(traces[c(2,3),'omega']==m.check$mcmc.list[[2]]$traces[c(1,19),'omega']))
	
	traces <- get.e0.parameter.traces(m$mcmc.list, burnin=10, thin=8)
	# original thin is 2, so here we thin additionally by 4 (60/4=15) 
	stopifnot(nrow(traces)==15)
	stopifnot(traces[2,'z']==m.check$mcmc.list[[1]]$traces[5,'z']) #(4+1)
	stopifnot(traces[14,'z']==m.check$mcmc.list[[2]]$traces[13,'z']) #(3*4 + 1)
	test.ok(test.name)
}

test.run.mcmc.simulation.auto <- function() {
	sim.dir <- tempfile()
	# run MCMC
	test.name <- 'running auto MCMC'
	start.test(test.name)
	m <- run.e0.mcmc(iter='auto', output.dir=sim.dir, thin=1,
					auto.conf=list(iter=10, iter.incr=5, max.loops=3, nr.chains=2, thin=1, burnin=5))
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 40)
	test.ok(test.name)

	test.name <- 'continuing auto MCMC'
	start.test(test.name)
	m <- continue.e0.mcmc(iter='auto', output.dir=sim.dir, auto.conf=list(max.loops=2))
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 60)
	test.ok(test.name)

	unlink(sim.dir, recursive=TRUE)
}
