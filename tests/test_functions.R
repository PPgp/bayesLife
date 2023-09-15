library(bayesLife)

start.test <- function(name, wpp.year = NULL) cat('\n<=== Starting test of', name, if(!is.null(wpp.year)) paste0('(WPP ', wpp.year, ')') else '', '====\n')
test.ok <- function(name) cat('\n==== Test of', name, 'OK.===>\n')

test.get.wpp.data <- function(wpp.year=2010) {
	test.name <- 'getting WPP data'
	start.test(test.name, wpp.year)
	ncountries <- list('2008'=158, '2010'=159, '2012'=162, '2015'=178)
	data <- bayesLife:::get.wpp.e0.data(wpp.year=wpp.year, present.year=if(wpp.year>2012) 2015 else 2010)
	stopifnot(length(dim(data$e0.matrix))==2)
	stopifnot(ncol(data$e0.matrix)==ncountries[[as.character(wpp.year)]])
	stopifnot(nrow(data$e0.matrix)==if(wpp.year>2012) 13 else 12)
	test.ok(test.name)
}

test.estimate.mcmc <- function(compression='None', wpp.year = 2019) {
	sim.dir <- tempfile()
    # run MCMC
    test.name <- 'estimating MCMC'
	start.test(test.name, wpp.year)
    m <- run.e0.mcmc(nr.chains = 1, iter = 10, thin = 1, output.dir = sim.dir, 
                     compression.type = compression, wpp.year = wpp.year,
                     mcmc.options = list(buffer.size = 5))
    stopifnot(m$mcmc.list[[1]]$finished.iter == 10)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 10)
	stopifnot(m$meta$mcmc.options$buffer.size == 5)
	ncountries <- nrow(get.countries.table(m))
	if(wpp.year == 2019)
	    stopifnot(ncountries == 180) # in include_2019 there are 180 default countries for e0 
	test.ok(test.name)
	
	# continue MCMC
	test.name <- 'continuing MCMC'
	start.test(test.name, wpp.year)
	m <- continue.e0.mcmc(iter=10, output.dir=sim.dir)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 20)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 20)
	stopifnot(!is.element(900, m$meta$regions$country_code)) # 'World' should not be included
	test.ok(test.name)
	
	# run MCMC for an aggregation
	test.name <- 'estimating MCMC for extra areas'
	start.test(test.name, wpp.year)
	data.dir <- file.path(find.package("bayesLife"), 'extdata')
	m <- run.e0.mcmc.extra(sim.dir=sim.dir, 
					my.e0.file=file.path(data.dir, 'my_e0_template.txt'), burnin=0)
	stopifnot(is.element(900, m$meta$regions$country_code)) # 'World' should be included
	test.ok(test.name)
	
	# run prediction
	test.name <- 'running projections'
	start.test(test.name, wpp.year)
	pred <- e0.predict(m, burnin=0, verbose=FALSE)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 20)
	stopifnot(!is.element(903, pred$mcmc.set$regions$country_code))
	stopifnot(all(dim(pred$joint.male$quantiles) == dim(pred$quantiles)))
	stopifnot(dim(pred$joint.male$quantiles)[3] == 17)
	npred <- dim(pred$e0.matrix.reconstructed)[2]
	t <- e0.trajectories.table(pred, "Australia", pi=80, both.sexes=TRUE)
	stopifnot(all(dim(t) == c(30, 3)))
	smpred <- summary(get.e0.jmale.prediction(pred))
	stopifnot(smpred$nr.traj == 20)
	test.ok(test.name)
	
	# run MCMC for another aggregation
	test.name <- 'running projections on extra areas'
	start.test(test.name, wpp.year)
	m <- run.e0.mcmc.extra(sim.dir=sim.dir, countries=903, burnin=0)
	# run prediction only for the area 903
	pred <- e0.predict.extra(sim.dir=sim.dir, verbose=FALSE)
	stopifnot(dim(pred$e0.matrix.reconstructed)[2] == npred+1)
	stopifnot(is.element(903, pred$mcmc.set$meta$regions$country_code))
	stopifnot(!is.null(bayesTFR:::get.trajectories(pred, 903)$trajectories))
	stopifnot(all(dim(pred$joint.male$quantiles) == dim(pred$quantiles)))
	stopifnot(dim(pred$joint.male$quantiles)[1] == (ncountries + 2))
	test.ok(test.name)
    
	test.name <- 'shifting the median'
	start.test(test.name, wpp.year)
    country <- 'Netherlands'
    country.idx <- get.country.object(country, m$meta)$index
    projs <- summary(pred, country=country)$projections
    e0.median.shift(sim.dir, country=country, shift=5.3, from=2051, to=2080)
	shifted.pred <- get.e0.prediction(sim.dir)
	shifted.projs <- summary(shifted.pred, country=country)$projections
	stopifnot(all(projs[8:13,c(1,3:dim(projs)[2])]+5.3 == shifted.projs[8:13,c(1,3:dim(projs)[2])]))
	stopifnot(all(projs[c(1:7, 14:17),c(1,3:dim(projs)[2])] == shifted.projs[c(1:7, 14:17),
								c(1,3:dim(projs)[2])]))
	test.ok(test.name)
	
	test.name <- 'resetting the median'
	start.test(test.name, wpp.year)
	shifted.pred <- e0.median.shift(sim.dir, country=country, reset = TRUE)
	shifted.projs <- summary(shifted.pred, country=country)$projections
	stopifnot(all(projs[,c(1,3:dim(projs)[2])] == shifted.projs[,c(1,3:dim(projs)[2])]))
	test.ok(test.name)
	
	test.name <- 'setting the median'
	start.test(test.name, wpp.year)
	expert.values <- c(90.5, 91, 93.8)
    shift <- expert.values - pred$quantiles[country.idx, '0.5',2:4] # Netherlands has index 106
	mod.pred <- e0.median.set(sim.dir, country=country, values=expert.values, years=2024)
	mod.projs <- summary(mod.pred, country=country)$projections
	stopifnot(all(mod.projs[2:4, c(1,3:dim(projs)[2])]==projs[2:4, c(1,3:dim(projs)[2])]+shift))
	stopifnot(all(mod.projs[c(1,5:17), c(1,3:dim(projs)[2])]==projs[c(1,5:17), c(1,3:dim(projs)[2])]))
	test.ok(test.name)
	
	test.name <- 'converting trajectories'
	start.test(test.name, wpp.year)
	convert.e0.trajectories(sim.dir, n=10)
	test.ok(test.name)
	
	test.name <- 'shifting medians to WPP'
	e0.shift.prediction.to.wpp(sim.dir)
	e0.shift.prediction.to.wpp(sim.dir, joint.male = TRUE)
	shifted.predF <- get.e0.prediction(sim.dir)
	shifted.predM <- get.e0.prediction(sim.dir, joint.male = TRUE)
	cntry <- 'Angola'
	shifted.projsF <- summary(shifted.predF, country = cntry)$projections
	shifted.projsM <- summary(shifted.predM, country = cntry)$projections
	shifted.projsF <- data.table::data.table(shifted.projsF)[, year := as.integer(rownames(shifted.projsF))]
	shifted.projsM <- data.table::data.table(shifted.projsM)[, year := as.integer(rownames(shifted.projsM))]
	e <- new.env()
	data("e0Fproj", "e0Mproj", package = paste0("wpp", wpp.year), envir = e)
	wppl <- merge(data.table::melt(data.table::data.table(e$e0Fproj), id.vars = c("country_code", "name"), 
	                            variable.name = "period", value.name = "e0F"),
	              data.table::melt(data.table::data.table(e$e0Mproj), id.vars = c("country_code", "name"), 
	                               variable.name = "period", value.name = "e0M"),
	              by = c("country_code", "name", "period"))
	wppl <- wppl[, year := as.integer(substr(period, 1,4)) + 3][name == cntry]
	datF <- merge(wppl, shifted.projsF[, c("year", "50%"), with = FALSE], by = "year")
	datM <- merge(wppl, shifted.projsM[, c("year", "50%"), with = FALSE], by = "year")
	stopifnot(all.equal(datF$e0F, datF[, `50%`]))
	stopifnot(all.equal(datM$e0M, datM[, `50%`]))
	stopifnot(!is.null(shifted.predF$median.shift) && !is.null(shifted.predM$median.shift))
	stopifnot(length(shifted.predF$median.shift) == nrow(get.countries.table(shifted.predF)) && 
	              length(shifted.predM$median.shift) == nrow(get.countries.table(shifted.predM)) )
	test.ok(test.name)
	
	test.name <- 'resetting all countries'
	e0.median.reset(sim.dir)
	new.predF <- get.e0.prediction(sim.dir)
	new.predM <- get.e0.prediction(sim.dir, joint.male = TRUE)
	stopifnot(is.null(new.predF$median.shift) && !is.null(new.predM$median.shift))
	e0.median.reset(sim.dir, joint.male = TRUE)
	new.predM2 <- get.e0.prediction(sim.dir, joint.male = TRUE)
	stopifnot(is.null(new.predM2$median.shift))
	test.ok(test.name)
	
	unlink(sim.dir, recursive=TRUE)
}


test.estimate.mcmc.with.suppl.data <- function(compression='None') {
	sim.dir <- tempfile()
    # run MCMC
    test.name <- 'estimating MCMC using supplemental data'
	start.test(test.name)
    m <- run.e0.mcmc(nr.chains = 1, iter = 30, thin = 1, output.dir = sim.dir, 
                     start.year = 1750, seed = 1, compression.type = compression,
                     mcmc.options = list(buffer.size = 10))
    stopifnot(length(m$meta$suppl.data$regions$country_code) == 29)
	stopifnot(all(dim(m$meta$suppl.data$e0.matrix) == c(40, 29)))
	test.ok(test.name)
	
	# continue MCMC
	test.name <- 'continuing MCMC with supplemental data'
	start.test(test.name)
	m <- continue.e0.mcmc(iter=10, output.dir=sim.dir)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 40)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 40)
	stopifnot(!is.element(900, m$meta$regions$country_code)) # 'World' should not be included
	test.ok(test.name)
	
	# run MCMC for an aggregation
	test.name <- 'estimating MCMC for extra areas with supplemental data'
	start.test(test.name)
	data.dir <- file.path(find.package("bayesLife"), 'extdata')
	m <- run.e0.mcmc.extra(sim.dir = sim.dir, 
					my.e0.file = file.path(data.dir, 'my_e0_template.txt'), 
					burnin = 0)
	stopifnot(is.element(900, m$meta$regions$country_code)) # 'World' should be included
	stopifnot(m$meta$mcmc.options$buffer.size == 10) # inherited from main run
	test.ok(test.name)
	
	# run prediction
	test.name <- 'running projections for simulation with supplemental data'
	start.test(test.name)
	pred <- e0.predict(m, burnin=10, verbose=FALSE, save.as.ascii=0)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 30)
	stopifnot(!is.element(903, pred$mcmc.set$regions$country_code))
	npred <- dim(pred$e0.matrix.reconstructed)[2]
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}


test.existing.simulation <- function() {
	test.name <- 'retrieving MCMC results'
	start.test(test.name)
	sim.dir <- file.path(find.package("bayesLife"), "ex-data", 'bayesLife.output')
	m <- get.e0.mcmc(sim.dir, low.memory=FALSE)
	stopifnot(length(m$mcmc.list)==1)
	stopifnot(dim(m$mcmc.list[[1]]$traces)[1]==60) 
	m <- get.e0.mcmc(sim.dir)
	#summary(m)
	#summary(e0.mcmc(m, 1), par.names.cs=NULL)
	stopifnot(bayesTFR:::get.total.iterations(m$mcmc.list) == 60)
	stopifnot(bayesTFR:::get.stored.mcmc.length(m$mcmc.list, burnin=35) == 25)
	test.ok(test.name)
	
	test.name <- 'retrieving projection results'
	start.test(test.name)
	pred <- get.e0.prediction(sim.dir)
	s <- summary(pred, country='Japan')
	stopifnot(s$nr.traj == 30)
	stopifnot(all(dim(s$projections)==c(17,11)))
	# comment out if thinned mcmcs are not included in the package
	#mb <- get.thinned.e0.mcmc(m, thin=2, burnin=30)
	#s <- summary(mb, meta.only=TRUE)
	#stopifnot(s$iters == 30)
	test.ok(test.name)
}

test.DLcurve <- function() {
	test.name <- 'plotting DL curves'
	start.test(test.name)
	sim.dir <- file.path(find.package("bayesLife"), "ex-data", 'bayesLife.output')
	m <- get.e0.mcmc(sim.dir)
	filename <- tempfile()
	png(filename=filename)
	e0.DLcurve.plot(m, 'Slovenia')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	test.name <- 'obtaining DL curves'
	start.test(test.name)
	e0 <- seq(40, 90, length=100)
	dl <- e0.country.dlcurves(e0, m, "Slovenia", burnin=10)
	stopifnot(all(dim(dl)==c(50,100)))
	# world distribution
	dlw <- e0.world.dlcurves(e0, m, burnin=10)
	stopifnot(all(dim(dlw)==c(50,100)))
	# median of the world DL in the e0 range of 50-60 is larger than the country-specific median in that range
	# check visually with:
	# e0.DLcurve.plot(m, 'Slovenia'); lines(e0, apply(dlw, 2, median), col="blue")
	stopifnot(all(apply(dlw[, e0 > 50 & e0 < 60], 2, median) > apply(dl[, e0 > 50 & e0 < 60], 2, median)))
	test.ok(test.name)
}

test.e0trajectories <- function() {
	test.name <- 'plotting e0 trajectories'
	start.test(test.name)
	sim.dir <- file.path(find.package("bayesLife"), "ex-data", 'bayesLife.output')
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

test.plot.all <- function() {
	test.name <- 'plotting e0 trajectories and DL curves for all countries'
	start.test(test.name)
	sim.dir <- file.path(find.package("bayesLife"), "ex-data", 'bayesLife.output')
	pred <- get.e0.prediction(sim.dir=sim.dir)
	mc <- get.e0.mcmc(sim.dir)
	dir <- tempfile()
	dir.create(dir)
	e0.trajectories.plot.all(pred, output.dir=dir, main='XXX trajs')
	trajf <- length(list.files(dir, pattern='png$', full.names=FALSE))
	e0.DLcurve.plot.all(mc, output.dir=dir, main='DL XXX', output.type='jpeg')
	dlf <- length(list.files(dir, pattern='jpeg$', full.names=FALSE))
	unlink(dir, recursive=TRUE)
	stopifnot(trajf == get.nr.countries(mc$meta))
	stopifnot(dlf == get.nr.countries(mc$meta))
	test.ok(test.name)
}


test.plot.density <- function() {
	test.name <- 'plotting parameter density'
	start.test(test.name)
	sim.dir <- file.path(find.package("bayesLife"), "ex-data", 'bayesLife.output')
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
	test.name <- 'creating e0 maps via rworldmap'
	start.test(test.name)
	sim.dir <- file.path(find.package("bayesLife"), "ex-data", 'bayesLife.output')
	pred <- get.e0.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	e0.map(pred, year=2098, device='png', device.args=list(filename=filename))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	test.name <- 'creating e0 maps of observed data'
	start.test(test.name)
	sim.dir <- file.path(find.package("bayesLife"), "ex-data", 'bayesLife.output')
	pred <- get.e0.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	e0.map(pred, year=1974, device='png', device.args=list(filename=filename))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	test.name <- 'creating parameter maps'
	start.test(test.name)
	sim.dir <- file.path(find.package("bayesLife"), "ex-data", 'bayesLife.output')
	pred <- get.e0.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	e0.map(pred, par.name='z.c', device='png', device.args=list(filename=filename))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	test.name <- 'creating e0 maps using ggplot2'
	filename <- paste0(tempfile(), ".png")
	e0.ggmap(pred, file.name = filename)
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	e0.ggmap(pred, same.scale = TRUE, file.name = filename)
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	e0.ggmap(pred, same.scale = TRUE, year = 2100, file.name = filename)
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
}
	
test.get.parameter.traces <- function() {
	test.name <- 'getting parameter traces'
	start.test(test.name)
	sim.dir <- file.path(find.package("bayesLife"), "ex-data", 'bayesLife.output')
	m <- get.e0.mcmc(sim.dir, low.memory=TRUE)
	traces <- get.e0.parameter.traces(m$mcmc.list, burnin=20, 
					thinning.index=c(4, 21, 39))
	stopifnot(nrow(traces)==3)
	m.check <- get.e0.mcmc(sim.dir, low.memory=FALSE, burnin=20)
	stopifnot(traces[1,'omega']==m.check$mcmc.list[[1]]$traces[4,'omega'])
	# the following is from the time when there were two chains in the example data
	# indices 21 and 39 in the collapsed traces correspond to indices 1 and 19, respectively, in chain 2
	# stopifnot(all(traces[c(2,3),'omega']==m.check$mcmc.list[[2]]$traces[c(1,19),'omega']))
	
	traces <- get.e0.parameter.traces(m$mcmc.list, burnin=20, thin=8) 
	stopifnot(nrow(traces)==5)
	#stopifnot(traces[2,'z']==m.check$mcmc.list[[1]]$traces[5,'z']) #(4+1)
	#stopifnot(traces[9,'z']==m.check$mcmc.list[[2]]$traces[13,'z']) #(3*4 + 1)
	test.ok(test.name)
}

test.run.mcmc.simulation.auto <- function(compression='None') {
	sim.dir <- tempfile()
	# run MCMC
	test.name <- 'running auto MCMC'
	start.test(test.name)
	m <- run.e0.mcmc(iter = 'auto', output.dir = sim.dir, thin = 1, compression.type = compression,
					mcmc.options = list(auto.conf = list(iter=10, iter.incr=5, max.loops=3, nr.chains=2, thin=1, burnin=5)))
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 40)
	test.ok(test.name)

	test.name <- 'continuing auto MCMC'
	start.test(test.name)
	m <- continue.e0.mcmc(iter='auto', output.dir=sim.dir, 
	                      auto.conf=list(max.loops=2))
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 60)
	test.ok(test.name)

	unlink(sim.dir, recursive=TRUE)
}

test.estimate.mcmc.with.overwrites <- function() {
	sim.dir <- tempfile()
	test.name <- 'estimating MCMC with country overwrites'
	start.test(test.name)
	overwrites <- data.frame(country_code=c(562, 686), #Niger, Senegal (index 17, 18)
							k.c.prior.up=c(5, 7),
							Triangle_3.c.prior.low=c(NA, 0),
							Triangle_3.c.prior.up=c(0, NA)
						)
	# run MCMC
    m <- run.e0.mcmc(nr.chains = 1, iter = 50, thin = 1, output.dir = sim.dir, 
                     mcmc.options = list(Triangle.c = list(prior.low =c(0, 0, -20, 0)),
                                         country.overwrites = overwrites),
    				seed = 10)
    iNiger <- get.country.object(562, m$meta)
    iSene <- get.country.object(686, m$meta)
    
    stopifnot((m$meta$country.bounds$k.c.prior.up[iNiger$index] == 5) && (m$meta$country.bounds$k.c.prior.up[iSene$index] == 7) && 
    			all(m$meta$country.bounds$k.c.prior.up[-c(iNiger$index,iSene$index)]==10))
    stopifnot((m$meta$country.bounds$Triangle_3.c.prior.low[iSene$index] == 0) &&  
    			all(m$meta$country.bounds$Triangle_3.c.prior.low[-iSene$index] == -20))
    #check traces		
    traces.Niger <- get.e0.parameter.traces.cs(m$mcmc.list, get.country.object('Niger', m$meta), 
    					par.names=c('Triangle.c', 'k.c'))
	stopifnot(all(traces.Niger[,'k.c_c562'] <= 5) && all(traces.Niger[,'Triangle.c_3_c562'] <= 0))
	traces.Sen <- get.e0.parameter.traces.cs(m$mcmc.list, get.country.object('Senegal', m$meta), 
    					par.names=c('Triangle.c', 'k.c'))
	stopifnot(all(traces.Sen[,'k.c_c686'] < 7) && all(traces.Sen[,'Triangle.c_3_c686'] >= 0))
	test.ok(test.name)
	
	test.name <- 'estimating MCMC for extra countries with country overwrites'
	start.test(test.name)
	overwrites <- data.frame(country_code = c(800, 900), #Uganda, World
							k.c.prior.up = c(3, 8),
							k.c.prior.low = c(2, NA))
	m <- run.e0.mcmc.extra(sim.dir = sim.dir, countries = c(800, 900), burnin = 0, 
	                       country.overwrites = overwrites)
	Ug <- get.country.object('Uganda', m$meta)
	Wrld <- get.country.object(900, m$meta)
	stopifnot((m$meta$country.bounds$k.c.prior.up[Ug$index] == 3) && (m$meta$country.bounds$k.c.prior.up[iSene$index] == 7) && 
    			(m$meta$country.bounds$k.c.prior.up[Wrld$index] == 8) && all(m$meta$country.bounds$k.c.prior.up[-c(iNiger$index, iSene$index, Ug$index,Wrld$index)]==10))
    stopifnot((m$meta$country.bounds$k.c.prior.low[Ug$index] == 2) && all(m$meta$country.bounds$k.c.prior.low[-Ug$index]==0))
    traces.Ug <- get.e0.parameter.traces.cs(m$mcmc.list, Ug, par.names=c('k.c'))
	stopifnot(all(traces.Ug[,'k.c_c800'] <= 3) && all(traces.Ug[,'k.c_c800'] > 2))
	traces.world <- get.e0.parameter.traces.cs(m$mcmc.list, Wrld, par.names=c('k.c'))
	stopifnot(all(traces.world[,'k.c_c900'] <= 8))
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}

test.run.mcmc.simulation.auto.parallel <- function() {
	sim.dir <- tempfile()
	# run MCMC
	test.name <- 'running auto MCMC in parallel'
	start.test(test.name)
	m <- run.e0.mcmc(iter='auto', output.dir=sim.dir, parallel=TRUE, cltype='SOCK', thin=1,
					mcmc.options = list(auto.conf=list(iter=10, iter.incr=5, max.loops=3, nr.chains=2, thin=1, burnin=5)))
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 40)
	test.ok(test.name)

	test.name <- 'continuing auto MCMC in parallel'
	start.test(test.name)
	m <- continue.e0.mcmc(iter='auto', output.dir=sim.dir, 
	                      auto.conf=list(max.loops=2), 
	                      parallel=TRUE, cltype='SOCK')
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 60)
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}

test.imputation <- function() {
	sim.dir <- tempfile()

	# run MCMC
	test.name <- 'running MCMC with missing values'
	start.test(test.name)
	my.e0.data <- data.frame(country_code=4, last.observed=1990)
	my.e0.file <- tempfile()
	write.table(my.e0.data, my.e0.file, sep='\t', row.names=FALSE)
	m <- run.e0.mcmc(iter=5, nr.chains=1, output.dir=sim.dir, my.e0.file=my.e0.file, thin=1)
	unlink(my.e0.file)
	cindex <- get.country.object(4, m$meta)$index
	stopifnot(m$mcmc.list[[1]]$finished.iter == 5)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 5)
	stopifnot(length(m$meta$Tc.index[[cindex]]) == 8)
	test.name <- 'running projections with imputation'
	start.test(test.name)
	pred <- e0.predict(m, burnin=0, verbose=FALSE)
	stopifnot(length(pred$joint.male$meta.changes$Tc.index[[cindex]])==14)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 5)	
	test.ok(test.name)
	
	test.name <- 'plotting imputed e0 trajectories'
	start.test(test.name)
	filename <- tempfile()
	png(filename=filename)
	e0.trajectories.plot(pred, 4, pi=c(90, 54))
	e0.trajectories.plot(pred, 4, pi=c(90, 54), both.sexes=TRUE)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	test.name <- 'running projections with male imputation'
	start.test(test.name)
	my.e0.data <- data.frame(country_code=4, last.observed=1979)
	my.e0.file <- tempfile()
	write.table(my.e0.data, my.e0.file, sep='\t', row.names=FALSE)
	pred <- e0.predict(m, burnin=0, my.e0.file=my.e0.file, verbose=FALSE, replace.output=TRUE)
	unlink(my.e0.file)
	stopifnot(length(pred$joint.male$meta.changes$Tc.index[[cindex]])==6)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 5)
	test.ok(test.name)

	test.name <- 'plotting imputed F&M e0 trajectories'
	start.test(test.name)
	filename <- tempfile()
	png(filename=filename)
	e0.trajectories.plot(pred, 4, pi=c(90, 54), both.sexes=TRUE)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	test.name <- 'running MCMC and prediction with missing values for extra country'
	start.test(test.name)
	my.e0.data <- data.frame(country_code=900, last.observed=1996)
	my.e0.file <- tempfile()
	write.table(my.e0.data, my.e0.file, sep='\t', row.names=FALSE)
	m <- run.e0.mcmc.extra(sim.dir=sim.dir, countries=900, my.e0.file=my.e0.file, burnin=0)
	my.e0M.data <- data.frame(country_code=900, last.observed=1989)
	write.table(my.e0M.data, my.e0.file, sep='\t', row.names=FALSE)
	pred <- e0.predict.extra(sim.dir, my.e0.file=my.e0.file, verbose=FALSE)
	unlink(my.e0.file)
	stopifnot(length(pred$joint.male$meta.changes$Tc.index[[cindex]])==6)
	stopifnot(length(pred$joint.male$meta.changes$Tc.index[[get.country.object(900, m$meta)$index]])==8)
	stopifnot(length(m$meta$Tc.index[[get.country.object(900, m$meta)$index]]) == 9)
	test.ok(test.name)
	
	test.name <- 'plotting imputed F&M e0 trajectories for extra countries'
	start.test(test.name)
	filename <- tempfile()
	png(filename=filename)
	e0.trajectories.plot(pred, 4, pi=c(90, 54), both.sexes=TRUE)
	e0.trajectories.plot(pred, 900, pi=c(90, 54), both.sexes=TRUE)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	unlink(sim.dir, recursive=TRUE)
}

test.my.locations.extra <- function() {
	sim.dir <- tempfile()
    # run MCMC
    test.name <- 'using my.locations.file in extra estimation'
	start.test(test.name)
    m <- run.e0.mcmc(nr.chains=1, iter=10, thin=1, output.dir=sim.dir)
    data.dim <- dim(m$meta$e0.matrix)
    # prepare data and locations
    my.e0.file <- file.path(find.package("bayesLife"), 'extdata', 'my_e0_template.txt')
    e0 <- bayesTFR:::read.tfr.file(file=my.e0.file)
    e0['country_code'] <- 9999
    fdata <- tempfile()
    write.table(e0, file=fdata, sep='\t', row.names=FALSE)
    new.locations <- data.frame(country_code=9999, name='my location', reg_code=9999, area_code=8888, reg_name='my region', area_name='my area', location_type=4)
    flocs <- tempfile()
    write.table(new.locations, file=flocs, sep='\t', row.names=FALSE)
    m <- run.e0.mcmc.extra(sim.dir=sim.dir, my.e0.file=fdata, my.locations.file=flocs, burnin=0)
    stopifnot(all(dim(m$meta$e0.matrix) == c(data.dim[1], data.dim[2]+1)))
	stopifnot(m$meta$regions$country_code[m$meta$nr.countries] == 9999)
	unlink(flocs)
	unlink(fdata)
	unlink(sim.dir, recursive=TRUE)
	test.ok(test.name)
}

test.reproduce.simulation <- function() {
    sim.dir <- tempfile()
    test.name <- 'reproducing simulation'
    start.test(test.name)
    seed <- 1234
    
    m1 <- run.e0.mcmc(iter=5, nr.chains=2, thin = 1, output.dir=sim.dir, start.year=1950, seed = seed, verbose = FALSE)
    res1 <- summary(m1)$results$statistics
    m2 <- run.e0.mcmc(iter=5, nr.chains=2, thin = 1, output.dir=sim.dir, start.year=1950, seed = seed, replace.output = TRUE, verbose = FALSE)
    res2 <- summary(m2)$results$statistics
    stopifnot(all(res1 == res2))
    
    m3 <- run.e0.mcmc(iter=5, nr.chains=2, thin = 1, output.dir=sim.dir, start.year=1950, replace.output = TRUE) # no seed
    res3 <- summary(m3)$results$statistics
    stopifnot(!all(res3 == res2))
    
    # parallel
    m1p <- run.e0.mcmc(iter=5, nr.chains=2, thin = 1, output.dir=sim.dir, start.year=1950, seed = seed, replace.output = TRUE, parallel = TRUE, ft_verbose = FALSE)
    res1p <- summary(m1p)$results$statistics
    m2p <- run.e0.mcmc(iter=5, nr.chains=2, thin = 1, output.dir=sim.dir, start.year=1950, seed = seed, replace.output = TRUE, parallel = TRUE, ft_verbose = FALSE)
    res2p <- summary(m2p)$results$statistics
    stopifnot(all(res1p == res2p))
    
    m3p <- run.e0.mcmc(iter=5, nr.chains=2, thin = 1, output.dir=sim.dir, start.year=1950, replace.output = TRUE, parallel = TRUE, ft_verbose = TRUE) # no seed
    res3p <- summary(m3p)$results$statistics
    stopifnot(!all(res3p == res2p))

    test.ok(test.name)
    unlink(sim.dir, recursive=TRUE)
}

test.subnational.predictions <- function() {
    sim.dir <- tempfile()
    test.name <- 'predicting subnational e0'
    start.test(test.name)
    
    my.sub.file <- file.path(find.package("bayesLife"), 'extdata', 'subnational_e0_template.txt')
    nat.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
    
    # Subnational female projections for Australia and Canada; start 1 time period before national projections
    preds <- e0.predict.subnat(c(36, 124), my.e0.file=my.sub.file,
                                sim.dir=nat.dir, output.dir=sim.dir, start.year=2018)
    stopifnot(all(names(preds) %in% c("36", "124")))
    stopifnot(nrow(get.countries.table(preds[["36"]]))==8)
    stopifnot(identical(preds[["124"]]$quantiles, get.rege0.prediction(sim.dir, 124)$quantiles))
    filename <- tempfile()
    pdf(file=filename)
    e0.trajectories.plot(preds[["36"]], "Queensland")
    e0.trajectories.plot(preds[["124"]], "Quebec")
    dev.off()
    size <- file.info(filename)['size']
    unlink(filename)
    stopifnot(size > 0)
    
    # generate predictions with different method
    preds.const <- e0.predict.subnat(c(36, 360), my.e0.file = my.sub.file, method = "shift",
                               sim.dir=nat.dir, output.dir=sim.dir, predict.jmale = TRUE, 
                               my.e0M.file = my.sub.file) # using the same e0 for male just for testing purposes
    stopifnot(identical(preds.const[["360"]]$quantiles, get.rege0.prediction(sim.dir, 360, method = "shift")$quantiles))
    stopifnot(!identical(preds.const[["360"]]$quantiles, preds[["360"]]$quantiles))
    
    stopifnot(!has.e0.jmale.prediction(preds[["36"]])) # default AR1 prediction does not have male predictions
    
    prF <- get.rege0.prediction(sim.dir, 360, method = "shift")
    prM1 <- get.e0.jmale.prediction(prF)
    prM2 <- get.rege0.prediction(sim.dir, 360, method = "shift", joint.male = TRUE)
    stopifnot(identical(prM1$quantiles, prM2$quantiles))
    stopifnot(!identical(prM1$quantiles, prF$quantiles))
    
    unlink(sim.dir, recursive=TRUE)
    test.ok(test.name)
    
    test.name <- 'predicting subnational e0 with imputation'
    start.test(test.name)
    sim.dir <- tempfile()
    preds <- e0.predict.subnat("Canada", my.e0.file = my.sub.file,
                                sim.dir = nat.dir, output.dir = sim.dir, start.year = 2021)
    stopifnot(all(names(preds) %in% c("124")))
    spred <- summary(preds[["124"]], "Quebec")
    stopifnot(min(spred$projection.years) == 2023)
    stopifnot(!has.e0.jmale.prediction(get.rege0.prediction(sim.dir, 124)))
    
    # generate male prediction
    predCan <- e0.jmale.predict.subnat(preds[["124"]], my.e0.file = my.sub.file)
    stopifnot(has.e0.jmale.prediction(get.rege0.prediction(sim.dir, 124)))
    
    filename <- tempfile()
    pdf(file=filename)
    e0.trajectories.plot(predCan, "Ontario", both.sexes = TRUE)
    dev.off()
    size <- file.info(filename)['size']
    unlink(filename)
    stopifnot(size > 0)
    
    unlink(sim.dir, recursive=TRUE)
    
    # Subnational projections for Canada with one region (Alberta) having imputed values
    sim.dir <- tempfile()
    mye0 <- read.delim(my.sub.file, check.names = FALSE)
    mye0[mye0$country_code == 124 & mye0$reg_code == 658, "last.observed"] <- 1999
    my.sub.file.mis <- tempfile()
    write.table(mye0, file = my.sub.file.mis, row.names = FALSE, quote = FALSE, sep = "\t")
    preds <- e0.predict.subnat("Canada", my.e0.file=my.sub.file.mis,
                                sim.dir=nat.dir, output.dir=sim.dir, start.year=2018, 
                               predict.jmale = TRUE, my.e0M.file = my.sub.file)
    obs <- preds[["124"]]$mcmc.set$meta$e0.matrix[,"658"]
    recon <- preds[["124"]]$e0.matrix.reconstructed[,"658"]
    stopifnot(all(is.na(obs[c("2003", "2008")]))) # observed is NA
    stopifnot(all(!is.na(obs["1998"]))) # last observed non-NA 
    stopifnot(all(!is.na(recon[c("2003", "2008")]))) # missing values were reconstructed
    unlink(sim.dir, recursive=TRUE)
    unlink(my.sub.file.mis)
    
    # Subnational projections for all countries in the file
    sim.dir <- tempfile()
    preds <- e0.predict.subnat(c(36, 360, 124), my.e0.file=my.sub.file, sim.dir=nat.dir, 
                                output.dir=sim.dir, end.year=2050)
    stopifnot(length(preds) == 3)
    
    # Retrieve results for all countries
    preds <- get.rege0.prediction(sim.dir)
    stopifnot(length(preds) == 3)
    t <- e0.trajectories.table(preds[["360"]], "Maluku")
    stopifnot(all(dim(t) == c(21, 5)))
    
    # Retrieve results for one country
    pred <- get.rege0.prediction(sim.dir, 124)
    spred <- summary(pred, "British Columbia")
    stopifnot(max(spred$projection.years) == 2048)
    
    # Retrieve trajectories
    trajs <- get.e0.trajectories(preds[["36"]], "Victoria")
    stopifnot(all(dim(trajs) == c(7, 30)))
    
    unlink(sim.dir, recursive=TRUE)
    
    # Annual subnational projections
    my.sub.file.annual.F <- tempfile()
    my.sub.file.annual.M <- tempfile()
    # female data
    datF <- data.table::fread(my.sub.file) # load the data and save only the last time period, pretending this is a single year data point
    datF <- datF[, .(country_code, country_name, reg_code, name, `2013` = `2010-2015`)]
    data.table::fwrite(datF, file = my.sub.file.annual.F, sep = "\t") # save it
    # male data
    datM <- data.table::copy(datF)[, `2013` := `2013`*0.95]
    data.table::fwrite(datM, file = my.sub.file.annual.M, sep = "\t")
    
    sim.dir <- tempfile()
    preds <- e0.predict.subnat(c(36, 124, 360), my.e0.file=my.sub.file.annual.F, sim.dir=nat.dir, 
                                output.dir=sim.dir, start.year = 2016, # will impute 2 years
                               end.year=2030, annual = TRUE, 
                               predict.jmale = TRUE, my.e0M.file = my.sub.file.annual.M)
    # Retrieve trajectories
    trajs <- get.e0.trajectories(preds[["124"]], "Alberta")
    stopifnot(all(dim(trajs) == c(16, 30))) # from 2015 to 2030
    stopifnot(all(c(2015, 2030) %in% as.integer(rownames(trajs))))
    
    pred <- get.rege0.prediction(sim.dir, 360)
    spred <- summary(pred, "Papua")
    stopifnot(max(spred$projection.years) == 2030)
    
    predM <- get.e0.jmale.prediction(pred)
    t <- tfr.trajectories.table(predM, "Bali")
    stopifnot(max(as.integer(rownames(t))) == 2030)
    stopifnot(all(c(2015, 2024, 2029) %in% as.integer(rownames(t))))
    
    unlink(sim.dir, recursive=TRUE)
    unlink(my.sub.file.annual.F)
    unlink(my.sub.file.annual.M)
    
    test.ok(test.name)
    
}

test.run.annual.simulation <- function(wpp.year = 2019) {
    sim.dir <- tempfile()
    
    test.name <- 'running MCMC with annual data'
    start.test(test.name, wpp.year)
    m <- run.e0.mcmc(iter = 5, nr.chains = 2, thin = 1, output.dir = sim.dir, 
                     annual = TRUE, present.year = 2018, wpp.year = wpp.year)
    stopifnot(get.total.iterations(m$mcmc.list, 0) == 10)
    stopifnot(all(1953:2018 %in% rownames(m$meta$e0.matrix)))
    test.ok(test.name)
    
    test.name <- 'running annual MCMC for extra country'
    start.test(test.name, wpp.year)
    countries <- c(900, 908)
    m <- run.e0.mcmc.extra(sim.dir = sim.dir, countries = countries, burnin = 0)
    stopifnot(countries %in% m$meta$regions$country_code)
    test.ok(test.name)
    
    test.name <- 'running annual projections'
    start.test(test.name, wpp.year)
    pred <- e0.predict(m, burnin=1, verbose = FALSE)
    spred <- summary(pred)
    tbl <- e0.trajectories.table(pred, "Japan")
    stopifnot(spred$nr.traj == 8)
    stopifnot(all(2019:2100 %in% spred$projection.years))
    stopifnot(all(c(908, 900) %in% get.countries.table(pred)$code))
    if(wpp.year > 2019)
        stopifnot('1900' %in% rownames(tbl)) # checks that supplemental data is used
    
    mpred <- get.e0.jmale.prediction(pred)
    smpred <- summary(mpred)
    mtbl <- e0.trajectories.table(mpred, "Japan")
    stopifnot(all(2019:2100 %in% smpred$projection.years))
    stopifnot(all(c(908, 900) %in% get.countries.table(mpred)$code))
    if(wpp.year > 2019)
        stopifnot('1900' %in% rownames(mtbl))
    
    test.ok(test.name)
    
    unlink(sim.dir, recursive=TRUE)
}

