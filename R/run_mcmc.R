if(getRversion() >= "2.15.1") utils::globalVariables(c("loess_sd"))
data(loess_sd, envir = environment())

update.ini.values <- function(nr.chains, annual = FALSE) {
    # set starting values
    #=====================
    get.init.values.between.low.and.up <- function(low, up)
        ifelse(rep(nr.chains==1, nr.chains), (low + up)/2, #seq(low, to=up, length=nr.chains)
               runif(nr.chains, low, up)
        )
    current <- e0mcmc.options(annual = annual)

    wp <- current$world.parameters
    for(par in names(wp)) {
        if(wp[[par]] == 4) {
            for (i in 1:4) {
                lname <- paste0("T", i)
                if(is.null(current[[par]][["ini"]][[lname]])) 
                    current[[par]][["ini"]][[lname]] <- get.init.values.between.low.and.up(
                                                    current[[par]][["ini.low"]][i], 
                                                    current[[par]][["ini.up"]][i])
            }
        } else {
            if(wp[[par]] == 1) {
                if(is.null(current[[par]][["ini"]]))
                    current[[par]][["ini"]] <- get.init.values.between.low.and.up(
                        current[[par]][["ini.low"]], current[[par]][["ini.up"]])
            }
        }
    }
    invisible(e0mcmc.options(current, annual = annual))
}

match.ini.to.chains <- function(nr.chains, annual = FALSE) {
    # propagate initial values for all chains if needed
    current <- e0mcmc.options(annual = annual)
    starting.values <- list()
    wp <- current$world.parameters
    for(par in names(wp)) {
        if(wp[[par]] == 4) {
            val <- as.list(current[[par]][["ini"]])
            starting.values[[par]] <- vector("list", length = 4)
            names(starting.values[[par]]) <- paste0("T", 1:4)
            for(i in 1:4) {
                lname <- paste0("T", i)
                starting.values[[par]][[lname]] <- .match.length.to.nr.chains(
                                                val[[i]], nr.chains, paste0(par, "[[", lname, "]]"))
            }
        } else {
            if(wp[[par]] == 1) {
                starting.values[[par]] <- .match.length.to.nr.chains(
                                                current[[par]][["ini"]], nr.chains, par)
            }
        }
    }
    return(starting.values)
}

.match.length.to.nr.chains <- function(val, nr.chains, var.name) {
    if (length(val) >= nr.chains) 
        return(val[1:nr.chains])
    if (length(val) == 1) 
        return(rep(val, nr.chains))
    warning(var.name, ' has the wrong length. Either 1 or ', nr.chains, 
                    ' is allowed.\nValue set to ', val[1], ' for all chains.')
    return(rep(val[1], nr.chains))
}

run.e0.mcmc <- function(sex=c("Female", "Male"), nr.chains = 3, iter = 160000, 
						output.dir = file.path(getwd(), 'bayesLife.output'), 
                        thin = 10, replace.output = FALSE, annual = FALSE,
                        start.year = 1873, present.year = 2020, wpp.year = 2019,
                        my.e0.file = NULL, my.locations.file = NULL, 
						constant.variance = FALSE, seed = NULL, parallel = FALSE, 
						nr.nodes = nr.chains, compression.type = 'None',
						verbose = FALSE, verbose.iter = 100, mcmc.options = NULL, ...) {
						 	
    dotargs <- list(...)
    if(any(names(dotargs) %in% (alloldargs <- .legacy.run.mcmc.args()))) {
        oldargs <- intersect(names(dotargs), alloldargs)
        oldagrs.nc <- intersect(names(dotargs), .legacy.run.mcmc.args.no.name.change())
        eg <- "E.g. mcmc.options = list("
        if(length(oldagrs.nc) > 0) 
            eg <- paste0(eg, oldagrs.nc[1], " = ", 
                         if(!is.list(dotargs[[oldagrs.nc[1]]]) && length(dotargs[[oldagrs.nc[1]]]) == 1) dotargs[[oldagrs.nc[1]]] else "...", ", ")
        eg <- paste0(eg, "z = list(ini.up = 0.7))\nSee ?e0mcmc.options")
        stop("Using arguments ", paste(oldargs, collapse = ","), " in bayesLife > 4.0 is obsolete. Use mcmc.options instead.\n", eg)
    }
    
	if(file.exists(output.dir)) {
		if(length(list.files(output.dir)) > 0 & !replace.output)
                        stop('Non-empty directory ', output.dir, 
                        ' already exists.\nSet replace.output=TRUE if you want to overwrite existing results.')
        unlink(output.dir, recursive=TRUE)
	}
    if(!is.null(seed)) set.seed(seed)
    dir.create(output.dir)
    old.opts <- e0mcmc.options(annual = annual)
    if(!is.null(mcmc.options))
        e0mcmc.options(mcmc.options, annual = annual)
    
    mcoptions <- update.ini.values(nr.chains, annual)
	auto.run <- FALSE
	auto.conf <- mcoptions$auto.conf
	if(iter == 'auto') { # defaults for auto-run (includes convergence diagnostics)
		iter <- auto.conf$iter
		nr.chains <- auto.conf$nr.chains
		auto.run <- TRUE		
	}      
	if (verbose) {
		cat('\nStarting Bayesian Hierarchical Model for Life Expectancy.\n')
		cat('=========================================================\n')
		cat('Using configuration for ', e0.options("use"))
		cat('Initialize simulation -', nr.chains, 'chain(s) in total.\n')
	}
	sex <- substr(match.arg(sex), 1, 1)
	bayesLife.mcmc.meta <- e0.mcmc.meta.ini(sex=sex, nr.chains = nr.chains,
                                   		start.year = start.year, present.year = present.year, 
                                        wpp.year = wpp.year, annual.simulation = annual, 
                                   		my.e0.file = my.e0.file, my.locations.file = my.locations.file,
                                        output.dir = output.dir, mcmc.options = mcoptions, 
                                        constant.variance = constant.variance, 
                                        compression.type = compression.type, verbose = verbose)
    store.bayesLife.meta.object(bayesLife.mcmc.meta, output.dir)
    starting.values <- match.ini.to.chains(nr.chains, annual = annual)
    iter <- .match.length.to.nr.chains(iter, nr.chains, "iter")

	if (parallel) { # run chains in parallel
		chain.set <- bayesTFR:::bDem.performParallel(nr.nodes, 1:nr.chains, mcmc.run.chain.e0, 
                                     initfun = mcoptions$parallel.init.function, seed = seed,
                                     meta = bayesLife.mcmc.meta, 
                                     thin = thin, iter = iter, 
                                     starting.values = starting.values,                                     
                                     verbose = verbose, verbose.iter = verbose.iter, ...)
	} else { # run chains sequentially
		chain.set <- list()
		for (chain in 1:nr.chains) {
			chain.set[[chain]] <- mcmc.run.chain.e0(chain, bayesLife.mcmc.meta, thin = thin, 
                                                iter = iter, starting.values = starting.values, 
                                                verbose = verbose, verbose.iter = verbose.iter)
		}
	}
	names(chain.set) <- 1:nr.chains
	mcmc.set <- structure(list(meta = bayesLife.mcmc.meta, mcmc.list = chain.set), 
	                      class='bayesLife.mcmc.set')
    cat('\nResults stored in', output.dir,'\n')
    
    if(auto.run) {
		diag <- try(e0.diagnose(sim.dir = output.dir, keep.thin.mcmc = TRUE, 
						thin = auto.conf$thin, burnin = auto.conf$burnin,
						verbose = verbose))
		if(auto.conf$max.loops > 1) {
			for(loop in 2:auto.conf$max.loops) {
				if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
				mcmc.set <- continue.e0.mcmc(iter = auto.conf$iter.incr, output.dir = output.dir, 
				                             nr.nodes = nr.nodes, parallel = parallel, 
				                             verbose = verbose, verbose.iter = verbose.iter, ...)
				diag <- try(e0.diagnose(sim.dir = output.dir, keep.thin.mcmc = TRUE, 
							thin = auto.conf$thin, burnin = auto.conf$burnin,
							verbose = verbose))
			}
		}
    }
    e0mcmc.options(old.opts, annual = annual)
    if (verbose) 
		cat('\nSimulation successfully finished!!!\n')
    invisible(mcmc.set)
}

.legacy.run.mcmc.args.no.name.change <- function()
    c("a", "delta", "tau", "outliers", "country.overwrites", "nu", 
      "dl.p1", "dl.p2", "sumTriangle.lim", "buffer.size", "auto.conf")

.legacy.run.mcmc.args <- function() {
    c(.legacy.run.mcmc.args.no.name.change(), 
      paste(c("Triangle", "k", "z", "lambda", "lambda.k", "lambda.z", "omega"), 
            rep(c("ini", "ini.low", "ini.up"), each = 7), sep = "."),
      paste(c("Triangle", "k", "z", "Triangle.c", "k.c", "z.c"),
            rep(c("prior.low", "prior.up"), each = 6)),
      paste(c("Triangle.c", "k.c", "z.c"),
            rep("ini.norm", 3))
    )
}

mcmc.run.chain.e0 <- function(chain.id, meta, thin = 1, iter = 100, starting.values = NULL, 
                              verbose = FALSE, verbose.iter = 10) {
	cat('\n\nChain nr.', chain.id, '\n')
    if (verbose) 
    	cat('************\n')
	this.sv <- list()
	for(var in names(starting.values)) {
		this.sv[[var]] <- if (is.list(starting.values[[var]])) sapply(starting.values[[var]], function(x) x[chain.id])
							else starting.values[[var]][chain.id]
	}
	if(iter[chain.id] < thin) {
	    warning("Argument iter is smaller than thin. Adjusted to ", thin, ".")
	    iter[chain.id] <- thin
	}
	mcmc <- e0.mcmc.ini(chain.id, meta, iter = iter[chain.id], ini.values = this.sv)

    if (verbose) {
        cat('Starting values:\n')
        print(unlist(mcmc[names(meta$mcmc.options$world.parameters)]))
        cat('Store initial values into ', mcmc$output.dir, '\n')
    }
	store.e0.mcmc(mcmc, append = FALSE, flush.buffer = TRUE, verbose = verbose)
	if (verbose) 
    	cat('Start sampling -', mcmc$iter, 'iterations in total.\n')
	mcmc <- do.call(meta$mcmc.options$estimation.function, 
	                list(mcmc, thin = thin, verbose = verbose, verbose.iter = verbose.iter))
    return(mcmc)
}
        
continue.e0.mcmc <- function(iter, chain.ids = NULL, 
                             output.dir = file.path(getwd(), 'bayesLife.output'), 
                             parallel = FALSE, nr.nodes = NULL, auto.conf = NULL, 
                             verbose = FALSE, verbose.iter=10, ...) {
        mcmc.set <- get.e0.mcmc(output.dir)
        opts <- mcmc.set$meta$mcmc.options
        auto.run <- FALSE
        if(iter == 'auto') { 
            # defaults for auto-run (includes convergence diagnostics)
            default.auto.conf <- opts$auto.conf
            if(is.null(auto.conf)) auto.conf <- list()
            # merge with new auto-conf
            for (par in names(default.auto.conf))
                if(is.null(auto.conf[[par]])) auto.conf[[par]] <- default.auto.conf[[par]]
			iter <- auto.conf$iter.incr
			auto.run <- TRUE
			fiter <- sapply(mcmc.set$mcmc.list, function(x) x$finished.iter)
			if (!all(fiter== fiter[1])) stop('All chains must be of the same length if the "auto" option is used.')
        }
        if (is.null(chain.ids) || auto.run) {
                chain.ids <- names(mcmc.set$mcmc.list)
        }
        if (parallel) { # run chains in parallel
                if(is.null(nr.nodes)) nr.nodes<-length(chain.ids)
                chain.list <- bayesTFR:::bDem.performParallel(nr.nodes, chain.ids, continue.e0.chain, 
                                                initfun=opts$parallel.init.function, 
                                                mcmc.list=mcmc.set$mcmc.list, iter=iter, 
                                                verbose=verbose, verbose.iter=verbose.iter, ...)
                for (i in 1:length(chain.ids))
                        mcmc.set$mcmc.list[[chain.ids[i]]] <- chain.list[[i]]
        } else { # run chains sequentially
                for (chain.id in chain.ids) {
                        mcmc.set$mcmc.list[[chain.id]] <- continue.e0.chain(chain.id, mcmc.set$mcmc.list, 
                                                        iter=iter, verbose=verbose, verbose.iter=verbose.iter)
                }
        }
        cat('\n')
        if(auto.run) {
        	diag <- try(e0.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
                                    thin=auto.conf$thin, burnin=auto.conf$burnin,
                                    verbose=verbose))
			if(auto.conf$max.loops>1) {
				for(loop in 2:auto.conf$max.loops) {
					if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
					mcmc.set <- continue.e0.mcmc(iter=auto.conf$iter.incr, output.dir=output.dir, nr.nodes=nr.nodes,
												 parallel=parallel, verbose=verbose, verbose.iter=verbose.iter, ...)
					diag <- try(e0.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
											thin=auto.conf$thin, burnin=auto.conf$burnin, verbose=verbose))
				}
			}
        }
        invisible(mcmc.set)
}

continue.e0.chain <- function(chain.id, mcmc.list, iter, verbose=FALSE, verbose.iter=10) {
        cat('\n\nChain nr.', chain.id, '\n')
        if (verbose)
                cat('************\n')
        mcmc <- mcmc.list[[chain.id]]
        mcmc$iter <- mcmc$finished.iter + iter
        if (verbose) 
                cat('Continue sampling -', iter, 'additional iterations,', mcmc$iter, 'iterations in total.\n')
        mcmc <- do.call(mcmc$meta$mcmc.options$estimation.function, 
                        list(mcmc, thin = mcmc$thin, start.iter=mcmc$finished.iter+1, 
                             verbose = verbose, verbose.iter = verbose.iter))
        return(mcmc)
}

run.e0.mcmc.extra <- function(sim.dir=file.path(getwd(), 'bayesLife.output'), 
								countries = NULL, my.e0.file = NULL, iter = NULL,
								thin = 1, burnin = 0, parallel = FALSE, nr.nodes = NULL, 
								my.locations.file = NULL, country.overwrites = NULL, 
								verbose = FALSE, verbose.iter = 100, ...) {
									
    mcmc.set <- get.e0.mcmc(sim.dir)
	Eini <- e0.mcmc.meta.ini.extra(mcmc.set, countries = countries, 
	                                my.e0.file = my.e0.file, 
							        my.locations.file = my.locations.file, 
							        burnin = burnin, country.overwrites = country.overwrites,
							        verbose = verbose)
	meta <- Eini$meta
	if(length(Eini$index) <= 0) {
		cat('\nNothing to be done.\n')
		return(invisible(mcmc.set))
	}
	chain.ids <- names(mcmc.set$mcmc.list)
	mcthin <- 1
	if(verbose) {
	    cat('\nMCMC for extra countries, using settings for ', e0.options("use"), "\n\n")
	}
	for (chain in chain.ids) { # update meta in each chain
		if(verbose) cat('Updating meta in chain', chain, '\n')
		mcmc.set$mcmc.list[[chain]]$meta <- meta
		mcmc.set$mcmc.list[[chain]] <- e0.mcmc.ini.extra(mcmc.set$mcmc.list[[chain]], countries=Eini$index,
												index.replace=Eini$index.replace)
		mcthin <- max(mcthin, mcmc.set$mcmc.list[[chain]]$thin)
	}
	mcthin <- mcmc.set$mcmc.list[[1]]$thin
	total.iter <- mcmc.set$mcmc.list[[1]]$length - bayesTFR:::get.thinned.burnin(mcmc.set$mcmc.list[[1]], burnin)
	thin <- max(thin, mcthin)
	post.idx <- if (thin > mcthin) unique(round(seq(thin, total.iter, by=thin/mcthin)))
				else 1:total.iter
	if (!is.null(mcmc.set$mcmc.list[[1]]$rng.state)) .Random.seed <- mcmc.set$mcmc.list[[1]]$rng.state
	
	if (parallel) { # run chains in parallel
		if(is.null(nr.nodes)) nr.nodes<-length(chain.ids)
		chain.list <- bayesTFR:::bDem.performParallel(nr.nodes, chain.ids, e0.mcmc.run.chain.extra, 
						initfun = mcmc.set$meta$mcmc.options$parallel.init.function, 
						mcmc.list=mcmc.set$mcmc.list, countries=Eini$index, 
						posterior.sample=post.idx, iter=iter, burnin=burnin, verbose=verbose, verbose.iter=verbose.iter, ...)
		for (i in 1:length(chain.ids))
			mcmc.set$mcmc.list[[chain.ids[i]]] <- chain.list[[i]]
	} else { # run chains sequentially
		for (chain.id in chain.ids) {
			mcmc.set$mcmc.list[[chain.id]] <- e0.mcmc.run.chain.extra(chain.id, mcmc.set$mcmc.list, 
												countries=Eini$index, posterior.sample=post.idx, iter=iter,  
												burnin=burnin, verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	store.bayesLife.meta.object(meta, meta$output.dir)
	mcmc.set$meta <- meta
	cat('\n')
	invisible(mcmc.set)
}
	
e0.mcmc.run.chain.extra <- function(chain.id, mcmc.list, countries, posterior.sample, 
												iter=NULL, burnin=2000, verbose=FALSE, verbose.iter=100) {
	cat('\n\nChain nr.', chain.id, '\n')
	if (verbose)
		cat('************\n')
	mcmc <- mcmc.list[[chain.id]]
		
	if (verbose) 
		cat('MCMC sampling for additional countries and regions.\n')

	mcmc <- e0.mcmc.sampling.extra(mcmc, mcmc.list=mcmc.list, countries=countries, 
									posterior.sample=posterior.sample, 
									iter=iter, burnin=burnin, 
									verbose=verbose, verbose.iter=verbose.iter)
	return(mcmc)
}


.get.Tcindex <- function(e0.matrix,  stop.if.less.than2=TRUE, cnames=NULL) {
	Tc.index <- list()
	for (country in 1:ncol(e0.matrix)) {
		Tc.index[[country]] <- which(!is.na(e0.matrix[,country]))
    	if(stop.if.less.than2 && length(Tc.index[[country]]) < 2) stop('Problem with ', cnames[country], 
    						". At least two data points must be observed.")
    }
	return(Tc.index)
}

.do.part.e0.mcmc.meta.ini <- function(data, meta) {
    opts <- meta$mcmc.options
	nr_countries <- ncol(data$e0.matrix)
    #T_end_c <- rep(NA, nr_countries)
    Tc.index <- .get.Tcindex(data$e0.matrix, cnames=data$regions$country_name)
	T <- nrow(data$e0.matrix)
	d.ct <- loessSD <- matrix(NA, nrow=T-1, ncol=nr_countries, 
							dimnames=list(rownames(data$e0.matrix)[1:(T-1)],
									  colnames(data$e0.matrix)))
	loessSD[,] <- 1
	for(i in 2:T) {
		nisna0 <- !is.na(data$e0.matrix[i-1,])
		nisna1 <- !is.na(data$e0.matrix[i,])
		nisna2 <- nisna1 & nisna0
		if (sum(nisna2) > 0) {
			d.ct[i-1,nisna2] <- data$e0.matrix[i,nisna2] - data$e0.matrix[i-1,nisna2]
			outliers <- nisna2 & ((d.ct[i-1,] < opts$outliers[1]) | (d.ct[i-1,] > opts$outliers[2]))
			d.ct[i-1,outliers] <- NA
		}
		if (sum(nisna0) > 0 && !meta$constant.variance)
			loessSD[i-1,nisna0]<- loess.lookup(data$e0.matrix[i-1,nisna0])
			#loessSD[i-1,nisna0]<-sapply(data$e0.matrix[i-1,nisna0],loess.lookup)
	}
	D.supp.ct <- loessSD.suppl <- NULL
	nr_countries.suppl <- 0
	suppl <- data$suppl.data
	if(!is.null(suppl$e0.matrix)) {
		nr_countries.suppl <- ncol(suppl$e0.matrix)
		suppl$Tc.index <- .get.Tcindex(suppl$e0.matrix, stop.if.less.than2=FALSE)
		# add first time point of the observed data to get the last increment of the supplemental data
		data.suppl <- rbind(suppl$e0.matrix, data$e0.matrix[1,suppl$index.to.all.countries])
		T <- nrow(data.suppl)
		d.suppl.ct <- loessSD.suppl <- matrix(NA, nrow=T-1, ncol=nr_countries.suppl)
		for(i in 2:T) {
			nisna0 <- !is.na(data.suppl[i-1,])
			nisna1 <- !is.na(data.suppl[i,])
			nisna2 <- nisna1 & nisna0
			if (sum(nisna2) > 0) {
				d.suppl.ct[i-1,nisna2] <- data.suppl[i,nisna2] - data.suppl[i-1,nisna2]
				outliers <- nisna2 & ((d.suppl.ct[i-1,] < opts$outliers[1]) | (d.suppl.ct[i-1,] > opts$outliers[2]))
				d.suppl.ct[i-1,outliers] <- NA
			}
			if (sum(nisna0) > 0)
				loessSD.suppl[i-1,nisna0]<- if(meta$constant.variance) 1 else loess.lookup(data.suppl[i-1,nisna0])
				#loessSD.suppl[i-1,nisna0]<- if(meta$constant.variance) 1 else sapply(data.suppl[i-1,nisna0],loess.lookup)
		}
		suppl$nr.countries <- nr_countries.suppl
		suppl$d.ct <- d.suppl.ct
		suppl$loessSD <- loessSD.suppl
	}
	data$nr.countries <- nr_countries
	data$Tc.index <- Tc.index
	data$d.ct <- d.ct
	data$loessSD <- loessSD
	data$suppl.data <- suppl
	bounds <- .do.country.specific.ini(nr_countries, c(data, meta))
	return(c(data, bounds))
}

.do.country.specific.ini <- function(nr_countries, meta) {
    opts <- meta$mcmc.options
	samplpars <- list()
    for (i in 1:4) {
    	samplpars[[paste0('Triangle_', i, '.c.prior.low')]] <- rep(opts$Triangle.c$prior.low[i], nr_countries)
    	samplpars[[paste0('Triangle_', i, '.c.prior.up')]] <- rep(opts$Triangle.c$prior.up[i], nr_countries)
    }
    samplpars$k.c.prior.low <- rep(opts$k.c$prior.low, nr_countries)
    samplpars$k.c.prior.up <- rep(opts$k.c$prior.up, nr_countries)
    samplpars$z.c.prior.low <- rep(opts$z.c$prior.low, nr_countries)
    samplpars$z.c.prior.up <- rep(opts$z.c$prior.up, nr_countries)
    for(parname in c(paste0('Triangle_', 1:4, '.c.prior.low'), 
    				paste0('Triangle_', 1:4, '.c.prior.up'), 'k.c.prior.low', 'k.c.prior.up', 'z.c.prior.low',
    						'z.c.prior.up'
    						))
    	names(samplpars[[parname]]) <- meta$regions$country_code

    if(!is.null(opts$country.overwrites)) {
    	for(row in 1:nrow(opts$country.overwrites)) {
    		country <- opts$country.overwrites[row, 'country_code']
    		for (col in colnames(opts$country.overwrites)) {
    			if(col == 'country_code') next
    			if(!is.element(col, names(samplpars))) {
    				warnings(col, ' is not a valid column name in country.overwrites. Column ignored.')
    				next
    			}
    			if(!is.element(as.character(country), names(samplpars[[col]]))) {
    				warnings(country, ' is not a valid country. Row ignored.')
    				break
    			}
 
    			if(!is.na(opts$country.overwrites[row,col]))
    				samplpars[[col]][as.character(country)] <- opts$country.overwrites[row,col]
    		}
    	}	
    }
	return(list(country.bounds = samplpars))
}


e0.mcmc.meta.ini <- function(sex = "F", nr.chains = 1, start.year = 1950, present.year = 2020, 
								wpp.year = 2019, my.e0.file = NULL, my.locations.file = NULL,
								annual.simulation = FALSE, output.dir = file.path(getwd(), 'bayesLife.output'),
								mcmc.options = NULL, ..., verbose=FALSE) {
	mcmc.input <- c(list(sex = sex, nr.chains = nr.chains,
						start.year = start.year, present.year = present.year, 
						wpp.year = wpp.year, my.e0.file = my.e0.file, annual.simulation = annual.simulation,
						output.dir = output.dir, mcmc.options = mcmc.options), list(...))

	if(present.year - 3 > wpp.year) 
	    warning("present.year is much larger then wpp.year. Make sure WPP data for present.year are available.")					
    data <- get.wpp.e0.data (sex, start.year = start.year, present.year = present.year, 
						wpp.year = wpp.year, my.e0.file = my.e0.file, 
						include.hiv = mcmc.options$include.hiv.countries,
						my.locations.file = my.locations.file, 
						annual = annual.simulation, verbose = verbose)
	part.ini <- .do.part.e0.mcmc.meta.ini(data, mcmc.input)
	new.meta <- c(mcmc.input, part.ini)
	if(!is.null(mcmc.options$meta.ini.fun))
	    new.meta <- do.call(mcmc.options$meta.ini.fun, list(new.meta))
	return(structure(new.meta, class = 'bayesLife.mcmc.meta'))
}

e0.mcmc.ini <- function(chain.id, mcmc.meta, iter = 100, ini.values = NULL,
                        verbose = FALSE) {
    opts <- mcmc.meta$mcmc.options
    Triangle.lim <- opts$sumTriangle.lim
    scale.Triangle <- function(Triangle) {
        # scale Triangle.ini if needed
        sTscale <- NULL
        sT <- sum(Triangle)
        if(sT > Triangle.lim[2]) sTscale <- Triangle.lim[2]
        if(sT < Triangle.lim[1]) sTscale <- Triangle.lim[1]
        if(!is.null(sTscale)) Triangle <- Triangle/sT * sTscale
        return(Triangle)
    }                            
	nr_countries <- mcmc.meta$nr.countries
    if (!exists(".Random.seed")) runif(1)
	
	ini.values[["Triangle"]] <- scale.Triangle(ini.values[["Triangle"]])
	mclist <- list(ini.values = ini.values)
	for(par in names(opts$world.parameters)) {
	    mclist[[par]] <- ini.values[[par]]
	    # legacy
	    mclist[[paste0(par, ".ini")]] <- ini.values[[par]]
	}
	mcmc <- structure(c(mclist, 
	                    list(output.dir = paste0('mc', chain.id), finished.iter = 1, length = 1,
        				    iter = iter, id = chain.id, traces = 0,
        				    traces.burnin = 0, rng.state = .Random.seed,
        				    compression.type = mcmc.meta$compression.type,
        				    meta = mcmc.meta)), 
        				class='bayesLife.mcmc')
	
    # starting values for each country
    samplpars <- mcmc.meta$country.bounds
    mcmc[['Triangle.c']] <- matrix(0, ncol = nr_countries, nrow = 4)
    for (i in 1:4)
		    mcmc[['Triangle.c']][i,] <- pmin(pmax(rnorm(nr_countries, 
		                                                mean = opts$Triangle.c$ini.norm[[1]][i], 
										                sd = opts$Triangle.c$ini.norm[[2]][i]), 
										        samplpars[[paste0('Triangle_', i, '.c.prior.low')]]), 
										    samplpars[[paste0('Triangle_', i, '.c.prior.up')]])
    mcmc[['Triangle.c']] <- apply(mcmc[['Triangle.c']], 2, scale.Triangle)
	mcmc[['k.c']] <- pmin(pmax(rnorm(nr_countries, 
	                                 opts$k.c$ini.norm[1], 
							         sd = opts$k.c$ini.norm[2]), 
							    samplpars$k.c.prior.low), 
							samplpars$k.c.prior.up)
	mcmc[['z.c']] <- pmin(pmax(rnorm(nr_countries, 
	                                 opts$z.c$ini.norm[1], 
							         sd = opts$z.c$ini.norm[2]), 
							    samplpars$z.c.prior.low), 
							samplpars$z.c.prior.up)
    return(mcmc) 
}

e0.mcmc.meta.ini.extra <- function(mcmc.set, countries = NULL, my.e0.file = NULL, 
                                   my.locations.file = NULL, burnin = 200, 
                                   country.overwrites = NULL, verbose = FALSE) {
	update.regions <- function(reg, ereg, id.replace, is.new, is.old) {
		nreg <- list()
		for(name in names(reg)) {
		    #if(!name %in% names(ereg)) next
		    if(is.character(reg[[name]]) || is.factor(reg[[name]])) {
		        reg[[name]][id.replace] <- as.character(ereg[[name]])[is.old]
		        nreg[[name]] <- c(as.character(reg[[name]]), 
		                          as.character(ereg[[name]])[is.new])
		    } else {
			    reg[[name]][id.replace] <- ereg[[name]][is.old]
			    nreg[[name]] <- c(reg[[name]], ereg[[name]][is.new])
		    }
		}
		return(nreg)
	}
	update.Tc.index <- function(Tci, eTci, id.replace, is.old) {
		nTci <- Tci
		j <- length(Tci) + 1
		old.counter <- 1
		for (i in 1:length(eTci)) {
			if (is.old[i]) {
				nTci[[id.replace[old.counter]]] <- eTci[[i]]
				old.counter <- old.counter + 1
			} else {
				nTci[[j]] <- eTci[[i]]
				j <- j+1
			}
		}
		return(nTci)
	}
	update.bounds <- function(bounds, ebounds, id.replace, is.new, is.old) {
		nbounds <- list()
		for (name in c(paste0('Triangle_', 1:4, '.c.prior.low'), 
    				paste0('Triangle_', 1:4, '.c.prior.up'), 
    				'k.c.prior.low', 'k.c.prior.up', 'z.c.prior.low','z.c.prior.up')) {
			bounds[[name]][id.replace] <- ebounds[[name]][is.old]
			nbounds[[name]] <- c(bounds[[name]], ebounds[[name]][is.new])
		}
		return(nbounds)
	}
	meta <- mcmc.set$meta
	#create e0 matrix only for the extra countries
	e0.with.regions <- set.e0.wpp.extra(meta, countries=countries, 
									  my.e0.file = my.e0.file, my.locations.file = my.locations.file, 
									  annual = meta$annual.simulation, verbose = verbose)
	if(is.null(e0.with.regions)) return(list(meta = meta, index = c()))
	# join old and new country.overwrites option; remove possible duplicates
	if(!is.null(country.overwrites)) { 
	    existing <- meta$mcmc.options$country.overwrites
	    if(!is.null(existing)) {
	        iexisting <- which(existing$country_code %in% country.overwrites$country_code)
	        if(length(iexisting) > 0)
	            existing <- existing[-iexisting,]
	    }
	    # rbind where there can be different columns in each dataset
	    meta$mcmc.options$country.overwrites <- as.data.frame(rbindlist(list(existing, country.overwrites), fill=TRUE, use.names = TRUE))
	}
	part.ini <- .do.part.e0.mcmc.meta.ini(e0.with.regions, meta)
	Emeta <- part.ini
						 		
	# join the new meta with the existing one
	is.old <- e0.with.regions$is_processed
	is.new <- !e0.with.regions$is_processed
	nold <- sum(is.old)
	nr_countries.all <- meta$nr.countries + Emeta$nr.countries - nold
	if (nold > 0) {
		codes.replace <- e0.with.regions$regions$country_code[is.old]
		id.replace <- unlist(sapply(codes.replace, get.country.object, meta=meta)['index',])
	} else {id.replace <- c()}
	new.meta <- list(nr.countries=nr_countries.all, mcmc.options = meta$mcmc.options)
					
	for (name in c('e0.matrix', 'e0.matrix.all', 'd.ct', 'loessSD')) {
		meta[[name]][,id.replace] <- Emeta[[name]][,is.old]
		new.meta[[name]] <- cbind(meta[[name]], Emeta[[name]][,is.new])
	}
	if(meta$mcmc.options$include.hiv.countries){
	    for(name in c("hiv.pred", "hiv.est")){
	        if(!is.null(meta$regions[[name]]) && is.null(Emeta$regions[[name]]))
	            Emeta$regions[[name]] <- rep(FALSE, length(Emeta$regions$country_code))
	    }
	}
	new.meta[['Tc.index']] <- update.Tc.index(meta$Tc.index, Emeta$Tc.index, id.replace, is.old)
	new.meta[['regions']] <- update.regions(meta$regions, Emeta$regions, id.replace, is.new, is.old)
	new.meta[['regionsDT']] <- create.regionsDT(new.meta[['regions']])
	if(is.null(meta$country.bounds)) { # simulation was created with previous versions of bayesLife
		meta$country.bounds <- .do.country.specific.ini(meta$nr.countries, meta)$country.bounds
	}
	new.meta[['country.bounds']] <- update.bounds(meta$country.bounds, Emeta$country.bounds, id.replace, is.new, is.old)

	if(!is.null(Emeta$suppl.data$e0.matrix)) {
		suppl.id.replace <- meta$suppl.data$index.from.all.countries[id.replace]
		suppl.id.replace <- suppl.id.replace[!is.na(suppl.id.replace)]
		suppl.is.old <- which(is.old)[which(is.element(meta$suppl.data$index.from.all.countries[id.replace], suppl.id.replace))]
		suppl.old <- Emeta$suppl.data$index.from.all.countries[suppl.is.old]
		suppl.is.new <- which(is.new & !is.na(Emeta$suppl.data$index.from.all.countries))
		suppl.new <- Emeta$suppl.data$index.from.all.countries[suppl.is.new]
		for (name in c('e0.matrix', 'd.ct', 'loessSD')) {
			meta$suppl.data[[name]][,suppl.id.replace] <- Emeta$suppl.data[[name]][,suppl.old]
			new.meta$suppl.data[[name]] <- cbind(meta$suppl.data[[name]], Emeta$suppl.data[[name]][,suppl.new])
		}
		suppl.is.old.tmp <- rep(FALSE, Emeta$suppl.data$nr.countries)
		suppl.is.old.tmp[suppl.is.old] <- TRUE
		new.meta$suppl.data$Tc.index <- update.Tc.index(meta$suppl.data$Tc.index, Emeta$suppl.data$Tc.index, 
											suppl.id.replace, suppl.is.old.tmp)
		new.meta$suppl.data$regions <- update.regions(meta$suppl.data$regions, Emeta$suppl.data$regions, 
												suppl.id.replace, suppl.new, suppl.old)
		new.meta$suppl.data$regionsDT <- create.regionsDT(new.meta$suppl.data$regions)
		n.new <- ncol(new.meta$suppl.data$e0.matrix) - ncol(meta$suppl.data$e0.matrix)
		new.meta$suppl.data$index.from.all.countries <- meta$suppl.data$index.from.all.countries
		new.meta$suppl.data$index.to.all.countries <- meta$suppl.data$index.to.all.countries
		new.meta$suppl.data$nr.countries <- ncol(new.meta$suppl.data$e0.matrix)
		if (n.new > 0) {
			new.meta$suppl.data$index.from.all.countries <- c(new.meta$suppl.data$index.from.all.countries, rep(NA, sum(is.new)))
			new.meta$suppl.data$index.from.all.countries[meta$nr.countries + suppl.is.new] <- seq(meta$suppl.data$nr.countries + 1, 
												length=n.new)
			new.meta$suppl.data$index.to.all.countries <- c(new.meta$suppl.data$index.to.all.countries, 
											seq(meta$nr.countries+1, new.meta$nr.countries)[suppl.is.new])
		}
	}
	index <- id.replace
	if (new.meta$nr.countries > meta$nr.countries) 
		index <- c(index, seq(meta$nr.countries+1, new.meta$nr.countries))
	for (item in names(new.meta)) {
		meta[[item]] <- new.meta[[item]]
	}

	return(list(meta=meta, index=index, index.replace=id.replace))
}

e0.mcmc.ini.extra <- function(mcmc, countries, index.replace = NULL) {
    opts <- mcmc$meta$mcmc.options
	nr.countries.extra <- length(countries)
	nreplace <- length(index.replace)
	if(nreplace > 0) {
    	for (i in 1:4)		
			mcmc$Triangle.c[i,index.replace] <- pmin(pmax(rnorm(nreplace, 
			                                        mean = opts$Triangle.c$ini.norm[[1]][i], 
										            sd=opts$Triangle.c$ini.norm[[2]][i]),
										        opts$Triangle.c$prior.low[i]), 
										    opts$Triangle.c$prior.up[i])
		mcmc$k.c[index.replace] <- pmin(pmax(rnorm(nreplace, opts$k.c$ini.norm[1], 
							    sd = opts$k.c$ini.norm[2]), opts$k.c$prior.low), 
							opts$k.c$prior.up)
		mcmc$z.c[index.replace] <- pmin(pmax(rnorm(nreplace, opts$z.c$ini.norm[1], 
							    sd = opts$z.c$ini.norm[2]), opts$z.c$prior.low), 
							opts$z.c$prior.up)
	}
	samplpars <- mcmc$meta$country.bounds
	if(nr.countries.extra > nreplace) {
		nextra <- nr.countries.extra-nreplace
		eidx <- (ncol(mcmc$Triangle.c)+1):(ncol(mcmc$Triangle.c)+nextra)
		mcmc$Triangle.c <- cbind(mcmc$Triangle.c, matrix(0, ncol=nextra, nrow=4))
		for (i in 1:4)		
			mcmc$Triangle.c[i,eidx] <- pmin(pmax(rnorm(nextra, 
			                              mean = opts$Triangle.c$ini.norm[[1]][i], 
										  sd = opts$Triangle.c$ini.norm[[2]][i]),
										samplpars[[paste0('Triangle_', i, '.c.prior.low')]][eidx]), 
										samplpars[[paste0('Triangle_', i, '.c.prior.up')]][eidx])
		mcmc$k.c <- c(mcmc$k.c, pmin(pmax(rnorm(nextra, opts$k.c$ini.norm[1], 
							sd = opts$k.c$ini.norm[2]), 
							samplpars$k.c.prior.low[eidx]), samplpars$k.c.prior.up[eidx]))
		mcmc$z.c <- c(mcmc$z.c, pmin(pmax(rnorm(nextra, opts$z.c$ini.norm[1], 
							sd = opts$z.c$ini.norm[2]), samplpars$z.c.prior.low[eidx]), samplpars$z.c.prior.up[eidx]))
	}
	return(mcmc)
}

e0.mcmc.meta.ini.subnat <- function(meta, country, start.year = 1950, present.year = 2020, 
                             my.e0.file = NULL, annual = NULL, verbose = FALSE) {
    meta$start.year <- start.year
    meta$present.year <- present.year
    if(is.null(annual)) annual <- meta$annual.simulation
    meta$annual.simulation <- annual
    
    data <- get.wpp.e0.subnat(country, start.year = start.year, present.year = present.year, 
                             my.e0.file = my.e0.file, annual = annual)
    data$nr.countries <- ncol(data$e0.matrix)
    data$Tc.index <- .get.Tcindex(data$e0.matrix, cnames=data$regions$country_name, stop.if.less.than2 = FALSE)
    data$subnat <- TRUE
    bounds <- .do.country.specific.ini(data$nr.countries, c(data, meta))
    this.meta <- c(data, bounds)
    for (item in names(meta))
        if(!(item %in% names(this.meta))) this.meta[[item]] <- meta[[item]]
    return(structure(this.meta, class = 'bayesLife.mcmc.meta'))
}

