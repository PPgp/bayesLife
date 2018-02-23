get.wpp.e0.data <- function(sex='M', start.year=1950, present.year=2015, wpp.year=2017, my.e0.file=NULL, 
							include.hiv=FALSE, my.locations.file=NULL, verbose=FALSE) {
	sex <- toupper(sex)
	if(sex != 'M' && sex != 'F')
		stop('Allowed values for argument "sex" are "M" and "F".')
	########################################
	# set data and match with areas
	########################################
	un.object <- read.UNe0(sex=sex, wpp.year=wpp.year, my.e0.file=my.e0.file, 
								present.year=present.year, verbose=verbose)
	data <- un.object$data.object$data
	# get region and area data
	locations <- bayesTFR:::read.UNlocations(data, wpp.year=wpp.year, my.locations.file=my.locations.file,
											package='bayesLife', verbose=verbose)
	loc_data <- locations$loc_data
	include <- locations$include
	prediction.only <- locations$prediction.only

	# include HIV/AIDS countries
	hiv.aids <- rep(FALSE, length(include))
	if(include.hiv) {
	    # find HIV/AIDS countries
	    hivincl <- merge(data[,c("country_code", "include_code")], 
	                     locations$loc_data[,c("country_code", "include_code")], 
	                     all.x=TRUE, by="country_code", sort=FALSE)
	    hivincl$include_code <- ifelse(hivincl$include_code.x >= 0, 
	                                   hivincl$include_code.x, hivincl$include_code.y)
	    hiv.aids <- hivincl$include_code == 3
	    include[hiv.aids] <- TRUE
	}
	data_incl <- data[include,]
	hiv.aids <- hiv.aids[include]
	nr_countries_estimation <- nrow(data_incl)
	if(any(!is.na(prediction.only))) { # move prediction countries at the end of data
		data_prediction <- data[prediction.only,]
		data_incl <- rbind(data_incl, data_prediction)
		hiv.aids <- c(hiv.aids, rep(FALSE, nrow(data_prediction)))
	}
	
	LEXmatrix.regions <- bayesTFR:::get.observed.time.matrix.and.regions(
							data_incl, loc_data, 
							start.year=start.year, 
							present.year=present.year)
	LEXmatrix.regions$regions$hiv.pred <- hiv.aids
	
	if (verbose) 
		cat('Dimension of the e0 matrix:', dim(LEXmatrix.regions$obs_matrix), '\n')

	LEXmatrixsuppl.regions <- bayesTFR:::.get.suppl.matrix.and.regions(un.object, LEXmatrix.regions, loc_data, 
									start.year, present.year)
	if(!is.null(un.object$suppl.data.object) && verbose) 
		cat('Dimension of the supplemental e0 matrix:', dim(LEXmatrixsuppl.regions$obs_matrix), '\n')
									
	return(list(e0.matrix=LEXmatrix.regions$obs_matrix, 
				e0.matrix.all=LEXmatrix.regions$obs_matrix_all, 
				regions=LEXmatrix.regions$regions, 
				nr.countries.estimation=nr_countries_estimation,
				suppl.data=bayesTFR:::.get.suppl.data.list(LEXmatrixsuppl.regions, matrix.name='e0.matrix')
				))
}

read.UNe0 <- function(sex, wpp.year, my.e0.file=NULL, ...) {
	un.dataset <- paste('e0', sex, sep='')
	un.suppl.dataset <- paste('e0', sex, '_supplemental', sep='')
	data <- bayesTFR:::do.read.un.file(un.dataset, wpp.year, my.file=my.e0.file, ...)
	suppl.data <- bayesTFR:::do.read.un.file(un.suppl.dataset, wpp.year, my.file=my.e0.file, ...)
	if(is.null(suppl.data$data)) suppl.data <- NULL
	return(list(data.object=data, suppl.data.object=suppl.data))
}

set.e0.wpp.extra <- function(meta, countries=NULL, my.e0.file=NULL, my.locations.file=NULL, verbose=FALSE) {
	#'countries' is a vector of country or region codes 
	un.object <- read.UNe0(sex=meta$sex, wpp.year=meta$wpp.year, my.e0.file=my.e0.file, 
							present.year=meta$present.year, verbose=verbose)
	data <- un.object$data.object
	extra.wpp <- bayesTFR:::.extra.matrix.regions(data=data, countries=countries, meta=meta, 
							package="bayesLife", my.locations.file=my.locations.file, verbose=verbose)
	if(!is.null(extra.wpp)) {
		extra.wpp <- list(e0.matrix=extra.wpp$tfr_matrix, 
						  e0.matrix.all=extra.wpp$tfr_matrix_all, 
						  regions=extra.wpp$regions, 
						  nr.countries.estimation=extra.wpp$nr_countries_estimation,
						  is_processed = extra.wpp$is_processed)
		locations <- bayesTFR:::read.UNlocations(data$data, wpp.year=meta$wpp.year, 
									my.locations.file=my.locations.file, package='bayesLife', verbose=verbose)
		suppl.wpp <- bayesTFR:::.get.suppl.matrix.and.regions(un.object, extra.wpp, locations$loc_data, 
									meta$start.year, meta$present.year)
		extra.wpp$suppl.data <- bayesTFR:::.get.suppl.data.list(suppl.wpp, matrix.name='e0.matrix')
	}
	return(extra.wpp)
}

get.wpp.e0.data.for.countries <- function(meta, sex='M', my.e0.file=NULL, my.locations.file=NULL, verbose=FALSE) {
	sex <- toupper(sex)
	if(sex != 'M' && sex != 'F')
		stop('Allowed values for argument "sex" are "M" and "F".')
	########################################
	# set data and match with areas
	########################################
	un.object <- read.UNe0(sex=sex, wpp.year=meta$wpp.year, present.year=meta$present.year, 
						my.e0.file=my.e0.file, verbose=verbose)
	data <- un.object$data.object$data
	# get region and area data
	locations <- bayesTFR:::read.UNlocations(data, wpp.year=meta$wpp.year, 
							my.locations.file=my.locations.file, package='bayesLife', verbose=verbose)
	loc_data <- locations$loc_data
	include <- c()
	for (i in 1:length(meta$regions$country_code)) { # put countries into the same order as in meta
		loc_index <- which(data$country_code == meta$regions$country_code[i])
		if(length(loc_index) <= 0) 
			stop('Country ', meta$regions$country_code[i], ' not found.')
		include <- c(include, loc_index)
	}
	data_incl <- data[include,]	
	LEXmatrix.regions <- bayesTFR:::get.observed.time.matrix.and.regions(
							data_incl, loc_data, 
							start.year=meta$start.year, 
							present.year=meta$present.year)
	if (verbose) 
		cat('Dimension of the e0 matrix:', dim(LEXmatrix.regions$obs_matrix), '\n')
	LEXmatrixsuppl.regions <- bayesTFR:::.get.suppl.matrix.and.regions(un.object, LEXmatrix.regions, loc_data, 
									meta$start.year, meta$present.year)
	if(!is.null(un.object$suppl.data.object) && verbose) 
		cat('Dimension of the supplemental e0 matrix:', dim(LEXmatrixsuppl.regions$obs_matrix), '\n')
											
	return(list(e0.matrix=LEXmatrix.regions$obs_matrix, 
				e0.matrix.all=LEXmatrix.regions$obs_matrix_all, 
				regions=LEXmatrix.regions$regions, 
				nr.countries.estimation=nrow(data_incl),
				suppl.data=bayesTFR:::.get.suppl.data.list(LEXmatrixsuppl.regions, matrix.name='e0.matrix')
				))
}

scale.hiv.trajectories <- function(trajectories = NULL, scale.to = NULL,
                                   logit.adjust = 0.001) {
    # Scale given trajectories to a data frame given by scale.to.
    # scale.to should have a column country_code. Other columns should 
    # be named by time periods, e.g. 2025-2030, 2030-2035 etc.
    # Argument trajectories is a data frame with columns country_code, Trajectory and 
    # time periods. The country codes and time periods are matched 
    # between the two datasets. There must be the same number of trajectories 
    # for all countries.
    # If trajectories is not given, dataset HIVprevTrajectories is used. 
    # Default for scale.to is dataset HIVprevalence.
    # The scaling is done using an adjusted logit.
    
    rep.row <- function(x,n) matrix(rep(x,each=n), nrow=n)

    expit <- function(x, adjust) { # inverse logit with adjustment
        a <- 1 - 2 * adjust
        (exp(x)- 1)/(2*a*(1+exp(x))) + 1/2
    }
    env <- new.env()
    if(is.null(scale.to)) {
        data("HIVprevalence", envir = env)
        scale.to <- env$HIVprevalence
    }
    if(is.null(trajectories)) {
        data("HIVprevTrajectories", envir = env)
        trajectories <- env$HIVprevTrajectories
    }
    if(! "country_code" %in% colnames(trajectories))
        stop("Column country_code is missing in dataset trajectories.")
    if(! "country_code" %in% colnames(scale.to))
        stop("Column country_code is missing in dataset scale.to.")
    
    cntries <- sort(intersect(scale.to$country_code, trajectories$country_code))
    cntries.char <- as.character(cntries)
    cols <- intersect(colnames(scale.to), colnames(trajectories))
    dont.include <- c("include_code", "name", "Trajectory")
    if(any(cols %in% dont.include))
        cols <- cols[-which(cols %in% dont.include)]
    if(length(cols) <= 1)
        stop("Datasets trajectories and scale.to do not have time periods in common. Check column names.")
    
    unhiv <- scale.to[scale.to$country_code %in% cntries, cols]
    unhivmtx <- unhiv[, -which(colnames(unhiv) == "country_code")]
    rownames(unhivmtx) <- unhiv$country_code
    # reorder so that countries are sorted
    unhivmtx <- unhivmtx[cntries.char, ]
    
    if(! "Trajectory" %in% colnames(trajectories))
        stop("Column Trajectory is missing in dataset trajectories.")
    
    trajs <- trajectories[order(trajectories$country_code, trajectories$Trajectory),]
    trajs <- trajs[trajs$country_code %in% cntries, ]
    trajs.red <- trajs[, cols]
    trajs.red <- trajs.red[,-which(colnames(trajs.red) == "country_code")]
    spl.trajs <- split(trajs.red, trajs$country_code)
    # compute medians for each country
    trajs.med.spl <- lapply(spl.trajs, function(m) apply(m, 2, median))
    trajs.med <- do.call(rbind, trajs.med.spl)
    # reorder so that countries are sorted
    trajs.med.mtx <- trajs.med[cntries.char,]
    
    # Put the median dataset an the scale.to dataset into the same shape 
    # as trajectories by repeating rows nr.traj times.
    nr.trajs <- length(unique(trajectories$Trajectory))
    trajs.med.big <- apply(trajs.med.mtx, 2, rep.row, nr.trajs)
    trajs.mtx <- as.matrix(trajs.red)
    unhiv.big <- apply(unhivmtx, 2, rep.row, nr.trajs)
    
    x <- (logit(trajs.mtx/100, adjust=logit.adjust) + 
              logit(unhiv.big/100, adjust=logit.adjust) - 
              logit(trajs.med.big/100, adjust=logit.adjust))
    scaled <- 100*pmax(expit(x, adjust=logit.adjust), 0)
    scaled <- cbind(data.frame(country_code = trajs$country_code,
                    Trajectory = trajs$Trajectory), 
                    scaled)
    return(scaled)
}
