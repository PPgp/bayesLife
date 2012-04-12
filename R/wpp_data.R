get.wpp.e0.data <- function(sex='M', start.year=1950, present.year=2010, wpp.year=2008, my.e0.file=NULL, 
							verbose=FALSE) {
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
	locations <- bayesTFR:::read.UNlocations(data, wpp.year=wpp.year, package='bayesLife', verbose=verbose)
	loc_data <- locations$loc_data
	include <- locations$include
	prediction.only <- locations$prediction.only

	data_incl <- data[include,]
	nr_countries_estimation <- nrow(data_incl)
	if(any(!is.na(prediction.only))) { # move prediction countries at the end of data
		data_prediction <- data[prediction.only,]
		data_incl <- rbind(data_incl, data_prediction)
	}
	
	LEXmatrix.regions <- bayesTFR:::get.observed.time.matrix.and.regions(
							data_incl, loc_data, 
							start.year=start.year, 
							present.year=present.year)
	if (verbose) 
		cat('Dimension of the e0 matrix:', dim(LEXmatrix.regions$obs_matrix), '\n')

	LEXmatrixsuppl.regions <- .get.suppl.matrix.and.regions(un.object, LEXmatrix.regions, loc_data, 
									start.year, present.year)
	if(!is.null(un.object$suppl.data.object) && verbose) 
		cat('Dimension of the supplemental e0 matrix:', dim(LEXmatrixsuppl.regions$obs_matrix), '\n')
									
	return(list(e0.matrix=LEXmatrix.regions$obs_matrix, 
				e0.matrix.all=LEXmatrix.regions$obs_matrix_all, 
				regions=LEXmatrix.regions$regions, 
				nr.countries.estimation=nr_countries_estimation,
				suppl.data=list(
					e0.matrix=if(!is.null(LEXmatrixsuppl.regions)) LEXmatrixsuppl.regions$obs_matrix else NULL,
					regions=if(!is.null(LEXmatrixsuppl.regions)) LEXmatrixsuppl.regions$regions else NULL,
					index.to.all.countries=if(!is.null(LEXmatrixsuppl.regions)) 
									LEXmatrixsuppl.regions$all.countries.index else NULL,
					index.from.all.countries=if(!is.null(LEXmatrixsuppl.regions)) 
									LEXmatrixsuppl.regions$index.from.all.countries else NULL)
				))
}

.get.suppl.matrix.and.regions <- function(un.object, LEXmatrix.regions, loc_data, start.year, present.year) {
	LEXmatrixsuppl.regions <- NULL
	if(is.null(un.object$suppl.data.object)) return(NULL)
	suppl.data <- un.object$suppl.data.object$data
	include <- which(is.element(suppl.data[,'country_code'], LEXmatrix.regions$regions$country_code))
	if(length(include) > 0)
		LEXmatrixsuppl.regions <- bayesTFR:::get.observed.time.matrix.and.regions(
							suppl.data[include,], loc_data, 
							start.year=start.year, 
							present.year=present.year)
	if(!is.null(LEXmatrixsuppl.regions)) {
		LEXmatrixsuppl.regions$all.countries.index <- c()
		index.from.all <- rep(NA, length(LEXmatrix.regions$regions$country_code))
		for(i in 1:length(suppl.data[include,'country_code'])) {
			incl.idx <- which(is.element(LEXmatrix.regions$regions$country_code,
												suppl.data[include,'country_code'][i]))
			LEXmatrixsuppl.regions$all.countries.index <- c(LEXmatrixsuppl.regions$all.countries.index, incl.idx)
			index.from.all[incl.idx] <- i
		}
		LEXmatrixsuppl.regions$index.from.all.countries <- index.from.all
	}			
	return(LEXmatrixsuppl.regions)
}


read.UNe0 <- function(sex, wpp.year, my.e0.file=NULL, ...) {
	un.file.name <- file.path(.find.package("bayesLife"), "data", paste('UN', wpp.year, 'e0', sex, '.txt', sep=''))
	un.suppl.file.name <- file.path(.find.package("bayesLife"), "data", paste('UN', wpp.year, 'e0', sex, '_supplemental.txt', sep=''))
	data <- bayesTFR:::do.read.un.file(un.file.name, wpp.year, my.file=my.e0.file, ...)
	suppl.data<- NULL
	if(file.exists(un.suppl.file.name)) 
		suppl.data <- bayesTFR:::do.read.un.file(un.suppl.file.name, wpp.year, ...)
	return(list(data.object=data, suppl.data.object=suppl.data))
}

set.e0.wpp.extra <- function(meta, countries=NULL, my.e0.file=NULL, verbose=FALSE) {
	#'countries' is a vector of country or region codes 
	un.object <- read.UNe0(sex=meta$sex, wpp.year=meta$wpp.year, my.e0.file=my.e0.file, 
							present.year=meta$present.year, verbose=verbose)
	data <- un.object$data.object
	extra.wpp <- bayesTFR:::.extra.matrix.regions(data=data, countries=countries, meta=meta, 
							package="bayesLife", verbose=verbose)
	if(!is.null(extra.wpp)) {
		extra.wpp <- list(e0.matrix=extra.wpp$tfr_matrix, 
						  e0.matrix.all=extra.wpp$tfr_matrix_all, 
						  regions=extra.wpp$regions, 
						  nr.countries.estimation=extra.wpp$nr_countries_estimation,
						  is_processed = extra.wpp$is_processed)
		locations <- bayesTFR:::read.UNlocations(data$data, wpp.year=meta$wpp.year, package='bayesLife', verbose=verbose)
		suppl.wpp <- .get.suppl.matrix.and.regions(un.object, extra.wpp, locations$loc_data, 
									meta$start.year, meta$present.year)
		extra.wpp$suppl.data=list(
			e0.matrix=if(!is.null(suppl.wpp)) suppl.wpp$obs_matrix else NULL,
			regions=if(!is.null(suppl.wpp)) suppl.wpp$regions else NULL,
			index.to.all.countries=if(!is.null(suppl.wpp)) suppl.wpp$all.countries.index else NULL,
			index.from.all.countries=if(!is.null(suppl.wpp)) suppl.wpp$index.from.all.countries else NULL)
	}
	return(extra.wpp)
}

get.wpp.e0.data.for.countries <- function(meta, sex='M', my.e0.file=NULL, verbose=FALSE) {
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
	locations <- bayesTFR:::read.UNlocations(data, wpp.year=meta$wpp.year, package='bayesLife', verbose=verbose)
	loc_data <- locations$loc_data
	include <- c()
	for (i in 1:length(meta$regions$country_code)) { # put countries into the same order as in meta
		loc_index <- which(data$country_code == meta$regions$country_code[i])
		if(length(loc_index) <= 0) 
			stop('Country ', data$country_code[i], ' not found.')
		include <- c(include, loc_index)
	}
	data_incl <- data[include,]	
	LEXmatrix.regions <- bayesTFR:::get.observed.time.matrix.and.regions(
							data_incl, loc_data, 
							start.year=meta$start.year, 
							present.year=meta$present.year)
	if (verbose) 
		cat('Dimension of the e0 matrix:', dim(LEXmatrix.regions$obs_matrix), '\n')
	LEXmatrixsuppl.regions <- .get.suppl.matrix.and.regions(un.object, LEXmatrix.regions, loc_data, 
									meta$start.year, meta$present.year)
	if(!is.null(un.object$suppl.data.object) && verbose) 
		cat('Dimension of the supplemental e0 matrix:', dim(LEXmatrixsuppl.regions$obs_matrix), '\n')
											
	return(list(e0.matrix=LEXmatrix.regions$obs_matrix, 
				e0.matrix.all=LEXmatrix.regions$obs_matrix_all, 
				regions=LEXmatrix.regions$regions, 
				nr.countries.estimation=nrow(data_incl),
				suppl.data=list(
					e0.matrix=if(!is.null(LEXmatrixsuppl.regions)) LEXmatrixsuppl.regions$obs_matrix else NULL,
					regions=if(!is.null(LEXmatrixsuppl.regions)) LEXmatrixsuppl.regions$regions else NULL,
					index.to.all.countries=if(!is.null(LEXmatrixsuppl.regions)) 
									LEXmatrixsuppl.regions$all.countries.index else NULL,
					index.from.all.countries=if(!is.null(LEXmatrixsuppl.regions)) 
									LEXmatrixsuppl.regions$index.from.all.countries else NULL)
				))
}
