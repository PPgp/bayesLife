get.wpp.e0.data <- function(sex='M', start.year=1950, present.year=2010, wpp.year=2008, my.e0.file=NULL, 
							verbose=FALSE) {
	sex <- toupper(sex)
	if(sex != 'M' && sex != 'F')
		stop('Allowed values for argument "sex" are "M" and "F".')
	########################################
	# set data and match with areas
	########################################
	data <- read.UNe0(sex=sex, wpp.year=wpp.year, my.e0.file=my.e0.file, 
								present.year=present.year, verbose=verbose)$data
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
	return(list(e0.matrix=LEXmatrix.regions$obs_matrix, 
				e0.matrix.all=LEXmatrix.regions$obs_matrix_all, 
				regions=LEXmatrix.regions$regions, 
				nr.countries.estimation=nr_countries_estimation))
}


read.UNe0 <- function(sex, wpp.year, my.e0.file=NULL, ...) {
	un.file.name <- file.path(.find.package("bayesLife"), "data", paste('UN', wpp.year, 'e0', sex, '.txt', sep=''))
	return(bayesTFR:::do.read.un.file(un.file.name, wpp.year, my.file=my.e0.file, ...))
}

set.e0.wpp.extra <- function(meta, countries=NULL, my.e0.file=NULL, verbose=FALSE) {
	#'countries' is a vector of country or region codes 
	data <- read.UNe0(sex=meta$sex, wpp.year=meta$wpp.year, my.e0.file=my.e0.file, 
							present.year=meta$present.year, verbose=verbose)
	extra.wpp <- bayesTFR:::.extra.matrix.regions(data=data, countries=countries, meta=meta, 
							package="bayesLife", verbose=verbose)
	if(!is.null(extra.wpp))
		extra.wpp <- list(e0.matrix=extra.wpp$tfr_matrix, 
						  e0.matrix.all=extra.wpp$tfr_matrix_all, 
						  regions=extra.wpp$regions, 
						  nr.countries.estimation=extra.wpp$nr_countries_estimation,
						  is_processed = extra.wpp$is_processed)
	return(extra.wpp)
}
