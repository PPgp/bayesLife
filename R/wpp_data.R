get.wpp.e0.data <- function(sex = 'M', start.year = 1950, present.year = 2015, 
                            wpp.year = 2017, my.e0.file = NULL, include.hiv = FALSE,
							my.locations.file = NULL, annual = FALSE, verbose = FALSE) {
	sex <- toupper(sex)
	if(sex != 'M' && sex != 'F')
		stop('Allowed values for argument "sex" are "M" and "F".')
	########################################
	# set data and match with areas
	########################################
	un.object <- read.UNe0(sex=sex, wpp.year=wpp.year, my.e0.file=my.e0.file, 
								present.year=present.year, annual = annual,
								verbose=verbose)
	data <- un.object$data.object$data
	# get region and area data
	locations <- bayesTFR:::read.UNlocations(data, wpp.year=wpp.year, my.locations.file=my.locations.file,
											package='bayesLife', verbose=verbose)
	loc_data <- locations$loc_data
	include <- locations$include
	prediction.only <- locations$prediction.only

	# include HIV/AIDS countries
	if(include.hiv) {
	    hiv.aids <- rep(FALSE, length(include))
	    # find HIV/AIDS countries
	    hivincl <- merge(data[,c("country_code", "include_code")], 
	                     unique(locations$loc_data[,c("country_code", "include_code")]), 
	                     all.x = TRUE, by = "country_code", sort = FALSE)
	    hivincl$include_code <- ifelse(hivincl$include_code.x >= 0, 
	                                   hivincl$include_code.x, hivincl$include_code.y)
	    hiv.aids <- hivincl$include_code == 3
	    include[hiv.aids] <- TRUE
	    hiv.aids <- hiv.aids[include]
	}
	
	data_incl <- data[include,]
	nr_countries_estimation <- nrow(data_incl)
	if(any(!is.na(prediction.only))) { # move prediction countries at the end of data
		data_prediction <- data[prediction.only,]
		data_incl <- rbind(data_incl, data_prediction)
		if(include.hiv)
		    hiv.aids <- c(hiv.aids, rep(FALSE, nrow(data_prediction)))
	}
	
	LEXmatrix.regions <- bayesTFR:::get.observed.time.matrix.and.regions(
							data_incl, loc_data, 
							start.year=start.year, 
							present.year=present.year, annual = annual, 
							interpolate = annual && is.null(my.e0.file))
	if(include.hiv)
    	LEXmatrix.regions$regions$hiv.pred <- hiv.aids
	
	if (verbose) 
		cat('Dimension of the e0 matrix:', dim(LEXmatrix.regions$obs_matrix), '\n')

	if(!annual) {
	    LEXmatrixsuppl.regions <- bayesTFR:::.get.suppl.matrix.and.regions(un.object, LEXmatrix.regions, loc_data, 
									start.year, present.year)
	    if(!is.null(un.object$suppl.data.object) && verbose) 
		    cat('Dimension of the supplemental e0 matrix:', dim(LEXmatrixsuppl.regions$obs_matrix), '\n')
	} else LEXmatrixsuppl.regions <- NULL
	suppl.data <- bayesTFR:::.get.suppl.data.list(LEXmatrixsuppl.regions, matrix.name='e0.matrix')
	suppl.data$regionsDT <- create.regionsDT(suppl.data$regions)
	return(list(e0.matrix = LEXmatrix.regions$obs_matrix, 
				e0.matrix.all = LEXmatrix.regions$obs_matrix_all, 
				regions = LEXmatrix.regions$regions, 
				regionsDT = create.regionsDT(LEXmatrix.regions$regions),
				nr.countries.estimation = nr_countries_estimation,
				suppl.data = suppl.data
				))
}

create.regionsDT <- function(reglist) {
    # convert regions from list to data.table
    if(is.null(reglist)) return(NULL)
    dt <- as.data.table(reglist)
    # convert factors to character
    fcols <- colnames(dt)[sapply(dt,class) == "factor"]
    if(length(fcols) > 0)
        dt[,(fcols):= lapply(.SD, as.character), .SDcols = fcols]
    return(dt)
}
    
read.UNe0 <- function(sex, wpp.year, my.e0.file=NULL, ...) {
	un.dataset <- paste('e0', sex, sep='')
	un.suppl.dataset <- paste('e0', sex, '_supplemental', sep='')
	data <- bayesTFR:::do.read.un.file(un.dataset, wpp.year, my.file=my.e0.file, ...)
	suppl.data <- bayesTFR:::do.read.un.file(un.suppl.dataset, wpp.year, my.file=my.e0.file, ...)
	if(is.null(suppl.data$data)) suppl.data <- NULL
	return(list(data.object=data, suppl.data.object=suppl.data))
}

set.e0.wpp.extra <- function(meta, countries=NULL, my.e0.file=NULL, my.locations.file=NULL, 
                             annual = FALSE, verbose=FALSE) {
	#'countries' is a vector of country or region codes 
	un.object <- read.UNe0(sex=meta$sex, wpp.year=meta$wpp.year, my.e0.file=my.e0.file, 
							present.year=meta$present.year, annual = annual, verbose=verbose)
	data <- un.object$data.object
	extra.wpp <- bayesTFR:::.extra.matrix.regions(data=data, countries=countries, meta=meta, 
							package="bayesLife", my.locations.file=my.locations.file, 
							annual = annual,
							interpolate = is.null(my.e0.file) && annual,
							verbose=verbose)
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
						my.e0.file=my.e0.file, annual = meta$annual.simulation, verbose=verbose)
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
							present.year=meta$present.year, annual = meta$annual, 
							interpolate = meta$annual && is.null(my.e0.file),)
	if (verbose) 
		cat('Dimension of the e0 matrix:', dim(LEXmatrix.regions$obs_matrix), '\n')
	if(is.null(meta$annual.simulation) || !meta$annual.simulation) {
	    LEXmatrixsuppl.regions <- bayesTFR:::.get.suppl.matrix.and.regions(un.object, LEXmatrix.regions, loc_data, 
									meta$start.year, meta$present.year)
	} else LEXmatrixsuppl.regions <- NULL
	if(!is.null(un.object$suppl.data.object) && verbose) 
		cat('Dimension of the supplemental e0 matrix:', dim(LEXmatrixsuppl.regions$obs_matrix), '\n')
											
	return(list(e0.matrix=LEXmatrix.regions$obs_matrix, 
				e0.matrix.all=LEXmatrix.regions$obs_matrix_all, 
				regions=LEXmatrix.regions$regions, 
				nr.countries.estimation=nrow(data_incl),
				suppl.data=bayesTFR:::.get.suppl.data.list(LEXmatrixsuppl.regions, matrix.name='e0.matrix')
				))
}

get.wpp.e0.subnat <- function(country, start.year=1950, present.year=2010, my.e0.file=NULL, annual = FALSE) {
    data <- bayesTFR:::do.read.subnat.file(my.e0.file, present.year = present.year)
    data <- data[data$country_code == country,]
    locations <- bayesTFR:::create.sublocation.dataset(data)
    loc_data <- locations$loc_data
    include <- locations$include
    prediction.only <- locations$prediction.only
    
    data_countries <- data[include | locations$prediction.only,]
    nr_countries_estimation <- sum(include)
    #stop("")
    LEXmatrix.regions <- bayesTFR:::get.observed.time.matrix.and.regions(
                                data_countries, loc_data, start.year = start.year, 
                                present.year = present.year, annual = annual, 
                                datacolnames=c(country.code='reg_code', country.name='name', reg.name='reg_name',
                                               reg.code='NA', area.name='country', area.code='country_code'))
    
    return(list(e0.matrix = LEXmatrix.regions$obs_matrix, 
                e0.matrix.all = LEXmatrix.regions$obs_matrix_all, 
                regions = LEXmatrix.regions$regions, 
                regionsDT = create.regionsDT(LEXmatrix.regions$regions),
                nr.countries.estimation = nr_countries_estimation,
                suppl.data = bayesTFR:::.get.suppl.data.list(NULL)
            ))
}

get.wpp.e0.subnat.joint <- function(country, meta, my.e0.file) {
    data <- bayesTFR:::do.read.subnat.file(my.e0.file, present.year = meta$present.year)
    data <- data[data$country_code == country,]
    locations <- bayesTFR:::create.sublocation.dataset(data)
    loc_data <- locations$loc_data
    include <- c()
    for (i in 1:length(meta$regions$country_code)) { # put regions into the same order as in meta
        loc_index <- which(data$reg_code == meta$regions$country_code[i])
        if(length(loc_index) <= 0) 
            stop('Region ', meta$regions$country_code[i], ' not found.')
        include <- c(include, loc_index)
    }
    data_countries <- data[include,]	
    
    LEXmatrix.regions <- bayesTFR:::get.observed.time.matrix.and.regions(
        data_countries, loc_data, start.year = meta$start.year, 
        present.year = meta$present.year, annual = meta$annual, 
        datacolnames=c(country.code='reg_code', country.name='name', reg.name='reg_name',
                       reg.code='NA', area.name='country', area.code='country_code'))
    
    return(list(e0.matrix = LEXmatrix.regions$obs_matrix, 
                e0.matrix.all = LEXmatrix.regions$obs_matrix_all, 
                regions = LEXmatrix.regions$regions, 
                regionsDT = create.regionsDT(LEXmatrix.regions$regions),
                nr.countries.estimation = nrow(data_countries),
                suppl.data = bayesTFR:::.get.suppl.data.list(NULL)
    ))
}