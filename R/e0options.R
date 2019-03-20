if(getRversion() >= "2.15.1") utils::globalVariables(c("loess_sd", ".e0options"))
data(loess_sd, envir = environment())

e0options <- function()
    .e0options

e0mcmc.options <- function(...) {
    e0.options("mcmc", ...)
}

e0pred.options <- function(...) {
    e0.options("pred", ...)
}


e0.options <- function(what, ...) {
    # this code was adapted from mclust.options 
    current <- .e0options
    if (nargs() == 1) return(current[[what]])
    args <- list(...)
    if (length(args) == 1 && is.null(names(args))) {
        arg <- args[[1]]
        switch(mode(arg), 
               list = args <- arg, 
               character = return(current[[what]][[arg]]), 
               stop("Invalid argument: ", dQuote(arg))
        )
    }
    if (length(args) == 0) return(current[[what]])
    if (is.null(names(args))) {
        if(is.list(args) && mode(unlist(args)) == "character") # retrieving multiple options
            return(current[[what]][unlist(args)])
        stop("Options must be given by name")
    }
    old.opts <- current[[what]]
    current[[what]] <- modifyList(old.opts, args)
    .e0options <<- current
    invisible(old.opts)
}

e0.options.default <- function() {
    structure(list(
        mcmc = e0.mcmc.options.default(),
        pred = e0.pred.options.default(),
        admin = list(package = "bayesLifeHIV")
    ))
}

e0.mcmc.options.default <- function() {
    pars <- list(
        a = c(13.215, 41.070, 9.235, 17.605, 2.84, 0.385),
        #a=c(15.7669391,40.9658241,0.2107961,19.8188061,2.9306625,0.400688628),
        delta = c(3.844, 4.035, 11.538, 5.639, 0.901, 0.4),
        #delta=c(1.887, 1.982, 1.99, 1.949, 0.995, 0.4), 
        tau = c(15.5976503,23.6500060,14.5056919,14.7185980,3.4514285,0.5667531),
        Triangle = structure(
            list(ini = list(T1 = NULL, T2 = NULL, T3 = NULL, T4 = NULL),
                 ini.low = c(10, 30, 0.1, 10),
                 ini.up  = c(30, 50, 10, 30),
                 prior.low = c(0, 0, -20, 0),
                 prior.up  = c(100, 100, 50, 100),
                 slice.width = c(10, 10, 10, 10)
            ), npar = 4),
        k = structure(
            list(ini = NULL, ini.low = 3, ini.up = 5, 
                 prior.low = 0, prior.up = 10
            ),
            npar = 1),
        z = structure(
            list(ini = NULL, ini.low = 0.0001, ini.up = 0.653, 
                 prior.low = 0, prior.up = 0.653, slice.width = 1),
            npar = 1),
        lambda = structure(
            list(ini = list(T1 = NULL, T2 = NULL, T3 = NULL, T4 = NULL),
                 ini.low = c(0.01, 0.01, 0.01, 0.01),
                 ini.up =  c(0.1, 0.1, 0.1, 0.1),
                 slice.width = c(0.1,0.1,0.1,0.1)
            ), npar = 4),
        lambda.k = structure(list(ini = NULL, ini.low = 0.3, ini.up = 1), npar = 1),
        lambda.z = structure(list(ini = NULL, ini.low = 1, ini.up = 40,
                                  slice.width = 10), npar = 1),
        omega = structure(list(ini = NULL, ini.low = 0.1, ini.up = 5,
                               slice.width = 1), npar = 1),
        Triangle.c = structure(
            list(ini.norm = list(mean = NULL, sd = c(2, 2, 2, 2)),
                 prior.low = c(0, 0, -20, 0), 
                 prior.up  = c(100, 100, 50, 100),
                 slice.width = c(10, 10, 10, 10)
            ), npar = 4),
        k.c = structure(list(ini.norm = c(mean = NA, sd = 2), 
                             prior.low = 0, 
                             prior.up = 10,
                             slice.width = 2), npar = 1),
        z.c = structure(list(ini.norm = c(mean = NA, sd = 0.2), 
                             prior.low = 0, prior.up = 0.653,
                             slice.width = 1), npar = 1),
        world.parameters = c(Triangle = 4, k = 1, z = 1, lambda = 4, 
                             lambda.k = 1, lambda.z = 1, omega = 1),
        country.parameters = c(Triangle.c = 4, k.c = 1, z.c = 1),
        country.overwrites = NULL,
        nu = 4, dl.p1 = 9, dl.p2 = 9, 
        sumTriangle.lim = c(30, 86),
        outliers = c(-5, 10),
        buffer.size = 100,
        auto.conf = list(max.loops = 5, iter = 160000, iter.incr = 20000, 
                         nr.chains = 3, thin = 225, burnin = 10000),
        estimation.function = "e0.mcmc.sampling",
        dlcurves.function = "e0.get.dlcurves",
        parallel.init.function = function(){library(bayesLife)},
        include.hiv.countries = FALSE
    )
    Triangle <- k <- z <- NULL # to avoid R check note "no visible binding ..."
    pars <- within(pars, {
        Triangle.c$ini.norm[["mean"]] <- round(Triangle$ini.low + (Triangle$ini.up - Triangle$ini.low)/2)
        k.c$ini.norm["mean"] <- round(k$ini.low + (k$ini.up - k$ini.low)/2)
        z.c$ini.norm["mean"] <- round(z$ini.low + (z$ini.up - z$ini.low)/2, 2)
    })
    pars
}

e0.pred.options.default <- function() {
    pars <- list(
        quantiles = c(0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 
                      0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 1)
    )
    pars
}

get.from.options <- function(name, opts, default = NULL) {
    # Return specific option given by name from a list of options (opts). 
    # If it's NULL, return the given default. 
    value <- opts[[name]]
    if(is.null(value)) value <- default
    return(value)
}

using.bayesLife <- function() {
    .e0options <<- e0.options.default()
}

.e0options <- e0.options.default()
