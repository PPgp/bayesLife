e0.predict.subnat <- function(countries, my.e0.file, sim.dir=file.path(getwd(), 'bayesLife.output'),
                              method = c("ar1", "shift", "scale"), predict.jmale = FALSE, my.e0M.file = NULL,
                               end.year=2100, start.year=NULL, output.dir = NULL, nr.traj=NULL, seed = NULL, 
                               ar.pars = c(rho = 0.95, U = 82.5, a = 0.0482, b = -0.0154), 
                               save.as.ascii = 0, verbose = TRUE, jmale.estimates = NULL, ...) {
    # Run subnational projections, using the Scale AR(1) model applied to a national bayesLife simulation 
    # sim.dir is the world-national simulation. Set output.dir to store results somewhere else.
    get.sd <- function(e0) {
        v <- if(e0 >= ar.pars["U"]) ar.pars["a"] else ar.pars["a"] + ar.pars["b"] * (e0 - ar.pars["U"])
        return(sqrt(v))
    }
    generate.trajectory.ar1 <- function(e0c, alpha0) {
        alpha <- c(alpha0, rep(NA, length(e0c)))
        if(is.na(alpha0)) return(alpha[-1])
        for(i in 2:(length(e0c)+1))
            alpha[i] <- ar.pars["rho"]*alpha[i-1] + rnorm(1, 0, get.sd(e0c[i-1]))
        return(alpha[-1] + e0c)
    }
    generate.trajectory.scale <- function(e0c, alpha0) {
        return(e0c * alpha0)
    }
    generate.trajectory.shift <- function(e0c, alpha0) {
        return(e0c + alpha0)
    }
    compute.alpha.ar1 <- function(e0, e0c)
        return(e0 - e0c)
    compute.alpha.scale <- function(e0, e0c)
        return(e0/e0c)
    compute.alpha.shift <- function(...)
        return(compute.alpha.ar1(...))
    
    joint.male.est <- subnat.gap.estimates()
    
    method <- match.arg(method)
    wpred <- get.e0.prediction(sim.dir) # contains national projections
    wdata <- wpred$e0.matrix.reconstructed
    wmeta <- wpred$mcmc.set$meta
    year.step <- if(wmeta$annual.simulation) 1 else 5
    if(!is.null(seed)) set.seed(seed)
    if (is.null(output.dir)) output.dir <- wmeta$output.dir
    quantiles.to.keep <- as.numeric(dimnames(wpred$quantiles)[[2]])
    e <- new.env()
    result <- list()
    orig.nr.traj <- nr.traj
    sex <- list(F = 'female', M = 'male')
    for (country in countries) {
        country.obj <- get.country.object(country, wmeta)
        if(is.null(country.obj$code)) {
            warning("Country ", country, " not found in the national projections. Use numerical codes or exact names.")
            next
        }
        if(verbose) cat('\n', country.obj$name, ': predicting e0 for', sex[[wmeta$sex]])
        starty <- if(is.null(start.year)) as.integer(dimnames(wdata)[[1]][wpred$present.year.index]) + year.step else start.year
        meta <- e0.mcmc.meta.ini.subnat(wmeta, country = country.obj$code, my.e0.file = my.e0.file, 
                                 start.year = 1900, present.year = starty - year.step, verbose = verbose)
        class(meta) <- "bayesLife.mcmc.meta"
        this.output.dir <- file.path(output.dir, 
                                     paste0('subnat_', method), paste0('c', country.obj$code))
        outdir <- file.path(this.output.dir, 'predictions')
        meta$output.dir <- this.output.dir
    
        wtrajs <- get.e0.trajectories(wpred, country.obj$code)
        nr.traj <- orig.nr.traj
        if(is.null(nr.traj)) nr.traj <- ncol(wtrajs)
        thinning.index <- round(seq(1, ncol(wtrajs), length = nr.traj))
        wtrajs <- wtrajs[as.integer(rownames(wtrajs)) <= end.year, thinning.index]
        nr.traj <- ncol(wtrajs)
        wyears <- as.integer(rownames(wtrajs))
        wend.year <- max(wyears)
        we0 <- wdata[, country.obj$index]
        we0obsy <- as.integer(names(we0))[-length(we0)]
        if(!wmeta$annual.simulation){
            seqy <- seq(min(we0obsy) - 3, max(wyears) + 2, by = 5)
            midy <- seqy + 3
        } else{
            seqy <- midy <-(min(we0obsy)-1):max(wyears)
        }
        presenty <- midy[cut(starty, seqy, labels=FALSE)-1]
        if(any(wyears < presenty))  # remove time periods from national trajectories before present year
            wtrajs <- wtrajs[-which(wyears < presenty),]
        if(presenty < min(wyears)) { # add observed data to national trajectories if present year is not there
            adde0 <- we0[-length(we0)][we0obsy >= presenty]
            adddata <- matrix(adde0, ncol=nr.traj, nrow=length(adde0), 
                        dimnames=list(names(adde0), colnames(wtrajs)))
            wtrajs <- rbind(adddata, wtrajs)
        }
        if(!presenty %in% rownames(meta$e0.matrix)) { # add NAs to e0 matrix
            addyears <- sort(seq(presenty, max(as.integer(rownames(meta$e0.matrix)) + year.step), by=-year.step))
            adddata <- matrix(NA, nrow=length(addyears), ncol=ncol(meta$e0.matrix),
                        dimnames=list(addyears, colnames(meta$e0.matrix)))
            for(e0name in c('e0.matrix', 'e0.matrix.all'))
                meta[[e0name]] <- rbind(meta[[e0name]], adddata)
        }
        # save meta to disk
        if(file.exists(this.output.dir)) unlink(this.output.dir, recursive=TRUE)
        dir.create(outdir, recursive=TRUE)
        bayesLife.mcmc.meta <- meta
        store.bayesLife.meta.object(bayesLife.mcmc.meta, this.output.dir)
    
        this.nr.project <- nrow(wtrajs) - 1
        nr.reg <- get.nr.countries(meta)
        PIs_cqp <- array(NA, c(nr.reg, length(quantiles.to.keep), nrow(wtrajs)),
                     dimnames=list(meta$regions$country_code, dimnames(wpred$quantiles)[[2]], dimnames(wtrajs)[[1]]))
        mean_sd <- array(NA, c(nr.reg, 2, nrow(wtrajs)))
        #meta$Tc.index <- .get.Tcindex(meta$e0.matrix, cnames = meta$regions$country_name)
        country.char <- as.character(country.obj$code)
        e0reconstructed <- meta$e0.matrix
        for(region in 1:nr.reg) {
            reg.obj <- get.country.object(region, meta, index=TRUE)
            regcode.char <- as.character(reg.obj$code)			
            rege0 <- get.observed.e0(region, meta, 'e0.matrix')
            rege0.last <- rep(rege0[length(rege0)], nr.traj)
            c.first <- do.call(paste0("compute.alpha.", method), 
                               list(rege0.last, we0[names(we0) %in% names(rege0.last)])) # set of initial scales
            #stop("")
            if(is.na(rege0.last[1])) { # impute where NA's at the end
                for(i in length(rege0):1) if(!is.na(rege0[i])) break
                widx <- which(names(we0) %in% names(rege0[i]))
                c.first <- rep(do.call(paste0("compute.alpha.", method), list(rege0[i], we0[widx])), 
                               nr.traj) # set of initial scales
                meta$Tc.index[region] <- i
                imptraj <- matrix(NA, nrow = length(rege0) - i, ncol = nr.traj) # trajectory matrix for imputation
                for(tr in 1:nr.traj) { # iterate over trajectories
                    imp.time <- i:(length(rege0)-1)
                    natval <- we0[widx + (imp.time - i + 1)]
                    impval <- do.call(paste0("generate.trajectory.", method), 
                                      list(natval, c.first[tr]))
                    imptraj[imp.time - i + 1, tr] <- impval
                    # set the starting scale of the projection to the end scale of the imputation (for this trajectory)
                    c.first[tr] <- do.call(paste0("compute.alpha.", method), 
                                       list(impval[length(impval)], natval[length(natval)]))
                }
                # impute using median
                rege0[(i+1):length(rege0)] <- e0reconstructed[(i+1):length(rege0),region] <- apply(imptraj, 1, median)
                rege0.last <- imptraj[nrow(imptraj),]
            } # end of imputation
        
            e0.pred <- matrix(NA, nrow = this.nr.project + 1, ncol = nr.traj, dimnames = list(rownames(wtrajs), NULL))
            e0.pred[1,] <- rege0.last
            
            # Projections
            for(s in 1:ncol(e0.pred)) { # iterate over trajectories
                e0.pred[-1, s] <- do.call(paste0("generate.trajectory.", method), 
                                                list(wtrajs[-1,s], c.first[s]))
            }
            trajectories <- e0.pred
            save(trajectories, file = file.path(outdir, paste0('traj_country', reg.obj$code, '.rda')))
      
            # compute quantiles
            PIs_cqp[region,,] <- apply(trajectories, 1, quantile, quantiles.to.keep, na.rm = TRUE)
            mean_sd[region,1,] <- apply(trajectories, 1, mean, na.rm = TRUE)
            mean_sd[region,2,] <- apply(trajectories, 1, sd, na.rm = TRUE) 	
        }
        present.year.index <- which(rownames(meta$e0.matrix.all) == rownames(wtrajs)[1])

        bayesLife.prediction <- structure(list(
            quantiles = PIs_cqp,
            traj.mean.sd = mean_sd,
            nr.traj=nr.traj,
            e0.matrix.reconstructed = e0reconstructed,
            output.directory = normalizePath(outdir),
            mcmc.set=list(meta=meta, mcmc.list=list()),
            nr.projections=this.nr.project,
            burnin=NA, end.year=wend.year, start.year=starty,
            method = method, ar.pars = ar.pars,
            present.year.index=present.year.index,
            present.year.index.all=present.year.index,
            country = country.obj),
            class='bayesLife.prediction')
        store.bayesLife.prediction(bayesLife.prediction, outdir)
        bayesTFR:::do.convert.trajectories(pred=bayesLife.prediction, n=save.as.ascii, output.dir=outdir, verbose=verbose)
        result[[as.character(country.obj$code)]] <- bayesLife.prediction
        
        if(predict.jmale && meta$sex == 'F') {
            if(verbose) cat(', male')
            if(is.null(my.e0M.file))
                stop("Argument my.e0M.file must be given if predict.male is TRUE.")
            result[[as.character(country.obj$code)]] <- e0.jmale.predict.subnat(result[[as.character(country.obj$code)]], 
                                                                                estimates = if(is.null(jmale.estimates)) joint.male.est else jmale.estimates,
                                                                         my.e0.file = my.e0M.file, 
                                                                         ..., save.as.ascii = save.as.ascii, verbose = verbose)
        }
    }
    cat('\nPrediction stored into', output.dir, '\n')
    invisible(result)
}

e0.jmale.predict.subnat <- function(e0.pred, estimates = NULL, gap.lim = c(0,18),	
                             max.e0.eq1.pred = 86, my.e0.file = NULL, save.as.ascii = 0, verbose = TRUE) {
    # Predicting subnational male e0 from female predictions. 
    # Argument estimates should be a list of the form given by subnat.gap.estimates()
    # If my.e0.file given, it should be a male e0 file. 
    
    meta <- e0.pred$mcmc.set$meta
    if (meta$sex != 'F') stop('The prediction object must be a result of FEMALE projections.')
    if(is.null(estimates))
        estimates <- subnat.gap.estimates()
    
    e0mwpp <- get.wpp.e0.subnat.joint(e0.pred$country$code, meta, my.e0.file = my.e0.file)
    if(nrow(meta$e0.matrix) > nrow(e0mwpp$e0.matrix)) { # add NAs to e0 matrix
        adddata <- matrix(NA, nrow=nrow(meta$e0.matrix) - nrow(e0mwpp$e0.matrix), 
                          ncol=ncol(e0mwpp$e0.matrix),
                          dimnames=list(setdiff(rownames(meta$e0.matrix), rownames(e0mwpp$e0.matrix)), 
                                        colnames(e0mwpp$e0.matrix)))
        for(e0name in c('e0.matrix', 'e0.matrix.all'))
            e0mwpp[[e0name]] <- rbind(e0mwpp[[e0name]], adddata)
    }
    e0m.data <- e0mwpp$e0.matrix
    meta.changes <- list(sex='M', e0.matrix=e0m.data, e0.matrix.all=e0mwpp$e0.matrix.all, suppl.data=e0mwpp$suppl.data)
    meta.changes$Tc.index <- .get.Tcindex(meta.changes$e0.matrix, cnames=meta$regions$country_name, stop.if.less.than2=FALSE)
    prediction.file <- file.path(e0.pred$output.directory, 'prediction.rda')
    joint.male <- e0.pred
    joint.male$output.directory <- file.path(e0.pred$output.directory, 'joint_male')
    joint.male$e0.matrix.reconstructed <- e0m.data
    joint.male$fit <- estimates
    joint.male$meta.changes <- meta.changes
    joint.male$mcmc.set <- NULL
    joint.male$joint.male <- NULL
    joint.male$pred.pars <- list(gap.lim=gap.lim, max.e0.eq1.pred=max.e0.eq1.pred)
    
    if(file.exists(joint.male$output.directory)) unlink(joint.male$output.directory, recursive=TRUE)
    dir.create(joint.male$output.directory, recursive=TRUE)
    bayesLife.prediction <- .do.jmale.predict(e0.pred, joint.male, 1:get.nr.countries(meta),  
                                              gap.lim=gap.lim, eq2.age.start=max.e0.eq1.pred, verbose = FALSE,
                                              supress.warnings = TRUE)
    save(bayesLife.prediction, file=prediction.file)
    bayesTFR:::do.convert.trajectories(pred=get.e0.jmale.prediction(bayesLife.prediction), n=save.as.ascii, 
                                       output.dir=joint.male$output.directory, verbose=verbose)
    invisible(bayesLife.prediction)
}

subnat.gap.estimates <- function()
    return(list(eq1 = list(coefficients = c(-1.076940553, e0.1953 = 0.002739663, Gprev = 0.942276337, e0 = 0.019435863, e0d75 = -0.099589171),
                                  sigma = 0.3902612, dof = 4.058482),
                       eq2 = list(coefficients = c(Gprev = 1), sigma = 0.4296343, dof = NULL))
            )
