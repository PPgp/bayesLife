.onLoad <- function (lib, pkg) {
    library.dynam("bayesLife", pkg, lib)
}

.onUnload <- function (libpath) {
  library.dynam.unload("bayesLife", libpath)
}

.onAttach <- function(lib, pkg)
{
    # unlock .e0options variable allowing its modification
    unlockBinding(".e0options", asNamespace("bayesLife")) 
    invisible()
}
