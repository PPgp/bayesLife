.First.lib <- function (lib, pkg) {
    library.dynam("bayesLife", pkg, lib)
}

.Last.lib <- function (libpath) {
  library.dynam.unload("bayesLife", libpath)
}
