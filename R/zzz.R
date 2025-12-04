.onLoad <- function(libname, pkgname) {
  # Register S7 methods for proper dispatch
  S7::methods_register()
}
