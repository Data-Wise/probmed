.onLoad <- function(libname, pkgname) {
  # Register S7 classes with S4 before method registration (needed for the
  # S7 method on the base `print` generic for GaugePmedResult).
  S7::S4_register(GaugePmedResult)
  # Register S7 methods for proper dispatch
  S7::methods_register()
}
