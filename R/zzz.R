.onLoad <- function(libname, pkgname) {
    # Register lavaan method if package is available
    if (requireNamespace("lavaan", quietly = TRUE)) {
        # Get the S4 class definition safely from the package namespace
        # This avoids issues if the class is not in the search path
        lavaan_class <- methods::getClass("lavaan", where = asNamespace("lavaan"))

        # Register the method
        # We use the internal function extract_mediation_lavaan defined in R/methods-extract-lavaan.R
        S7::method(extract_mediation, lavaan_class) <- extract_mediation_lavaan
    }
}
