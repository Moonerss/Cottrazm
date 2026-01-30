check_pkgs <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = T)) {
      stop(sprintf('This function requires the %s package.', pkg))
    }
  }
}
