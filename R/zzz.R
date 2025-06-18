.onAttach <- function(libname, pkgname) {
  msg <- paste0(
    "Attaching package: ", pkgname, "\n\n",
    "New features available. Type\n",
    " - BTSR.models() to see the list of available models.\n",
    " - help('BTSR-package') for a full description of the\n",
    "   models structure.\n"
  )
  packageStartupMessage(msg)
}
