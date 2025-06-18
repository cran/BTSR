#***********************************************************************
# internal function.
# Not supposed to be called by the user.
# Format a message to be passed to stop, warning, cat, etc
#***********************************************************************
# Version 1.0.0:
#  - added this function
#***********************************************************************
.format.message <- function(msg) {
  # Split message into lines
  lines <- strsplit(msg, "\n")[[1]]
  # Get the longest line's length
  max_length <- max(nchar(lines))
  # Create dashes
  line <- paste(rep("-", max_length), collapse = "")
  # Ensure msg prints correctly
  formatted_msg <- paste(lines, collapse = "\n")
  paste0("\n", line, "\n", formatted_msg, "\n", line, "\n")
}

#***********************************************************************
# internal function.
# Not supposed to be called by the user.
# Issues a message and stop the code
#***********************************************************************
# Version 1.0.0:
#  - added this function
#***********************************************************************
.stop.with.message <- function(msg) {
  stop(.format.message(msg), call. = FALSE)
}

#***********************************************************************
# internal function.
# Not supposed to be called by the user.
# Issues a warning message. Do not stop the code
#***********************************************************************
# Version 1.0.0:
#  - added this function
#***********************************************************************
.warning.with.message <- function(msg) {
  warning(.format.message(msg), immediate. = TRUE, call. = FALSE)
}

#***********************************************************************
# internal function.
# Not supposed to be called by the user.
# Error messages from FORTRAN
#***********************************************************************
# Version 1.0.0:
#  - added this function
#***********************************************************************
.code.to.message <- function(code, method = NULL) {
  # check if is the output from a fitting algorithm
  # check if the method used is Nelder-Mead
  is.algo <- !is.null(method)
  is.nelder <- ifelse(is.algo, startsWith(tolower(method), "n"), FALSE)
  algo <- is.algo && is.nelder

  switch(as.character(code),
    "0" = {
      msg <- "SUCCESSFUL TERMINATION"
    },
    "1" = {
      msg <- ifelse(algo,
        "MAXIMUM NO. OF FUNCTION EVALUATIONS EXCEEDED",
        "Unknown error.\n  - Set `iprint = 1` in the `control` list to print details"
      )
      code2 <- "fail"
    },
    "2" = {
      msg <- ifelse(algo,
        "NOP < 1 OR STOPCR <= 0",
        ""
      )
      code2 <- "fail"
    },
    "11" = {
      msg <- "Fail to evaluate g11(y): non-finite value computed"
      code2 <- "link"
    },
    "12" = {
      msg <- "Fail to evaluate g12(y): non-finite value computed"
      code2 <- "link"
    },
    "1112" = {
      msg <- "Fail to evaluate g11(y) and g12(y): non-finite value computed"
      code2 <- "link"
    },
    "11" = {
      msg <- "mu smaller than the lower limit."
      code2 <- "linkpar"
    },
    "12" = {
      msg <- "nu smaller than the lower limit."
      code2 <- "linkpar"
    },
    "21" = {
      msg <- "mu bigger than the upper limit."
      code2 <- "linkpar"
    },
    "22" = {
      msg <- "nu bigger than the upper limit."
      code2 <- "linkpar"
    },
    "123" = {
      msg <- "yt out-of-bounds."
      code2 <- "linkpar"
    },
    "91" = {
      msg <- "Fail to compute g12(y.start): y.start is out of bounds"
      code2 <- "start"
    },
    "92" = {
      msg <- "Fail to compute g22(vt.start): vt.start is out of bounds"
      code2 <- "start"
    },
    "9991" = {
      msg <- "Fail to compute g12(y.start): non-finite value"
      code2 <- "startlink"
    },
    "9992" = {
      msg <- "Fail to compute g22(vt.start): non-finite value"
      code2 <- "startlink"
    }
  )

  if (code != 0) {
    msg2 <- switch(code2,
      "link" = " - Try changing the link function",
      "linkpar" = " - Try changing the link function or some parameter values",
      "start" = " - Please, select a different starting value",
      "startlink" = " - Try, changing the link function or the starting value",
      "fail" = "FAIL / FUNCTION DID NOT CONVERGE!"
    )

    .warning.with.message(
      paste0(
        if (code2 != "fail") "Revision Required. ",
        if (is.algo) paste0(ifelse(is.nelder, "Nelder-Mead", "L-BFGS-B")),
        "\n", msg, "\n", msg2,
        if (code2 != "fail") "\nReturning FORTRAN output without further processing"
      )
    )
  }
  invisible(msg)
}
