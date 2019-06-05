#' @importFrom stats optim
#
# Repeats call to optim if maximum iterations is reached prior to convergence.
# Default maximum number of attempts is 10
#
# @param argList list; input values for function optim
# @param maxit integer; maximum number of attempts to reach convergence
#
# @return NULL if optim() was not succesful; list of optim results otherwise
#
.myOptim <- function(argList, maxit = 10L) {

  wll <- list()

  # ensures that at least one iteration is completed
  maxit <- max(maxit, 1L)

  for (i in 1L:maxit) {

    # call optim()
    tmp <- do.call(what = optim, args = argList)

    if (tmp$convergence == 1L) {
      # if maximum iterations hit, reset par input to current value
      # store result and attempt again if not at maxit limit
      argList[[ "par" ]] <- tmp$par
      wll[[ i ]] <- tmp
      if (length(x = wll) == maxit) {
        message(paste("optim() did not converge after", maxit, "iterations"))
        return( NULL )
      }
    } else if (tmp$convergence == 0L) {
      # if convergence attained, store results and break out of loop
      wll[[ i ]] <- tmp
      break
    } else if (tmp$convergence > 1L & !is.null(x = tmp$message)) {
      # if convergence flag is failure of another type, print message and
      # return NULL
      message(paste("optim() failed with flag", tmp$convergence, "and message\n",
                    tmp$message))
      return( NULL )
    } else if (tmp$convergence > 1L) {
      # if convergence flag is failure of another type without a message;
      # print failure and return NULL
      message("optim() failed with flag", tmp$convergence, "\n")
      return( NULL )
    }

  }

  return( wll )

}
