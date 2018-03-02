coerce.multivariate <- function(x, outcome.target) {
  ## Warning:  This functon assumes that get.univariate.target has been called first, to
  ## verify the coherency of the target.  This means that the target exists in the forest object, and that
  ## it contains outcome statistics.
  ## If this is a multivarate family, we coerce the object, based on outcome.target
  ## into a univaritate regression or classification object.
  x$univariate <- TRUE
  if (x$family == "regr+" | x$family == "class+" | x$family == "mix+") {
    ##  Coerce the mulitvariate object into a univariate object.
    x.coerced <- unlist(list(x$classOutput, x$regrOutput), recursive = FALSE)[[outcome.target]]
    x$univariate <- FALSE
    x$yvar <- x$yvar[, outcome.target]
    ## Test for factors.  Ordered factors are treated as factors!
    if (is.factor(x$yvar) || is.ordered(x$yvar)) {
      x$family <- "class"
    }
    else {
      x$family <- "regr"
    }
    ## Make various assignments to the coerced object.
    x$predicted <- x.coerced$predicted
    x$predicted.oob <- x.coerced$predicted.oob
    x$class <- x.coerced$class
    x$class.oob <- x.coerced$class.oob
    x$err.rate <- x.coerced$err.rate
    x$importance <- x.coerced$importance
    x$yvar.names <- outcome.target
  }
  x$outcome.target <- outcome.target
  x
}
get.bootstrap <- function (bootstrap) {
  ## Convert bootstrap option into native code parameter.
  if (bootstrap == "by.root") {
    bootstrap <- 0
  }
    else if (bootstrap == "by.node") {
      bootstrap <- 2^19
    }
      else if (bootstrap == "none") {
        bootstrap <- 2^20
      }
        else if (bootstrap == "by.user") {
          bootstrap <- 2^19 + 2^20
        }
        else {
          stop("Invalid choice for 'bootstrap' option:  ", bootstrap)
        }
  return (bootstrap)
}
get.samptype <- function (samptype) {
  ## Convert samptype option into native code parameter.
  if (samptype == "swr") {
    bits <- 0
  }
    else if (samptype == "swor") {
      bits <- 2^12
    }
      else {
        stop("Invalid choice for 'samptype' option:  ", samptype)
      }
  return (bits)
}
get.cr.bits <- function (fmly) {
  if (fmly == "surv-CR") {
    return(2^21)
  } else {
    return(0)
  }
}
get.na.action <- function (na.action) {
  if (na.action == "na.omit") {
    ## This is the high byte!
    na.action <- 0
  }
    else if (na.action == "na.impute") {
      ## This is the high byte!
      na.action <- 2^4
      ## To recover the original functionality in which the split
      ## statistic uses missing in-node imputed values, uncomment 
      ## the following statement:
      ## na.action <- 0
    }
    else {
        stop("Invalid choice for 'na.action' option:  ", na.action)
      }
  return (na.action)
}
get.forest <- function (forest) {
  ## Convert forest option into native code parameter.
  if (!is.null(forest)) {
    if (forest == TRUE) {
      forest <- 2^5
    }
      else if (forest == FALSE) {
        forest <- 0
      }
        else {
          stop("Invalid choice for 'forest' option:  ", forest)
        }
  }
    else {
      stop("Invalid choice for 'forest' option:  ", forest)
    }
  return (forest)
}
get.forest.wt <- function (grow.equivalent, bootstrap, weight) {
    ## Convert weight option into native code parameter.
    if (!is.null(weight)) {
        if (weight == FALSE) {
            wght.bits <- 0
        }
        else if (grow.equivalent == TRUE) {
            if (weight == TRUE) {
                wght.bits <- 2^0 + 2^1
            }
            else if (weight == "inbag") {
                wght.bits <- 2^0 + 2^1
            }
            else if (weight == "oob") {
                wght.bits <- 2^0 + 2^2
            }
            else if (weight == "all") {
                wght.bits <- 2^0 + 2^1 + 2^2
            }
            else {
                stop("Invalid choice for 'weight' option:  ", weight)
            }
        }
        else if (grow.equivalent == FALSE) {
            if (weight == TRUE) {
                wght.bits <- 2^0 + 2^1 + 2^2
            }
            else if (weight == "all") {
                wght.bits <- 2^0 + 2^1 + 2^2
            }
            else {
                stop("Invalid choice for 'weight' option:  ", weight)
            }
        }
        else {
            stop("Invalid choice for 'grow.equivalent' in weight:  ", grow.equivalent)
        }
    }
    else {
        stop("Invalid choice for 'weight' option:  ", weight)
    }
    return (wght.bits)
}
get.importance <-  function (importance) {
  ## Convert importance option into native code parameter.
  if (!is.null(importance)) {
    ## Override lazy values.
    if (importance == TRUE) {
      importance <- "permute"
    }
      else if (importance == FALSE) {
        importance <- "none"
      }
    if (importance == "none") {
      importance <- 0
    }
      else if (importance == "anti.ensemble") {
        importance <- 2^25 + 0
      }
        else if (importance == "permute.ensemble") {
          importance <- 2^25 + 2^8
        }
          else if (importance == "random.ensemble") {
            importance <- 2^25 + 2^9
          }
            else if (importance == "anti.joint.ensemble") {
              importance <- 2^25 + 2^10 + 0
            }
              else if (importance == "permute.joint.ensemble") {
                importance <- 2^25 + 2^10 + 2^8
              }
                else if (importance == "random.joint.ensemble") {
                  importance <- 2^25 + 2^10 + 2^9
                }
                  else if (importance == "anti") {
                    importance <- 2^25 + 2^24 + 0
                  }
                    else if (importance == "permute") {
                      importance <- 2^25 + 2^24 + 2^8
                    }
                      else if (importance == "random") {
                        importance <- 2^25 + 2^24 + 2^9
                      }
                        else if (importance == "anti.joint") {
                          importance <- 2^25 + 2^24 + 2^10 + 0
                        }
                          else if (importance == "permute.joint") {
                            importance <- 2^25 + 2^24 + 2^10 + 2^8
                          }
                            else if (importance == "random.joint") {
                              importance <- 2^25 + 2^24 + 2^10 + 2^9
                            }
                              else {
                                stop("Invalid choice for 'importance' option:  ", importance)
                              }
  }
    else {
      stop("Invalid choice for 'importance' option:  ", importance)
    }
  return (importance)
}
get.impute.only <-  function (impute.only, nMiss) {
  if (impute.only) {
    if (nMiss > 0) {
      return (2^16)
    }
      else {
        stop("Data has no missing values, using 'impute' makes no sense.")
      }
  }
    else {
      return (0)
    }
}
get.outcome <- function (outcome) {
  ## Convert outcome option into native code parameter.
  if (outcome == "train") {
    outcome <- 0
  }
    else if (outcome == "test") {
      outcome <- 2^17
    }
      else {
        stop("Invalid choice for 'outcome' option:  ", outcome)
      }
  return (outcome)
}
get.perf <-  function (perf, impute.only, family) {
  ## first deal with impute.only where there is no performance
  if (impute.only == TRUE) {
    return("none")
  }
  ## now deal with non-classification
  if (family != "class") {
    if (is.null(perf)) {
      return("default")
    }
    perf <- match.arg(perf, c("none", "default", "standard"))##only allowed values
    if (perf == "standard") {
      perf <- "default"
    }
    return(perf)
  }
  ## now deal with classification
  if (is.null(perf)) {
    return("default")
  }
  perf <- match.arg(perf, c("none", "default", "standard", "misclass", "brier", "g.mean"))##only allowed values
  if (perf == "standard" || perf == "misclass") {
    perf <- "default"
  }
  perf
}
get.perf.bits <- function (perf) {
  if (perf == "default") {
    return (2^2)
  }
  else if (perf == "g.mean") {
    return (2^2 + 2^14)
  }
  else if (perf == "brier") {
    return (2^2 + 2^3)
  }
  else {#everything else becomes "none"
    return (0)
  }
}
get.rfq <- function(rfq) {
  if (is.null(rfq)) {
    rfq <- FALSE
  }
  rfq
}
get.rfq.bits <- function (rfq, family) {
    result <- 0
    if (family == "class") {
        if (rfq) {
            result <- 2^15
        }
    }
    return (result)
}
get.proximity <- function (grow.equivalent, proximity) {
  ## Convert proximity option into native code parameter.
    if (!is.null(proximity)) {
      if (proximity == FALSE) {
        prox.bits <- 0
      }
        else if (grow.equivalent == TRUE) {
          if (proximity == TRUE) {
            prox.bits <- 2^28 + 2^29
          }
            else if (proximity == "inbag") {
              prox.bits <- 2^28 + 2^29
            }
              else if (proximity == "oob") {
                prox.bits <- 2^28 + 2^30
              }
                else if (proximity == "all") {
                  prox.bits <- 2^28 + 2^29 + 2^30
                }
                  else {
                    stop("Invalid choice for 'proximity' option:  ", proximity)
                  }
        }
          else if (grow.equivalent == FALSE) {
            if (proximity == TRUE) {
              prox.bits <- 2^28 + 2^29 + 2^30
            }
              else if (proximity == "all") {
                prox.bits <- 2^28 + 2^29 + 2^30
              }
                else {
                  stop("Invalid choice for 'proximity' option:  ", proximity)
                }
          }
            else {
              stop("Invalid choice for 'grow.equivalent' in proximity:  ", grow.equivalent)
            }
    }
      else {
        stop("Invalid choice for 'proximity' option:  ", proximity)
      }
    return (prox.bits)
  }
 
  get.membership <- function (membership) {
    ## Convert option into native code parameter.
    bits <- 0
    if (!is.null(membership)) {
      if (membership == TRUE) {
        bits <- 2^6
      }
        else if (membership != FALSE) {
          stop("Invalid choice for 'membership' option:  ", membership)
        }
    }
      else {
        stop("Invalid choice for 'membership' option:  ", membership)
      }
    return (bits)
  }
  get.rf.cores <- function () {
    if (is.null(getOption("rf.cores"))) {
      if(!is.na(as.numeric(Sys.getenv("RF_CORES")))) {
        options(rf.cores = as.integer(Sys.getenv("RF_CORES")))
      }
    }
    return (getOption("rf.cores", -1L))
  }
  get.seed <- function (seed) {
    if ((is.null(seed)) || (abs(seed) < 1)) {
      seed <- runif(1,1,1e6)
    }
    seed <- -round(abs(seed))
    return (seed)
  }
  get.split.depth <- function (split.depth) {
    ## Convert split.depth option into native code parameter.
    if (!is.null(split.depth)) {
      if (split.depth == "all.trees") {
        split.depth <- 2^22
      }
        else if (split.depth == "by.tree") {
          split.depth <- 2^23
        }
          else if (split.depth == FALSE) {
            split.depth <- 0
          }
            else {
              stop("Invalid choice for 'split.depth' option:  ", split.depth)
            }
    }
      else {
        stop("Invalid choice for 'split.depth' option:  ", split.depth)
      }
    return (split.depth)
  }
  get.split.null <- function (split.null) {
    ## Convert split.null option into native code parameter.
    if (!is.null(split.null)) {
      if (split.null == TRUE) {
        split.null <- 2^18
      }
      else if (split.null == FALSE) {
        split.null <- 0
      }
      else {
        stop("Invalid choice for 'split.null' option:  ", split.null)
      }
    }
    else {
      stop("Invalid choice for 'split.null' option:  ", split.null)
    }
    return (split.null)
  }
  get.split.cust <- function (split.cust) {
    ## Convert split.cust option into native code parameter.
    if (!is.null(split.cust)) {
      if ((split.cust >= 1) && (split.cust <= 16)) {
        ## Bit shift eight left.
        split.cust <- 256 * (split.cust - 1)
      }
        else {
          stop("Invalid choice for 'split.cust' option:  ", split.cust)
        }
    }
      else {
        split.cust <- 0
      }
    return (split.cust)
  }
  get.statistics <- function (statistics) {
    ## Convert statistics option into native code parameter.
    if (!is.null(statistics)) {
      if (statistics == TRUE) {
        statistics <- 2^27
      }
        else if (statistics == FALSE) {
          statistics <- 0
        }
          else {
            stop("Invalid choice for 'statistics' option:  ", statistics)
          }
    }
      else {
        stop("Invalid choice for 'statistics' option:  ", statistics)
      }
    return (statistics)
  }
  get.terminal.qualts <- function (terminal.qualts, incoming.flag) {
    ## Convert option into native code parameter.  This 
    ## is sensitive to incoming and outgoing data 
    ## (from the native code perspective).
    bits <- 0
    if (is.null(incoming.flag)) {
      ## Do nothing.  This ensures backwards compatibility with
      ## versions prior to these bits being flagged in the grow object.
    }
      else if (incoming.flag) {
        bits <- bits + 2^17
      }
    if (!is.null(terminal.qualts)) {
      if (terminal.qualts == TRUE) {
        bits <- bits + 2^16
      }
        else if (terminal.qualts != FALSE) {
          stop("Invalid choice for 'terminal.qualts' option:  ", terminal.qualts)
        }
    }
      else {
        stop("Invalid choice for 'terminal.qualts' option:  ", terminal.qualts)
      }
    return (bits)
  }
  get.terminal.quants <- function (terminal.quants, incoming.flag) {
    ## Convert option into native code parameter.  This 
    ## is sensitive to incoming and outgoing data 
    ## (from the native code perspective).
    bits <- 0
    if (is.null(incoming.flag)) {
      ## Do nothing.  This ensures backwards compatibility with
      ## versions prior to these bits being flagged in the grow object.
    }
      else if (incoming.flag) {
        bits <- bits + 2^19
      }
    if (!is.null(terminal.quants)) {
      if (terminal.quants == TRUE) {
        bits <- bits + 2^18
      }
        else if (terminal.quants != FALSE) {
          stop("Invalid choice for 'terminal.quants' option:  ", terminal.quants)
        }
    }
      else {
        stop("Invalid choice for 'terminal.quants' option:  ", terminal.quants)
      }
    return (bits)
  }
  get.trace <- function (do.trace) {
    ## Convert trace into native code parameter.
    if (!is.logical(do.trace)) {
      if (do.trace >= 1) {
        do.trace <- round(do.trace)
      }
        else {
          do.trace <- 0
        }
    }
      else {
        do.trace <- 1 * do.trace
      }
    return (do.trace)
  }
  get.tree.err <- function (tree.err) {
    ## Convert tree.err option into native code parameter.
    if (!is.null(tree.err)) {
      if (tree.err == FALSE) {
          tree.err <- 0
      }
        else if (tree.err == TRUE) {
            tree.err <- 2^13
        }
          else {
            stop("Invalid choice for 'tree.err' option:  ", tree.err)
          }
    }
      else {
        stop("Invalid choice for 'tree.err' option:  ", tree.err)
      }
    return (tree.err)
  }
  get.var.used <- function (var.used) {
    ## Convert var.used option into native code parameter.
    if (!is.null(var.used)) {
      if (var.used == "all.trees") {
        var.used <- 2^12
      }
        else if (var.used == "by.tree") {
          var.used <- 2^13
        }
          else if (var.used == FALSE) {
            var.used <- 0
          }
            else {
              stop("Invalid choice for 'var.used' option:  ", var.used)
            }
    }
      else {
        stop("Invalid choice for 'var.used' option:  ", var.used)
      }
    return (var.used)
  }
  get.vimp.only <-  function (vimp.only) {
    ## Convert vimp.only option into native code parameter.
    if (!is.null(vimp.only)) {
      if (vimp.only) {
        return (2^27)
      }
        else if (!vimp.only) {
          return (0)
        }
          else {
            stop("Invalid choice for 'vimp.only' option:  ", vimp.only)
          }
    }
      else {
        stop("Invalid choice for 'vimp.only' option:  ", vimp.only)
      }
  }
  ## HIDDEN VARIABLES FOLLOW:
  is.hidden.impute.only <-  function (user.option) {
    if (is.null(user.option$impute.only)) {
      FALSE
    }
      else {
        as.logical(as.character(user.option$impute.only))
      }
  }
  is.hidden.terminal.qualts <-  function (user.option) {
    ## Default value is !FALSE
    if (is.null(user.option$terminal.qualts)) {
      !FALSE
    }
      else {
        as.logical(as.character(user.option$terminal.qualts))
      }
  }
  is.hidden.terminal.quants <-  function (user.option) {
    ## Default value is FALSE
    if (is.null(user.option$terminal.quants)) {
      FALSE
    }
      else {
        as.logical(as.character(user.option$terminal.quants))
      }
  }
  is.hidden.perf.type <-  function (user.option) {
    ## Default value is NULL
    if (is.null(user.option$perf.type)) {
      NULL
    }
      else {
        as.character(user.option$perf.type)
      }
  }
  is.hidden.rfq <-  function (user.option) {
    if (is.null(user.option$rfq)) {
      NULL
    }
      else {
        as.logical(as.character(user.option$rfq))
      }
  }
is.hidden.ytry <-  function (user.option) {
    if (is.null(user.option$ytry)) {
        NULL
    }
    else {
        as.integer(user.option$ytry)
    }
}
is.hidden.htry <-  function (user.option) {
    if (is.null(user.option$htry)) {
        htry = 0
    }
    else {
        htry = as.integer(user.option$htry)
        if (htry < 0) {
            stop("Invalid choice for 'htry' option:  ", user.option$htry)
        }
    }
    return (htry)
}
