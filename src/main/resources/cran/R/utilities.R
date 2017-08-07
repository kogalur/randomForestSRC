get.bootstrap <- function (bootstrap) {
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
    na.action <- 0
  }
    else if (na.action == "na.impute") {
      na.action <- 2^4
    }
    else {
        stop("Invalid choice for 'na.action' option:  ", na.action)
      }
  return (na.action)
}
get.forest <- function (forest) {
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
  if (!is.null(weight)) {
    if (weight == FALSE) {
      weight <- 0
    }
      else if (grow.equivalent == TRUE) {
        if ((bootstrap != "by.root") && (bootstrap != "by.user")) {
          weight <- 2^0 + 2^2
        }
          else {
            if (weight == TRUE) {
              weight <- 2^0
            }
              else if (weight == "inbag") {
                weight <- 2^0
              }
                else if (weight == "oob") {
                  weight <- 2^0 + 2^1
                }
                  else if (weight == "all") {
                    weight <- 2^0 + 2^2
                  }
                    else {
                      stop("Invalid choice for 'weight' option:  ", weight)
                    }
          }
      }
        else if (grow.equivalent == FALSE) {
            if (weight == TRUE) {
                weight <- 2^0 + 2^2
            }
            else if (weight == "all") {
                weight <- 2^0 + 2^2
            }
            else {
                stop("Invalid choice for 'weight' option:  ", weight)
            }
        }
        else {
            stop("Invalid choice for 'weight' option:  ", weight)
        }
    }
  else {
      stop("Invalid choice for 'weight' option:  ", weight)
  }
  return (weight)
}
get.importance <-  function (importance) {
  if (!is.null(importance)) {
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
get.perf <-  function (perf, impute.only, family, perf.type) {
    if (impute.only == TRUE) {
        result = "none"
    }
    else {
        if (is.null(perf)) {
            result = "default"
        }
        else {
            result = perf
        }
    }
    if (result == "default") {
        if (family == "class") {
            if (!is.null(perf.type)) {
                if (perf.type == "g.mean") {
                    result = "g.mean"
                }
                if (perf.type == "g.mean.rfq") {
                    result = "g.mean.rfq"
                }
                if (perf.type == "brier") {
                    result = "brier"
                }
            }
        }
    }
    else {
        result = "none"
    }
    return (result)
}
get.perf.bits <- function (perf) {
    if (perf == "default") {
        return (2^2)
    }
    else if (perf == "g.mean") {
        return (2^2 + 2^14)
    }
    else if (perf == "g.mean.rfq") {
        return (2^2 + 2^15)
    }
    else if (perf == "brier") {
        return (2^2 + 2^3)
    }
    else {
        return (0)
    }
}
get.proximity <- function (grow.equivalent, proximity) {
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
    if (!is.null(split.cust)) {
      if ((split.cust >= 1) && (split.cust <= 16)) {
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
  get.trace <- function (do.trace) {
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
  get.var.used <- function (var.used) {
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
  get.membership <- function (membership) {
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
  get.terminal.qualts <- function (terminal.qualts, incoming.flag) {
    bits <- 0
    if (is.null(incoming.flag)) {
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
    bits <- 0
    if (is.null(incoming.flag)) {
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
  get.tree.err <- function (tree.err) {
    if (!is.null(tree.err)) {
      if (tree.err == FALSE) {
        tree.err <- 2^13
      }
        else if (tree.err == TRUE) {
          tree.err <- 0
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
  is.hidden.impute.only <-  function (user.option) {
    if (is.null(user.option$impute.only)) {
      FALSE
    }
      else {
        as.logical(as.character(user.option$impute.only))
      }
  }
  is.hidden.terminal.qualts <-  function (user.option) {
    if (is.null(user.option$terminal.qualts)) {
      !FALSE
    }
      else {
        as.logical(as.character(user.option$terminal.qualts))
      }
  }
  is.hidden.terminal.quants <-  function (user.option) {
    if (is.null(user.option$terminal.quants)) {
      FALSE
    }
      else {
        as.logical(as.character(user.option$terminal.quants))
      }
  }
  is.hidden.perf.type <-  function (user.option) {
    if (is.null(user.option$perf.type)) {
      NULL
    }
      else {
        as.character(user.option$perf.type)
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
  get.univariate.target <- function(x, outcome.target = NULL) {
    if (x$family == "regr+" | x$family == "class+" | x$family == "mix+") {
      if (is.null(outcome.target)) {
        target <- match(c("regrOutput", "classOutput"), names(x))
        target <- target[!is.na(target)]
        if(length(target) > 0) {
          do.break <- FALSE
          for (i in target) {
            for (j in 1:length(x[[i]])) {
              if (length(x[[i]][[j]]) > 0) {
                outcome.target <- names(x[[i]][j])
                do.break <- TRUE
                break
              }
            }
            if (do.break == TRUE) {
              break
            }
          }
        }
          else {
            stop("No outcomes found in object.  Please contact technical support.")
          }
      }
      else {
        if (sum(is.element(outcome.target, x$yvar.names)) != 1) {
          stop("User must specify one and only one outcome.target for multivariate families.")
        }
        target <- match(c("regrOutput", "classOutput"), names(x))
        target <- target[!is.na(target)]
        found = FALSE
        if(length(target) > 0) {
          do.break <- FALSE
          for (i in target) {
            for (j in 1:length(x[[i]])) {
              if (length(x[[i]][[j]]) > 0) {
                if (outcome.target == names(x[[i]][j])) {
                  found = TRUE
                  do.break <- TRUE
                  break
                }
              }
            }
            if (do.break == TRUE) {
              break
            }
          }     
        }
        if (!found) {
          stop("Target outcome not found in object.  Re-run analysis with target outcome specified.")
        }
      }
    }
    outcome.target
  }
  coerce.multivariate <- function(x, outcome.target) {
    x$univariate <- TRUE
    if (x$family == "regr+" | x$family == "class+" | x$family == "mix+") {
      x.coerced <- unlist(list(x$classOutput, x$regrOutput), recursive = FALSE)[[outcome.target]]
      x$univariate <- FALSE
      x$yvar <- x$yvar[, outcome.target]
      if (is.factor(x$yvar) || is.ordered(x$yvar)) {
        x$family <- "class"
      }
      else {
        x$family <- "regr"
      }
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
