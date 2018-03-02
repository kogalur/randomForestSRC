predict.rfsrc <-
  function(object,
           newdata,
           m.target=NULL,
           importance = c(FALSE, TRUE, "none", "permute", "random", "anti", "permute.ensemble", "random.ensemble", "anti.ensemble")[1],
           na.action = c("na.omit", "na.impute"),
           outcome = c("train", "test"),
           proximity = FALSE,
           forest.wt = FALSE,
           ptn.count = 0,
            
           var.used = c(FALSE, "all.trees", "by.tree"),
           split.depth = c(FALSE, "all.trees", "by.tree"),
           seed = NULL,
           do.trace = FALSE,
           membership = FALSE,
           statistics = FALSE,
           ...)
{
  result.predict <- generic.predict.rfsrc(object,
                                          newdata,
                                          m.target = m.target,
                                          importance = importance,
                                          na.action = na.action,
                                          outcome = outcome,
                                          proximity = proximity,
                                          forest.wt = forest.wt,                                          
                                          ptn.count = ptn.count,
                                           
                                          var.used = var.used,
                                          split.depth = split.depth,
                                          seed = seed,
                                          do.trace = do.trace,
                                          membership = membership,
                                          statistics = statistics,
                                          ...)
  return(result.predict)
}
