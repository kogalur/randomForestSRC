get.base.learner <- function(interact.depth = 2,
                             interact.rule = c("multiplication", "division", "addition", "subtraction"),
                             synthetic.depth = 0,
                             dim.reduce = TRUE) {
    ## No verificiation is done here, or before conveyance to the
    ## native code.  So be careful.
    ## When augmentation is of the interaction type, interact.depth is
    ## the number of sub-terminal nodes grown to
    ## determine interactions.  Thus, a value of 2 or greater is
    ## necessary for interactions. The number of interactions considered is interact.depth - 1.
    ## This does not guarantee
    ## that there will be interactions.  It is just a pre-requisite.
    ## When augmentation is of the synthetic type, synthetic.depth is the
    ## number of sub-terminal nodes in the sub-tree that is grown to
    ## determine the synthetic x-variable.  Thus, a value of 2 or
    ## great is necessary for synthetic interactions.
  ## Don't mess with the order here.  The matching of position is relevant to the C-side code.
  valid.rule <- c("multiplication", "division", "addition", "subtraction")
  rule <- match.arg(interact.rule, valid.rule)
  base.learner <- list(as.integer(interact.depth), as.integer(which(valid.rule == rule)), as.integer(synthetic.depth), as.integer(dim.reduce))
    names(base.learner) = c("interact.depth", "learner.rule", "synthetic.depth", "dim.reduce")
    class(base.learner) = "base.learner"
  return (base.learner)
}
get.lot <- function(hdim = 5, treesize = function(x){min(50, x * .25)}, lag = 8, strikeout = 3) {
    ## The size of tree can be specified as a function or an integer.
    ## If a function is specified, it MUST be processed downstream and
    ## converted to an integer before passing the lot object into the
    ## native code.
    if (!is.function(treesize) && !is.numeric(treesize)) {
        stop("treesize must be a function or number specifying size of tree")
    }
    else {
        if (is.function(treesize)) {
            lot = list(as.integer(hdim), treesize, as.integer(lag), as.integer(strikeout))
        }
        if (is.numeric(treesize)) {
            lot = list(as.integer(hdim), as.integer(treesize), as.integer(lag), as.integer(strikeout))
        }
    }
    names(lot) = c("hdim", "treesize", "lag", "strikeout")
    class(lot) = "lot"
    return (lot)
}
