#'@title Identify the outliers
#'@description Identify the outliers
#'@param x  Vector
#'@param thr significant level
#'@export
#'@aliases
grubbs.flag <- function(x,thr) {
  require(outliers)
  outliers <- NULL
  test <- x
  grubbs.result <- grubbs.test(test)
  pv <- grubbs.result$p.value
  # throw an error if there are too few values for the Grubb's test
  if (length(test) < 3 ) stop("Grubb's test requires > 2 input values")
  while(pv < thr) {
    outliers <- c(outliers,as.numeric(strsplit(grubbs.result$alternative," ")[[1]][3]))
    test <- x[!x %in% outliers]
    # stop if all but two values are flagged as outliers
    if (length(test) < 3 ) {
      warning("All but two values flagged as outliers")
      break
    }
    grubbs.result <- grubbs.test(test)
    pv <- grubbs.result$p.value
  }
  return(data.frame(X=x,Outlier=(x %in% outliers)))
}
