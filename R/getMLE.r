#' Calculate the weighted MLE
#'
#' This function calculate the weighted MLE for the input covariate matrix x, response vector y, and weight vector w.
#' It returns a list with three elements: par, the weighted MLE; msg, the fitting message; iter, the number of itterations used.

#' @param x the input covariate matrix
#' @param y the input response vector
#' @param w the wight vector
#' @keywords getMLE
#' @export
#' @examples
#' library(OSMAC)
#' dat <- adult.train
#' X <- as.matrix(dat[,c(1,3,5,12:13)])
#' X <- t(t(X) / apply(X, 2, sd))
#' X <- cbind(1, X)
#' Y <- as.numeric(dat[,15]) - 1
#' n <- dim(X)[1]
#' d <- dim(X)[2]
#' getMLE(X, Y, 1)
#' $par
#'                           Age         Fnlwgt  Education.num   Loss.captial
#'     -8.6366072      0.6374174      0.0648296      0.8780786      0.2342951
#' Hours.per.week
#'      0.5249214
#'
#' $message
#' [1] "Successful convergence"
#'
#' $iter
#' [1] 5

getMLE <- function(x, y, w) {
    beta <- rep(0, d)
    loop  <- 1
    Loop  <- 100
    msg <- "NA"
    while (loop <= Loop) {
        pr <- c(1 - 1 / (1 + exp(x %*% beta)))
        H <- t(x) %*% (pr * (1 - pr) * w * x)
        S <- colSums((y - pr) * w * x)
        tryCatch(
            {shs <- NA
             shs <- solve(H, S) },
            error=function(e){
                cat("\n ERROR :", loop, conditionMessage(e), "\n")})
        if (is.na(shs[1])) {
            msg <- "Not converge"
            beta <- loop <- NA
            break
        }
        beta.new <- beta + shs
        tlr  <- sum((beta.new - beta)^2)
        beta  <- beta.new
        if(tlr < 0.000001) {
            msg <- "Successful convergence"
            break
        }
        if (loop == Loop)
            warning("Maximum iteration reached")
        loop  <- loop + 1
    }
    list(par=beta, message=msg, iter=loop)
}

