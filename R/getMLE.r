#'This function calculate the weighted MLE for the input 
#'   covariate matrix x, 
#'   response vector y, 
#'   and weight vector w.
#'This function returns a list with three elements: 
#'   par, the weighted MLE; 
#'   msg, the fitting message; 
#'   iter, the number of itterations used. 
#'
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

