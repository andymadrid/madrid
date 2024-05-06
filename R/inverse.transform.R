#' @export
inverse.transform <- function(tAge, adult.age=20) {
	if (tAge < 0) {
		iAge <- (exp(tAge + log(adult.age + 1)))-1
	}
	else {
		iAge <- (tAge*(adult.age+1))+(adult.age)
	}
	return(iAge)
}
