#' @export
transform_age <- function(age, adult.age=20) {
	if (age <= adult.age) {
		tAge <- log(age+1)-log(adult.age+1)
	}
	else {
		tAge <- (age-adult.age)/(adult.age+1)
	}
	return(tAge)
}
