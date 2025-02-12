
test_clr = function(){

	library(crumblr)
	library(RUnit)
	library(compositions)

	# set probability of each category
	prob <- c(0.1, 0.2, 0.3, 0.5)

	# number of total counts
	countsTotal <- 300

	# number of samples
	n_samples <- 100

	# simulate info for each sample
	info <- data.frame(Age = rgamma(n_samples, 50, 1))
	rownames(info) <- paste0("sample_", 1:n_samples)

	# simulate counts from multinomial
	counts <- t(rmultinom(n_samples, size = countsTotal, prob = prob))
	colnames(counts) <- paste0("cat_", 1:length(prob))
	rownames(counts) <- paste0("sample_", 1:n_samples)


	# check CLR
	r1 = compositions::clr(counts + 0.5)
	r2 = crumblr::clr(counts)
	checkEqualsNumeric(r1, r2)

	# check clrInv	 
	i1 = compositions::clrInv(r2)
	i2 = crumblr::clrInv(r2)
	checkEqualsNumeric(i1, i2)

	# check crumblr	 
	r1 = compositions::clr(counts + 0.5)
	cobj <- crumblr(counts)
	checkEqualsNumeric(t(cobj$E), r1)


}