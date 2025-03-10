
generate_data = function(){
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

	list(counts = counts, info = info)
}


test_clr = function(){

	library(crumblr)
	library(RUnit)
	library(compositions)

	data = generate_data()

	# check CLR
	r1 = compositions::clr(data$counts + 0.5)
	r2 = crumblr::clr(data$counts)
	checkEqualsNumeric(r1, r2)
}


test_clrInv = function(){

	library(crumblr)
	library(RUnit)
	library(compositions)

	data = generate_data()
	
	# check clrInv	 
	r1 = compositions::clr(data$counts + 0.5)
	r2 = crumblr::clr(data$counts)

	i1 = compositions::clrInv(r2)
	i2 = crumblr::clrInv(r2)
	checkEqualsNumeric(i1, i2)
}


test_clrInv = function(){

	library(crumblr)
	library(RUnit)
	library(compositions)

	data = generate_data()
	
	# check crumblr	 
	r1 = compositions::clr(data$counts + 0.5)
	cobj <- crumblr(data$counts)
	checkEqualsNumeric(t(cobj$E), r1)
}




