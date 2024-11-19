setwd(paste0('/mnt/ExtraDrive1/Work/desktop_data/2024_papers/',
	'akhs-abundance-trends/akhs-abundance-trends/data/'))
load('akpv_datacube.rda')
setwd(paste0('/mnt/ExtraDrive1/Work/desktop_data/2024_papers/',	
	'akhs-abundance-trends/akhs-abundance-trends/data-raw/'))
load('dTerr.rda')
load('dGlac.rda')
load('akpvpolys_sf.rda')

dstk = rbind(dTerr[,c('polyid', 'stockid', 'yr')],
	dGlac[,c('polyid', 'stockid', 'yr')])

tabpath = paste0('/mnt/ExtraDrive1/Work/desktop_data/2024_papers/',
	'akhs-abundance-trends/akhs-abundance-trends/tables/')

stockid_names = data.frame(stockid = 1:12, stocknames = 
	c('Aleutians',
	'Pribilofs',
	'Bristol Bay',
	'North Kodiak',
	'South Kodiak',
	'Prince William Sound',
	'Cook Inlet/Shelikof Strait',
	'Glacier Bay/Icy Strait',
	'Lynn Canal/Stephens Passage',
	'Sitka/Chatham Strait',
	'Dixon/Cape Decision',
	'Clarence Strait')
)
	

#-------------------------------------------------------------------------------
#                 Tab_sample_sizes
#-------------------------------------------------------------------------------


num_sampled = matrix(NA, nrow = 12, ncol = 4)
for(i in 1:12) { 
num_sampled[i,1] = i
num_sampled[i,2] = sum(akpvpolys_sf$stockid == i)
num_sampled[i,3] = length(unique(dTerr[dTerr$stockid == i, 'polyid'])) + 
	length(unique(dGlac[dGlac$stockid == i, 'polyid']))
num_sampled[i,4] = length(unique(
		dTerr[dTerr$stockid == i & dTerr$count > 0, 'polyid'])) + 
	length(unique(
		dGlac[dGlac$stockid == i & dGlac$count > 0, 'polyid']))
}
num_sampled

stkSSUsamples_meanmedian = NULL
for(i in 1:12) {
	dtmp = dstk[dstk$stockid == i,]
	dtmp[,'polyid'] = as.factor(as.character(dtmp$polyid))
	stkSSUsamples_meanmedian = rbind(stkSSUsamples_meanmedian,
		c(mean(table(dtmp$polyid)), median(table(dtmp$polyid))))
}
c(mean(table(dstk$polyid)), median(table(dstk$polyid)))

num_sampled = cbind(num_sampled, stkSSUsamples_meanmedian)
num_sampled = rbind(num_sampled, c(rep(NA, times = 4), c(mean(table(dstk$polyid)), median(table(dstk$polyid)))))
colnames(num_sampled) = c('stockid', 'n_SSU', 'n_sampled_onceplus','n_nonzero',
	'times_sampled_mean', 'times_sampled_median')

write.table(num_sampled, file = paste0(tabpath, 'num_sampled.csv'), quote = FALSE,
	sep = ',', row.names = FALSE)


#-------------------------------------------------------------------------------
#                 Tab_statewide_and_stocks
#-------------------------------------------------------------------------------

# State-wide Abundance
	pop = matrix(NA, nrow = 1000, ncol = ncol(akpv_datacube[[1]]))
	for(i in 1:1000) pop[i,] = 
		apply(akpv_datacube[[i]], 2, sum)
	est = apply(pop, 2, mean)
	se = sqrt(apply(pop, 2, var))
	bot = apply(pop, 2, quantile, prob = .025)
	top = apply(pop, 2, quantile, prob = .975)
	cv = se/est
	current_est_table = data.frame(stock = 'statewide', est = est[28], 
		se = se[28], bot = bot[28], top = top[28], cv = cv[28])
	for(j in 1:12) {
		stock_id = j
		pop = matrix(NA, nrow = 1000, ncol = ncol(akpv_datacube[[1]]))
		for(i in 1:1000) pop[i,] = 
			apply(akpv_datacube[[i]][attr(akpv_datacube[[i]], 'stockid') == 
				stock_id,], 2, sum)
		est = apply(pop, 2, mean)
		se = sqrt(apply(pop, 2, var))
		bot = apply(pop, 2, quantile, prob = .025)
		top = apply(pop, 2, quantile, prob = .975)
		cv = se/est
		current_est_table = rbind(current_est_table,
			data.frame(stock = stock_id, est = est[28], 
			se = se[28], bot = bot[28], top = top[28], cv = cv[28]))
	}
	


#-------------------------------------------------------------------------------
#                 Probability of decreasing
#-------------------------------------------------------------------------------

	pDecrease = rep(NA, times = 12)
	for(j in 1:12) {
		stock_id = j
		pop = matrix(NA, nrow = 1000, ncol = ncol(akpv_datacube[[1]]))
		for(i in 1:1000) pop[i,] = 
			apply(akpv_datacube[[i]][
				attr(akpv_datacube[[i]], 'stockid') == stock_id,], 2, sum)
		maxi = ncol(akpv_datacube[[1]])
		trendlen = 8
		propTrendMat = NULL
		for(i in 1:(maxi - trendlen + 1))
			propTrendMat = cbind(propTrendMat,
				100*(exp(apply(pop,1,function(v){coef(lm(I(log(y))~x, 
					data.frame(x=1:8,y = v[i:(i + trendlen - 1)])))[2]}))-1))
		pDecrease[j] = sum(propTrendMat[,21] < 0)/1000
	}

current_est_table = cbind(current_est_table, c(NA, pDecrease))
colnames(current_est_table)[7] = 'pDecrease'

# ----------- Table of Results

write.csv(current_est_table, file = paste0(tabpath,'current_est_table.csv'),
	quote = FALSE, row.names = FALSE)




