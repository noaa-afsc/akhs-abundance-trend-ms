load(paste0('/mnt/ExtraDrive1/Work/desktop_data/2022_papers/',
	'HSsurv/data/', 'akpv_datacube.rda'))
load(paste0('/mnt/ExtraDrive1/Work/desktop_data/2022_papers/',
	'HSsurv/data/','dTerr.rda'))
load(paste0('/mnt/ExtraDrive1/Work/desktop_data/2022_papers/',
	'HSsurv/data/','dGlac.rda'))


figpath = paste0('/mnt/ExtraDrive1/Work/desktop_data/2022_papers/',
	'HSsurv/inst/scripts/DataCube/figures/')

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
#                 Fig_abundance_by_stock
#-------------------------------------------------------------------------------

dTerr = dTerr[dTerr$yr > 1995,]
dTerr$polyid = as.factor(as.character(dTerr$polyid))
dTerr$stockid = as.factor(as.character(dTerr$stockid))
dGlac = dGlac[dGlac$yr > 1995,]
dGlac$polyid = as.factor(as.character(dGlac$polyid))
dGlac$stockid = as.factor(as.character(dGlac$stockid))
# rename column so consistent with terrestrial datasets
names(dGlac)[5] = 'count'
# remove records with missing count values
dGlac = dGlac[!is.na(dGlac$count),]
dstk = rbind(dTerr[,c('polyid', 'stockid', 'yr')],
	dGlac[,c('polyid', 'stockid', 'yr')])

#sitebyyearmedians = matrix(NA, nrow = dim(akpv_datacube[[1]])[1], 
#	ncol = dim(akpv_datacube[[1]])[2])
#for(i in 1:dim(akpv_datacube[[1]])[1]) {
#	for(j in 1:dim(akpv_datacube[[1]])[2]) {
#		sitebyyearmedians[i,j] =
#			median(unlist(lapply(akpv_datacube,function(x){x[i,j]})))
#	}
#}
#rownames(sitebyyearmedians) = rownames(akpv_datacube[[1]])
#colnames(sitebyyearmedians) = colnames(akpv_datacube[[1]])
#attr(sitebyyearmedians, 'stockid') = attr(akpv_datacube[[1]], 'stockid')
#attr(sitebyyearmedians, 'stockid') == stock_id

mode = function(x) {	
	dout = density(x)
	dout$x[dout$y == max(dout$y)]
}


plot_abundance = function(stock_id){
	
#	years = as.character(1996:2023)
#	summedians = rep(NA, times = 28)
#	for(i in 1:28)
#		summedians[i] = sum(sitebyyearmedians[attr(sitebyyearmedians, 
#			'stockid') == stock_id, years[i]])

	cexLab =  2.0
	cexAxis = 1.5
	cexYax2 = 1.5
	err_bar_lwd = 1.7
	abu_pts_cex = 1.7
	marTop = c(0,5,3,1)
	marBot = c(4,5,0,1)
	cex_main = 1.4
	pop = matrix(NA, nrow = 1000, ncol = ncol(akpv_datacube[[1]]))
	for(i in 1:1000) pop[i,] = 
		apply(akpv_datacube[[i]][attr(akpv_datacube[[i]], 'stockid') == 
			stock_id,], 2, sum)
	bot = apply(pop, 2, quantile, prob = .025)
	top = apply(pop, 2, quantile, prob = .975)
	# get average abundance by polyid and year by averaging over MCMC samples
	abu_poly_by_year = akpv_datacube[[1]][attr(akpv_datacube[[1]], 'stockid') == 
				stock_id,]
	for(i in 2:length(akpv_datacube)) 
		abu_poly_by_year = abu_poly_by_year +
			akpv_datacube[[i]][attr(akpv_datacube[[i]], 'stockid') == 
					stock_id,]
	abu_poly_by_year = abu_poly_by_year/length(akpv_datacube)
	rownames(abu_poly_by_year)
	# which polyids were sampled in year i
	fraci = rep(0, times = ncol(akpv_datacube[[1]]))
	for(i in 1:ncol(akpv_datacube[[1]])) {
		ind = rownames(abu_poly_by_year) %in% 
			levels(as.factor(as.character(dstk[dstk$stockid == stock_id & 
				dstk$yr == 1995 + i,'polyid'])))
		fraci[i] = sum(abu_poly_by_year[ind,i])/sum(abu_poly_by_year[,i])
	}
	par(mar = marTop)
	plot(c(1,ncol(akpv_datacube[[1]])), c(min(bot),max(top)), type = 'n',
		xaxt = 'n', ylab = 'Abundance', cex.main = cex_main,
		xlab = '', cex.lab = cexLab, cex.axis = cexAxis,
		main = paste0('Stock ', stock_id, ': ', 
			stockid_names[stock_id,'stocknames']))
	pchtype = rep(1, times = ncol(akpv_datacube[[1]]))
	pchtype[fraci > 0] = 19
	points(apply(pop, 2, mean), pch = pchtype, cex = 2, lwd = 2)
	for(i in 1:ncol(akpv_datacube[[1]]))
		lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)
	points(apply(pop, 2, mean), pch = pchtype, cex = 2, lwd = 2)
#	points(summedians, pch = 3, cex = 2, lwd = 2)
	
	# plot the fractions sampled per year
	par(mar = marBot)
	plot(c(1,ncol(akpv_datacube[[1]])), c(0,1), type = 'n',
		xaxt = 'n', yaxt = 'n', cex.lab = cexLab, cex.axis = cexAxis,
		xlab = '', ylab = '')
	rect(xleft = (1:ncol(akpv_datacube[[1]]))-.5, ybottom = rep(0, times = ncol(akpv_datacube[[1]])),
		xright = (1:ncol(akpv_datacube[[1]]))+.5, ytop = fraci, col = rep('blue', times = 22))
	axis(1, at = c(5,10,15,20,25), labels = c(2000, 2005, 2010, 2015, 2020),
		cex.axis = cexAxis)
	axis(2, at = c(0,.5), 
		labels = c(0, .5), cex.axis = cexYax2)
	mtext('Year', side = 1, padj = 2.2, cex = .65*cexLab)
	mtext('Effort', side = 2, padj = -2.2, cex = .65*cexLab)
}


pdf(paste0(figpath, 'Fig_abundance_by_stock.pdf'), width = 8.5, height = 11)

	layout(matrix(1:24, ncol = 3), heights = rep(c(2,1.2), times= 12))

	plot_abundance(stock_id = 1)
	plot_abundance(stock_id = 4)
	plot_abundance(stock_id = 7)
	plot_abundance(stock_id = 10)
	plot_abundance(stock_id = 2)
	plot_abundance(stock_id = 5)
	plot_abundance(stock_id = 8)
	plot_abundance(stock_id = 11)
	plot_abundance(stock_id = 3)
	plot_abundance(stock_id = 6)
	plot_abundance(stock_id = 9)
	plot_abundance(stock_id = 12)
	
	layout(1)
	
dev.off()

png(paste0(figpath, 'Fig_abundance_by_stock.png'), width = 1000, height = 1100)

	layout(matrix(1:24, ncol = 3), heights = rep(c(2,1.2), times= 12))

	plot_abundance(stock_id = 1)
	plot_abundance(stock_id = 4)
	plot_abundance(stock_id = 7)
	plot_abundance(stock_id = 10)
	plot_abundance(stock_id = 2)
	plot_abundance(stock_id = 5)
	plot_abundance(stock_id = 8)
	plot_abundance(stock_id = 11)
	plot_abundance(stock_id = 3)
	plot_abundance(stock_id = 6)
	plot_abundance(stock_id = 9)
	plot_abundance(stock_id = 12)
	
	layout(1)
	
dev.off()



#-------------------------------------------------------------------------------
#                 Fig_8yr_trend_by_stock
#-------------------------------------------------------------------------------

plot_trend_absolute = function(stock_id){
	pop = matrix(NA, nrow = 1000, ncol = ncol(akpv_datacube[[1]]))
	for(i in 1:1000) pop[i,] = 
		apply(akpv_datacube[[i]][
			attr(akpv_datacube[[i]], 'stockid') == stock_id,], 2, sum)
  maxi = ncol(akpv_datacube[[1]])
  trendlen = 8
  linTrendMat = NULL
  for(i in 1:(maxi - trendlen + 1))
    linTrendMat = cbind(linTrendMat,
      apply(pop,1,function(v){coef(lm(y~x, 
        data.frame(x=1:8,y = v[i:(i + trendlen - 1)])))[2]}))
  bot = apply(linTrendMat, 2, quantile, prob = .025)
  top = apply(linTrendMat, 2, quantile, prob = .975)
  par(mar = c(5,5,5,1))
    plot(c(1,length(top)), c(min(bot),max(top)), type = 'n',
      xaxt = 'n', ylab = 'Trailing 8-Year Trend (Seals/Year)', 
      cex.main = 1.5, xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
			main = paste0('Stock ', stock_id, ': ', 
				stockid_names[stock_id,'stocknames']))
    axis(1, at = c(1, 5, 10, 15, 20), 
			labels = c(2003, 2007, 2012, 2017, 2022),
			cex.axis = 1.5)
    lines(c(1, 22), c(0,0), lty = 2, lwd = 3, col ='red')
    points(apply(linTrendMat, 2, mean), pch = 19, cex = 2)
    for(i in 1:length(top))
      lines(c(i,i), c(bot[i], top[i]), lty = 1, lwd = 2)
}


pdf(paste0(figpath, 'Fig_8yr_trend_by_stock.pdf'), width = 13, height = 17)

	layout(matrix(1:12, nrow = 4, ncol = 3, byrow = TRUE))
	
	plot_trend_absolute(stock_id = 1)
	plot_trend_absolute(stock_id = 2)
	plot_trend_absolute(stock_id = 3)
	plot_trend_absolute(stock_id = 4)
	plot_trend_absolute(stock_id = 5)
	plot_trend_absolute(stock_id = 6)
	plot_trend_absolute(stock_id = 7)
	plot_trend_absolute(stock_id = 8)
	plot_trend_absolute(stock_id = 9)
	plot_trend_absolute(stock_id = 10)
	plot_trend_absolute(stock_id = 11)
	plot_trend_absolute(stock_id = 12)

	layout(1)
	
dev.off()

png(paste0(figpath, 'Fig_8yr_trend_by_stock.png'), width = 960, height = 1240)

	layout(matrix(1:12, nrow = 4, ncol = 3, byrow = TRUE))
	
	plot_trend_absolute(stock_id = 1)
	plot_trend_absolute(stock_id = 2)
	plot_trend_absolute(stock_id = 3)
	plot_trend_absolute(stock_id = 4)
	plot_trend_absolute(stock_id = 5)
	plot_trend_absolute(stock_id = 6)
	plot_trend_absolute(stock_id = 7)
	plot_trend_absolute(stock_id = 8)
	plot_trend_absolute(stock_id = 9)
	plot_trend_absolute(stock_id = 10)
	plot_trend_absolute(stock_id = 11)
	plot_trend_absolute(stock_id = 12)

	layout(1)
	
dev.off()


#-------------------------------------------------------------------------------
#                 Fig_8yr_precent_trend_by_stock
#-------------------------------------------------------------------------------

plot_trend_percent = function(stock_id){
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
  bot = apply(propTrendMat,2,quantile, prob = .025)
  top = apply(propTrendMat,2,quantile, prob = .975)
  par(mar = c(5,5,5,1))
    plot(c(1,length(top)), c(min(bot),max(top)), type = 'n',
      xaxt = 'n', ylab = 'Trend (%/Year)', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
			main = paste0('Stock ', stock_id, ': ', 
				stockid_names[stock_id,'stocknames']), 
				ylim = c(-6,9))
    axis(1, at = c(1,5,10,15,20), labels = c(2003, 2007, 2012, 2017, 2022),
     cex.axis = 1.5)
    lines(c(1,22),c(0,0), lty = 2, lwd = 3, col ='red')
    points(apply(propTrendMat,2,mean), pch = 19, cex = 2)
    for(i in 1:length(top))
      lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)
}

pdf(paste0(figpath, 'Fig_8yr_precent_trend_by_stock.pdf'), 
	width = 13, height = 17)

	layout(matrix(1:12, nrow = 4, ncol = 3, byrow = TRUE))

	plot_trend_percent(stock_id = 1)
	plot_trend_percent(stock_id = 2)
	plot_trend_percent(stock_id = 3)
	plot_trend_percent(stock_id = 4)
	plot_trend_percent(stock_id = 5)
	plot_trend_percent(stock_id = 6)
	plot_trend_percent(stock_id = 7)
	plot_trend_percent(stock_id = 8)
	plot_trend_percent(stock_id = 9)
	plot_trend_percent(stock_id = 10)
	plot_trend_percent(stock_id = 11)
	plot_trend_percent(stock_id = 12)

	layout(1)
	
dev.off()

png(paste0(figpath, 'Fig_8yr_precent_trend_by_stock.png'), width = 960, height = 1240)

	layout(matrix(1:12, nrow = 4, ncol = 3, byrow = TRUE))

	plot_trend_percent(stock_id = 1)
	plot_trend_percent(stock_id = 2)
	plot_trend_percent(stock_id = 3)
	plot_trend_percent(stock_id = 4)
	plot_trend_percent(stock_id = 5)
	plot_trend_percent(stock_id = 6)
	plot_trend_percent(stock_id = 7)
	plot_trend_percent(stock_id = 8)
	plot_trend_percent(stock_id = 9)
	plot_trend_percent(stock_id = 10)
	plot_trend_percent(stock_id = 11)
	plot_trend_percent(stock_id = 12)

	layout(1)
	
dev.off()

#-------------------------------------------------------------------------------
#                 Fig_CV_by_stock
#-------------------------------------------------------------------------------

plot_CV = function(stock_id){
	pop = matrix(NA, nrow = 1000, ncol = ncol(akpv_datacube[[1]]))
	for(i in 1:1000) pop[i,] = 
		apply(akpv_datacube[[i]][
			attr(akpv_datacube[[i]], 'stockid') == stock_id,], 2, sum)
  CV = sqrt(apply(pop,2,var))/apply(pop,2,mean)
  plot(1:length(CV), CV, 
      xaxt = 'n', ylab = 'CV', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2, cex.axis = 1.5, type = 'l', lwd = 3,
			main = paste0('Stock ', stock_id, ': ', 
				stockid_names[stock_id,'stocknames']))
  axis(1, at = c(5,10,15,20, 25), labels = c(2000, 2005, 2010, 2015, 2020),
   cex.axis = 1.5)
}

pdf(paste0(figpath, 'Fig_CV_by_stock.pdf'), 
	width = 13, height = 17)

	layout(matrix(1:12, nrow = 4, ncol = 3, byrow = TRUE))

	plot_CV(stock_id = 1)
	plot_CV(stock_id = 2)
	plot_CV(stock_id = 3)
	plot_CV(stock_id = 4)
	plot_CV(stock_id = 5)
	plot_CV(stock_id = 6)
	plot_CV(stock_id = 7)
	plot_CV(stock_id = 8)
	plot_CV(stock_id = 9)
	plot_CV(stock_id = 10)
	plot_CV(stock_id = 11)
	plot_CV(stock_id = 12)

	layout(1)
	
dev.off()

#-------------------------------------------------------------------------------
#                 Fig_statewide
#-------------------------------------------------------------------------------

pdf(paste0(figpath, 'Fig_statewide.pdf'), width = 11, height = 8.5)

	layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))

# State-wide Abundance
	pop = matrix(NA, nrow = 1000, ncol = ncol(akpv_datacube[[1]]))
	for(i in 1:1000) pop[i,] = 
		apply(akpv_datacube[[i]], 2, sum)
	bot = apply(pop, 2, quantile, prob = .025)
	top = apply(pop, 2, quantile, prob = .975)
	par(mar = c(5,5,5,1))
	plot(c(1,ncol(akpv_datacube[[1]])), c(min(bot),max(top)), type = 'n',
		xaxt = 'n', ylab = 'Estimated Abundance', cex.main = 1.5,
		xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
		main = 'State-wide Abundance')
	axis(1, at = c(5, 10, 15, 20, 25), labels = c(2000, 2005, 2010, 2015, 2020),
		cex.axis = 1.5)
	points(apply(pop,2,mean), pch = 19, cex = 2)
	for(i in 1:ncol(akpv_datacube[[1]]))
		lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)

# State-wide CV
  CV = sqrt(apply(pop,2,var))/apply(pop,2,mean)
  plot(1:length(CV), CV, 
      xaxt = 'n', ylab = 'CV', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
      main = 'CV by Year', type = 'l', lwd = 3)
  axis(1, at = c(5,10,15,20, 25), labels = c(2000, 2005, 2010, 2015, 2020),
   cex.axis = 1.5)

# State-wide Trailing 8-year Trend
 
  maxi = ncol(akpv_datacube[[1]])
  trendlen = 8
  linTrendMat = NULL
  for(i in 1:(maxi - trendlen + 1))
    linTrendMat = cbind(linTrendMat,
      apply(pop,1,function(v){coef(lm(y~x, 
        data.frame(x=1:8,y = v[i:(i + trendlen - 1)])))[2]}))
  bot = apply(linTrendMat,2,quantile, prob = .025)
  top = apply(linTrendMat,2,quantile, prob = .975)
  par(mar = c(5,5,5,1))
    plot(c(1,length(top)), c(min(bot),max(top)), type = 'n',
      xaxt = 'n', ylab = 'Trend (Seals/Year)', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
      main = paste0('Trailing ', trendlen, '-Year Trend by Year'))
    axis(1, at = c(1,5,10,15,20), labels = c(2003, 2007, 2012, 2017, 2022),
     cex.axis = 1.5)
    lines(c(1,22),c(0,0), lty = 2, lwd = 3, col ='red')
    points(apply(linTrendMat,2,mean), pch = 19, cex = 2)
    for(i in 1:length(top))
      lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)

# moving trailing 8-yr trend estimates multiplicative        
  maxi = ncol(akpv_datacube[[1]])
  trendlen = 8
  propTrendMat = NULL
  for(i in 1:(maxi - trendlen + 1))
    propTrendMat = cbind(propTrendMat,
      100*(exp(apply(pop,1,function(v){coef(lm(I(log(y))~x, 
        data.frame(x=1:8,y = v[i:(i + trendlen - 1)])))[2]}))-1))
  bot = apply(propTrendMat,2,quantile, prob = .025)
  top = apply(propTrendMat,2,quantile, prob = .975)
  par(mar = c(5,5,5,1))
    plot(c(1,length(top)), c(min(bot),max(top)), type = 'n',
      xaxt = 'n', ylab = 'Trend (%/Year)', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
      main = paste0('Trailing ', trendlen, '-Year Trend by Year'))
    axis(1, at = c(1,5,10,15,20), labels = c(2003, 2007, 2012, 2017, 2022),
     cex.axis = 1.5)
    lines(c(1,22),c(0,0), lty = 2, lwd = 3, col ='red')
    points(apply(propTrendMat,2,mean), pch = 19, cex = 2)
    for(i in 1:length(top))
      lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)

	layout(1)

dev.off()


cbind(2003:2023,apply(propTrendMat,2,mean))
cbind(2003:2023,apply(linTrendMat,2,mean))

png(paste0(figpath, 'Fig_statewide.png'), width = 620, height = 620)

	layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))

# State-wide Abundance
	pop = matrix(NA, nrow = 1000, ncol = ncol(akpv_datacube[[1]]))
	for(i in 1:1000) pop[i,] = 
		apply(akpv_datacube[[i]], 2, sum)
	bot = apply(pop, 2, quantile, prob = .025)
	top = apply(pop, 2, quantile, prob = .975)
	par(mar = c(5,5,5,1))
	plot(c(1,ncol(akpv_datacube[[1]])), c(min(bot),max(top)), type = 'n',
		xaxt = 'n', ylab = 'Estimated Abundance', cex.main = 1.5,
		xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
		main = 'State-wide Abundance')
	axis(1, at = c(5, 10, 15, 20, 25), labels = c(2000, 2005, 2010, 2015, 2020),
		cex.axis = 1.5)
	points(apply(pop,2,mean), pch = 19, cex = 2)
	for(i in 1:ncol(akpv_datacube[[1]]))
		lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)

# State-wide CV
  CV = sqrt(apply(pop,2,var))/apply(pop,2,mean)
  plot(1:length(CV), CV, 
      xaxt = 'n', ylab = 'CV', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
      main = 'CV by Year', type = 'l', lwd = 3)
  axis(1, at = c(5,10,15,20, 25), labels = c(2000, 2005, 2010, 2015, 2020),
   cex.axis = 1.5)

# State-wide Trailing 8-year Trend
 
  maxi = ncol(akpv_datacube[[1]])
  trendlen = 8
  linTrendMat = NULL
  for(i in 1:(maxi - trendlen + 1))
    linTrendMat = cbind(linTrendMat,
      apply(pop,1,function(v){coef(lm(y~x, 
        data.frame(x=1:8,y = v[i:(i + trendlen - 1)])))[2]}))
  bot = apply(linTrendMat,2,quantile, prob = .025)
  top = apply(linTrendMat,2,quantile, prob = .975)
  par(mar = c(5,5,5,1))
    plot(c(1,length(top)), c(min(bot),max(top)), type = 'n',
      xaxt = 'n', ylab = 'Trend (Seals/Year)', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
      main = paste0('Trailing ', trendlen, '-Year Trend by Year'))
    axis(1, at = c(1,5,10,15,20), labels = c(2003, 2007, 2012, 2017, 2022),
     cex.axis = 1.5)
    lines(c(1,22),c(0,0), lty = 2, lwd = 3, col ='red')
    points(apply(linTrendMat,2,mean), pch = 19, cex = 2)
    for(i in 1:length(top))
      lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)

# moving trailing 8-yr trend estimates multiplicative        
  maxi = ncol(akpv_datacube[[1]])
  trendlen = 8
  propTrendMat = NULL
  for(i in 1:(maxi - trendlen + 1))
    propTrendMat = cbind(propTrendMat,
      100*(exp(apply(pop,1,function(v){coef(lm(I(log(y))~x, 
        data.frame(x=1:8,y = v[i:(i + trendlen - 1)])))[2]}))-1))
  bot = apply(propTrendMat,2,quantile, prob = .025)
  top = apply(propTrendMat,2,quantile, prob = .975)
  par(mar = c(5,5,5,1))
    plot(c(1,length(top)), c(min(bot),max(top)), type = 'n',
      xaxt = 'n', ylab = 'Trend (%/Year)', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
      main = paste0('Trailing ', trendlen, '-Year Trend by Year'))
    axis(1, at = c(1,5,10,15,20), labels = c(2003, 2007, 2012, 2017, 2022),
     cex.axis = 1.5)
    lines(c(1,22),c(0,0), lty = 2, lwd = 3, col ='red')
    points(apply(propTrendMat,2,mean), pch = 19, cex = 2)
    for(i in 1:length(top))
      lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)

	layout(1)

dev.off()


#-------------------------------------------------------------------------------
#                 Fig_statewide 3-panel
#-------------------------------------------------------------------------------

png(paste0(figpath, 'Fig_statewide_3panel.png'), width = 620, height = 720)

	layout(matrix(1:2, nrow = 2, ncol = 1), heights = c(1.05, 1.2))

# State-wide Abundance

	pop = matrix(NA, nrow = 1000, ncol = ncol(akpv_datacube[[1]]))
	for(i in 1:1000) pop[i,] = 
		apply(akpv_datacube[[i]], 2, sum)
	bot = apply(pop, 2, quantile, prob = .025)
	top = apply(pop, 2, quantile, prob = .975)
	par(mar = c(0,5,1,1))
	plot(c(1.1,ncol(akpv_datacube[[1]])), c(min(bot),max(top)), type = 'n',
		xaxt = 'n', ylab = 'Estimated Abundance', cex.main = 1.5,
		xlab = 'Year', cex.lab = 2.5, cex.axis = 1.5,
#		main = 'State-wide Abundance'
	)
	grid(nx = NULL,
     ny = NA,
     lty = 2, col = "gray", lwd = 2)
#	axis(1, at = c(5, 10, 15, 20, 25), labels = c(2000, 2005, 2010, 2015, 2020),
#		cex.axis = 1.5)
	points(apply(pop,2,mean), pch = 19, cex = 2)
	for(i in 1:ncol(akpv_datacube[[1]]))
		lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)


# State-wide Trailing 8-year Trend
 
#  maxi = ncol(akpv_datacube[[1]])
#  trendlen = 8
#  linTrendMat = NULL
#  for(i in 1:(maxi - trendlen + 1))
#    linTrendMat = cbind(linTrendMat,
#      apply(pop,1,function(v){coef(lm(y~x, 
#        data.frame(x=1:8,y = v[i:(i + trendlen - 1)])))[2]}))
#  bot = apply(linTrendMat,2,quantile, prob = .025)
#  top = apply(linTrendMat,2,quantile, prob = .975)
#  par(mar = c(0,5,0,1))
#	plot(c(1.1,ncol(akpv_datacube[[1]])), c(min(bot),max(top)), type = 'n',
#      xaxt = 'n', yaxt = 'n', ylab = 'Trend (Seals/Year)', cex.main = 1.5,
#      xlab = 'Year', cex.lab = 2.5, cex.axis = 1.5,
#      main = paste0('Trailing ', trendlen, '-Year Trend by Year')
#      )
#    axis(2, at = c(-2000,0,2000,4000), labels = c(-2000, 0, 2000, 4000),
#     cex.axis = 1.5)
#	grid(nx = NULL,
#     ny = NA,
#     lty = 2, col = "gray", lwd = 2)
#    lines(c(1,28),c(0,0), lty = 2, lwd = 3, col ='red')
#    points(8:28, apply(linTrendMat,2,mean), pch = 19, cex = 2)
#    for(i in 1:length(top))
#      lines(c(i,i) + 7,c(bot[i],top[i]), lty = 1, lwd = 2)

# moving trailing 8-yr trend estimates multiplicative        
  maxi = ncol(akpv_datacube[[1]])
  trendlen = 8
  propTrendMat = NULL
  for(i in 1:(maxi - trendlen + 1))
    propTrendMat = cbind(propTrendMat,
      100*(exp(apply(pop,1,function(v){coef(lm(I(log(y))~x, 
        data.frame(x=1:8,y = v[i:(i + trendlen - 1)])))[2]}))-1))
  bot = apply(propTrendMat,2,quantile, prob = .025)
  top = apply(propTrendMat,2,quantile, prob = .975)
  par(mar = c(5,5,0,1))
	plot(c(1.1,ncol(akpv_datacube[[1]])), c(min(bot),max(top)), type = 'n',
      xaxt = 'n', ylab = 'Trend (%/Year)', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2.5, cex.axis = 1.5,
#      main = paste0('Trailing ', trendlen, '-Year Trend by Year')
      )
		axis(1, at = c(5, 10, 15, 20, 25), labels = c(2000, 2005, 2010, 2015, 2020),
			cex.axis = 1.5)
	grid(nx = NULL,
     ny = NA,
     lty = 2, col = "gray", lwd = 2)
    lines(c(1,28),c(0,0), lty = 2, lwd = 3, col ='red')
    points(8:28, apply(propTrendMat,2,mean), pch = 19, cex = 2)
    for(i in 1:length(top))
      lines(c(i,i) + 7, c(bot[i],top[i]), lty = 1, lwd = 2)

	layout(1)

dev.off()

#-------------------------------------------------------------------------------
#                 Tab_sample_sizes
#-------------------------------------------------------------------------------

load(paste0('/mnt/ExtraDrive1/Work/desktop_data/2022_papers/',
	'HSsurv/data/', 'akpvpolys_sf.rda'))

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

colnames(num_sampled) = c('stockid', 'n_SSU', 'n_sampled','n_nonzero')

write.table(num_sampled, file = paste0(figpath, 'num_sampled.csv'), quote = FALSE,
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
write.csv(current_est_table, file = paste0(		
	'/mnt/ExtraDrive1/Work/desktop_data/2022_papers/',
	'HSsurv/inst/scripts/DataCube/current_est_table.csv'),
	quote = FALSE, row.names = FALSE)
