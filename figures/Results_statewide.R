load(paste0('/mnt/ExtraDrive1/Work/desktop_data/2024_papers/',
	'akhs-abundance-trends/akhs-abundance-trends/data/', 'akpv_datacube.rda'))
load(paste0('/mnt/ExtraDrive1/Work/desktop_data/2024_papers/',	
	'akhs-abundance-trends/akhs-abundance-trends/data-raw/','dTerr.rda'))
load(paste0('/mnt/ExtraDrive1/Work/desktop_data/2024_papers/',
	'akhs-abundance-trends/akhs-abundance-trends/data-raw/','dGlac.rda'))


figpath = paste0('/mnt/ExtraDrive1/Work/desktop_data/2024_papers/',
	'akhs-abundance-trends/akhs-abundance-trends/figures/png/')

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
#                 Fig_statewide 3-panel
#-------------------------------------------------------------------------------

png(paste0(figpath, 'Fig_statewide_abu_trend.png'), width = 620, height = 720)

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
	)
	grid(nx = NULL,
     ny = NA,
     lty = 2, col = "gray", lwd = 2)
	polygon(x = c(1:28,28:1),y = c(top,rev(bot)), col = 'grey', border = 'grey')
	lines(1:28, apply(pop,2,mean), pch = 19, lwd = 2)
	points(apply(pop,2,mean), pch = 19, cex = 2)
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
      )
		axis(1, at = c(5, 10, 15, 20, 25), labels = c(2000, 2005, 2010, 2015, 2020),
			cex.axis = 1.5)
	grid(nx = NULL,
     ny = NA,
     lty = 2, col = "gray", lwd = 2)
  polygon(x = c(8:28,28:8),y = c(top,rev(bot)), col = 'grey', border = 'grey')
	lines(8:28, apply(propTrendMat,2,mean), lwd = 2)
	points(8:28, apply(propTrendMat,2,mean), pch = 19, cex = 2)
  lines(c(1,28),c(0,0), lty = 2, lwd = 3, col ='black')

	layout(1)

dev.off()

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


mode = function(x) {	
	dout = density(x)
	dout$x[dout$y == max(dout$y)]
}

stock_id = 1
plot_abundance = function(stock_id){
	
	cexLab =  2.0
	cexAxis = 1.5
	cexYax2 = 1.5
	err_bar_lwd = 1.7
	abu_pts_cex = 1.7
	marTop = c(0,5,3,1)
	marBot = c(4,5,0,1)
	cex_main = 1.7
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
	grid(nx = NULL,
     ny = NA,
     lty = 2, col = "gray", lwd = 2)
	pchtype = rep(3, times = ncol(akpv_datacube[[1]]))
	pchtype[fraci > 0] = 19
	polygon(x = c(1:28,28:1),y = c(top,rev(bot)), col = 'grey', border = 'grey')
	vals = 	apply(pop, 2, mean)
	lines(vals, lwd = 1.5)
	valsNA = vals
	valsNA[fraci == 0] = NA
	points(vals, pch = 19, cex = 1.9, lwd = 1.5)
	points(valsNA, pch = 19, cex = 1.9, lwd = 1.5)
	
	# plot the fractions sampled per year
	par(mar = marBot)
	plot(c(1,ncol(akpv_datacube[[1]])), c(0,1), type = 'n',
		xaxt = 'n', yaxt = 'n', cex.lab = cexLab, cex.axis = cexAxis,
		xlab = '', ylab = '')
	grid(nx = NULL,
     ny = NA,
     lty = 2, col = "gray", lwd = 2)
	rect(xleft = (1:ncol(akpv_datacube[[1]]))-.5, ybottom = rep(0, times = ncol(akpv_datacube[[1]])),
		xright = (1:ncol(akpv_datacube[[1]]))+.5, ytop = fraci, col = rep('grey', times = 22))
	axis(1, at = c(5,10,15,20,25), labels = c(2000, 2005, 2010, 2015, 2020),
		cex.axis = cexAxis)
	axis(2, at = c(0,.5), 
		labels = c(0, .5), cex.axis = cexYax2)
	mtext('Year', side = 1, padj = 2.2, cex = .65*cexLab)
	mtext('Effort', side = 2, padj = -2.2, cex = .65*cexLab)
}

#  --------------------- PDF

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

#  --------------------- PNG

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
				ylim = c(-7,9))
    axis(1, at = c(1,5,10,15,20), labels = c(2003, 2007, 2012, 2017, 2022),
     cex.axis = 1.5)
	grid(nx = NULL,
     ny = NA,
     lty = 2, col = "gray", lwd = 2)
	polygon(x = c(1:21,21:1),y = c(top,rev(bot)), col = 'grey', border = 'grey')
  lines(c(1, 22), c(0,0), lty = 2, lwd = 3)
  lines(apply(propTrendMat, 2, mean), pch = 19, lwd = 3)
  points(apply(propTrendMat, 2, mean), pch = 19, cex = 2)
}

#  --------------------- PDF

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

#  --------------------- PNG

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

################################################################################
################################################################################
################################################################################
################################################################################
#                                EXTRAS
################################################################################
################################################################################
################################################################################
################################################################################

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
	polygon(x = c(1:21,21:1),y = c(top,rev(bot)), col = 'grey', border = 'grey')
  lines(c(1, 22), c(0,0), lty = 2, lwd = 3)
  lines(apply(linTrendMat, 2, mean), pch = 19, lwd = 3)
  points(apply(linTrendMat, 2, mean), pch = 19, cex = 2)
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

stkSSUsamples_meanmedian = NULL
for(i in 1:12) {
	dtmp = dstk[dstk$stockid == i,]
	dtmp[,'polyid'] = as.factor(as.character(dtmp$polyid))
	stkSSUsamples_meanmedian = rbind(stkSSUsamples_meanmedian,
		c(mean(table(dtmp$polyid)), median(table(dtmp$polyid))))
}
c(mean(table(dstk$polyid)), median(table(dstk$polyid)))

num_sampled = cbind(num_sampled, stkSSUsamples_meanmedian)
colnames(num_sampled) = c('stockid', 'n_SSU', 'n_sampled_onceplus','n_nonzero',
	'times_sampled_mean', 'times_sampled_median')

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



#-------------------------------------------------------------------------------
#                 Haulout Covariates
#-------------------------------------------------------------------------------

in_path = paste0('/mnt/ExtraDrive1/Work/desktop_data/2024_papers',
	'/HSsurv/HSsurv/data/')

# colors '#D81B60','#1E88E5', '#FFC107','#004D40'
# red, blue, yellow, green
# rgb proportion: (.847, .106, .376), (.118, .533, .898),(1, .757, .027)
load(paste0(in_path,'dHOterr.rda'))
load(paste0(in_path,'dTerr.rda'))
	
# get the aerial survey data for only stocks 3 through 7
dstk = dTerr[dTerr$stockid == 3 | dTerr$stockid == 4 | dTerr$stockid == 5 | 
	dTerr$stockid == 6 | dTerr$stockid == 7,]
# get the haul-out data for only stocks 3 through 7
dHO = dHOterr[dHOterr$stockid == 3 | dHOterr$stockid == 4 | 
	dHOterr$stockid == 5 | dHOterr$stockid == 6 | dHOterr$stockid == 7,]

# only use data after 1995
dstk = dstk[dstk$yr > 1995,]
# make sure polyid and stockid are factors in count data
dstk$polyid = as.factor(as.character(dstk$polyid))
dstk$stockid = as.factor(as.character(dstk$stockid))
IDs = levels(dstk$polyid)
#order by datadatetime within speno, and remove duplicate times
dstk = dstk[order(dstk$polyid, dstk$survey_dt),]

# make sure speno is factor in haulout data
dHO$speno = as.factor(as.character(dHO$speno))
HOIDs = levels(dHO$speno)
#order by datadatetime within speno
dHO = dHO[order(dHO$speno, dHO$date_time),]
# check for any duplicated times (must be increasing within animal)
any(dHO$yrhr0[2:dim(dHO)[1]] - dHO$yrhr0[1:(dim(dHO)[1] - 1)] == 0)

# inspect average haulout proportions by seal ID
y = dHO$y
# create binary data from binned data
y[dHO$y < 0.5] = 0
y[y != 0] = 1
# create design matrix
X <- model.matrix(y ~ dystd + I(dystd^2) + tdstd + I(tdstd^2)
  + hrstd + I(hrstd^2), data = dHO)
# load the MCMC results from fitting haul-out model


#pdf(file = paste0(out_path,'HO_betas.pdf'), width = 15, height = 6)
png(paste0(figpath, 'Fig_HO_expl_var.png'), width = 1000, height = 1300)

layout(matrix(1:12, nrow = 4, byrow = TRUE))

cex_lab = 2.1
cex_axis = 2.5
adj_1 = 0.5
padj_1 = 2.3
adj_2 = 0.5
padj_2 = -2.1
cex_main = 5

load(paste0(in_path, 'MHO_1to2.rda'))
MHO = MHO_1to2

lenbeta = length(MHO[['beta']][[1]])

# fixed effects
beta0 = unlist(lapply(MHO[['beta']], function(x){x[[1]]}))
beta1 = unlist(lapply(MHO[['beta']], function(x){x[[2]]}))
beta2 = unlist(lapply(MHO[['beta']], function(x){x[[3]]}))
beta3 = unlist(lapply(MHO[['beta']], function(x){x[[4]]}))
beta4 = unlist(lapply(MHO[['beta']], function(x){x[[5]]}))
if(lenbeta > 5) {
beta5 = unlist(lapply(MHO[['beta']], function(x){x[[6]]}))
beta6 = unlist(lapply(MHO[['beta']], function(x){x[[7]]}))
}
gam = unlist(lapply(MHO[['gam']], function(x){x[[1]]}))
	
  par(mar = c(6,7,8,3))
  plot(c(1,80),c(0,1), type = 'n',
    ylab = '', xlab = '', xaxt = 'n',
    cex.lab = cex_lab, cex.axis = cex_axis)
  x = ((1:75) - 30)/30
  HOvec = matrix(0, nrow = length(beta1), ncol = 75)
  for(k in 1:length(beta1)) {
    Xb = beta0[k] + mean(unlist(MHO[['gam']][k])) + beta1[k]*x + beta2[k]*x^2
    HOvec[k,] = exp(Xb)/(1+exp(Xb))
  }
  bot = apply(HOvec,2,quantile, probs = .05)
  top = apply(HOvec,2,quantile, probs = .95)
	polygon(x = c(1:75,75:1),y = c(top,rev(bot)), col = 'grey', border = 'grey')
  lines(1:75, apply(HOvec,2,mean), lwd = 4)
	axis(1, at = (0:4)*20, labels = (0:4)*20, cex.axis = 2.5, padj = 0.4)
	mtext('Days since 15 July', side = 1, cex = cex_lab, 
		adj = adj_1, padj = padj_1)
	mtext('Haul-out Probability', side = 2, cex = cex_lab, 
		adj = adj_2, padj = padj_2)

	
# beta - hour-of-day effect

	plot(c(-12,12),c(0,1), type = 'n',
			ylab = '', xaxt = 'n',
			xlab = '',
			cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
			main = 'Stocks 1-2')
	x = (-60:60)/30
	HOvec = matrix(0, nrow = length(beta1), ncol = length(x))
	for(k in 1:length(beta1)) {
			Xb = beta0[k] + mean(unlist(MHO[['gam']][k])) + 
				beta5[k]*x + beta6[k]*x^2
			HOvec[k,] = exp(Xb)/(1+exp(Xb))
	}
  bot = apply(HOvec,2,quantile, probs = .05)
  top = apply(HOvec,2,quantile, probs = .95)
	polygon(x = c((-60:60)/5,rev((-60:60)/5)),y = c(top,rev(bot)), 
		col = 'grey', border = 'grey')
	lines((-60:60)/5, apply(HOvec,2,mean), lwd = 4)
	axis(1, at = (-2:2)*5, labels = (-2:2)*5, cex.axis = 2.5, padj = 0.4)
	mtext('Hours from Solar Noon', side = 1, cex = cex_lab, 
			adj = adj_1, padj = padj_1)
	mtext('Haul-out Probability', side = 2, cex = cex_lab, 
			adj = adj_2, padj = padj_2)

# beta - tide effects (or, hours from solar noon for glacial haulout)

  plot(c(-5,5), c(0,1), type = 'n',
    ylab = '', xlab = '', xaxt = 'n',
    cex.lab = cex_lab, cex.axis = cex_axis)
  x = (-50:50)/25
  xtran = 2.5*x
  HOvec = matrix(0, nrow = length(beta3), ncol = length(x))
  for(k in 1:length(beta3)) {
    Xb = beta0[k] + mean(unlist(MHO[['gam']][k])) + beta3[k]*x + beta4[k]*x^2
    HOvec[k,] = exp(Xb)/(1+exp(Xb))
  }
  bot = apply(HOvec,2,quantile, probs = .05)
  top = apply(HOvec,2,quantile, probs = .95)
	polygon(x = c(xtran,rev(xtran)),y = c(top,rev(bot)), 
		col = 'grey', border = 'grey')
  lines(xtran, apply(HOvec,2,mean), lwd = 4)
	axis(1, at = (-2:2)*2, labels = (-2:2)*2, cex.axis = 2.5, padj = 0.4)
	mtext('Hours from Low Tide', side = 1, cex = cex_lab, 
		adj = adj_1, padj = padj_1)
	mtext('Haul-out Probability', side = 2, cex = cex_lab, 
		adj = adj_2, padj = padj_2)

# beta - date effects

load(paste0(in_path, 'MHO_3to7.rda'))
MHO = MHO_3to7

lenbeta = length(MHO[['beta']][[1]])

# fixed effects
beta0 = unlist(lapply(MHO[['beta']], function(x){x[[1]]}))
beta1 = unlist(lapply(MHO[['beta']], function(x){x[[2]]}))
beta2 = unlist(lapply(MHO[['beta']], function(x){x[[3]]}))
beta3 = unlist(lapply(MHO[['beta']], function(x){x[[4]]}))
beta4 = unlist(lapply(MHO[['beta']], function(x){x[[5]]}))
if(lenbeta > 5) {
beta5 = unlist(lapply(MHO[['beta']], function(x){x[[6]]}))
beta6 = unlist(lapply(MHO[['beta']], function(x){x[[7]]}))
}
gam = unlist(lapply(MHO[['gam']], function(x){x[[1]]}))

 par(mar = c(6,7,8,3))
  plot(c(1,80),c(0,1), type = 'n',
    ylab = '', xlab = '', xaxt = 'n',
    cex.lab = cex_lab, cex.axis = cex_axis)
  x = ((1:75) - 30)/30
  HOvec = matrix(0, nrow = length(beta1), ncol = 75)
  for(k in 1:length(beta1)) {
    Xb = beta0[k] + mean(unlist(MHO[['gam']][k])) + beta1[k]*x + beta2[k]*x^2
    HOvec[k,] = exp(Xb)/(1+exp(Xb))
  }
  bot = apply(HOvec,2,quantile, probs = .05)
  top = apply(HOvec,2,quantile, probs = .95)
	polygon(x = c(1:75,75:1),y = c(top,rev(bot)), col = 'grey', border = 'grey')
  lines(1:75, apply(HOvec,2,mean), lwd = 4)
	axis(1, at = (0:4)*20, labels = (0:4)*20, cex.axis = 2.5, padj = 0.4)
	mtext('Days since 15 July', side = 1, cex = cex_lab, 
		adj = adj_1, padj = padj_1)
	mtext('Haul-out Probability', side = 2, cex = cex_lab, 
		adj = adj_2, padj = padj_2)

# beta - hour-of-day effect

	plot(c(-12,12),c(0,1), type = 'n',
			ylab = '', xaxt = 'n',
			xlab = '',
			cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
			main = 'Stocks 3-7')
	x = (-60:60)/30
	HOvec = matrix(0, nrow = length(beta1), ncol = length(x))
	for(k in 1:length(beta1)) {
			Xb = beta0[k] + mean(unlist(MHO[['gam']][k])) + 
				beta5[k]*x + beta6[k]*x^2
			HOvec[k,] = exp(Xb)/(1+exp(Xb))
	}
  bot = apply(HOvec,2,quantile, probs = .05)
  top = apply(HOvec,2,quantile, probs = .95)
	polygon(x = c((-60:60)/5,rev((-60:60)/5)),y = c(top,rev(bot)), 
		col = 'grey', border = 'grey')
	lines((-60:60)/5, apply(HOvec,2,mean), lwd = 4)
	axis(1, at = (-2:2)*5, labels = (-2:2)*5, cex.axis = 2.5, padj = 0.4)
	mtext('Hours from Solar Noon', side = 1, cex = cex_lab, 
			adj = adj_1, padj = padj_1)
	mtext('Haul-out Probability', side = 2, cex = cex_lab, 
			adj = adj_2, padj = padj_2)

# beta - tide effects (or, hours from solar noon for glacial haulout)

  plot(c(-5,5), c(0,1), type = 'n',
    ylab = '', xlab = '', xaxt = 'n',
    cex.lab = cex_lab, cex.axis = cex_axis)
  x = (-50:50)/25
  xtran = 2.5*x
  HOvec = matrix(0, nrow = length(beta3), ncol = length(x))
  for(k in 1:length(beta3)) {
    Xb = beta0[k] + mean(unlist(MHO[['gam']][k])) + beta3[k]*x + beta4[k]*x^2
    HOvec[k,] = exp(Xb)/(1+exp(Xb))
  }
  bot = apply(HOvec,2,quantile, probs = .05)
  top = apply(HOvec,2,quantile, probs = .95)
	polygon(x = c(xtran,rev(xtran)),y = c(top,rev(bot)), 
		col = 'grey', border = 'grey')
  lines(xtran, apply(HOvec,2,mean), lwd = 4)
	axis(1, at = (-2:2)*2, labels = (-2:2)*2, cex.axis = 2.5, padj = 0.4)
	mtext('Hours from Low Tide', side = 1, cex = cex_lab, 
		adj = adj_1, padj = padj_1)
	mtext('Haul-out Probability', side = 2, cex = cex_lab, 
		adj = adj_2, padj = padj_2)
	

# beta - date effects

load(paste0(in_path, 'MHO_8to12.rda'))
MHO = MHO_8to12

lenbeta = length(MHO[['beta']][[1]])

# fixed effects
beta0 = unlist(lapply(MHO[['beta']], function(x){x[[1]]}))
beta1 = unlist(lapply(MHO[['beta']], function(x){x[[2]]}))
beta2 = unlist(lapply(MHO[['beta']], function(x){x[[3]]}))
beta3 = unlist(lapply(MHO[['beta']], function(x){x[[4]]}))
beta4 = unlist(lapply(MHO[['beta']], function(x){x[[5]]}))
if(lenbeta > 5) {
beta5 = unlist(lapply(MHO[['beta']], function(x){x[[6]]}))
beta6 = unlist(lapply(MHO[['beta']], function(x){x[[7]]}))
}
gam = unlist(lapply(MHO[['gam']], function(x){x[[1]]}))

 par(mar = c(6,7,8,3))
  plot(c(1,80),c(0,1), type = 'n',
    ylab = '', xlab = '', xaxt = 'n',
    cex.lab = cex_lab, cex.axis = cex_axis)
  x = ((1:75) - 30)/30
  HOvec = matrix(0, nrow = length(beta1), ncol = 75)
  for(k in 1:length(beta1)) {
    Xb = beta0[k] + mean(unlist(MHO[['gam']][k])) + beta1[k]*x + beta2[k]*x^2
    HOvec[k,] = exp(Xb)/(1+exp(Xb))
  }
  bot = apply(HOvec,2,quantile, probs = .05)
  top = apply(HOvec,2,quantile, probs = .95)
	polygon(x = c(1:75,75:1),y = c(top,rev(bot)), col = 'grey', border = 'grey')
  lines(1:75, apply(HOvec,2,mean), lwd = 4)
	axis(1, at = (0:4)*20, labels = (0:4)*20, cex.axis = 2.5, padj = 0.4)
	mtext('Days since 15 July', side = 1, cex = cex_lab, 
		adj = adj_1, padj = padj_1)
	mtext('Haul-out Probability', side = 2, cex = cex_lab, 
		adj = adj_2, padj = padj_2)

# beta - hour-of-day effect

	plot(c(-12,12),c(0,1), type = 'n',
			ylab = '', xaxt = 'n',
			xlab = '',
			cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
			main = 'Stocks 8-12')
	x = (-60:60)/30
	HOvec = matrix(0, nrow = length(beta1), ncol = length(x))
	for(k in 1:length(beta1)) {
			Xb = beta0[k] + mean(unlist(MHO[['gam']][k])) + 
				beta5[k]*x + beta6[k]*x^2
			HOvec[k,] = exp(Xb)/(1+exp(Xb))
	}
  bot = apply(HOvec,2,quantile, probs = .05)
  top = apply(HOvec,2,quantile, probs = .95)
	polygon(x = c((-60:60)/5,rev((-60:60)/5)),y = c(top,rev(bot)), 
		col = 'grey', border = 'grey')
	lines((-60:60)/5, apply(HOvec,2,mean), lwd = 4)
	axis(1, at = (-2:2)*5, labels = (-2:2)*5, cex.axis = 2.5, padj = 0.4)
	mtext('Hours from Solar Noon', side = 1, cex = cex_lab, 
			adj = adj_1, padj = padj_1)
	mtext('Haul-out Probability', side = 2, cex = cex_lab, 
			adj = adj_2, padj = padj_2)

# beta - tide effects (or, hours from solar noon for glacial haulout)

  plot(c(-5,5), c(0,1), type = 'n',
    ylab = '', xlab = '', xaxt = 'n',
    cex.lab = cex_lab, cex.axis = cex_axis)
  x = (-50:50)/25
  xtran = 2.5*x
  HOvec = matrix(0, nrow = length(beta3), ncol = length(x))
  for(k in 1:length(beta3)) {
    Xb = beta0[k] + mean(unlist(MHO[['gam']][k])) + beta3[k]*x + beta4[k]*x^2
    HOvec[k,] = exp(Xb)/(1+exp(Xb))
  }
  bot = apply(HOvec,2,quantile, probs = .05)
  top = apply(HOvec,2,quantile, probs = .95)
	polygon(x = c(xtran,rev(xtran)),y = c(top,rev(bot)), 
		col = 'grey', border = 'grey')
  lines(xtran, apply(HOvec,2,mean), lwd = 4)
	axis(1, at = (-2:2)*2, labels = (-2:2)*2, cex.axis = 2.5, padj = 0.4)
	mtext('Hours from Low Tide', side = 1, cex = cex_lab, 
		adj = adj_1, padj = padj_1)
	mtext('Haul-out Probability', side = 2, cex = cex_lab, 
		adj = adj_2, padj = padj_2)
	
###    GLACIAL DATA
	
load(paste0(in_path, 'MHO_glac.rda'))
MHO = MHO_glac

lenbeta = length(MHO[['beta']][[1]])

# fixed effects
beta0 = unlist(lapply(MHO[['beta']], function(x){x[[1]]}))
beta1 = unlist(lapply(MHO[['beta']], function(x){x[[2]]}))
beta2 = unlist(lapply(MHO[['beta']], function(x){x[[3]]}))
beta3 = unlist(lapply(MHO[['beta']], function(x){x[[4]]}))
beta4 = unlist(lapply(MHO[['beta']], function(x){x[[5]]}))

# beta - date effects

 par(mar = c(6,7,8,3))
  plot(c(1,80),c(0,1), type = 'n',
    ylab = '', xlab = '', xaxt = 'n',
    cex.lab = cex_lab, cex.axis = cex_axis)
  x = ((1:75) - 30)/30
  HOvec = matrix(0, nrow = length(beta1), ncol = 75)
  for(k in 1:length(beta1)) {
    Xb = beta0[k] + mean(unlist(MHO[['gam']][k])) + beta1[k]*x + beta2[k]*x^2
    HOvec[k,] = exp(Xb)/(1+exp(Xb))
  }
  bot = apply(HOvec,2,quantile, probs = .05)
  top = apply(HOvec,2,quantile, probs = .95)
	polygon(x = c(1:75,75:1),y = c(top,rev(bot)), col = 'grey', border = 'grey')
  lines(1:75, apply(HOvec,2,mean), lwd = 4)
	axis(1, at = (0:4)*20, labels = (0:4)*20, cex.axis = 2.5, padj = 0.4)
	mtext('Days since 15 July', side = 1, cex = cex_lab, 
		adj = adj_1, padj = padj_1)
	mtext('Haul-out Probability', side = 2, cex = cex_lab, 
		adj = adj_2, padj = padj_2)

# beta - hour-of-day effect

	plot(c(-12,12),c(0,1), type = 'n',
			ylab = '', xaxt = 'n',
			xlab = '',
			cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
			main = 'Glacial Data')
	x = (-60:60)/30
	HOvec = matrix(0, nrow = length(beta1), ncol = length(x))
	for(k in 1:length(beta1)) {
			Xb = beta0[k] + mean(unlist(MHO[['gam']][k])) + 
				beta3[k]*x + beta4[k]*x^2
			HOvec[k,] = exp(Xb)/(1+exp(Xb))
	}
  bot = apply(HOvec,2,quantile, probs = .05)
  top = apply(HOvec,2,quantile, probs = .95)
	polygon(x = c((-60:60)/5,rev((-60:60)/5)),y = c(top,rev(bot)), 
		col = 'grey', border = 'grey')
	lines((-60:60)/5, apply(HOvec,2,mean), lwd = 4)
	axis(1, at = (-2:2)*5, labels = (-2:2)*5, cex.axis = 2.5, padj = 0.4)
	mtext('Hours from Solar Noon', side = 1, cex = cex_lab, 
			adj = adj_1, padj = padj_1)
	mtext('Haul-out Probability', side = 2, cex = cex_lab, 
			adj = adj_2, padj = padj_2)

dev.off()


#-------------------------------------------------------------------------------
#                 Plain Abundance
#-------------------------------------------------------------------------------


stock_id = 1
plot_simple_abundance = function(stock_id){
	
#	years = as.character(1996:2023)
#	summedians = rep(NA, times = 28)
#	for(i in 1:28)
#		summedians[i] = sum(sitebyyearmedians[attr(sitebyyearmedians, 
#			'stockid') == stock_id, years[i]])

	cexLab =  3.0
	cexAxis = 2.5
	abu_pts_cex = 1.7
	marTop = c(6,6,0,0)
	cex_main = 1.7
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
	par(mar = marTop)
	plot(c(1,ncol(akpv_datacube[[1]])), c(min(bot),max(top)), type = 'n',
		xaxt = 'n', ylab = 'Abundance', cex.main = cex_main,
		xlab = '', cex.lab = cexLab, cex.axis = cexAxis)
	grid(nx = NULL,
     ny = NA,
     lty = 2, col = "gray", lwd = 3)
	vals = 	apply(pop, 2, mean)
	lines(vals, lwd = 3.5)
	points(vals, pch = 19, cex = 3.0)
	axis(1, at = c(5,10,15,20,25), labels = c(2000, 2005, 2010, 2015, 2020),
		cex.axis = cexAxis)
	mtext('Year', side = 1, padj = 2.2, cex = cexLab)
}

for(i in 1:12) {
	stock_id = i
	png(paste0(figpath, 'Fig_AbuSimple_stock_',stock_id,'.png'), 
		width = 600, height = 600)
		plot_simple_abundance(i)
	dev.off()
}
