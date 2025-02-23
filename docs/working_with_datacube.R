# note the path to the datacube
library(here)
load(here::here('data/akpv_datacube.rda'))

# check out the structure of the datacube
class(akpv_datacube) # it is a list
str(akpv_datacube[[1]]) # the first element is a 194 x 27 matrix
length(akpv_datacube)  # there are 1000 list items, each element of which is 
#  a 194 x 27 matrix
attributes(akpv_datacube[[1]]) # notice the attributes of each matrix
rownames(akpv_datacube[[1]]) # the polyid's are the row names
colnames(akpv_datacube[[1]]) # the years are the column names
attr(akpv_datacube[[1]], 'stockid') # there is an extra attribute for 
#  stockid, which is equal to length as rownames, and in the same order

# the datacube is composed of 1000 matrices, where each matrix is an MCMC sample
# if we want to see the MCMC samples for the first year for the first site
unlist(lapply(akpv_datacube,function(x){x[1,1]}))

# we can also use names rather than numbers.  Here is polyid "OA00" for 1997
unlist(lapply(akpv_datacube,function(x){x["OA00","1997"]}))

# find all polyids for a given stock in any matrix using stockid attribute
attr(akpv_datacube[[1]], 'stockid') == 6
# make a subset of those for the first MCMC sample in year 1996
akpv_datacube[[1]][attr(akpv_datacube[[1]], 'stockid') == 6, "1996"]
# sum those to get the first MCMC sample for total abundance of stock 6 in 1996
sum(akpv_datacube[[1]][attr(akpv_datacube[[1]], 'stockid') == 6, "1996"])
# get the sum for the 2nd MCMC sample 
i = 2 # MCMC sample
sum(akpv_datacube[[i]][attr(akpv_datacube[[i]], 'stockid') == 6, "1996"])

mode = function(x) {	
	dout = density(x)
	mean(dout$x[dout$y == max(dout$y)])
}

sitebyyearmedians = matrix(NA, nrow = dim(akpv_datacube[[1]])[1], 
	ncol = dim(akpv_datacube[[1]])[2])
for(i in 1:dim(akpv_datacube[[1]])[1]) {
	for(j in 1:dim(akpv_datacube[[1]])[2]) {
		sitebyyearmedians[i,j] =
			median(unlist(lapply(akpv_datacube,function(x){x[i,j]})))
	}
}
rownames(sitebyyearmedians) = rownames(akpv_datacube[[1]])
colnames(sitebyyearmedians) = colnames(akpv_datacube[[1]])
attr(sitebyyearmedians, 'stockid') = attr(akpv_datacube[[1]], 'stockid')
stock_id = 5
attr(sitebyyearmedians, 'stockid') == stock_id
years = as.character(1996:2023)
summedians = rep(NA, times = 28)
for(i in 1:28)
	summedians[i] = sum(sitebyyearmedians[attr(sitebyyearmedians, 
		'stockid') == stock_id, years[i]])
plot(years, summedians)

#Nmin

nmin_list <- vector("list", length=12L)

for (stock_id in 1:12) {
pop = matrix(NA, nrow = 1000, ncol = ncol(akpv_datacube[[1]]))
for(i in 1:1000) pop[i,] =
apply(akpv_datacube[[i]][attr(akpv_datacube[[i]], 'stockid') ==
stock_id,], 2, sum)
Nmin_allyears = apply(pop,2,quantile, prob = .2)
Nmin_2023 = Nmin_allyears[28]
nmin_list[[stock_id]] <- Nmin_2023
}

library(datapasta)
tibble(Nmin = unlist(nmin_list), stock_id = 1:12) |> 
	datapasta::tribble_format()
