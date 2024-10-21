# set a path to the datacube
mypath = paste0('/mnt/ExtraDrive1/Work/desktop_data/2022_papers/',
	'HSsurv2023/data/')
# load the datacube
load(paste0(mypath, 'akpv_datacube.rda'))

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

