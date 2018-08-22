source('alra.R')

# The B and NK cells from the purified of the manuscript.
if (!file.exists('data/b_nk_example.RDS')) {
    download.file('https://gauss.math.yale.edu/~gcl22/alra/b_nk_example.RDS', 'data/b_nk_example.RDS')
    download.file('https://gauss.math.yale.edu/~gcl22/alra/labels_example.RDS', 'data/labels_example.RDS')
}

A <- readRDS('data/b_nk_example.RDS')
labels <- readRDS('data/labels_example.RDS')

# Library and log normalize the data
A_norm <- normalize_data(A)

A_norm_completed <- alra(A_norm)[[3]]

# Check the completion. Note that the results improve when using the entire dataset (see the ALRA-paper repo for those codes), as opposed to this subset.
#NCAM1 should be expressed in all NK cells, but is only expressed in 4% of of NK cells in the original data.  
#CR2 should be expressed in all B cells, but is only expressed in 1% of of B cells in the original data.  
print(aggregate(A_norm[,c('NCAM1','CR2')], by=list(" "=labels),FUN=function(x) round(c(percent=100*sum(x>0)/length(x)),1)))
print(aggregate(A_norm_completed[,c('NCAM1','CR2')], by=list(" "=labels),FUN=function(x) round(c(percent=100*sum(x>0)/length(x)),1)))
