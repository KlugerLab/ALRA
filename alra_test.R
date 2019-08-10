# Demo for running ALRA on real data. Runtime <10 minutes on standard laptop.

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


# Choose k. 
k_choice <- choose_k(A_norm)

# For the results in the paper, automatically chosen k worked quite well, but in
# some cases you might want to take a closer look, as we do here. The k is
# chosen based on the spacings between the singular values, as it can be quite
# hard to identify the ``beginning of noise'' from just looking at the spectrum
# itself. Uncomment the code below to plot them 

#library(ggplot2)
#library(gridExtra)
#df <- data.frame(x=1:100,y=k_choice$d)
#g1<-ggplot(df,aes(x=x,y=y),) + geom_point(size=1)  + geom_line(size=0.5)+ geom_vline(xintercept=k_choice$k)   + theme( axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_i') + ggtitle('Singular values')
#df <- data.frame(x=2:100,y=diff(k_choice$d))[3:99,]
#g2<-ggplot(df,aes(x=x,y=y),) + geom_point(size=1)  + geom_line(size=0.5)+ geom_vline(xintercept=k_choice$k+1)   + theme(axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_{i} - s_{i-1}') + ggtitle('Singular value spacings')
#grid.arrange(g1,g2,nrow=1)

A_norm_completed <- alra(A_norm,k=k_choice$k)[[3]]

# Check the completion. Note that the results improve when using the entire dataset (see the ALRA-paper repo for those codes), as opposed to this subset.
#NCAM1 should be expressed in all NK cells, but is only expressed in 4% of of NK cells in the original data.  
#CR2 should be expressed in all B cells, but is only expressed in 1% of of B cells in the original data.  
print(aggregate(A_norm[,c('NCAM1','CR2')], by=list(" "=labels),FUN=function(x) round(c(percent=100*sum(x>0)/length(x)),1)))
print(aggregate(A_norm_completed[,c('NCAM1','CR2')], by=list(" "=labels),FUN=function(x) round(c(percent=100*sum(x>0)/length(x)),1)))
