library(rsvd)

normalize_data <- function (A) {
    #  Simple convenience function to library and log normalize a matrix

    totalUMIPerCell <- rowSums(A);
    A_norm <- sweep(A, 1, totalUMIPerCell, '/');
    A_norm <- A_norm * 10E3
    A_norm <- log(A_norm +1);
}

choose_k <- function (A_norm,K=100, pval_thresh=1E-10, noise_start=80) {
    #  Heuristic for choosing rank k for the low rank approximation based on
    #  statistics of the spacings between consecutive singular values. Finds
    #  the smallest singular value \sigma_i such that $\sigma_i - \sigma_{i-1}
    #  is significantly different than spacings in the tail of the singular values.
    # 
    #
    # Args:
    #   A_norm: The log-transformed expression matrix of cells (rows) vs. genes (columns)
    #   pval_thresh : The threshold for ``significance''
    #   noise_start : Index for which all smaller singular values are considered noise
    #
    # Returns:
    #   A list with three items
    #       1) Chosen k
    #       2) P values of each possible k 
    #       3) Singular values of the matrix A_norm
    noise_svals <- noise_start:K
    rsvd_out <- rsvd(A_norm,K)
    diffs <- diff(rsvd_out$d)
    pvals <- pnorm(diffs,mean(diffs[noise_svals-1]),sd(diffs[noise_svals-1]))
    k <- max(which( pvals  <pval_thresh))
    return (list( k=k,pvals =pvals,d=rsvd_out$d))
}

alra <- function( A_norm, k=0,q=10) {
    # Computes the k-rank approximation to A_norm and adjusts it according to the
    # error distribution learned from the negative values.
    #
    # Args:
    #   A_norm: The log-transformed expression matrix of cells (rows) vs. genes (columns)
    #   k : the rank of the rank-k approximation. Set to 0 for automated choice of k.
    #   q : the number of additional power iterations in randomized SVD
    #
    # Returns:
    #   A list with three items
    #       1) The rank k approximation of A_norm.
    #       2) The rank k approximation of A_norm, adaptively thresholded
    #       3) The rank k approximation of A_norm, adaptively thresholded and with the first two moments of the non-zero values matched to the first two moments of the non-zeros of A_norm
    # Example:
    #     result.completed <- adjusted_svd(A_norm,15)
    #     A_norm_rank15 <- result.completed[[1]]     # The low rank approximation for reference purposes...not suggested for matrix completion
    #     A_norm_rank15_cor <- result.completed[[3]] # The actual adjusted, completed matrix
    if (k ==0 ) {
        k_choice <- choose_k(A_norm)
        k <-  k_choice$k
        print(sprintf("Chose k=%d",k))
    }

    print("Getting nonzeros")
    originally_nonzero <- A_norm >0 

    print("Randomized SVD")
    fastDecomp_noc <- rsvd(A_norm, k,q=q);
    A_norm_rank_k <- fastDecomp_noc$u[,1:k]%*%diag(fastDecomp_noc$d[1:k])%*% t(fastDecomp_noc$v[,1:k])


    print("Find mins")
    A_norm_rank_k_mins <- abs(apply(A_norm_rank_k,2,min))
    print("Sweep")
    A_norm_rank_k_cor <- replace(A_norm_rank_k, A_norm_rank_k <= A_norm_rank_k_mins[col(A_norm_rank_k)], 0)


    sd_nonzero <- function(x) sd(x[!x == 0])
    sigma_1 <- apply(A_norm_rank_k_cor, 2, sd_nonzero)
    sigma_2 <- apply(A_norm, 2, sd_nonzero)
    mu_1 <- colSums(A_norm_rank_k_cor)/colSums(!!A_norm_rank_k_cor)
    mu_2 <- colSums(A_norm)/colSums(!!A_norm)

    toscale <- !is.na(sigma_1) & !is.na(sigma_2)

    print(sprintf("Scaling all except for %d columns", sum(!toscale)))

    sigma_1_2 <- sigma_2/sigma_1
    toadd  <- -1*mu_1*sigma_2/sigma_1 + mu_2

    A_norm_rank_k_temp <- A_norm_rank_k_cor[,toscale]
    A_norm_rank_k_temp <- sweep(A_norm_rank_k_temp,2, sigma_1_2[toscale],FUN = "*")
    A_norm_rank_k_temp <- sweep(A_norm_rank_k_temp,2, toadd[toscale],FUN = "+")

    A_norm_rank_k_cor_sc <- A_norm_rank_k_cor
    A_norm_rank_k_cor_sc[,toscale] <- A_norm_rank_k_temp
    A_norm_rank_k_cor_sc[A_norm_rank_k_cor==0] = 0

    lt0 <- A_norm_rank_k_cor_sc  <0
    A_norm_rank_k_cor_sc[lt0] <- 0 
    print(sprintf("%.2f%% of the values became negative in the scaling process and were set to zero", 100*sum(lt0)/(nrow(A_norm)*ncol(A_norm))))

    A_norm_rank_k_cor_sc[originally_nonzero & A_norm_rank_k_cor_sc ==0] <- A_norm[originally_nonzero & A_norm_rank_k_cor_sc ==0]

    colnames(A_norm_rank_k_cor) <- colnames(A_norm)
    colnames(A_norm_rank_k_cor_sc) <- colnames(A_norm)
    colnames(A_norm_rank_k) <- colnames(A_norm)

    original_nz <- sum(A_norm>0)/(nrow(A_norm)*ncol(A_norm))
    completed_nz <- sum(A_norm_rank_k_cor_sc>0)/(nrow(A_norm)*ncol(A_norm))
    print(sprintf("The matrix went from %.2f%% nonzero to %.2f%% nonzero", 100*original_nz, 100*completed_nz))

    list(A_norm_rank_k=A_norm_rank_k,A_norm_rank_k_cor =A_norm_rank_k_cor, A_norm_rank_k_cor_sc=A_norm_rank_k_cor_sc)
}
