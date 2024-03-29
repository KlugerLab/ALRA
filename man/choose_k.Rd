% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alra.R
\name{choose_k}
\alias{choose_k}
\title{Heuristic for choosing rank k for the low rank approximation based on
statistics of the spacings between consecutive singular values. Finds
the smallest singular value \eqn{\sigma_i} such that \eqn{\sigma_i - \sigma_{i-1}}
is significantly different than spacings in the tail of the singular values.}
\usage{
choose_k(
  A_norm,
  K = 100,
  thresh = 6,
  noise_start = 80,
  q = 2,
  use.mkl = F,
  mkl.seed = -1
)
}
\arguments{
\item{A_norm}{The log-transformed expression matrix of cells (rows) vs. genes (columns)}

\item{K}{Number of singular values to compute. Must be less than the
smallest dimension of the matrix.}

\item{thresh}{Number of standard deviations away from the ``noise'' singular
values which you consider to be signal}

\item{noise_start}{Index for which all smaller singular values are
considered noise}

\item{q}{Number of additional power iterations}

\item{use.mkl}{Use the Intel MKL based implementation of SVD. Needs to be
installed from https://github.com/KlugerLab/rpca-mkl.}

\item{mkl.seed}{Only relevant if use.mkl=T. Set the seed for the random
generator for the Intel MKL implementation of SVD. Any number <0 will
use the current timestamp. If use.mkl=F, set the seed using
set.seed() function as usual.}
}
\value{
A list with three items
\enumerate{
\item Chosen k
\item P values of each possible k
\item Singular values of the matrix A_norm
}
}
\description{
Heuristic for choosing rank k for the low rank approximation based on
statistics of the spacings between consecutive singular values. Finds
the smallest singular value \eqn{\sigma_i} such that \eqn{\sigma_i - \sigma_{i-1}}
is significantly different than spacings in the tail of the singular values.
}
