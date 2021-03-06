% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pls_func.R
\name{plot_explained_variance}
\alias{plot_explained_variance}
\title{Plot explained variance of plsr object}
\usage{
plot_explained_variance(plsr_obj)
}
\arguments{
\item{plsr_obj}{The plsr object.}
}
\description{
Calculates and plots the variance explained in the original data X and Y by each additional latent variable.
}
\examples{
plsr_obj = pls(rating_data, tracking_data,10,10)
plot_explained_variance(plsr_obj)
}
