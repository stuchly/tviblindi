% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hompath.R
\name{make_snowman}
\alias{make_snowman}
\title{Create a synthetic snowman dataset}
\usage{
make_snowman(N)
}
\arguments{
\item{N}{number of points in a single segment of the snowman. (Since the bottom part of the snowman is cut off, total number of points will be around 2.5 times \code{N}).}
}
\value{
A matrix of coordinates. Rows correspond to points.
}
\description{
\code{make_snowman} outputs a matrix of coordinates in a space of dimension 2. When plotted, the points approximate the shape of a snowman, consisting of three circular segments. This is intended as test data for topological analysis.
See vignette of package \code{tviblindi} for a demonstration of use.
}
