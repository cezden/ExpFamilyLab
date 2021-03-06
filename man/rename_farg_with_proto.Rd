% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplification.R
\name{rename_farg_with_proto}
\alias{rename_farg_with_proto}
\title{Change the name of 1st formal argument}
\usage{
rename_farg_with_proto(f, res.proto)
}
\arguments{
\item{f}{(function)}

\item{res.proto}{"prototype" function with desired formal name, e.g. \code{function(x) NULL}}
}
\value{
function \code{f} with resubstituted first argument
}
\description{
Change the name of 1st formal argument
}

