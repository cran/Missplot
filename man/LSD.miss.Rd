% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LSD.miss.R
\name{LSD.miss}
\alias{LSD.miss}
\title{Analysis of LSD when there is one missing observation}
\usage{
LSD.miss(data, row.miss, col.miss)
}
\arguments{
\item{data}{A data-frame containing response in the first column ,row number in the second column,
column number in the third column , treatment number in the fourth column corresponding to the response value.}

\item{row.miss}{Row number corresponding to the missing observation.}

\item{col.miss}{Column number corresponding to the missing observation.}
}
\value{
A data-frame containing x.hat , SSE x.hat , x_double.hat , SSE x_double.hat,
F statistics , p-value.
}
\description{
Analysis of LSD when there is one missing observation
}
\details{
In design of experiments in LSD setup if there is one missing observation
present in the design , we can use the function LSD.miss to estimate the missing observation for testing
the differential effects for the treatments. Here, we estimate the missing obsevation by
minimizing the SSE of the design.
}
\section{Credits}{
 Credits to Professor Surupa Chakraborty for building the theoritical concepts of Design of Experiment
and Professor Madhura Dasgupta for basic concepts for R.
}

\section{Remark}{
 Information on row number and column number corresponding to the missing observation
is to be known.
}

\examples{
#Observation corresponding to the second row and third column is missing in the data
data=data.frame(res=rnorm(16,35,20),row_no=rep(1:4,each=4),col_no=rep(1:4,times=4),
    treat=c(1,2,3,4,2,3,4,1,3,4,1,2,4,1,2,3))
LSD.miss(data,2,3)
}
\author{
Saheli Datta , Shantanu Nayek
}
