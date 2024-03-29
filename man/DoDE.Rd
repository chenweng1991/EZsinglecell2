% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/General_functions.R
\name{DoDE}
\alias{DoDE}
\title{DoDE}
\usage{
DoDE(tri.dummy, cpcol, onlyoneSample = F, cpus = 16)
}
\arguments{
\item{tri.dummy}{this is generated from  DE.gettripple}

\item{cpcol}{the column in tri.dummy$info, the contents of which are used for iteratively compare with one another}

\item{onlyoneSample}{If true, regress out batch effect. Notice, there should be a "Sample" column in in tri.dummy$info that indicate sample or donor or batch}

\item{cpus}{a number of cpus being used for calculation, default is 16}
}
\value{
return a list that includes all DE result iteratively
}
\description{
This is the main function for calculating differentially expressed genes
}
\examples{
ROCKvsnorock.endo.de<-DoDE(ROCKvsnorock.endo.tri.dummy,"name",onlyoneSample=T,cpus=16)
}
