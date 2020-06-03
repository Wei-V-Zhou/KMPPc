
\name{dataprocess}
\alias{dataprocess}
\title{dataprocess}
\usage{
dataprocess(clinicaldata = NULL, exprSet = NULL, x_axis = NULL, y_axis = NULL)
}
\arguments{
\item{clinicaldata}{The data can be loaded by the results of dataexplore. You can
load("survival_inputdata.Rdata") to get both of datasets.}

\item{exprSet}{The data can be loaded by the results of dataexplore. You can
load("survival_inputdata.Rdata") to get both of datasets.}

\item{x_axis}{This argument is the x axis for time.}

\item{y_axis}{This argument is the y axis for status.}
}

\description{
dataprocess for TCGA
}
\details{
By default, this function you can filter some dirty data to prepare for kmplot.
}
\examples{
## Process the TCGA data automatically
dataprocess()
## Plot kmfit using "DFS_MONTHS" and "DFS_STATUS"
dataprocess(x_axis = "DFS_MONTHS", y_axis = "DFS_STATUS")
}
\references{
Giuseppe Cardillo. (1993). Kmplot plot the kaplan-meier estimation of the survival function. Geotectonica Et Metallogenia.
}
\author{
Wei Zhou <247328181@qq.com>
}
