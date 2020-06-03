
\name{dataexplore}
\alias{dataexplore}
\title{dataexplore}
\usage{
dataexplore(cancerType, studyId = NULL, dataType = NULL)
}
\arguments{
\item{cancerType}{The cancertype you'd like to explore, including "BreastCancer",
"PancreaticCancer", "Glioma"}

\item{studyId}{This argument you can explore by yourself, because there is too
many researches or database you can take advantage of}

\item{dataType}{It provides "mrna", "CNA", "mutation" for more information to dig.}
}

\description{
dataexplore for TCGA
}
\details{
By default, this function first input cancerType such as "BreastCancer" to
explore the dataset whatever you like.
}
\examples{
## Explore the interested dataset by input the specified cancertype
   dataexplore(cancerType = "BreastCancer", dataType = "mrna")
   
## If you have interested datatype, you can input, such as mrna, CNA, mutation, etc
   dataexplore(cancerType = "BreastCancer", dataType = "mrna")
}
\references{
Giuseppe Cardillo. (1993). Kmplot plot the kaplan-meier estimation of the survival function. Geotectonica Et Metallogenia.
}
\author{
Wei Zhou <247328181@qq.com>
}
