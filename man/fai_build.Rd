\name{fai_build}
\alias{fai_build}
\title{
Build a .fai-index for the given FASTA file.
}
\description{
Build a .fai-index for the given FASTA file.
}
\usage{
fai_build( filename )
}
\arguments{
\item{filename}{Name of the FASTA file for which an index file should be built.}
}
\details{
Use .Call("FAI_build", filename ) to eliminate the overhead of using the R wrapper function.
}
\value{
TRUE if call succeeds, FALSE if it fails.
}
\author{
Ulrich Wittelsbuerger
}
\seealso{
fai_open
}
