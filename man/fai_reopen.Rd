\name{fai_reopen}
\alias{fai_reopen}
\title{
Reopen a FAIhandle that has become stale.
}
\description{
Reopen a FAIhandle that has become stale, e.g. by restarting R or loading a workspace containing a FAIhandle variable.
}
\usage{
fai_reopen( faifh )
}
\arguments{
  \item{faifh}{A FAIhandle to a .fai-indexed FASTA file}
}
\details{
Use .Call("FAI_reopen", faifh ) to eliminate the slight overhead of using the R wrapper function.
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
