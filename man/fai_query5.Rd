\name{fai_query5}
\alias{fai_query5}
\title{
Extract a part of a FASTA sequence.
}
\description{
Return a part of the a FASTA sequence.
}
\usage{
fai_query5( faifh, sequencename, beginpos, endpos, resultstring )
}
\arguments{
\item{faifh}{FAIhandle as returned by fai_open}
\item{sequencename}{Identifier of a sequence in the fasta file}
\item{beginpos}{Start position of the subsequence to extract}
\item{endpos}{String variable into which to store the subsequence}
\item{resultstring}{End position of the subsequence to extract}
}
\details{
Use .Call("FAI_query5", faifh, sequencename, beginpos, endpos, resultstring ) to eliminate the overhead of using the R wrapper function.
}
\value{
TRUE if call succeeds, FALSE if it fails.
}
%XXXXXXXXXXXXXXX
\examples{

##
##	Example :
##
faifile  <- system.file("extdata", "ex.fasta", package = "WhopGenome")
faifh <- fai_open( faifile )
stopifnot( faifh != NULL )
result <- ""
fai_query5( faifh , "1", 9 , 20, result )
print( result )
fai_close( faifh )

}
\author{
Ulrich Wittelsbuerger
}
\seealso{
fai_open
}
