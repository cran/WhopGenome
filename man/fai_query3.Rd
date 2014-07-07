\name{fai_query3}
\alias{fai_query3}
\title{
Extract a part of a FASTA sequence.
}
\description{
Return a part of a FASTA sequence.
}
\usage{
fai_query3( faifh, regionstring, resultstring )
}
\arguments{
\item{faifh}{FAIhandle as returned by fai_open}
\item{resultstring}{End position of the subsequence to extract}
\item{regionstring}{String of the form sequencename:beginpos-endpos e.g. "MTRR1mouse:20-40" specifying the sequence and region}
}
\details{
Use .Call("FAI_query3", faifh, regionstring, resultstring ) to eliminate the overhead of using the R wrapper function.
NOTE: The numbers in the function names fai_query3 and fai_query5 functions are not related to 3' / 5' strands but to the number of arguments.
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
fai_query3( faifh , "1:100-200", result )
print( result )
fai_close( faifh )

}
\author{
Ulrich Wittelsbuerger
}
\seealso{
fai_open
}
