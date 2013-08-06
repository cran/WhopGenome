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
The query3/query5 functions are not related to 3' / 5' strands.
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
