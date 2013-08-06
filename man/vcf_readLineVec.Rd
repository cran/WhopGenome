\name{vcf_readLineVec}
\alias{vcf_readLineVec}
\alias{vcf_readLineVecFiltered}
\title{
Read a line of data from the given VCF file and return the fields as vector elements
}
\description{
Read a line of data from the given VCF file and return the fields as vector elements
}
\usage{
vcf_readLineVec(vcffh)
vcf_readLineVecFiltered(vcffh)
}
\arguments{
  \item{vcffh}{VCF file handle}
}
\details{
The latter version applies filtering set up with vcf_addfilter.
Use .Call("VCF_readLineTSV", vcffh ) or .Call("VCF_readLineTSVFiltered", vcffh ) respectively
to eliminate the overhead of using the R wrapper function.
}
\value{
A vector where each element is a field from a line of data in the VCF
}
\author{
Ulrich Wittelsbuerger
}
\seealso{
vcf_addfilter,vcf_describefilters
}
