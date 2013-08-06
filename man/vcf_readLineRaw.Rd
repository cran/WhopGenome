\name{vcf_readLineRaw}
\alias{vcf_readLineRaw}
\alias{vcf_readLineRawFiltered}
\title{
Read a line of data from the given VCF file and return it as a string without postprocessing.
}
\description{
Read a line of data from the given VCF file and return it as a string without postprocessing.
}
\usage{
vcf_readLineRaw(vcffh,stri)
}
\arguments{
  \item{vcffh}{VCF file handle}
  \item{stri}{String variable into which to read the line}
}
\details{
Use .Call("VCF_readLineRaw", vcffh, stri ) and .Call("VCF_readLineRawFiltered", vcffh, stri ) respectively,
to eliminate the overhead of using the R wrapper function.
}
\value{
A raw string representing a line of data from the file
}
\author{
Ulrich Wittelsbuerger
}
