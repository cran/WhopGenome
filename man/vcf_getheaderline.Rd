\name{vcf_getheaderline}
\alias{vcf_getheaderline}
\title{
Return one of the header lines of the VCF file
}
\description{
Return one of the header lines of the VCF file
}
\usage{
vcf_getheaderline(vcff, whichnum)
}
\arguments{
  \item{vcff}{VCF file handle}
  \item{whichnum}{Number of header line to retrieve}
}
\details{
Use .Call("VCF_getHeaderLine", vcff, whichnum ) to eliminate the overhead of using the R wrapper function.
}
\value{
A string containing the full header line.
}
\author{
Ulrich Wittelsbuerger
}
