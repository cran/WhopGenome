\name{vcf_clearfilters}
\alias{vcf_clearfilters}
\title{
Removes all filter steps.
}
\description{
Removes all active filters, no pre-filtering of returned lines will take place.
There is no function to undo this step.
}
\usage{
vcf_clearfilters(vcffh)
}
\arguments{
  \item{vcffh}{VCF file handle}
}
\details{
Use .Call("VCF_clearFilters", vcffh ) to eliminate the overhead of using the R wrapper function.
}
\value{
None.
}
\author{
Ulrich Wittelsbuerger
}
