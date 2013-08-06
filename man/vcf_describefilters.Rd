\name{vcf_describefilters}
\alias{vcf_describefilters}
\title{
Prints description of current filter rules
}
\description{
Prints a better understandable description of the filter rules currently active for the given VCF file.
}
\usage{
vcf_describefilters(vcffh)
}
\arguments{
  \item{vcffh}{VCF file handle}
}
\details{
Use .Call("VCF_describeFilterConfig", filename ) to eliminate the overhead of using the R wrapper function.
Note the different naming of the library function!
}
\value{
None.
}
\author{
Ulrich Wittelsbuerger
}
