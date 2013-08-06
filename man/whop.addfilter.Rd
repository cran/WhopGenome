\name{vcf_addfilter}
\alias{vcf_addfilter}
\title{
Add a condition for SNP filtering from VCF files.
}
\description{
Add a condition for filtering SNPs based on the INFO column in a certain VCF file.
}
\usage{
vcf_addfilter(vcf, fieldnam, cmptype, cmpvalue1, cmpvalue2 = 0, action)
}
\arguments{
  \item{vcf}{VCF file handle}
  \item{fieldnam}{INFO field name}
  \item{cmptype}{Type of comparison to perform. See Details}
  \item{cmpvalue1}{Comparison reference value 1 or lower bound}
  \item{cmpvalue2}{Comparison reference value 2 or upper bound}
  \item{action}{Action to take if comparison matches : NOP, SKIP, KEEP or fails: SKIP_NOT, KEEP_NOT}
}
\details{
Comparison types:
	DOES_EXIST		cmpvalue1 is a key of a 'key=value' pair and exists in the INFO field's data. E.g. used to test for H2 meaning Hapmap2 membership of a SNP
	INT_CMP
	INT_CMP_OO
	INT_CMP_OC
	INT_CMP_CO
	INT_CMP_CC
	FLT_CMP
	FLT_CMP_OO
	FLT_CMP_OC
	FLT_CMP_CO
	FLT_CMP_CC
}
\value{
Success status
}
\author{
Ulrich Wittelsbuerger
}
\examples{

## The function is currently defined as
function (vcf, fieldnam, cmptype, cmpvalue1, cmpvalue2 = 0, action)
{
    fieldnam = toupper(fieldnam)
    cmptype = switch(EXPR = cmptype, DOES_EXIST = 0, INT_CMP = 1,
        INT_CMP_OO = 2, INT_CMP_OC = 3, INT_CMP_CO = 4, INT_CMP_CC = 5,
        FLT_CMP = 6, FLT_CMP_OO = 7, FLT_CMP_OC = 8, FLT_CMP_CO = 9,
        FLT_CMP_CC = 10, -1)
    if (cmptype == -1)
        stop("Illegal comparison type!")
    if (cmptype >= 1 && cmptype <= 5) {
        cmp1 = as.integer(cmpvalue1)
        cmp2 = as.integer(cmpvalue2)
    }
    else if (cmptype >= 7 && cmptype <= 10) {
        cmp1 = as.double(cmpvalue1)
        cmp2 = as.double(cmpvalue2)
    }
    if (cmptype >= 2 && cmptype <= 5) {
        if (is.na(cmp1) || is.na(cmp2))
            stop("cmpvalue1 and cmpvalue2 cannot be NA with this comparisong type!")
    }
    else if (cmptype >= 7 && cmptype <= 10) {
        if (is.na(cmp1) || is.na(cmp2))
            stop("cmpvalue1 and cmpvalue2 cannot be NA with this comparisong type!")
    }
    cmpaction = switch(action, NOP = 0, SKIP = 1, KEEP = 2, SKIP_NOT = 129,
        KEEP_NOT = 130, -1)
    if (cmpaction == -1)
        stop("Unknown action selected!")
    cat("CmpType=", cmptype, "\n")
    cat("ActionType=", cmpaction, "\n")
    cat("CmpValues=", cmp1, ",", cmp2, "\n")
  }
}
