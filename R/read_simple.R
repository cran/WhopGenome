###
###
###		Reading
###
###
###


#
#	returns a data-frame where the rownames are the field names of the VCF
#
vcf_readLineDF <- function( vcffh )	return(	data.frame( .Call("VCF_readLineTSV",vcffh,PACKAGE="WhopGenome"), row.names=.Call("VCF_getFieldNames",vcffh,PACKAGE="WhopGenome") )	)


#
#	reads a line and returns a vector with the fields of that line
#
vcf_readLineVec <- function( vcffh )	return( .Call("VCF_readLineTSV",vcffh,PACKAGE="WhopGenome")	)
vcf_readLineVecFiltered <- function( vcffh )	.Call("VCF_readLineTSVFiltered",vcffh,PACKAGE="WhopGenome")


#
#	reads a line and stores in 'str'
#
vcf_readLineRaw <- function( vcffh, stri )	.Call("VCF_readLineRaw",vcffh,stri, PACKAGE="WhopGenome")
vcf_readLineRawFiltered <- function( vcffh, stri )	.Call("VCF_readLineRawFiltered",vcffh,stri,PACKAGE="WhopGenome")

#
#
#
vcf_countBiallelicSNPs <- function( vcffh )	.Call("VCF_countBiallelicSNPs",vcffh,PACKAGE="WhopGenome")
vcf_countSNPs <- function( vcffh )	.Call("VCF_countSNPs",vcffh,PACKAGE="WhopGenome")
vcf_parseNextSNP <- function( vcffh )	.Call("VCF_parseNextSNP",vcffh,PACKAGE="WhopGenome")
vcf_parseNextLine <- function( vcffh )	.Call("VCF_parseNextLine",vcffh,PACKAGE="WhopGenome")
vcf_getChrom <- function( vcffh )	.Call("VCF_getChrom",vcffh,PACKAGE="WhopGenome")
vcf_getPos <- function( vcffh )	.Call("VCF_getPos",vcffh,PACKAGE="WhopGenome")
vcf_getID <- function( vcffh )	.Call("VCF_getID",vcffh,PACKAGE="WhopGenome")
vcf_getRef <- function( vcffh )	.Call("VCF_getRef",vcffh,PACKAGE="WhopGenome")
vcf_getAlt <- function( vcffh )	.Call("VCF_getAlt",vcffh,PACKAGE="WhopGenome")
vcf_getQual <- function( vcffh )	.Call("VCF_getQual",vcffh,PACKAGE="WhopGenome")
vcf_getFilter <- function( vcffh )	.Call("VCF_getFilter",vcffh,PACKAGE="WhopGenome")
vcf_getInfo <- function( vcffh )	.Call("VCF_getInfo",vcffh,PACKAGE="WhopGenome")
vcf_getInfoField <- function( vcffh, fieldnam )	.Call("VCF_getInfoField",vcffh,fieldnam,PACKAGE="WhopGenome")
vcf_getSample <- function( vcffh, stridx )	.Call("VCF_getSample",vcffh,stridx,PACKAGE="WhopGenome")

#
#
#
VCF_read_snp_diplo_bial_int_altpresence <- function( vcffh, mat ) .Call("", vcffh, mat, PACKAGE="WhopGenome" )
VCF_read_snp_diplo_bial_str_allelechars <- function( vcffh, mat ) .Call("", vcffh, mat, PACKAGE="WhopGenome" )
VCF_read_snp_diplo_bial_str_01 <- function( vcffh, mat ) .Call("", vcffh, mat, PACKAGE="WhopGenome" )
VCF_read_snp_diplo_bial_int_nuclcodes <- function( vcffh, mat ) .Call("", vcffh, mat, PACKAGE="WhopGenome" )
VCF_read_snp_diplo_bial_str_nuclcodes <- function( vcffh, mat ) .Call("", vcffh, mat, PACKAGE="WhopGenome" )



#
#
#

