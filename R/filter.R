###
###
###		Data filtering setup and management
###
###

#TODO : clearfilters, removefilter, describefilters, replacefilter

#
#	Add a filter
#
vcf_addfilter <- function( vcf, fieldnam, cmptype, cmpvalue1, cmpvalue2=0, action )
{
	#if( fieldnam == "FILTERPASS" )
	
	#1. add fieldnam to ruleset of vcf
	#	error if ruleset full; name might already be in there
	fieldnam = toupper(fieldnam)
	
	#2. set cmptype
	#
	cmptype = switch( EXPR = cmptype, DOES_EXIST = 0,
						INT_CMP=1,
						INT_CMP_OO=2,
						INT_CMP_OC=3,
						INT_CMP_CO=4,
						INT_CMP_CC=5,
						FLT_CMP=6,
						FLT_CMP_OO=7,
						FLT_CMP_OC=8,
						FLT_CMP_CO=9,
						FLT_CMP_CC=10,
						-1
					);
	if( cmptype == -1 )
		stop("Illegal comparison type!");

	#3. set cmp parameters
	#

		#
		#
	if( cmptype >= 1 && cmptype <= 5 )
	{
		cmp1 = as.integer( cmpvalue1 )
		cmp2 = as.integer( cmpvalue2 )
	}
	else if( cmptype >= 7 && cmptype <= 10 )
	{
		cmp1 = as.double( cmpvalue1 )
		cmp2 = as.double( cmpvalue2 )
	}

		#
		#
	if( cmptype >= 2 && cmptype <= 5 )
	{
		if( is.na(cmp1) || is.na( cmp2 ) )
			stop("cmpvalue1 and cmpvalue2 cannot be NA with this comparisong type!");
	}
	else if( cmptype >= 7 && cmptype <= 10 )
	{
		if( is.na(cmp1) || is.na( cmp2 ) )
			stop("cmpvalue1 and cmpvalue2 cannot be NA with this comparisong type!");
	}
	
	#4. set action on cmp-result
	#
	cmpaction = switch( action , NOP=0,SKIP=1,KEEP=2,SKIP_NOT=0x81,KEEP_NOT=0x82, -1 )
	if( cmpaction == -1 )
		stop("Unknown action selected!");

	
	#
	#
	cat("CmpType=",cmptype,"\n")
	cat("ActionType=",cmpaction,"\n")
	cat("CmpValues=",cmp1,",",cmp2,"\n")
	
	#	add rule parameters
	#
	return(
		.Call("VCF_addFilter",vcf,fieldnam,cmptype,cmpaction,cmp1,cmp2,PACKAGE="WhopGenome")
		)

	#
}#...addFilter

##
##
##
##
##
vcf_clearfilters <- function( vcffh )	.Call("VCF_clearFilters" , vcffh , PACKAGE="WhopGenome")

##
##
##
##
##
vcf_describefilters <- function( vcffh )	.Call("VCF_describeFilterConfig" , vcffh , PACKAGE="WhopGenome")


#
#TODO:
#
#	removefilter, replacefilter
#
#removefilter: vcffh, filterentrynum
#replacefilter: = addfilter + filternum-parameter
#
#
