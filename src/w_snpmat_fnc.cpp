/*
**
**		WhopGen
**
**		WHOle genome population genetics with PopGENome
**
**
**		Whopgen - Reading data into matrices
**
**
**

	TODO:



**
*/

//*
//*			INCLUDES
//*

#include	"w_common.h"
#include	"functor_common.h"




	//

using namespace std;


//*
//*			DEFINES
//*



//*
//*			STRUCTS
//*


//*
//*			CLASSES
//*



//*
//*			DATA
//*




//*
//*			EXTERNS
//*

extern	SEXP	nucleotide_A;
extern	SEXP	nucleotide_C;
extern	SEXP	nucleotide_G;
extern	SEXP	nucleotide_T;
extern	SEXP	nucleotide_N;


//*
//*			CODE
//*



/*! Reads SNP lines from the VCF and processes them using the given functor
**
**	- 	Writes data directly into the given R-matrix in <mat>
**
**
*/
EXPORT SEXP helper_read_intmatrix_functored( SEXP vcfptr, MatrixLoaderBaseClass &fn, SEXP mat )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("helper_read_intmatrix_functored :: Parameter not a VCFhandle EXTPTR!\n");
		RBool::False();
	}
	fn.vcf = f;
	
	df2("A\n");
	
	//	any samples selected ?
	//
	if( f->num_wanted_samples < 1 )
	{
		Rprintf("helper_read_intmatrix_functored :: No samples selected!\n");
		RBool::False();
	}
	
	df2("B\n");

	//	without a FORMAT field, it is not clear how to discover the genotype subfield per sample!
	//
	unsigned int	samplefieldindex = f->getFirstSampleFieldIndex();
	if( samplefieldindex <= FORMAT )
	{
		Rprintf("helper_read_intmatrix_functored :: VCF does not have a FORMAT field!\n");
		RBool::False();
	}
	
	df2("C\n");
	
	//	is mat a matrix?
	//
	RMatrix m(mat);
	if( false == m.isValid() )
	{
		Rprintf("helper_read_intmatrix_functored :: Second parameter is not an integer matrix!\n");
		RBool::False();
	}
	
	df2("D\n");

	//	biallelic matrices are integer-typed
	//		save 4 bytes per entry over standard 'double' datatype
	//
	if( m.getType() != INTSXP )
	{
		df0("helper_read_intmatrix_functored :: Parameter not a integer matrix!\n");
		return RBool::False();
	}
	
	df2("E\n");
	
	
	//	enough rows to hold all samples?
	//
	unsigned int	nrow = m.numRows();
	if( f->num_wanted_samples > (unsigned)nrow )
	{
		Rprintf("helper_read_intmatrix_functored :: %d samples selected but matrix offers only rows for %d samples!\n",f->num_wanted_samples,nrow);
		return  RBool::False();
	}
	
	//
	//
	//
	int				nonbialcols=0,
					bialcols=0
						;
	unsigned int	ncol = m.numCols();
	int*		 	ptr = m.getIntPtr();
	
	//
	//
	if( ptr == 0 )
	{
		Rprintf("helper_read_intmatrix_functored :: ERROR : Could not get access to the matrix in form of an int*!\n");
		return  RBool::False();
	}
	
	//
	char			textbuffer[256];		//for sprintf()-ing SNP positions for matrix-column names
	char			*fieldptr=0;
	unsigned int	per_column = 0;		//vars here to find out what to clear
	unsigned int	per_row = 0;		//	when too little data exists
	unsigned int	column_stepsize = nrow;
	int				snppos=-1;
	
	SEXP colnamvec = m.getColNames();
	if( R_NilValue == colnamvec )
		Rprintf("helper_read_intmatrix_functored :: WARNING : matrix has no column names vector! not setting SNP positions in matrix!\n");

	//
	//
	for( ; per_column < ncol ; per_column ++ )
	{
		
		//-
		//
		//	Get a valid biallelic SNP line from the file
		//
		//-

		//char	*refptr=0;
		//char	*altptr=0;
		fn.refptr=0;
		fn.altptr=0;
		if( false == fn.findNextLine( f, &fn.refptr, &fn.altptr ) )
		{
			//df1("No more lines!\n");
			//break;
			return RBool::False();
		}

		//--------------
		//
		//	POSTCONDITION : got a line from VCF and it is a valid SNP
		//
		//--------------
		
		//
		fieldptr = (char*)f->getFieldPtr( POS );
		if( 0 == fieldptr )
		{
			df1("POS field not found!\n");
			//break;
			return RBool::False();
		}
		
		//
		snppos = atoi( fieldptr );
		if( 0 >= snppos )
		{
			df1("POS field specifies 0 or negative value!\n");
			//break;
			return RBool::False();
		}

		//--------------
		//
		//	POSTCONDITIONS :
		//		- got a line from VCF and it is a valid SNP
		//		- successfully identified the GT subfield
		//
		//--------------
		
		
		//	check whether filtering needs to be done
		//
		if( /*bFilterActivated && */(false == filterLine( f )) )
		{
			df1("helper_read_intmatrix_functored :: Line (pos %d)has been filtered away\n",snppos);
			continue;
		}
		
		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-
		fn.bHadRef = false;
		fn.bHadAlt = false;
		
		//
		for( per_row = 0; per_row < f->num_wanted_samples ; per_row ++ )
		{
			//	- get field of sample
			//
			fieldptr = (char*)f->getFieldPtr( f->wanted_samples[per_row] );
			if( fieldptr == 0 )
			{
				Rprintf("helper_read_intmatrix_functored :: ERROR when trying to get sample %d (matrix row %d) in file!\n",f->wanted_samples[per_row],per_row);
				Rprintf("	per_row =%d\nwanted_sample[per_row]=%d\n",per_row, f->wanted_samples[per_row] );
				Rprintf("	baseindex=%d, field = %d\n",samplefieldindex, (samplefieldindex + f->wanted_samples[per_row]) );
				Rprintf("	numparsedfields=%d\n",f->numParsedFields());
				return RBool::False();
			}
			
			/*
			**	- parse the GT subfield ( regexp: [0-9]+[\/\|][0-9]+ )
			*/
			if( false == fn.processSampleGTi( fieldptr, &ptr[per_row] ) )
			{
				Rprintf("(!!) Error during sample processing (pos=%d,matrix col=%d,row=%d): '%s'\n",snppos,per_column,per_row, &fn.errormessage[0] );
				//break;
				return RBool::False();
			}

			//
		}//...for( all rows == all samples )
		
		
		//----------------------------------------
		//
		//	clear any unused matrix rows of this column, if it had both alleles
		//
		if( true==fn.bHadAlt && true==fn.bHadRef )
		{
			
			//	fill unused matrix elements in current column
			//
			for( ; per_row < nrow ; per_row ++ )
			{
				ptr[per_row] = -2;
			}
			
			//
			ptr += column_stepsize;
			
			bialcols++;
			
			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != colnamvec )
			{
				snprintf(textbuffer,sizeof(textbuffer)-2,"%d",snppos);
				SET_STRING_ELT( colnamvec, per_column, mkChar(textbuffer) );
			}
		
			//FIXME any rows beyond the ones we need : zero them out here
			//	or as a separate step after this double for-loop
		}
		//
		//	if the last VCF-line's individuals did not provide both alleles
		//		reuse current column and try the next line
		//
		else
		{

			if( R_NilValue != colnamvec )
			{
				SET_STRING_ELT( colnamvec, per_column, fn.EmptyMatrixColumnNameSEXP );
			}
			
			per_column--;	//re-use the matrix column for the next SNP
			nonbialcols++;	//count how many non-biallelic columns we found
		}
		//
		//
		//
		//----------------------------------------
		
		//
	}//...for( each column in the matrix )

	return ( (0<per_column)?RBool::True():RBool::False() );
}


/*! Reads SNP lines from the VCF and processes them using the given functor
**
**	- 	Writes data directly into the given R-matrix in <mat>
**
**
*/
EXPORT SEXP helper_read_strmatrix_functored( SEXP vcfptr, MatrixLoaderBaseClass &fn, SEXP mat )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("helper_read_strmatrix_functored :: Parameter not a VCFhandle EXTPTR!\n");
		RBool::False();//return R_NilValue;
	}
	fn.vcf = f;
	
	df2("A\n");
	
	//	any samples selected ?
	//
	if( f->num_wanted_samples < 1 )
	{
		Rprintf("helper_read_strmatrix_functored :: No samples selected!\n");
		return RBool::False();
	}
	
	df2("B\n");

	//	without a FORMAT field, it is not clear how to discover the genotype subfield per sample!
	//
	unsigned int	samplefieldindex = f->getFirstSampleFieldIndex();
	if( samplefieldindex <= FORMAT )
	{
		Rprintf("helper_read_strmatrix_functored :: VCF does not have a FORMAT field!\n");
		return RBool::False();
	}
	
	df2("C\n");
	
	//	is mat a matrix?
	//
	RMatrix m(mat);
	if( false == m.isValid() )
	{
		Rprintf("helper_read_strmatrix_functored :: Second parameter is not an integer matrix!\n");
		return RBool::False();
	}
	
	df2("D\n");

	//	biallelic matrices are integer-typed
	//		save 4 bytes per entry over standard 'double' datatype
	//
	if( m.getType() != STRSXP )
	{
		df0("helper_read_strmatrix_functored :: Parameter not a string matrix!\n");
		return RBool::False();
	}
	
	df2("E\n");
	
	
	//	enough rows to hold all samples?
	//
	unsigned int	nrow = m.numRows();
	if( f->num_wanted_samples > (unsigned)nrow )
	{
		Rprintf("helper_read_strmatrix_functored :: %d samples selected but matrix offers only rows for %d samples!\n",f->num_wanted_samples,nrow);
		return RBool::False();
	}
	
	//
	//
	//
	int				nonbialcols=0,
					bialcols=0
						;
	unsigned int	ncol = m.numCols();
	SEXP*		 	ptr = m.getStrPtr();
	
	//
	//
	if( ptr == 0 )
	{
		Rprintf("helper_read_strmatrix_functored :: ERROR : Could not get access to the matrix in form of an SEXP!\n");
		return  RBool::False();
	}
	
	//
	char			textbuffer[256];		//for sprintf()-ing SNP positions for matrix-column names
	char			*fieldptr=0;
	unsigned int	per_column = 0;		//vars here to find out what to clear
	unsigned int	per_row = 0;		//	when too little data exists
	unsigned int	column_stepsize = nrow;
	int				snppos=-1;
	
	SEXP colnamvec = m.getColNames();
//	if( R_NilValue == colnamvec )
//		Rprintf("helper_read_intmatrix_functored :: WARNING : matrix has no column names vector! not setting SNP positions in matrix!\n");

	//
	//
	for( ; per_column < ncol ; per_column ++ )
	{
		
		//-
		//
		//	Get a valid biallelic SNP line from the file
		//
		//-

		//char	*refptr=0;
		//char	*altptr=0;
		fn.refptr=0;
		fn.altptr=0;
		if( false == fn.findNextLine( f, &fn.refptr, &fn.altptr ) )
		{
			//df1("No more lines!\n");
			//break;
			return RBool::False();
		}

		//--------------
		//
		//	POSTCONDITION : got a line from VCF and it is a valid SNP
		//
		//--------------
		
		//
		fieldptr = (char*)f->getFieldPtr( POS );
		if( 0 == fieldptr )
		{
			df1("POS field not found!\n");
			//break;
			return RBool::False();
		}
		
		//
		snppos = atoi( fieldptr );
		if( 0 >= snppos )
		{
			df1("POS field specifies 0 or negative value!\n");
			//break;
			return RBool::False();
		}

		//--------------
		//
		//	POSTCONDITIONS :
		//		- got a line from VCF and it is a valid SNP
		//		- successfully identified the GT subfield
		//
		//--------------
		
		
		//	check whether filtering needs to be done
		//
		if( /*bFilterActivated && */(false == filterLine( f )) )
		{
			df1("helper_read_strmatrix_functored :: Line (pos %d)has been filtered away\n",snppos);
			continue;
		}
		
		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-
		fn.bHadRef = false;
		fn.bHadAlt = false;
		
		//
		for( per_row = 0; per_row < f->num_wanted_samples ; per_row ++ )
		{
			//	- get field of sample
			//
			fieldptr = (char*)f->getFieldPtr( f->wanted_samples[per_row] );
			if( fieldptr == 0 )
			{
				Rprintf("helper_read_strmatrix_functored :: ERROR when trying to get sample %d (matrix row %d) in file!\n",f->wanted_samples[per_row],per_row);
				Rprintf("	per_row =%d\nwanted_sample[per_row]=%d\n",per_row, f->wanted_samples[per_row] );
				Rprintf("	baseindex=%d, field = %d\n",samplefieldindex, (samplefieldindex + f->wanted_samples[per_row]) );
				Rprintf("	numparsedfields=%d\n",f->numParsedFields());
				return RBool::False();
			}
			
			/*
			**	- parse the GT subfield ( regexp: [0-9]+[\/\|][0-9]+ )
			*/
			if( false == fn.processSampleGTs( fieldptr, &ptr[per_row] ) )
			{
				Rprintf("(!!) Error during sample processing (pos=%d,matrix col=%d,row=%d): '%s'\n",snppos,per_column,per_row, &fn.errormessage[0] );
				//break;
				return RBool::False();
			}

			//
		}//...for( all rows == all samples )
		
		
		//----------------------------------------
		//
		//	clear any unused matrix rows of this column, if it had both alleles
		//
		if( true==fn.bHadAlt && true==fn.bHadRef )
		{
			
			//
			for( ; per_row < nrow ; per_row ++ )
			{
				ptr[per_row] = fn.EmptyMatrixCellSEXP;
			}
			
			//
			ptr += column_stepsize;
			
			bialcols++;
			
			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != colnamvec )
			{
				snprintf(textbuffer,sizeof(textbuffer)-2,"%d",snppos);
				SET_STRING_ELT( colnamvec, per_column, mkChar(textbuffer) );
			}
		
			//FIXME any rows beyond the ones we need : zero them out here
			//	or as a separate step after this double for-loop
		}
		//
		//	if the last VCF-line's individuals did not provide both alleles
		//		reuse current column and try the next line
		//
		else
		{

			if( R_NilValue != colnamvec )
			{
				SET_STRING_ELT( colnamvec, per_column, fn.EmptyMatrixColumnNameSEXP );
			}
			
			per_column--;	//re-use the matrix column for the next SNP
			nonbialcols++;	//count how many non-biallelic columns we found
		}
		//
		//
		//
		//----------------------------------------
		
		//
	}//...for( each column in the matrix )

	return ( (0<per_column)?RBool::True():RBool::False() );
}







/*!	write into a INT matrix biallelic SNPS as 0 (=only REF/REF present) or 1 (=ALT/ALT or REF/ALT or ALT/REF present)
**
**
**
*/
EXPORT SEXP read_snp_diplo_bial_int_altpresence( SEXP vcfptr, SEXP mat )
{
	//
	//
	class MatrixLoader_DiploBialSNP : public MatrixLoaderBaseClass
	{
	public:

		//-
		//-		VIRTUAL METHODS
		//-
		
		//!
		virtual		bool	findNextLine( vcff * f, char** refptr, char** altptr )
		{
			return findNextLine_BiallelicSNP( f, refptr, altptr );
		}
		
		
		virtual		bool	processSampleGTs( const char * samplePtr, SEXP*ptr )
		{
			return false;
		}
		
		//!
		virtual		bool	processSampleGTi( const char * samplePtr, int*ptr )
		{
						//
			//
			const char* fcopy = samplePtr;
	//			Rprintf("%d: %c,%c,%c //",snppos,fcopy[0],fcopy[1],fcopy[2]);
			int left_allele = fcopy[0]-'0';//parseDecInt( fcopy );
			
			//make sure there are at least two alleles, regardless of whether phased or unphased
			//
			if( fcopy[1] != '/' && fcopy[1] != '|' )
			{
				//Rprintf("ERROR : unexpected character '%c' in Genotype field at position %d\n",fcopy[1],snppos);
				snprintf( errormessage, sizeof(errormessage)/sizeof(errormessage[0]), "ERROR : unexpected character '%c' in Genotype (GT) field\n", fcopy[1] );
				df0("	=> Syntax error in GT field (%s)!\n",samplePtr);
				return false;
			}
			
			//get second allele
			//
			int right_allele = fcopy[2]-'0';//parseDecInt( fcopy );
			
			//make sure we're handling a truly diploid species,
			//	i.e. exactly two alleles, no more | or / which divide the alleles in the genotype field)
			//
				/*\t separates samples*/
				/*0-byte ends last sample in the line (no \n ?? ) */
				/*: separates genotype from other per-sample data*/
			if( fcopy[3] != '\t' && fcopy[3] != 0 && fcopy[3] != ':' )
			{
				//Rprintf(
				snprintf( errormessage, sizeof(errormessage)/sizeof(errormessage[0]), "Syntax error in GT field (%s)!\n",samplePtr);
				return false;
			}

	//			Rprintf("left=%d, right=%d\n",left_allele,right_allele);

			//
			//
			if( (left_allele==1)||(right_allele==1) )
			{
				bHadAlt = true;	//memorise that we got samples with alternate alleles
				*ptr = 1;//ptr[per_row] = 1;
			}
			else
			{
				bHadRef = true;	//memorise that we got samples with reference alleles
				*ptr = 0;//ptr[per_row] = 0;
			}
			
			return true;

		}//...processSampleGT


	} f;//...class MatrixLoader_DiploBialSNP
	
	
	return helper_read_intmatrix_functored( vcfptr, f, mat );
}














/*!	Write into a SEXP matrix biallelic SNPS as TC,CT,TT,CC.... allele pairs
**
**
**
*/
EXPORT SEXP read_snp_diplo_bial_str_allelechars( SEXP vcfptr, SEXP mat )
{
	
	class MatrixLoader_DiploBialSNPStrAlleles : public MatrixLoaderBaseClass
	{
	public:
	
		MatrixLoader_DiploBialSNPStrAlleles()
		{
			EmptyMatrixCellSEXP = mkChar("EMPTY");
		}

		//-
		//-		VIRTUAL METHODS
		//-
		
		//!
		virtual		bool	findNextLine( vcff * f, char** refptr, char** altptr )
		{
			return findNextLine_BiallelicSNP( f, refptr, altptr );
		}
		
		
		//!
		virtual		bool	processSampleGTi( const char * samplePtr, int*ptr )
		{
			return false;
		}
		
		//!
		virtual		bool	processSampleGTs( const char * samplePtr, SEXP*ptr )
		{
						//
			//
			const char* fcopy = samplePtr;
	//			Rprintf("%d: %c,%c,%c //",snppos,fcopy[0],fcopy[1],fcopy[2]);
			int left_allele = fcopy[0]-'0';
			
			//make sure there are at least two alleles, regardless of whether phased or unphased
			//
			if( fcopy[1] != '/' && fcopy[1] != '|' )
			{
				//Rprintf("ERROR : unexpected character '%c' in Genotype field at position %d\n",fcopy[1],snppos);
				snprintf( errormessage, sizeof(errormessage)/sizeof(errormessage[0]), "ERROR : unexpected character '%c' in Genotype (GT) field\n", fcopy[1] );
				df0("	=> Syntax error in GT field (%s)!\n",samplePtr);
				return false;
			}
			
			//get second allele
			//
			int right_allele = fcopy[2]-'0';
			
			//make sure we're handling a truly diploid species,
			//	i.e. exactly two alleles, no more | or / which divide the alleles in the genotype field)
			//
				/*\t separates samples*/
				/*0-byte ends last sample in the line (no \n ?? ) */
				/*: separates genotype from other per-sample data*/
			if( fcopy[3] != '\t' && fcopy[3] != 0 && fcopy[3] != ':' )
			{
				//Rprintf(
				snprintf( errormessage, sizeof(errormessage)/sizeof(errormessage[0]), "Syntax error in GT field (%s)!\n",samplePtr);
				return false;
			}

	//			Rprintf("left=%d, right=%d\n",left_allele,right_allele);
			char	all[3];
			all[0] = left_allele ? *refptr : *altptr;//fcopy[0];
			all[1] = right_allele ? *refptr : *altptr;//fcopy[2];
			all[2] = 0;
			
			*ptr = mkChar( (const char*)&all );

			//
			//
			if( (left_allele==1)||(right_allele==1) )
			{
				bHadAlt = true;	//memorise that we got samples with alternate alleles
			}
			else
			{
				bHadRef = true;	//memorise that we got samples with reference alleles
			}
			
			return true;

		}//...processSampleGT


	} f ;
	//...class MatrixLoader_DiploBialSNPStrAlleles
	
	return helper_read_strmatrix_functored( vcfptr, f, mat );
}
















/*!	Write into a SEXP matrix biallelic SNPS as 00,01,10,11.. allele pairs (0=ref, 1=alt)
**
**
**
*/
EXPORT SEXP read_snp_diplo_bial_str_01( SEXP vcfptr, SEXP mat )
{
	
	class MatrixLoader_DiploBialSNPStr01 : public MatrixLoaderBaseClass
	{
	public:
	
		MatrixLoader_DiploBialSNPStr01()
		{
			EmptyMatrixCellSEXP = mkChar("EMPTY");
		}

		//-
		//-		VIRTUAL METHODS
		//-
		
		//!
		virtual		bool	findNextLine( vcff * f, char** refptr, char** altptr )
		{
			return findNextLine_BiallelicSNP( f, refptr, altptr );
		}
		
		
		//!
		virtual		bool	processSampleGTi( const char * samplePtr, int*ptr )
		{
			return false;
		}
		
		//!
		virtual		bool	processSampleGTs( const char * samplePtr, SEXP*ptr )
		{
						//
			//
			const char* fcopy = samplePtr;
	//			Rprintf("%d: %c,%c,%c //",snppos,fcopy[0],fcopy[1],fcopy[2]);
			int left_allele = fcopy[0]-'0';
			
			//make sure there are at least two alleles, regardless of whether phased or unphased
			//
			if( fcopy[1] != '/' && fcopy[1] != '|' )
			{
				//Rprintf("ERROR : unexpected character '%c' in Genotype field at position %d\n",fcopy[1],snppos);
				snprintf( errormessage, sizeof(errormessage)/sizeof(errormessage[0]), "ERROR : unexpected character '%c' in Genotype (GT) field\n", fcopy[1] );
				df0("	=> Syntax error in GT field (%s)!\n",samplePtr);
				return false;
			}
			
			//get second allele
			//
			int right_allele = fcopy[2]-'0';
			
			//make sure we're handling a truly diploid species,
			//	i.e. exactly two alleles, no more | or / which divide the alleles in the genotype field)
			//
				/*\t separates samples*/
				/*0-byte ends last sample in the line (no \n ?? ) */
				/*: separates genotype from other per-sample data*/
			if( fcopy[3] != '\t' && fcopy[3] != 0 && fcopy[3] != ':' )
			{
				//Rprintf(
				snprintf( errormessage, sizeof(errormessage)/sizeof(errormessage[0]), "Syntax error in GT field (%s)!\n",samplePtr);
				return false;
			}

	//			Rprintf("left=%d, right=%d\n",left_allele,right_allele);
			char	all[3];
			all[0] = fcopy[0];
			all[1] = fcopy[2];
			all[2] = 0;
			
			*ptr = mkChar( (const char*)&all );

			//
			//
			if( (left_allele==1)||(right_allele==1) )
			{
				bHadAlt = true;	//memorise that we got samples with alternate alleles
			}
			else
			{
				bHadRef = true;	//memorise that we got samples with reference alleles
			}
			
			return true;

		}//...processSampleGT


	} f ;
	//...class MatrixLoader_DiploBialSNPStr01
	
	return helper_read_strmatrix_functored( vcfptr, f, mat );

}

//nucleotide_mapping














extern	char	nucleotide_mapping[];


/*!	Write into a SEXP matrix biallelic SNPS as 00,01,10,11.. allele pairs (0=ref, 1=alt)
**
**
**
*/
EXPORT SEXP read_snp_diplo_bial_int_nuclcodes( SEXP vcfptr, SEXP mat )
{
	
	class MatrixLoader_DiploBialSNP_Int_NucleotideCodes : public MatrixLoaderBaseClass
	{
	public:
	
		MatrixLoader_DiploBialSNP_Int_NucleotideCodes()
		{
			EmptyMatrixCellSEXP = R_NilValue;//mkChar("55");
			EmptyMatrixCellInteger = 66;
		}

		//-
		//-		VIRTUAL METHODS
		//-
		
		//!
		virtual		bool	findNextLine( vcff * f, char** refptr, char** altptr )
		{
			return findNextLine_BiallelicSNP( f, refptr, altptr );
		}
		
		
		//!
		virtual		bool	processSampleGTs( const char * samplePtr, SEXP*ptr )
		{
			return false;
		}
		
		//!
		virtual		bool	processSampleGTi( const char * samplePtr, int*ptr )
		{
			//
			//
			const char* fcopy = samplePtr;
	//			Rprintf("%d: %c,%c,%c //",snppos,fcopy[0],fcopy[1],fcopy[2]);
			int left_allele = fcopy[0]-'0';
			
			//make sure there are at least two alleles, regardless of whether phased or unphased
			//
			if( fcopy[1] != '/' && fcopy[1] != '|' )	//FIXME : VCFv3.3 allows '\' also - obsolete since VCFv4.0
			{
				//Rprintf("ERROR : unexpected character '%c' in Genotype field at position %d\n",fcopy[1],snppos);
				snprintf( errormessage, sizeof(errormessage)/sizeof(errormessage[0]), "ERROR : unexpected character '%c' in Genotype (GT) field\n", fcopy[1] );
				df0("	=> Syntax error in GT field (%s)!\n",samplePtr);
				return false;
			}
			
			//get second allele
			//
			int right_allele = fcopy[2]-'0';
			
			//make sure we're handling a truly diploid species,
			//	i.e. exactly two alleles, no more | or / which divide the alleles in the genotype field)
			//
				/*\t separates samples*/
				/*0-byte ends last sample in the line (no \n ?? ) */
				/*: separates genotype from other per-sample data*/
			if( fcopy[3] != '\t' && fcopy[3] != 0 && fcopy[3] != ':' )
			{
				//Rprintf(
				snprintf( errormessage, sizeof(errormessage)/sizeof(errormessage[0]), "Syntax error in GT field (%s)!\n",samplePtr);
				return false;
			}

	//			Rprintf("left=%d, right=%d\n",left_allele,right_allele);
			*ptr = 
						( nucleotide_mapping[ int((0==left_allele)?*refptr:*altptr)&0xFF ] * 10 )
					+	nucleotide_mapping[ int((0==right_allele)?*refptr:*altptr)&0xFF ];

			//
			//
			if( (left_allele==1)||(right_allele==1) )
			{
				bHadAlt = true;	//memorise that we got samples with alternate alleles
			}
			else
			{
				bHadRef = true;	//memorise that we got samples with reference alleles
			}
			
			return true;

		}//...processSampleGT


	} f ;
	//...class MatrixLoader_DiploBialSNPStr01
	
	return helper_read_intmatrix_functored( vcfptr, f, mat );

}


















/*!	Write into a SEXP matrix diploid, biallelic SNP alleles as pairs of nucleotide-codes 1=T or U,2=C,3=G,4=A,5=N : 15, 22, 31, 24...
**
**
**
*/
EXPORT SEXP read_snp_diplo_bial_str_nuclcodes( SEXP vcfptr, SEXP mat )
{
	
	class MatrixLoader_DiploBialSNPStrNucleotideCodes : public MatrixLoaderBaseClass
	{
	public:
	
		MatrixLoader_DiploBialSNPStrNucleotideCodes()
		{
			EmptyMatrixCellSEXP = mkChar("55");
		}

		//-
		//-		VIRTUAL METHODS
		//-
		
		//!
		virtual		bool	findNextLine( vcff * f, char** refptr, char** altptr )
		{
			return findNextLine_BiallelicSNP( f, refptr, altptr );
		}
		
		
		//!
		virtual		bool	processSampleGTi( const char * samplePtr, int*ptr )
		{
			return false;
		}
		
		//!
		virtual		bool	processSampleGTs( const char * samplePtr, SEXP*ptr )
		{
						//
			//
			const char* fcopy = samplePtr;
	//			Rprintf("%d: %c,%c,%c //",snppos,fcopy[0],fcopy[1],fcopy[2]);
			int left_allele = fcopy[0]-'0';
			
			//make sure there are at least two alleles, regardless of whether phased or unphased
			//
			if( fcopy[1] != '/' && fcopy[1] != '|' )
			{
				//Rprintf("ERROR : unexpected character '%c' in Genotype field at position %d\n",fcopy[1],snppos);
				snprintf( errormessage, sizeof(errormessage)/sizeof(errormessage[0]), "ERROR : unexpected character '%c' in Genotype (GT) field\n", fcopy[1] );
				df0("	=> Syntax error in GT field (%s)!\n",samplePtr);
				return false;
			}
			
			//get second allele
			//
			int right_allele = fcopy[2]-'0';
			
			//make sure we're handling a truly diploid species,
			//	i.e. exactly two alleles, no more | or / which divide the alleles in the genotype field)
			//
				/*\t separates samples*/
				/*0-byte ends last sample in the line (no \n ?? ) */
				/*: separates genotype from other per-sample data*/
			if( fcopy[3] != '\t' && fcopy[3] != 0 && fcopy[3] != ':' )
			{
				snprintf( errormessage, sizeof(errormessage)/sizeof(errormessage[0]), "Syntax error in GT field (%s)!\n",samplePtr);
				return false;
			}

	//		Rprintf("left=%d, right=%d\n",left_allele,right_allele);
			char	all[3];
			all[0] = nucleotide_mapping[ int(left_allele?*refptr:*altptr)&0xFF ] + '0';
			all[1] = nucleotide_mapping[ int(right_allele?*refptr:*altptr)&0xFF ] + '0';
			all[2] = 0;
			
			*ptr = mkChar( (const char*)&all );	//TODO optimise this : use 5x5 matrix with once-initialised mkChar() entries: "11","12", "13",... / "21","22","23... ...

			//
			//
			if( (left_allele==1)||(right_allele==1) )
			{
				bHadAlt = true;	//memorise that we got samples with alternate alleles
			}
			else
			{
				bHadRef = true;	//memorise that we got samples with reference alleles
			}
			
			//snprintf( errormessage , "No error\n" );
			return true;

		}//...processSampleGT


	} f ;
	//...class MatrixLoader_DiploBialSNPStr01
	
	return helper_read_strmatrix_functored( vcfptr, f, mat );

}
