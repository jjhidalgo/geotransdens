       SUBROUTINE ASSEMBLE_INTO_SHUFFLED_BAND_NO_SYM
     ;(IA_COLS    ,IA_ROWS   ,IDSC_COLS    ,IDSC_ROWS   ,I_TYPE_A   
     ;,A          ,A_DSC     ,FACTOR       ,KXX         ,LBLOCK_NR   
     &,LMXNDL     ,LNNDEL    ,NUMEL        ,NBAND) 
  

********************************************************************************
*
* PURPOSE  To assemble a matrix into another one.
*
*
* DESCRIPTION This is a decision routine. It call the appropriate routine
*			according to the type of origin matrix being the destiny matrix a
*			banded non symmetric one.
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A						Origin matrix
*  A_DSC					Destiny matrix.
*  iad					Columns. Index array of WatSolve.
*  iadd					Diagonal. Index array of WatSolve.
*  iadn					Number of columns. Index array of WatSolve.
*
* EXTERNAL VARIABLES: SCALARS
*
*  FACTOR					Factor that multiplies the matrix to be assembled.
*						For time-weighted schemes pourposes (theta-scheme).
*
*  IOMET					Solving method.
*							0 --> Picard iterations or lineal problem.
*							1 --> Newton's method.
*  I_TYPE_A               Type of A matrix.
*                         The possible types of matrices are:
*                             1 -->   nodewise Vector.
*                             2 -->   elementwise vector
*                             3 -->   derivative type matrix
*                             4 -->   Full matrix (elementwise).
*                             5 -->   Symmetric matrix (elementwise).
*                             6 -->   Symmetric matrix without diagonal (elementwise).
*                                     (It is supposed tha the sum of the terms
*                                     out of the diagonal is equal to the diagonal with
*                                     the opposite sign).
*                             7 -->   Symmetric banded matrix.
*                             8 -->   Non symmetric banded matrix.
*                             9 -->   Sparse matrix (as the one used by WatSolve).
*  I_TYPE_A_DSC			Type of A_DSC matrix.
*  IA_COLS				Number of columns (first dimension) of A matrix.
*  IA_ROWS				Number of rows (second dimension) of A matrix.
*  IDSC_COLS				Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS				Number of rowss (second dimension) of A_DSC matrix.
*  maxnb					Maximun adjacents nodes. Watsolve parameter.
*  maxnn					Maximun number of unknowns. Watsolve parameter.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*	ASSEMBLE_VECTOR_INTO_SHUFFLED_BAND
*	ASSEMBLE_FULL_INTO_SHUFFLED_BAND_NO_SYM
*
* HISTORY: First coding: JHG (7-2003)
*
********************************************************************************



		IMPLICIT NONE
	
		INTEGER*4:: I_TYPE_A,IA_COLS, IA_ROWS, IDSC_COLS
     &			   ,IDSC_ROWS,LMXNDL,LBLOCK_NR,NBAND,NUMEL
     &              

          INTEGER*4:: KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

		REAL*8::  A(IA_ROWS,IA_COLS),
     &			  A_DSC(IDSC_ROWS,IDSC_COLS), 
     &			  FACTOR




		SELECT CASE (I_TYPE_A)

			CASE(1)

                   CALL ASSEMBLE_VECTOR_NOD_INTO_SHUFFLED_BAND
     &                 (FACTOR     ,IA_ROWS     ,IDSC_COLS   ,IDSC_ROWS
     &                 ,LBLOCK_NR  ,NBAND       ,A           ,A_DSC)

			CASE (2)

				CALL ASSEMBLE_VEC_ELEM_INTO_SHUFFLED_BAND
     &                (FACTOR    ,IA_COLS     ,IA_ROWS    ,IDSC_COLS
     &                ,IDSC_ROWS ,LBLOCK_NR   ,LMXNDL     ,NBAND
     &                ,A         ,A_DSC       ,KXX        ,LNNDEL)

			CASE (3)

				CALL ASSEMBLE_DERIVATIVE_INTO_SHUFFLED_BAND_NO_SYM
     &          (A         ,A_DSC     ,FACTOR   ,IA_COLS   ,IA_ROWS
     &          ,IDSC_COLS ,IDSC_ROWS ,KXX      ,LBLOCK_NR ,LNNDEL
     &          ,NBAND)

			CASE (4)

				CALL ASSEMBLE_FULL_INTO_SHUFFLED_BAND_NO_SYM
     ;(FACTOR    ,IA_COLS   ,IA_ROWS    ,IDSC_COLS  ,IDSC_ROWS
     ;,LMXNDL    ,NBAND     ,LBLOCK_NR  ,A          ,A_DSC
     ;,LNNDEL    ,KXX)

              CASE (5)
	            
	            CALL ASSEMBLE_SYM_INTO_SHUFFLED_BAND_NO_SYM
     ;(FACTOR    ,IA_COLS   ,IA_ROWS    ,IDSC_COLS  ,IDSC_ROWS
     ;,LMXNDL    ,NBAND     ,LBLOCK_NR  ,A          ,A_DSC
     ;,LNNDEL    ,KXX)


              CASE (6)

				CALL ASSEMBLE_SYM_NO_DIAG_INTO_SHUFFLED_BAND_NO_SYM
     &                (A        ,A_DSC    ,FACTOR   ,IA_COLS  ,IA_ROWS
     &                ,IDSC_COLS,IDSC_ROWS,LNNDEL   ,KXX      ,LMXNDL
     &                ,NBAND    ,LBLOCK_NR)

              CASE (10)

				CALL ASSEMBLE_CROSSTERM_INTO_SHUFFLED_BAND_NO_SYM
     ;(FACTOR    ,IA_COLS   ,IA_ROWS    ,IDSC_COLS  ,IDSC_ROWS
     ;,LMXNDL    ,NBAND     ,LBLOCK_NR  ,A          ,A_DSC
     ;,LNNDEL    ,KXX)
			CASE DEFAULT

				!CALL ERROR 
				!Nonsense cases.

		END SELECT

       RETURN
       END
