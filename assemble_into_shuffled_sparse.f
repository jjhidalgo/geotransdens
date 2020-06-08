       SUBROUTINE ASSEMBLE_INTO_SHUFFLED_SPARSE
     ;(FACTOR   ,IA_COLS   ,IA_ROWS    ,I_TYPE_A   ,LBLOCK_NR  
     ;,LMXNDL   ,NN        ,NB         ,NUMEL      ,A  
     ;,A_DSC    ,IAD_D     ,IADD_D     ,IADN_D     ,KXX
     ;,LNNDEL)


********************************************************************************
*
* PURPOSE  To assemble a matrix into another one.
*
*
* DESCRIPTION This is a decision routine. It call the appropriate routine
*			according to the type of origin matrix being the destiny matrix a
*			sparse one.
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
*  I_TYPE_A				Type of A matrix.
*						The possible types of matrices are:
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
*  IA_COLS				Number of columns (first dimension) of A matrix.
*  IA_ROWS				Number of rows (second dimension) of A matrix.
*  IDSC_COLS				Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS				Number of rowss (second dimension) of A_DSC matrix.
*  maxnb					Maximun adjacents nodes. Watsolve parameter.
*  maxnn					Maximun number of unknowns. Watsolve parameter.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*	ASSEMBLE_FULL_INTO_SHUFFLED_SPARSE
*	ASSEMBLE_SYM_INTO_SHUFFLED_SPARSE
*	ASSEMBLE_SYM_NO_DIAG_INTO_SHUFFLED_SPARSE
*
* HISTORY: First coding: JHG (7-2003)
*
********************************************************************************



	IMPLICIT NONE
! EXTERNAL VARIABLES, SCALARS
      INTEGER*4 I_TYPE_A,NN,NB, LBLOCK_NR
     ;         ,LMXNDL , NUMEL, IA_COLS,IA_ROWS
      REAL*8 FACTOR

! EXTERNAL VARIABLES, ARRAYS
      INTEGER*4 IAD_D(NB,NN),IADD_D(NN),IADN_D(NN),KXX(LMXNDL,NUMEL)
     ;          ,LNNDEL(NUMEL)
      REAL*8 A(IA_ROWS,IA_COLS),A_DSC(NB,NN) 

! INTERNAL VARIABLES, SCALARS

! INTERNAL VARIABLES, ARRAYS	





		SELECT CASE (I_TYPE_A)

			CASE (1)

				CALL ASSEMBLE_VECTOR_NODE_INTO_SHUFFLED_SPARSE
     ;(FACTOR      ,IA_ROWS     ,LBLOCK_NR      ,NN        ,NB
     ;,A           ,A_DSC       ,IAD_D          ,IADN_D)


	        CASE (2)

	            CALL ASSEMBLE_VECTOR_ELEM_INTO_SHUFFLED_SPARSE 
     ;(FACTOR    ,IA_COLS  ,IA_ROWS       ,LBLOCK_NR        ,NN        
     ;,NB        ,A        ,A_DSC         ,IAD_D            ,IADN_D
     ;,LNNDEL    ,KXX) 


          

	        CASE (3)

	             CALL ASSEMBLE_DERIVATIVE_INTO_SHUFFLED_SPARSE
     ;(FACTOR     ,IA_COLS    ,IA_ROWS    ,LMXNDL     ,LBLOCK_NR  
     ;,NB         ,NN         ,NUMEL      ,A          ,A_DSC      
     ;,IAD_D      ,IADN_D     ,KXX        ,LNNDEL)



			CASE (4)

				CALL ASSEMBLE_FULL_INTO_SHUFFLED_SPARSE
     ;(FACTOR    ,IA_COLS    ,IA_ROWS   ,NUMEL     ,LMXNDL     
     ;,LBLOCK_NR ,NB         ,NN        ,IAD_D     ,IADN_D     
     ;,KXX       ,LNNDEL     ,A         ,A_DSC)      

			CASE (5)

				CALL ASSEMBLE_SYM_INTO_SHUFFLED_SPARSE
     ;(FACTOR     ,IA_COLS   ,IA_ROWS    ,LBLOCK_NR  ,LMXNDL    
     ;,NB         ,NN        ,NUMEL      ,A          ,A_DSC     
     ;,IAD_D      ,IADN_D    ,KXX        ,LNNDEL)


			CASE (6)

				CALL ASSEMBLE_SYM_ND_INTO_SHUFFLED_SPARSE
     ;(FACTOR    ,IA_COLS   ,IA_ROWS    ,LBLOCK_NR ,LMXNDL   
     ;,NB        ,NN        ,NUMEL      ,A         ,A_DSC    
     ;,KXX       ,LNNDEL    ,IAD_D      ,IADN_D    ,IADD_D) 


              CASE(10)

	           CALL ASSEMBLE_CROSSTERM_INTO_SHUFFLED_SPARSE
     ;(FACTOR    ,IA_COLS    ,IA_ROWS   ,NUMEL     ,LMXNDL     
     ;,LBLOCK_NR ,NB         ,NN        ,IAD_D     ,IADD_D
     ;,IADN_D    ,KXX        ,LNNDEL    ,A         ,A_DSC) 


			CASE DEFAULT

				!CALL ERROR 

		END SELECT

       RETURN
       END
