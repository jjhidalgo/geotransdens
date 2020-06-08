	SUBROUTINE ASSEMBLE_VECTOR_NOD_INTO_SHUFFLED_BAND
     ;(FACTOR     ,IA_ROWS     ,IDSC_COLS   ,IDSC_ROWS
     ;,LBLOCK_NR  ,NBAND       ,A           ,A_DSC)
   	   

********************************************************************************
*
* PURPOSE  To assemble a vector into a shuffled banded matrix symmetric or not.
*
*	
* DESCRIPTION 
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A						Origin vector.
*  A_DSC					Destiny matrix (banded).
*
* EXTERNAL VARIABLES: SCALARS
*
*  FACTOR					Factor that multiplies the matrix to be assembled.
*						For time-weighted schemes pourposes (theta-scheme).
*  IA_ROWS				Number of rows (second dimension) of A matrix.
*  IDSC_COLS				Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS				Number of rows (second dimension) of A_DSC matrix.
*  LBLOCK_NR				The block of the non-shufled matrix where the element-wise matrix
*						should go. 
*							if 1, the block of derivatives of flow eqn to head
*							if 2, the block of derivatives of flow eqn to mass fraction
*							if 3, the block of derivatives of tpt eqn to head
*							if 4, the block of derivatives of tpt eqn to mass fraction
*  NBAND					Bandwith of A_DSC.
*
*
* INTERNAL VARIABLES: SCALARS
*
*  IDIAG					Position of the diagonal in A_DSC.
*  I						Counter for nodes in DO... END DO statment.
*  INODE					Number of node.
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* OBSERVATIONS
*
* IA_ROWS = NUMNP
*
* HISTORY: First coding: JHG (10-2003)
*
********************************************************************************



		IMPLICIT NONE
	
		INTEGER*4:: IA_ROWS, IDSC_COLS, IDSC_ROWS,
     &			    NBAND, LBLOCK_NR,IDIAG,INODE,NODE,LMDIAG


		REAL*8::  A(IA_ROWS), A_DSC(IDSC_ROWS,IDSC_COLS)

		REAL*8::FACTOR

			

C------------------------- the components of the vector are added to the
C------------------------- diagonal of the matrix..

          LMDIAG = 2*NBAND + 2
		DO INODE = 1, IA_ROWS
			
			SELECT CASE (LBLOCK_NR)

				CASE (1)

					NODE = (INODE - 1)*2 + 1
					IDIAG = LMDIAG

				CASE (2)

					NODE = (INODE - 1)*2 + 1
					IDIAG = LMDIAG +1
			
				CASE(3)

					NODE = INODE*2
					IDIAG = LMDIAG - 1

				CASE(4)

					NODE  = INODE*2
					IDIAG = LMDIAG

			END SELECT

			

			A_DSC(NODE,IDIAG) = A_DSC(NODE,IDIAG) + A(INODE)*FACTOR			


		END DO ! INODE = 1, IA_ROWS


	END SUBROUTINE ASSEMBLE_VECTOR_NOD_INTO_SHUFFLED_BAND
