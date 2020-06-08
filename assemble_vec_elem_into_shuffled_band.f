	SUBROUTINE ASSEMBLE_VEC_ELEM_INTO_SHUFFLED_BAND
     ;(FACTOR     ,IA_COLS      ,IA_ROWS     ,IDSC_COLS   ,IDSC_ROWS
     ;,LBLOCK_NR  ,LMXNDL       ,NBAND       ,A           ,A_DSC
     ;,KXX        ,LNNDEL)
   	   

********************************************************************************
*
* PURPOSE  To assemble an elementwise vector into a shuffled banded matrix symmetric or not.
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
*  IA_COLS				Number of columns (first dimension) of A matrix.
*  IA_ROWS				Number of rows (second dimension) of A matrix.
*  IDSC_COLS				Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS				Number of rows (second dimension) of A_DSC matrix.
*  LMXNDL					Maximum number of nodes in an element.
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
*  L						Counter for elements in DO... END DO statment.
*  NNUD					Number of nodes of current element.
*  NODE					Number of node I in element L.
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY: First coding: JHG (10-2003)
*
********************************************************************************



		IMPLICIT NONE
	
		INTEGER*4:: IA_COLS, IA_ROWS, IDSC_COLS, IDSC_ROWS, LMXNDL,
     &			    NBAND, LBLOCK_NR

		INTEGER*4::IDIAG,I,L,NNUD,NODE

          INTEGER*4  KXX(LMXNDL,IA_ROWS), LNNDEL(IA_ROWS)

		REAL*8::  A(IA_ROWS,IA_COLS), A_DSC(IDSC_ROWS,IDSC_COLS)
     &			 
		REAL*8::FACTOR


		


C------------------------- For each element,

	 DO L = 1,IA_ROWS !IA_ROWS = NUMEL

		NNUD = LNNDEL(L)

C------------------------- the components of the vector are added to the
C------------------------- diagonal of the matrix..
		DO I = 1, NNUD
			
			SELECT CASE (LBLOCK_NR)

				CASE(1)

					NODE = (KXX(I,L) - 1)*2 + 1
			        IDIAG = 2*NBAND + 2

				CASE(2)

					NODE = (KXX(I,L) - 1)*2 + 1
			        IDIAG = 2*NBAND +3
			
				CASE(3)

					NODE = KXX(I,L)*2
					IDIAG = 2*NBAND + 1

				CASE(4)

					NODE = KXX(I,L)*2
					IDIAG = 2*NBAND + 2

			END SELECT

			

			A_DSC(NODE,IDIAG) = A_DSC(NODE,IDIAG) + A(L,I)*FACTOR			

		END DO ! I = 1, NNUD

	END DO ! L = 1, IA_ROWS



	END SUBROUTINE ASSEMBLE_VEC_ELEM_INTO_SHUFFLED_BAND
