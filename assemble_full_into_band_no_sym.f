	SUBROUTINE ASSEMBLE_FULL_INTO_BAND_NO_SYM
     &		  (A        ,A_DSC    ,FACTOR   ,IA_COLS  ,IDSC_COLS
     &		  ,IDSC_ROWS,KXX      ,LMXNDL   ,LNNDEL   ,NUMEL
     &          ,NBAND)


********************************************************************************
*
* PURPOSE  To assemble a full matrix into a non-symmetric banded matrix.
*
*
* DESCRIPTION The subroutine assembles the elements (multiplied by a factor)
*			of a symmetric matrix into symmetric banded one.
*
*
*			The full matrix must be stored in a vector in the following
*			fashion:
*
*					 | 11  12  13 |
*					 | 21  22  23 |   =>  |11, 12, 13, 21, 22, 23,31, 32, 33|
*					 | 31  32  33 |
*
*			Then, The element (I,J) is located the position 
*			(I - 1)*NNUD + J in the vector.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A						Origin matrix (full)
*  A_DSC					Destiny matrix (non-symmetric banded).
*
* EXTERNAL VARIABLES: SCALARS
*
*  FACTOR					Factor that multiplies the matrix to be assembled.
*						For time-weighted schemes pourposes (theta-scheme).
*
*  IA_COLS				Number of columns (first dimension) of A matrix.
*  NUMEL				Number of rows (second dimension) of A matrix.
*  IDSC_COLS				Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS				Number of rowss (second dimension) of A_DSC matrix.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY: First coding: JHG (7-2003)
*
********************************************************************************



		IMPLICIT NONE
	
		INTEGER:: IA_COLS, NUMEL, IDSC_COLS,IDSC_ROWS

		INTEGER::NNUD,I,J,IJCOL,JICOL,NOD_IJ,NOD_JI,L,LMXNDL,NODE1
     &			,NODE2,NBAND,DIAG


		REAL*8::  A(NUMEL,IA_COLS),
     &			  A_DSC(IDSC_ROWS,IDSC_COLS), 
     &			  FACTOR

		INTEGER::LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)


	 DIAG = NBAND + 1

	 DO L=1,NUMEL

		NNUD=LNNDEL(L)

		DO I=1,NNUD-1

			NODE1=KXX(I,L)

			DO J=I+1,NNUD

				NODE2 = KXX(J,L)

                                          ! first we sum element ij....				
				NOD_IJ = (I - 1)*NNUD + J 
	            
				IF (NODE1.LT.NODE2) THEN
				    IJCOL = DIAG + IABS(NODE1-NODE2)
	            ELSE
				    IJCOL = DIAG - IABS(NODE1-NODE2)
                  ENDIF 	          
				     
				A_DSC(NODE1,IJCOL) = A_DSC(NODE1,IJCOL)
     &                              + A(L,NOD_IJ)*FACTOR

                                          ! and then  we sum element ji
				NOD_JI = (J - 1)*NNUD + I 
				IF (NODE2.LT.NODE1) THEN
				    JICOL = DIAG + IABS(NODE1-NODE2)
	            ELSE
				    JICOL = DIAG - IABS(NODE1-NODE2)
                  ENDIF 	          
				A_DSC(NODE2,JICOL) = A_DSC(NODE2,JICOL)
     &                              + A(L,NOD_JI)*FACTOR

			END DO ! J=I+1,NNUD

		END DO ! I=1,NNUD-1


C------------------------- Diagonal.		
		DO I=1,NNUD

			NODE1 = KXX(I,L)
			NOD_IJ = (I - 1)*NNUD + I
			A_DSC(NODE1,DIAG) = A_DSC(NODE1,DIAG)+A(L,NOD_IJ)*FACTOR

		END DO !I=1,NNUD

	END DO ! L=1,NUMEL

	END SUBROUTINE ASSEMBLE_FULL_INTO_BAND_NO_SYM