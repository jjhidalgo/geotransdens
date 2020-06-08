      SUBROUTINE ASSEMBLE_SYM_NO_DIAG_INTO_BAND_NO_SYM
     &          (A, A_DSC, FACTOR, IA_COLS, IA_ROWS, IDSC_COLS,IDSC_ROWS
     &          ,KXX,LNNDEL,LMXNDL,NBAND)

********************************************************************************
*
* PURPOSE  To assemble a symmetric (without diagonal) matrix into a  non-symetric banded matrix.
*
*
* DESCRIPTION The subroutine assembles the elements (multiplied by a factor)
*             of a symmetric (without diagonal) matrix into non-symmetric banded one.
*             The matrix is supposed to be symmetric and the diagonal is supposed
*             to be equal to the sum of the elements outside of the diagonal with the
*             opposite sign.
*
*             The symmetric (without diagonal) matrix must be stored in a vector
*             in the following fashion:
*
*                      | --  12  13 |
*                      | --  --  23 |   =>  |12, 13, 23|
*                      | --  --  -- |
*
*             Then, The element (I,J) is located the position 
*             (I - 1)*NNUD + J - I*(I+1)/2 in the vector.
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A                      Origin matrix (full, symmetric)
*  A_DSC                  Destiny matrix (symmetric banded).
*
* EXTERNAL VARIABLES: SCALARS
*
*  FACTOR                 Factor that multiplies the matrix to be assembled.
*                         For time-weighted schemes pourposes (theta-scheme).
*
*  IA_COLS                Number of columns (first dimension) of A matrix.
*  NUMEL              Number of rows (second dimension) of A matrix.
*  IDSC_COLS              Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS              Number of rowss (second dimension) of A_DSC matrix.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY: First coding: JHG (7-2003)
*
********************************************************************************



          IMPLICIT NONE
      
          INTEGER*4:: IA_COLS, IA_ROWS, IDSC_COLS,IDSC_ROWS

          INTEGER*4::NNUD,I,J,IJCOL,JICOL,IPOS_AIJ,L,LMXNDL,NODE1
     &             ,NODE2,NBAND,IDIAG

         INTEGER*4::LNNDEL(IA_ROWS),KXX(LMXNDL,IA_ROWS)

          REAL*8::  A(IA_ROWS,IA_COLS),
     &              A_DSC(IDSC_ROWS,IDSC_COLS), 
     &              FACTOR



       IDIAG = NBAND + 1

       DO L=1,IA_ROWS

          NNUD=LNNDEL(L)

          DO I=1,NNUD-1

          
              DO J=I+1,NNUD

                  NODE1 = KXX(I,L)
                  NODE2 = KXX(J,L) 

                  IF (NODE1.LT.NODE2) THEN
				    IJCOL = IDIAG + IABS(NODE1-NODE2)
	                JICOL = IDIAG - IABS(NODE1-NODE2)
	            ELSE
				    IJCOL = IDIAG - IABS(NODE1-NODE2)
	                JICOL = IDIAG + IABS(NODE1-NODE2)
                  ENDIF               
                  

                  
                  IPOS_AIJ = (I - 1)*NNUD + J - I*(I+1)/2

                  A_DSC(NODE1,IJCOL)= A_DSC(NODE1,IJCOL) 
     &                              + A(L,IPOS_AIJ)*FACTOR
				A_DSC(NODE2,JICOL)= A_DSC(NODE2,JICOL)
     &                              + A(L,IPOS_AIJ)*FACTOR
C------------------------- Diagonal.          
                  A_DSC(NODE1,IDIAG) = A_DSC(NODE1,IDIAG)
     &                                - A(L,IPOS_AIJ) * FACTOR

                  A_DSC(NODE2,IDIAG) = A_DSC(NODE2,IDIAG)
     &                                - A(L,IPOS_AIJ) * FACTOR

              END DO ! J=I+1,NNUD

          END DO ! I=1,NNUD-1


      END DO ! L=1,IA_ROWS
      


	END SUBROUTINE ASSEMBLE_SYM_NO_DIAG_INTO_BAND_NO_SYM