      SUBROUTINE ASSEMBLE_DERIVATIVE_INTO_SHUFFLED_BAND_NO_SYM
     &          (A         ,A_DSC     ,FACTOR   ,IA_COLS   ,IA_ROWS
     &          ,IDSC_COLS ,IDSC_ROWS ,KXX      ,LBLOCK_NR ,LNNDEL
     &          ,NBAND)


********************************************************************************
*
* PURPOSE  To assemble a matrix of derivatives into a non-symmetric banded matrix.
*
*
* DESCRIPTION The subroutine assembles the elements (multiplied by a factor)
*             of a matrix of derivatives into symmetric banded one.
*
*
*             The matrix of derivatives must be stored elementwise
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A                      Origin matrix (elementwise derivative)
*  A_DSC                  Destiny matrix (non-symmetric banded).
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
* HISTORY: First coding: JHG (2-2004)
*
********************************************************************************





              
          IMPLICIT NONE
         
          INTEGER*4::I,IA_COLS, IA_ROWS, IDSC_COLS, IDSC_ROWS
     &              ,IJCOL,INODE, IROW,J, JNODE,JROW,L
     &              ,LBLOCK_NR,NBAND, NNUD,LMDIAG

          INTEGER*4::LNNDEL(IA_ROWS),KXX(IA_COLS,IA_ROWS)

          REAL*8:: DI_DJ,FACTOR

          REAL*8::  A(IA_ROWS,IA_COLS),
     &              A_DSC(IDSC_ROWS,IDSC_COLS)

          
          LMDIAG = NBAND*2 +2
         
          DO L=1,IA_ROWS ! = NUMEL

              NNUD = LNNDEL(L)

              DO I=1,NNUD

                  INODE = KXX(I,L)

                  DI_DJ = FACTOR * A(L,I)


                  DO J=1,NNUD

                      JNODE = KXX(J,L)

                      

                      SELECT CASE(LBLOCK_NR)

                          CASE(1) !block dfludflu

                              IROW = (INODE - 1) * 2 + 1
                              JROW = (JNODE - 1) * 2 + 1

                          CASE(2)  !block dfludtra

                              IROW = (INODE - 1) * 2 + 1
                              JROW = JNODE * 2 

                          CASE(3)  !block dtradflu

                              IROW = INODE * 2 
                              JROW = (JNODE - 1) * 2 + 1

                          CASE(4)   !block dtradtra

                              IROW = INODE * 2
                              JROW = JNODE * 2

                      END SELECT !LBLOCK_NR
                      
					IF (IROW .GT. JROW) THEN 
						IJCOL = LMDIAG - IABS(IROW-JROW)
	                ELSE
	                    IJCOL = LMDIAG + IABS (IROW-JROW)
	                ENDIF
                      
                      A_DSC(IROW,IJCOL) =A_DSC(IROW,IJCOL) + DI_DJ


                  END DO !J=1,NNUD

              END DO    !I=1,NNUD

          END DO      !L=1,NUMEL

      END SUBROUTINE ASSEMBLE_DERIVATIVE_INTO_SHUFFLED_BAND_NO_SYM