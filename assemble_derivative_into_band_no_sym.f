      SUBROUTINE ASSEMBLE_DERIVATIVE_INTO_BAND_NO_SYM
     &          (A         ,A_DSC     ,FACTOR   ,IA_COLS   ,IA_ROWS
     &          ,IDSC_COLS ,IDSC_ROWS ,KXX      ,LNNDEL
     &          ,NBAND)


********************************************************************************
*
* PURPOSE  To assemble a matrix of derivatives into a non-symmetric banded matrix.
*
*
* DESCRIPTION The subroutine assembles the elements (multiplied by a factor)
*             of a matrix of derivatives into non-symmetric banded one.
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
         
          INTEGER*4::I,IA_COLS, IA_ROWS, IDIAG, IDSC_COLS, IDSC_ROWS
     &              ,IJCOL,INODE, J,JICOL,JNODE,L
     &              ,NBAND, NNUD

          INTEGER*4::LNNDEL(IA_ROWS),KXX(IA_COLS,IA_ROWS)

          REAL*8:: DI_DJ,FACTOR

          REAL*8::  A(IA_ROWS,IA_COLS),
     &              A_DSC(IDSC_ROWS,IDSC_COLS)

          

          IDIAG = NBAND + 1

          DO L=1,IA_ROWS ! = NUMEL

              NNUD = LNNDEL(L)

              DO I=1,NNUD

                  INODE = KXX(I,L)
				DI_DJ = FACTOR * A(L,I)

                  DO J=I+1,NNUD

                      JNODE = KXX(J,L)

                      IF (INODE.GT.JNODE) THEN
						IJCOL = IDIAG - IABS(INODE-JNODE)
						JICOL = IDIAG + IABS(INODE-JNODE)
	                ELSE
						IJCOL = IDIAG + IABS(INODE-JNODE)
						JICOL = IDIAG - IABS(INODE-JNODE)
					ENDIF
                      
                      A_DSC(INODE,IJCOL) =A_DSC(INODE,IJCOL) + DI_DJ

                      A_DSC(JNODE,JICOL) =A_DSC(JNODE,JICOL) + DI_DJ

         
                  END DO !J=I+1,NNUD
              END DO    !I=1,NNUD

C------------------------- Diagonal.      
              DO I=1,NNUD

				INODE = KXX(I,L)
                  DI_DJ = FACTOR * A(L,I)

                  A_DSC(INODE,IDIAG) = A_DSC(INODE,IDIAG) + DI_DJ

              END DO !I=1,NNUD


          END DO      !L=1,NUMEL


      END SUBROUTINE ASSEMBLE_DERIVATIVE_INTO_BAND_NO_SYM