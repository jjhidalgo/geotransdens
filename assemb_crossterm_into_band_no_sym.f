      SUBROUTINE ASSEMBLE_CROSSTERM_INTO_BAND_NO_SYM
     &          (A        ,A_DSC    ,FACTOR   ,IA_COLS  ,IDSC_COLS
     &          ,IDSC_ROWS,KXX      ,LMXNDL   ,LNNDEL   ,IA_ROWS
     &          ,NBAND)


********************************************************************************
*
* PURPOSE  To assemble a matrix whose rows are equal and the first element is
*          equal to the sum of the rest of th elements with the opposite sing
*          into a shuffled non-symmetric banded matrix.
*
*
* DESCRIPTION The subroutine assembles the elements (multiplied by a factor)
*             of a 'row-wise' matrix into shuffled symmetric banded one.
*
*
*             The 'row-wise' matrix is a matrix whiose rows are all equal.
*             Furthermore, the elemetn A1J = -SUMJ(AIJ), J>2.
*             It must be stored in a vector in the following
*             fashion:
*
*                      | 11  12  13 |
*                      | 11  12  13 |   =>  |12, 13|
*                      | 11  12  13 |
*
*             Then, The element (I,J) is located the position 
*             (J - 1) in the vector (I > 1)
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A                      Origin matrix (full)
*  A_DSC                  Destiny matrix (non-symmetric banded).
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  FACTOR                 Factor that multiplies the matrix to be assembled.
*                         For time-weighted schemes pourposes (theta-scheme).
*
*  IA_COLS                Number of columns (first dimension) of A matrix.
*  IA_ROWS                Number of rows (second dimension) of A matrix.
*  IDSC_COLS              Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS              Number of rowss (second dimension) of A_DSC matrix.
*  LMXNDL                 Maximum number of nodes in an element
*  NBAND                  Matrix bandwidth
*
* INTERNAL VARIABLES: ARRAYS
*
*
*  ROWVECTOR             Auxiliar vecto containing a row of A
*  I                     Counter
*  IDIAG                 Position of the diagonal in A_DSC
*  IJCOL                 Counter
*  J                     Counter
*  JICOL                 Counter
*  L                     Element counter
*  NODE1                 Node number
*  NODE2                 Node number
*  NNUD                  Numbers of nodes in an element
*
* HISTORY: First coding: JHG (11-2004)
*
********************************************************************************


      IMPLICIT NONE

C-------------------- External

      INTEGER*4:: IA_COLS, IA_ROWS, IDSC_COLS,IDSC_ROWS

      INTEGER*4::LMXNDL, NBAND

      REAL*8:: FACTOR

      INTEGER*4::LNNDEL(IA_ROWS),KXX(LMXNDL,IA_ROWS)

      REAL*8::  A(IA_ROWS,IA_COLS), A_DSC(IDSC_ROWS,IDSC_COLS)

C-------------------- Internal

      INTEGER*4::I,IDIAG,IJCOL,J,JICOL,L,NODE1,NODE2,NNUD
      
      REAL*8,ALLOCATABLE:: ROWVECTOR(:)

C--------------------
C--------------------

      ALLOCATE(ROWVECTOR(LMXNDL))

      IDIAG = NBAND + 1

      DO L=1,IA_ROWS
          NNUD=LNNDEL(L)

C--------------------Fills local vector with row and
C--------------------calculate missing element

          ROWVECTOR=0D0

          DO I=2,NNUD
              ROWVECTOR(1) = ROWVECTOR(1)-A(L,I-1)
              ROWVECTOR (I) = A(L,I-1)
          ENDDO

          DO I=1,NNUD-1

              NODE1=KXX(I,L)

              DO J=I+1,NNUD

                  NODE2 = KXX(J,L)

C--------------------first we sum element ij....


                  IF (NODE1.LT.NODE2) THEN
                      IJCOL = IDIAG + IABS(NODE1-NODE2)
                  ELSE
                      IJCOL = IDIAG - IABS(NODE1-NODE2)
                  ENDIF

                  A_DSC(NODE1,IJCOL) = A_DSC(NODE1,IJCOL)
     &                               + ROWVECTOR(J)*FACTOR

C--------------------and then  we sum element ji
                  IF (NODE2.LT.NODE1) THEN
                      JICOL = IDIAG + IABS(NODE1-NODE2)
                  ELSE
                      JICOL = IDIAG - IABS(NODE1-NODE2)
                  ENDIF
                  A_DSC(NODE2,JICOL) = A_DSC(NODE2,JICOL)
     &            + ROWVECTOR(I)*FACTOR

              END DO ! J=I+1,NNUD

          END DO ! I=1,NNUD-1


C------------------------- Diagonal.

          DO I=1,NNUD

              NODE1 = KXX(I,L)

              A_DSC(NODE1,IDIAG)=A_DSC(NODE1,IDIAG)+ROWVECTOR(I)*FACTOR

          END DO !I=1,NNUD

      END DO ! L=1,IA_ROWS

      DEALLOCATE(ROWVECTOR)

      END SUBROUTINE ASSEMBLE_CROSSTERM_INTO_BAND_NO_SYM
