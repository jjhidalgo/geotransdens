      SUBROUTINE ASSEMBLE_VEC_ELEM_INTO_BAND
     &           (A        ,A_DSC    ,FACTOR   ,IDSC_COLS,IDSC_ROWS
     &           ,KXX      ,LMXNDL   ,LNNDEL   ,NBAND    ,NUMEL)




********************************************************************************
*
* PURPOSE  To assemble a elementwise vector into a banded matrix symmetric or not.
*
*     
* DESCRIPTION The subroutine assembles the elements (multiplied by a factor)
*             of a vector into a banded one. The vector is added to
*             the diagonal of the matrix.
*             The vector is supposed to be stored nodewise. Then, the number
*             of rows must be equal to the number of nodes in the finite
*             element grid.
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A                      Origin vector. (elementwise)
*  A_DSC                  Destiny matrix (banded).
*  LNNDEL                 Number of nodes at every element
*
* EXTERNAL VARIABLES: SCALARS
*
*  FACTOR                 Factor that multiplies the matrix to be assembled.
*                         For time-weighted schemes pourposes (theta-scheme).
*  IDSC_COLS              Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS              Number of rows (second dimension) of A_DSC matrix.
*  LMXNDL                 Maximum number of nodes in an element.
*  NBAND                  Bandwith of A_DSC.
*  NUMEL                  Number of elements.
*
*
* INTERNAL VARIABLES: SCALARS
*
*  IDIAG                  Position of the diagonal in A_DSC.
*  I                      Counter for nodes in DO... END DO statment.
*  L                      Counter for elements in DO... END DO statment.
*  NNUD                   Number of nodes of current element.
*  NODE                   Number of node I in element L.
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY: First coding: JHG (7-2003)
*
********************************************************************************



          IMPLICIT NONE
      
          INTEGER*4:: IDIAG,IDSC_COLS, IDSC_ROWS, INODE,L,
     &                LMXNDL,NBAND,NNUD,NUMEL,I

          INTEGER*4::LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)


          REAL*8::  A(NUMEL,LMXNDL), A_DSC(IDSC_ROWS,IDSC_COLS)

          REAL*8::FACTOR

C------------------------- First, the position of the diagonal is calculated.

       IDIAG = NBAND + 1



C------------------------- The components of the vector are added to the
C------------------------- diagonal of the matrix..
       DO L = 1, NUMEL

          NNUD=LNNDEL(L)

          DO I = 1, NNUD
              
			INODE = KXX(I,L) 
              
              A_DSC(INODE,IDIAG) = A_DSC(INODE,IDIAG) + 
     &                             A(L,I)*FACTOR           

          END DO ! INODE = 1, NNUD

       END DO !L = 1, NUMEL


      END SUBROUTINE ASSEMBLE_VEC_ELEM_INTO_BAND