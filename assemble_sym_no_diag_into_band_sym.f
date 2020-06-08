      SUBROUTINE ASSEMBLE_SYM_NO_DIAG_INTO_BAND_SYM
     &           (A        ,A_DSC    ,FACTOR   ,IA_COLS  ,IDSC_COLS
     &           ,IDSC_ROWS,KXX      ,LMXNDL   ,LNNDEL   ,NUMEL    
     &           ,NBAND)


********************************************************************************
*
* PURPOSE  To assemble a symmetric (without diagonal) matrix into a symetric banded matrix.
*
*
* DESCRIPTION The subroutine assembles the elements (multiplied by a factor)
*             of a symmetric (without diagonal) matrix into symmetric banded one.
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
      
          INTEGER:: IA_COLS, IDSC_COLS,IDSC_ROWS,NUMEL

          INTEGER::NNUD,I,J,IJCOL,KEXT,IJ_A,L,LMXNDL,INODE,JNODE,
     &             NBAND,NBAND1


          REAL*8::  A(NUMEL,IA_COLS),
     &              A_DSC(IDSC_ROWS,IDSC_COLS), 
     &              FACTOR

          INTEGER*4::LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)


       NBAND1=NBAND+1

       DO L=1,NUMEL

          NNUD=LNNDEL(L)

          DO I=1,NNUD-1

              INODE=KXX(I,L)

              DO J=I+1,NNUD

                  JNODE = KXX(J,L)
                  IJCOL = NBAND1 - IABS(INODE-JNODE)
                  KEXT = MAX0(INODE,JNODE)
                  IJ_A = (I - 1)*NNUD + J - I*(I+1)/2
                  A_DSC(KEXT,IJCOL) = A_DSC(KEXT,IJCOL)
     &                              + A(L,IJ_A)*FACTOR

C------------------------- Diagonal.          
                  A_DSC(INODE,NBAND1) = A_DSC(INODE,NBAND1)
     &                                - A(L,IJ_A) * FACTOR

                  A_DSC(JNODE,NBAND1) = A_DSC(JNODE,NBAND1)
     &                                - A(L,IJ_A) * FACTOR

              END DO ! J=I+1,NNUD


          END DO ! I=1,NNUD-1


      END DO ! L=1,NUMEL

       
          
      END SUBROUTINE ASSEMBLE_SYM_NO_DIAG_INTO_BAND_SYM