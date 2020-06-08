      SUBROUTINE ASSEMBLE_FULL_INTO_BAND_SYM
     &           (A        ,A_DSC    ,FACTOR   ,IA_COLS  ,IDSC_COLS,
     &           IDSC_ROWS ,KXX      ,LMXNDL   ,LNNDEL   ,NBAND    
     &           ,NUMEL)


********************************************************************************
*
* PURPOSE  To assemble a full matrix into a symetric banded matrix.
*
*
* DESCRIPTION The subroutine assembles the elements (multiplied by a factor)
*             of a full matrix into symmetric banded one.
*             The full matrix is supposed to be symmetric, otherwise it would
*             be nonsense to store it in a band symmetric matrix.
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
*  ASSEMBLE_SYM_INTO_BAND_SYM
*
* HISTORY: First coding: JHG (7-2003)
*
********************************************************************************

      IMPLICIT NONE

      INTEGER*4:: IA_COLS, NUMEL,IDSC_COLS,NBAND,
     &            IDSC_ROWS,LMXNDL
	  INTEGER*4::LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)
      REAL*8::A(IA_COLS, NUMEL), A_DSC(IDSC_COLS,IDSC_ROWS)
      REAL*8::FACTOR

      CALL ASSEMBLE_SYM_INTO_BAND_SYM
     &          (A        ,A_DSC    ,FACTOR   ,IA_COLS  ,IDSC_COLS
     &          ,IDSC_ROWS,KXX      ,LMXNDL   ,LNNDEL   ,NUMEL    
     &          ,NBAND)

          
      END SUBROUTINE ASSEMBLE_FULL_INTO_BAND_SYM
