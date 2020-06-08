       SUBROUTINE ASSEMBLE
     ;(FACTOR    ,IA_COLS      ,IA_ROWS   ,IDSC_COLS  ,IDSC_ROWS  
     ;,I_TYPE_A  ,I_TYPE_A_DSC ,LMXNDL    ,NB         ,NBAND  
     ;,NN        ,NUMEL        ,A         ,A_DSC      ,IAD_S      
     ;,IADD_S    ,IADN_S    ,KXX          ,LNNDEL)

********************************************************************************
*
* PURPOSE  To assemble a matrix into another one.
*
*
* DESCRIPTION This is a decision routine. It call the appropriate routine
*             according to the type of destiny matrix.
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A                      Origin matrix
*  A_DSC                  Destiny matrix.
*  iad                    Columns. Index array of WatSolve.
*  iadd                   Diagonal. Index array of WatSolve.
*  iadn                   Number of columns. Index array of WatSolve.
*
* EXTERNAL VARIABLES: SCALARS
*
*  FACTOR                 Factor that multiplies the matrix to be assembled.
*                         For time-weighted schemes pourposes (theta-scheme).
*
*  IOMET                  Solving method.
*                             0 --> Picard iterations or lineal problem.
*                             1 --> Newton's method.
*  I_TYPE_A_DSC           Type of descomposed matrix.
*  I_TYPE_A               Type of A matrix.
*                         The possible types of matrices are:
*                             1 -->   nodewise Vector.
*                             2 -->   elementwise vector
*                             3 -->   derivative type matrix
*                             4 -->   Full matrix (elementwise).
*                             5 -->   Symmetric matrix (elementwise).
*                             6 -->   Symmetric matrix without diagonal (elementwise).
*                                     (It is supposed tha the sum of the terms
*                                     out of the diagonal is equal to the diagonal with
*                                     the opposite sign).
*                             7 -->   Symmetric banded matrix.
*                             8 -->   Non symmetric banded matrix.
*                             9 -->   Sparse matrix (as the one used by WatSolve).
*  A_DSC matrix.
*  IA_COLS                Number of columns (first dimension) of A matrix.
*  IA_ROWS                Number of rows (second dimension) of A matrix.
*  IDSC_COLS              Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS              Number of rowss (second dimension) of A_DSC matrix.
*  maxnb                  Maximun adjacents nodes. Watsolve parameter.
*  maxnn                  Maximun number of unknowns. Watsolve parameter.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ASSEMBLE_INTO_BAND_SYM
*  ASSEMBLE_INTO_BAND_NO_SYM
*  ASSEMBLE_INTO_SPARSE
*
* HISTORY: First coding: JHG (7-2003)
*
********************************************************************************



          IMPLICIT NONE
      
          INTEGER*4:: I_TYPE_A, I_TYPE_A_DSC,
     &                IA_COLS, IA_ROWS, IDSC_COLS,IDSC_ROWS,
     &                LMXNDL,NBAND,NUMEL,NB,NN


          INTEGER*4:: KXX(LMXNDL,NUMEL),LNNDEL(NUMEL),IAD_S(NB,NN)
     &                ,IADD_S(NN),IADN_S(NN)


          REAL*8::  A(IA_ROWS,IA_COLS), A_DSC(IDSC_ROWS,IDSC_COLS), 
     &              FACTOR

      


          SELECT CASE (I_TYPE_A_DSC)

              CASE (7)

                  CALL ASSEMBLE_INTO_BAND_SYM
     &            (A        ,A_DSC    ,FACTOR   ,I_TYPE_A 
     &            ,IA_COLS  ,IA_ROWS  ,IDSC_COLS,IDSC_ROWS
     &            ,KXX      ,LMXNDL   ,LNNDEL   ,NBAND)


              CASE (8)
              
                  CALL ASSEMBLE_INTO_BAND_NO_SYM
     &           (A        ,A_DSC    ,FACTOR   ,I_TYPE_A 
     &           ,IA_COLS  ,IA_ROWS  ,IDSC_COLS,IDSC_ROWS
     &           ,KXX      ,LMXNDL   ,LNNDEL   ,NBAND)



              CASE (9,10)

                  CALL ASSEMBLE_INTO_SPARSE
     &           (A          ,A_DSC      ,FACTOR     ,KXX      ,I_TYPE_A
     &           ,IA_COLS    ,IA_ROWS    ,IAD_S      ,IADD_S   ,IADN_S
     &           ,LMXNDL     ,LNNDEL     ,NB         ,NN       ,NUMEL)


              CASE DEFAULT

                  !CALL ERROR 

          END SELECT

       RETURN
       END
