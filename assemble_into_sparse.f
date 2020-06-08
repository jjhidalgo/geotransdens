       SUBROUTINE ASSEMBLE_INTO_SPARSE 
     & (A          ,A_DSC      ,FACTOR     ,KXX        ,I_TYPE_A
     & ,IA_COLS    ,IA_ROWS    ,IAD_S      ,IADD_S     ,IADN_S
     &,LMXNDL      ,LNNDEL     ,NB         ,NN         ,NUMEL)

********************************************************************************
*
* PURPOSE  To assemble a matrix into another one.
*
*
* DESCRIPTION This is a decision routine. It call the appropriate routine
*             according to the type of origin matrix being the destiny matrix a
*             sparse one.
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
*  I_TYPE_A_DSC           Type of A_DSC matrix.
*  IA_COLS                Number of columns (first dimension) of A matrix.
*  IA_ROWS                Number of rows (second dimension) of A matrix.
*  IDSC_COLS              Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS              Number of rowss (second dimension) of A_DSC matrix.
*  maxnb                  Maximun adjacents nodes. Watsolve parameter.
*  maxnn                  Maximun number of unknowns. Watsolve parameter.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*     ASSEMBLE_FULL_INTO_SPARSE
*     ASSEMBLE_SYM_INTO_SPARSE
*     ASSEMBLE_SYM_NO_DIAG_INTO_SPARSE
*
* HISTORY: First coding: JHG (7-2003)
*
********************************************************************************



          IMPLICIT NONE


*  EXTERNAL VARIABLES: SCALARS
      INTEGER*4 I_TYPE_A ,LMXNDL ,IA_COLS,IA_ROWS,NB, NN,NUMEL        
      REAL*8 FACTOR

*  EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 KXX(LMXNDL,NUMEL),IAD_S(NB,NN),IADD_S(NN)
     ;         ,IADN_S(NN),LNNDEL(NUMEL)
      REAL*8   A_DSC(NB,NN), A(IA_ROWS,IA_COLS)

*  INTERNAL VARIABLES: SCALARS

*  INTERNAL VARIABLES: ARRAYS





          SELECT CASE (I_TYPE_A)

              CASE (1)
                  CALL ASSEMBLE_VECTOR_NODE_INTO_SPARSE
     &(A_DSC    ,FACTOR   ,IADD_S     ,NB    ,NN    ,A(1,1))

              CASE (2)

	            CALL ASSEMBLE_VECTOR_ELEM_INTO_SPARSE
     ;(FACTOR   ,LMXNDL    ,NB       ,NN       ,IA_ROWS    ,A
     ;,A_DSC    ,IADD_S    ,KXX      ,LNNDEL)


              CASE (3)

	            CALL ASSEMBLE_DERIVATIVE_INTO_SPARSE
     ;(FACTOR     ,IA_COLS    ,NB         ,IA_ROWS    ,NN      
     ;,LMXNDL     ,A          ,A_DSC      ,IAD_S      ,IADN_S     
     ;,KXX        ,LNNDEL)

              CASE(4)

                  CALL ASSEMBLE_FULL_INTO_SPARSE
     ;(A           ,A_DSC       ,FACTOR     ,IA_COLS    ,IAD_S
     ;,IADN_S      ,KXX         ,LNNDEL     ,LMXNDL     ,NB       
     ;,NN          ,IA_ROWS)
               
              CASE (5)

                  CALL ASSEMBLE_SYM_INTO_SPARSE
     &(A            ,A_DSC        ,FACTOR       ,IAD_S       ,IADN_S         
     &,KXX          ,IA_COLS      ,LMXNDL       ,LNNDEL      ,NB        
     &,NN           ,IA_ROWS)

              CASE (6)

                  CALL ASSEMBLE_SYM_NO_DIAG_INTO_SPARSE
     ;(A           ,A_DSC       ,FACTOR       ,KXX          ,IA_COLS   
     ;,IAD_S       ,IADD_S      ,IADN_S       ,LMXNDL       ,LNNDEL   
     ;,NB          ,NN          ,IA_ROWS)      
     
     
              CASE (10) 
			
				CALL ASSEMBLE_CROSSTERM_INTO_SPARSE
     ;(FACTOR     ,IA_COLS    ,NB         ,NUMEL      ,NN      
     ;,LMXNDL     ,A          ,A_DSC      ,IAD_S      ,IADD_S
     ;,IADN_S     ,KXX        ,LNNDEL)

              CASE DEFAULT

                  !CALL ERROR 

          END SELECT

       RETURN
       END
