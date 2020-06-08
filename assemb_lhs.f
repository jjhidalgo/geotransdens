       SUBROUTINE ASSEMB_LHS
     &(IA_COLS    ,IA_ROWS      ,IDERA_COLS ,IDERA_ROWS  ,IDERB_COLS
     &,IDERB_ROWS ,IDSC_COLS    ,IDSC_ROWS  ,ID_COLS     ,ID_ROWS
     &,INDSSTR    ,IONEWT
     &,I_TYPE_A   ,I_TYPE_A_DSC ,I_TYPE_D   ,I_TYP_DERA  ,I_TYP_DERB
     &,LMXNDL       ,MAXNB      ,NBAND       ,NUMEL
     &,NUMNP      ,THETA        ,TINC       ,A           ,A_DSC
     &,D          ,DER_A        ,DER_B      ,IAD_S      
     &,IADD_S     ,IADN_S       ,KXX        ,LNNDEL      ,EPS)

********************************************************************************
*
* PURPOSE To compute the Left Hand Side term of the flow or transport equations
*
*
* DESCRIPTION This subroutine is a generic one, computes flow and transport,
*             linear and non-linear coefficient matrices.
*
*            Depending on the used solving method, the derivative of model
*            equation or the stiffness and storage matrix are stored.
*
*            Different storages are available.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A                        Stiffnes matrix
*  A_DSC                    System matrix.
*  D                        Storage matrix.
*  DER_A                    Matrix containing the derivatives of equation to solve.
*  iad                      Columns. Index array of WatSolve.
*  iadd                     Diagonal. Index array of WatSolve.
*  iadn                      Number of columns. Index array of WatSolve.
*
* EXTERNAL VARIABLES: SCALARS
*
*  FACTOR                    Factor that multiplies the matrix to be assembled.
*                        For time-weighted schemes pourposes (theta-scheme).
*
*  IOMET                    Solving method.
*                            0 --> Picard iterations or lineal problem.
*                            1 --> Newton's method.
*  I_TYPE_A                Type of A matrix.
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
*  I_TYPE_A_DSC             Type of A_DSC matrix.
*  I_TYPE_D                 Type of D matrix.
*  I_TYPE_DER_A             Type of DER_A matrix.
*  IA_COLS                  Number of columns (first dimension) of A matrix.
*  IA_ROWS                  Number of rows (second dimension) of A matrix.
*  ID_COLS                  Number of columns (first dimension) of D matrix.
*  ID_ROWS                  Number of rows (second dimension) of A matrix.
*  IDER_COLS                Number of columns (first dimension) of DER_A array.
*  IDER_ROWS                Number of rows (second dimension) of DER_A array.
*  IDSC_COLS                Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS                Number of rowss (second dimension) of A_DSC matrix.
*  maxnb                    Maximun adjacents nodes. Watsolve parameter.
*  NUMNP                    Maximun number of unknowns. Watsolve parameter.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*    ASSEMBLE
*
* HISTORY: First coding: JHG (7-2003)
*
********************************************************************************

*#############################################################
*
*    We need to multiply by a factor the matrices.
*    It cannot be done like this:
*
*        A = THETA * A
*
*    since we need A matrix for inverse problem and for next
*    time step. Then it should be done inside ASSEMBLE
*    routine.
*
*#############################################################

      IMPLICIT NONE
C-------------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4::I_TYP_DERA         ,I_TYP_DERB         ,I_TYPE_A
     &          ,I_TYPE_A_DSC       ,I_TYPE_D ,IA_COLS  ,IA_ROWS
     &          ,ID_COLS  ,ID_ROWS  ,IDERA_COLS         ,IDERA_ROWS
     &          ,IDERB_COLS         ,IDERB_ROWS         ,IDSC_COLS
     &          ,IDSC_ROWS,INDSSTR  ,IONEWT   ,LMXNDL   ,MAXNB
     &          ,NBAND    ,NUMEL    ,NUMNP

      REAL*8::EPS      ,THETA    ,TINC

C-------------------------- EXTERNAL VARIABLES: ARRAYS

      INTEGER*4::IAD_S(MAXNB,NUMNP) ,IADD_S(NUMNP)  ,IADN_S(NUMNP)
     &          ,KXX(LMXNDL,NUMEL)  ,LNNDEL(NUMEL)
 
      REAL*8::A(IA_ROWS,IA_COLS)           ,A_DSC(IDSC_ROWS,IDSC_COLS)
     &       ,D(ID_ROWS,ID_COLS)           ,DER_A(IDERA_ROWS,IDERA_COLS)
     &       ,DER_B(IDERB_ROWS,IDERB_COLS)

C-------------------------- INTERNAL VARIABLES: SCALARS    

      REAL*8::TINCINV

c------------------------- Initializes system matrix

      A_DSC = 0D0


C------------------------- Storage matrix and stiffness matrix are assembled
C------------------------- always.


      CALL ASSEMBLE
     &    (THETA     ,IA_COLS      ,IA_ROWS   ,IDSC_COLS  ,IDSC_ROWS
     &    ,I_TYPE_A  ,I_TYPE_A_DSC ,LMXNDL    ,MAXNB      ,NBAND
     &    ,NUMNP     ,NUMEL        ,A         ,A_DSC      ,IAD_S
     &    ,IADD_S    ,IADN_S       ,KXX        ,LNNDEL)



      IF (INDSSTR .GT. 0) THEN

          TINCINV = 1D0/TINC

          CALL ASSEMBLE
     &        (TINCINV   ,ID_COLS      ,ID_ROWS   ,IDSC_COLS  ,IDSC_ROWS
     &        ,I_TYPE_D  ,I_TYPE_A_DSC ,LMXNDL    ,MAXNB      ,NBAND
     &        ,NUMNP     ,NUMEL        ,D         ,A_DSC      ,IAD_S
     &        ,IADD_S    ,IADN_S       ,KXX       ,LNNDEL)

      ENDIF

C------------------------- If the solving method is Newton's one,

      IF (IONEWT.EQ.1)  THEN

C------------------------- the derivatives are assembled.
          CALL ASSEMBLE
     &        (EPS       ,IDERA_COLS   ,IDERA_ROWS,IDSC_COLS  ,IDSC_ROWS
     &        ,I_TYP_DERA,I_TYPE_A_DSC ,LMXNDL    ,MAXNB      ,NBAND
     &        ,NUMNP     ,NUMEL        ,DER_A     ,A_DSC      ,IAD_S
     &        ,IADD_S    ,IADN_S       ,KXX       ,LNNDEL)

C-------------------------- The vector of derivatives to the independent term 

          CALL ASSEMBLE
     &        (-EPS      ,IDERB_COLS   ,IDERB_ROWS,IDSC_COLS  ,IDSC_ROWS
     &        ,I_TYP_DERB,I_TYPE_A_DSC ,LMXNDL    ,MAXNB      ,NBAND
     &        ,NUMNP     ,NUMEL        ,DER_B     ,A_DSC      ,IAD_S
     &        ,IADD_S    ,IADN_S       ,KXX       ,LNNDEL)

      END IF !IONEWT.EQ.1

      END SUBROUTINE ASSEMB_LHS
