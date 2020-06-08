      SUBROUTINE FUN_PAR
     ;(NPAR     ,OBJPAR   ,COVPAR   ,PARC     ,PARM
     ;,RESID    ,WGT_UNK)

********************************************************************************
*
* PURPOSE Computes the objective function of unknown parameters , defined as
*
*                                    * t     -1     *
*            OBJPAR=SUM LAMBDA  (p -p )    V   (p -p )     j=1,NTYPAR
*                    j        j   j  j      j    j  j
*
* DESCRIPTION Very simple, due to COVPAR structure.
*
*    - Step 0: Declaration of variables
*    - step 1: Calculates residuals of parameters
*    - Step 2: Calculates product of COVPAR*RESID
*    - Step 3: Calculates the product WEIGHT*RESID*COVPAR*RESID, updating
*              the objective function of parameters
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVPAR                 Inverse of the a priori covariance matrix of param.
*  PARC                   Vector containing calculated values for all
*                         unknown parameters
*  PARM                   Vector containing measured/expected values for all
*                         unknown parameters
*  RESID                  Auxiliar array for residuals/weighted residuals
*  WGT_UNK                Array containing parameters objectivefunctions weights 
*                         of problems unknowns
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  NPAR                   Total number of parameters to be estimated
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  OBJPAR                 Parameters contribution to objective function
*
* INTERNAL VARIABLES: SCALARS
*
*  IPAR                   Dummy counter of estimated parameters
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  MUL_SYMMAT_VEC         Product of symmetric matrix (vector storage)*vector
*  ZERO_ARRAY             Sets to 0 a given array of real records
*
* HISTORY: First coding (AAR) Jan-2003
*
********************************************************************************

C_______________________ Step 0: Declaration of variables. Initialization

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NPAR
                                                                 ! Real external
      REAL*8 OBJPAR,COVPAR(NPAR*(NPAR+1)/2),PARC(NPAR)
     ;      ,PARM(NPAR),RESID(NPAR,2),WGT_UNK(NPAR)
                                                              ! Integer internal
      INTEGER*4 IPAR
                                                                 ! Real internal
      OBJPAR=0.D0
      CALL ZERO_ARRAY(RESID,2*NPAR)

C_________________ Step 1: Calculates residuals of parameters

      DO IPAR=1,NPAR
        RESID(IPAR,1)=PARC(IPAR)-PARM(IPAR)
      END DO

C____________________________________ Step 2: Calculates product of COVPAR*RESID

      CALL MUL_SYMMAT_VEC (NPAR,COVPAR,RESID(1,2),RESID(1,1))

C_______________________ Step 3: Calculates final product and updates obj. func.

      DO IPAR=1,NPAR
           OBJPAR=OBJPAR
     ;           +WGT_UNK(IPAR)*RESID(IPAR,1)*RESID(IPAR,2)
      END DO

      RETURN
      END
