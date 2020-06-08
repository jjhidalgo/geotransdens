       SUBROUTINE GRADIENT
     ;(GNORM    ,IDIMCOV  ,MAINF    ,NBANDCOV ,NDEVS    ,NFLAGS   
     ;,NPAR     ,NSTAT    ,NUMTOBS  ,COVINV   ,COVPAR   ,FOBJ_WGT 
     ;,GRAD     ,IFLAGS   ,PARC     ,PARM     ,RESID    ,VJAC     
     ;,VOBS     ,VOBSC    ,WGT_UNK  ,MEASTYP)  

********************************************************************************
*
* PURPOSE Computes gradient vector, defined as
*
*                        t  -1       *                       -1       *
*     g=2*SUM (LAMBDA * J *V  * (m -m ) ) + 2*SUM (LAMBDA * V  * (p -p ) )
*          k         k   k  k     k  k         j         j   j     j  j
*
*         where k is type of measurement and j is parameter type. J\k, V\k are
*         the related parts of jacobian and inverse of the a priori covariance
*         matrix, respectively. m\k and p\j are calculated values of measurem.
*         and parameters, respect. m* and p* are measurements of state variab.
*         and parameters, respectively
*
* DESCRIPTION
*
*  - Step 0: Declaration of variables
*  - Step 1: Initialize Gradient
*  - Step 2: Computes components of gradient related to observations
*  - Step 3: Computes components of gradient related to parameters
*    - Step 3.1: Initializes array of residuals
*    - Step 3.2: Calculates residuals of parameters
*    - Step 3.3: Calculates weighted residuals
*    - Step 3.4: Applies final product, with the objective function weight
*  - Step 4: Computes Gradient Euclidian Gnorm
*  - Step 5: Echoes Gradient and norm to MAINF
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix of measurements
*  COVPAR                 Inverse of the covariance matrix of parameters
*  FOBJ_WGT               Array containing all objective function weights for
*                         state variables (heads, concentrations, etc)
*  GRAD                   Vector containing objective function gradient         
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  PARC                   Vector containing calculated values for unknown
*                         parameters                                            
*  PARM                   Vector containing measured values for all             
*                         parameters                                            
*  RESID                  Auxiliar array for residuals of parameters
*  VJAC                   Jacobian matrix                                       
*  VOBS                   Observation value                                     
*  VOBSC                  Value of simulated value corresponding to observation 
*  WGT_UNK                Array containing weights related to unknown parameters
*                         (objective function calculation purposes)
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  GNORM                  Gradient norm
*  IDIMCOV                Used to dimension COVINV
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NBANDCOV               Band width of the meas. covariance matrix           
*  NDEVS                  Number of devices                              
*  NFLAGS                 Maximum number of allowed flags                       
*  NPAR                   Total number of parameters to be estimated            
*  NSTAT                  Maximum number of state variables whose data is used
*                         for calibration (used for dimensioning)
*  NUMTOBS                Total number of observations                          
*
* INTERNAL VARIABLES: SCALARS
*
*  IPAR                   Dummy counter of estimated parameters
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  COMP_GRAD                                                                    
*  IO_SUB                                                                       
*  MUL_SYMMAT_VEC                                                               
*  ZERO_ARRAY                                                                   
*
* HISTORY: AAR (First coding), Jan 2003
*          AAR (Revision)    , Aug 2003
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

       IMPLICIT NONE
                                                              ! Integer external
       INTEGER*4 NPAR,NFLAGS,IDIMCOV,NBANDCOV,NDEVS,NUMTOBS,MAINF,NSTAT
     ;          ,IFLAGS(NFLAGS)
                                                                 ! Real external
       REAL*8 GNORM
     ;       ,GRAD(NPAR),PARC(NPAR),PARM(NPAR),COVINV(IDIMCOV)
     ;       ,VOBS(NUMTOBS),VOBSC(NUMTOBS+NDEVS),RESID(NPAR,2)
     ;       ,COVPAR(NPAR*(NPAR+1)/2),WGT_UNK(NPAR)
     ;       ,VJAC(NUMTOBS,NPAR),FOBJ_WGT(NSTAT),MEASTYP(NUMTOBS)
                                                              ! Integer internal
       INTEGER*4 IPAR

       IF(IFLAGS(3).EQ.1) CALL IO_SUB('GRADIENT',0)

C_______________________ Step 1: Initialize Gradient

       CALL ZERO_ARRAY(GRAD,NPAR)

C_______________________ Step 2: Computes components of gradient related to 
C_______________________         observations

       CALL COMP_GRAD
     ; (IDIMCOV  ,NBANDCOV ,NDEVS    ,NPAR     ,NSTAT    ,NUMTOBS
     ; ,COVINV   ,FOBJ_WGT ,GRAD     ,MEASTYP  ,VJAC     ,VOBS
     ; ,VOBSC)

C_______________________ Step 3: Computes components of gradient related to 
C_______________________         parameters

C_______________________ Step 3.1: Initializes array of residuals

       CALL ZERO_ARRAY(RESID,2*NPAR)

C_______________________ Step 3.2: Calculates residuals of parameters

       DO IPAR=1,NPAR
         RESID(IPAR,1)=PARC(IPAR)-PARM(IPAR)
       END DO

C_______________________ Step 3.3: Calculates weighted residuals

       CALL MUL_SYMMAT_VEC (NPAR,COVPAR,RESID(1,2),RESID(1,1))

C_______________________ Step 3.4: Applies final product, with the objective 
C_______________________           function weight

       DO IPAR=1,NPAR
         GRAD(IPAR)=GRAD(IPAR)+2.D0*WGT_UNK(IPAR)*RESID(IPAR,2)
       END DO

C_______________________  Step 4: Computes Gradient Euclidian Gnorm

       GNORM=0.D0
       DO IPAR=1,NPAR
         GNORM=GNORM+GRAD(IPAR)*GRAD(IPAR)
       END DO
       GNORM=DSQRT(GNORM)

C_______________________  Step 5: Echoes Gradient and norm to MAINF

       IF (IFLAGS(7).NE.0) THEN
         WRITE(MAINF,*) ' COMPONENTES GRADIENTE SON'
         DO IPAR=1,NPAR
           WRITE(MAINF,*)' NP ES=',IPAR,' GRAD(NP)=',GRAD(IPAR)
         END DO
         WRITE(MAINF,*) ' LA NORMA DEL GRADIENTE ES: ',GNORM
       END IF

       IF(IFLAGS(3).EQ.1) CALL IO_SUB('GRADIENT',1)

       RETURN
       END
