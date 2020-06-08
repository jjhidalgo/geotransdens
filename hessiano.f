       SUBROUTINE HESSIANO
     ;(IDIMCOV  ,IOSCALING  ,MAINF    ,NBANDCOV  ,NFLAGS   ,NPAR     
     ;,NSTAT    ,NUMTOBS    ,COVINV   ,COVPAR    ,FOBJ_WGT ,HESS     
     ;,IFLAGS   ,VJAC       ,WGT_UNK  ,MEASTYP)

********************************************************************************
*
* PURPOSE Compute and scale Hessian matrix, defined as
*
*                        t  -1                       -1
*     H=  SUM (LAMBDA * J *V  * J ) + SUM (LAMBDA * V   )
*          k         k   k  k    k               j    j 
*
*
* DESCRIPTION
*
*  - Step 0: Declaration of variables
*  - Step 1: Initializes Hessian matrix
*  - Step 2: Computes components related to observations
*  - Step 3: Computes components related to parameters
*  - Step 4: Echoes to MAINF Hessian matrix without scaling
*  - Step 5: Scales Hessian (does not include the diagonal)
*                  Hij\esc=Hij/SQRT(Hii)/SQRT(Hjj)
*  - Step 6: Echoes to MAINF the scaled Hessian matrix 
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the observations covariance matrix
*  COVPAR                 Inverse of the parameters covariance matrix         
*  FOBJ_WGT               Array containing all objective function weights for   
*                         state variables (heads, concentrations, etc)          
*  HESS                   Hessian matrix of objective function.                 
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  VJAC                   Jacobian matrix                                       
*  WGT_UNK                Array containing weights related to unknown parameters
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCOV                Used to dimension COVINV                            
*  IOSCALING              If 1, hessian matrix is scaled
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NFLAGS                 Maximum number of allowed flags                       
*  NPAR                   Total number of parameters to be estimated            
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)               
*  NUMTOBS                Total number of observations                          
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Dummy counter     
*  IDIAG                  Dummy counter                                         
*  IDIAG1                 Dummy counter                                         
*  IPAR                   Dummy counter                                         
*  IPOSHESS               Position (I,J) at HESS, COVPAR
*  JPAR                   Dummy counter                        
*  NP                     Dummy counter
*  NP1                    Dummy counter                                         
*  SQR1                   Counter
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  COMP_HESS                                                                    
*  IO_SUB                                                                       
*  ZERO_ARRAY                                                                   
*
* HISTORY: AAR (First coding), Jan-2003
*          AAR (Revision)    , Aug-2003
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

       IMPLICIT NONE
                                                              ! Integer external
       INTEGER*4 NPAR,NFLAGS,IDIMCOV,NBANDCOV,NUMTOBS,MAINF,NSTAT
     ;          ,IOSCALING
     ;          ,IFLAGS(NFLAGS),MEASTYP(NUMTOBS)
                                                                 ! Real external
       REAL*8 HESS(NPAR*(NPAR+1)/2),COVINV(IDIMCOV),FOBJ_WGT(NSTAT)
     ;       ,COVPAR(NPAR*(NPAR+1)/2),VJAC(NUMTOBS,NPAR),WGT_UNK(NPAR)
                                                              ! Integer internal
       INTEGER*4 IPAR,NP,NP1,IDIAG1,I,IDIAG,IPOSHESS,JPAR
                                                                       ! Real internal
       REAL*8 SQR1

       IF(IFLAGS(3).EQ.1) CALL IO_SUB('HESSIANO',0)

C_______________________ Step 1: Initializes Hessian matrix

       CALL ZERO_ARRAY(HESS,NPAR*(NPAR+1)/2)

C_______________________ Step 2: Computes components related to observations

       CALL COMP_HESS
     ;(IDIMCOV  ,NBANDCOV ,NPAR     ,NSTAT    ,NUMTOBS    ,COVINV   
     ;,FOBJ_WGT ,HESS     ,MEASTYP  ,VJAC)

C_______________________ Step 3: Computes components related to parameters

       DO IPAR=1,NPAR
         
         DO JPAR=1,IPAR
           IPOSHESS=IPAR*(IPAR-1)/2+JPAR
           HESS(IPOSHESS)=HESS(IPOSHESS)+COVPAR(IPOSHESS)*WGT_UNK(IPAR)
         END DO

       END DO

C_______________________ Step 4: Echoes to MAINF Hessian matrix without scaling

       IF (IFLAGS(5).LT.0) WRITE(MAINF,*) ' HESS_NO ESC',HESS

       IF (IOSCALING.NE.1) RETURN            ! Does not scale the hessian matrix

C_______________________ Step 5: Scales Hessian (does not include the diagonal)
C_______________________               Hij\esc=Hij/SQRT(Hii)/SQRT(Hjj)

       I=1
       IDIAG=0
       DO NP1=1,NPAR
         IDIAG=IDIAG+NP1
         SQR1=DSQRT(HESS(IDIAG))
         IDIAG1=0                                              
         DO NP=1,NP1-1
           IDIAG1=IDIAG1+NP
           HESS(I)=HESS(I)/(SQR1*DSQRT(HESS(IDIAG1)))
           I=I+1
         END DO
         I=I+1
       END DO

C_______________________ Step 6: Echoes to MAINF the scaled Hessian matrix 

       IF (IFLAGS(6).LT.0) WRITE(MAINF,*) ' HES  ESC',HESS

       IF(IFLAGS(3).EQ.1) CALL IO_SUB('HESSIANO',1)

       RETURN
       END 
