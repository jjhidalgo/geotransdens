       SUBROUTINE CHECK_PAR_ANGLE 
     ;(COSMIN  ,GNORM  ,NPAR  ,NUMIN  ,XMARQ  ,DLT_PAR  ,GRAD)

********************************************************************************
*
* PURPOSE Increase Marquardt's parameters if gradient vector and param. incr.
*         vector are parallel
*
* DESCRIPTION Flow chart:
*
*  - Step 0: Declaration of variables
*  - Step 1: Computes Cosin of angle between gradient and parameters increment 
*            vector
*  - Step 2: Updates XMARQ in order to increase angle
*
* EXTERNAL VARIABLES: ARRAYS
*
*  DLT_PAR                Vector containing increments of unknown parameters
*  GRAD                   Vector containing objective function gradient
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  COSMIN                 XMARQ is multiplied by NUMIN if the cosinus of the    
*                         angle between the angle and the parameters increment  
*                         is less then COSMIN during MAXICOS iterations.        
*  GNORM                  Gradient norm
*  NPAR                   Number of inverse problem unknowns
*  NUMIN                  Value to divide XMARQ in apropiate iterations         
*  XMARQ                  Marquardt's parameter
*
* INTERNAL VARIABLES: SCALARS

*  PNORM                  Norm of parameters increment array
*  VCOS                   Angle gradient-increment of parameters
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY:  AAR   Revision (August-2003)
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

       IMPLICIT NONE
                                                              ! Integer external
       INTEGER*4 NUMIN,NPAR
                                                                 ! Real external
       REAL*8 COSMIN,GNORM,XMARQ,DLT_PAR(NPAR),GRAD(NPAR)
                                                              ! Integer internal
       INTEGER*4 IPAR
                                                                 ! Real internal
       REAL*8 PNORM,VCOS
        
C_______________________ Step 1: Computes Cosin of angle between gradient and 
C_______________________        parameters increment vector

       VCOS=0D0
       PNORM=0D0

       DO IPAR=1,NPAR
          VCOS=VCOS-GRAD(IPAR)*DLT_PAR(IPAR)
          PNORM=PNORM+(DLT_PAR(IPAR))**2
       END DO

       PNORM=DSQRT(PNORM)
       VCOS=VCOS/(GNORM*PNORM)

C_______________________ Step 2: Updates XMARQ in order to increase angle

       IF (VCOS.LT.COSMIN)THEN
         XMARQ=XMARQ*NUMIN+1D-3
       ENDIF

       RETURN
       END
