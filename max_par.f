       SUBROUTINE MAX_PAR 
     ;(ALF     ,MAINF    ,NFLAGS  ,NPAR    ,PERMXARIT  ,PERMXLOG  
     ;,XMAXIM  ,DLT_PAR  ,INDPAR  ,IFLAGS  ,PARAUX)

*****************************************************************************
*
* PURPOSE  Calculates the maximum increment of unknown parameters and reduces 
*          it if it is necessary.
*
* DESCRIPTION
*
*  - Step 0: Declaration of variables
*  - Step 1: Initializes maximum increments to be calc.
*  - Step 2: Looks for maximum increment in the components of the unknown 
*            parameters vector. Aritmethic increments are standardized to the 
*            value of the last good value of corresponding unknown parameter
*  -Step 3: Looks for maximum increment, arit. or log
*  -Step 4: Calculates the reduction factor
*  -Step 5: If reduction factor is (0,1), applies the correction, updating 
*           array DLT_PAR. Otherwise, flags reduction factor with 1.0
*  -Step 6: Echoes final set of increments
*
* EXTERNAL VARIABLES: ARRAYS
*
*  DLT_PAR                Array containing param's increment at every iteration.             
*  INDPAR                 Array of 0's and 1's. 0 means that the parameter is   
*                         estimated aritmethically, otherwise logarithmically.  
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  PARAUX                 Auxiliary variable to store the best computed         
*                         parameters                                            
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  ALF                    Reduction factor of the unknown parameters                                                                    
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFLAGS                 Maximum number of allowed flags                       
*  NPAR                   Total number of parameters to be estimated            
*  PERMXARIT              Maximum relative change per iteration for the rest of 
*                         parameters                                            
*  PERMXLOG               Maximum change per iteration of log-transformed       
*                         variables                                             
*  XMAXIM                 Maximum change in the parameters                          
*
* INTERNAL VARIABLES: SCALARS
*
*  IPAR                   Dummy counter of unknown parameters
*  XMAXARIT               Maximum increment of param. estimated arithmetically
*  XMAXLOG                Maximum increment of param. estimated logarithmically
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY: AAR (Revision)    , Aug 2003
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

       IMPLICIT NONE
                                                              ! Integer external
       INTEGER*4 MAINF,NPAR,NFLAGS,INDPAR(NPAR),IFLAGS(NFLAGS)
                                                                 ! Real external
       REAL*8 PERMXLOG,PERMXARIT,ALF,XMAXIM,DLT_PAR(NPAR),PARAUX(NPAR)
                                                              ! Integer internal
       INTEGER*4 IPAR
                                                                 ! Real internal
       REAL*8 XMAXARIT,XMAXLOG

C_______________________ Step 1: Initializes maximum increments to be calc.

       XMAXARIT=0.D0
       XMAXLOG=0.D0

C_______________________ Step 2: Looks for maximum increment in the components 
C_______________________         of the unknown parmeters vector. Aritmethic 
C_______________________         increments are standardized to the value of 
C_______________________         the last good value of corresponding unknown
C_______________________         parameter

       DO IPAR=1,NPAR

          IF (INDPAR(IPAR).EQ.1) THEN               ! Log-transform increment

             XMAXLOG=DMAX1(XMAXLOG,DABS(DLT_PAR(IPAR)))

          ELSE                                      ! Arithmetic increment

             IF (PARAUX(IPAR).NE.0D0)               ! Avoids null
     ;         XMAXARIT=DMAX1(XMAXARIT,DABS(DLT_PAR(IPAR)/PARAUX(IPAR)))

          END IF ! INDPAR(IPAR).EQ.1

       END DO ! IPAR=1,NPAR

C_______________________ Step 3: Looks for maximum increment, arit. or log

       XMAXIM=DMAX1(XMAXARIT,XMAXLOG)

C_______________________ Step 4: Calculates the reduction factor

       IF (XMAXLOG.NE.0D0.AND.XMAXARIT.NE.0D0) THEN
          ALF=DMIN1(PERMXLOG/XMAXLOG,PERMXARIT/XMAXARIT)
       ELSE IF (XMAXLOG.EQ.0D0) THEN
          ALF=PERMXARIT/XMAXARIT
       ELSE IF (XMAXARIT.EQ.0D0) THEN
          ALF=PERMXLOG/XMAXLOG
       END IF

C_______________________ Step 5: If reduction factor is (0,1), applies the
C_______________________         correction, updating array DLT_PAR. Otherwise,
C_______________________         flags reduction factor with 1.0

       IF (ALF.LT.1D0.AND.ALF.GT.0D0) THEN

          DO IPAR=1,NPAR
            DLT_PAR(IPAR)=ALF*DLT_PAR(IPAR)
          END DO

       ELSE

          ALF=1.D0

       END IF

C_______________________ Step 6: Echoes final set of increments

       IF (IFLAGS(9).NE.0)
     ;   WRITE(MAINF,*)'INCREM. DE PARAM. REAL (ya corregidos) ES'
     ;                                                     ,DLT_PAR

       RETURN
       END
