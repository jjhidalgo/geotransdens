       SUBROUTINE CHECK_STOP
     ;(DMINF     ,EPS        ,FNEW      ,FOLD      ,GMNOR
     ;,GMNOR1    ,GNORM      ,GNORM1    ,IOWPI     ,ISUMFO     
     ;,MAINF     ,MAXITER    ,MIN_STOP  ,NITERF1   ,NITERF2    
     ;,NMTERF1   ,NUMITER    ,XMAXIM)        

********************************************************************************
*
* PURPOSE This subroutine checks and writes the reason of the minimization 
*         process to stop.
*
* DESCRIPTION Main output variable is MIN_STOP
*
*             MIN_STOP=0 : Minimization process has to continue
*             MIN_STOP=1 : Minimization process has converged due to:
*                          1) Small change of parameters
*                          2) Small change in the objective function
*                          3) Small gradient norm
*                          4) Gradient norm is smaller than the first one
*                          5) Maximum number of iterations
*             MIN_STOP=2 : Minimization process has not converged due to:
*                          1) Maximum number of iterations
*                          2) Maximum number of bad iterations
*                          In this case, one more simulation will be done
*                          with last good set of unknown parameters
*
* EXTERNAL VARIABLES: ARRAYS
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  DMINF                  Convergence criterion                                 
*  EPS                    Convergence criterion. Algorithm stops if maximum     
*                         relative change in one parameter is smaller than EPS  
*  FNEW                   Objective function value in the current iteration     
*  FOLD                   Objective function computed value in last iteration   
*  GMNOR                  Algorithm stops if gradient norm becomes smaller      
*                         than GMNOR                                            
*  GMNOR1                 Convergence criterion                                 
*  GNORM                  Gradient norm                                         
*  GNORM1                 Gradient norm at first iteration
*  IOWPI                  Allows writing minimization process information
*  ISUMFO                 Control variable. If ISUMFO=0 obj. func. diminishes
*                         (good iteration) else obj. func. increases (bad iter)
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MAXITER                Maximum number of iterations                          
*  MIN_STOP               Control variable for stopping minimization process    
*  NITERF1                Fall iterations since last good iteration         
*  NITERF2                Number of total fall iterations
*  NMTERF1                Maximum number of failed iterations                   
*  NUMITER                Current iteration in inverse problem process          
*  XMAXIM                 Maximum change in the parameters                          
*
* INTERNAL VARIABLES: SCALARS
*
*  ITER1                  Convergence iteration
*  INDEX                  Convergence criterion flag
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  MIN_PROCESS_INFO       Writes minimization process information
*
*
* HISTORY: AAR (Revision)    , Aug 2003
*
********************************************************************************

C______________________________________________ Step 0: Declaration of variables

       IMPLICIT NONE
                                                              ! Integer external
       INTEGER*4 NUMITER,MAXITER,MIN_STOP,NITERF1,NITERF2,NMTERF1
     ;          ,IOWPI,MAINF,ISUMFO
                                                                 ! Real external
       REAL*8 XMAXIM,DMINF,EPS,FNEW,FOLD,GMNOR,GMNOR1,GNORM,GNORM1
                                                              ! Integer internal
       INTEGER*4 ITER1,INDEX
                                                                    ! Characters
       CHARACTER*50 MESS
  
       IF (ISUMFO.EQ.0) THEN    ! GOOD ITERATION
           
         ITER1=NUMITER
         MIN_STOP=1

C_______________________ Small increment of parameters

         IF (NUMITER.GT.1. AND. XMAXIM.LT.EPS) THEN
            MESS=' SMALL INCREMENT OF PARAMETERS                  '
            INDEX=3

C_______________________ Small gradient norm

         ELSE IF (GNORM.LT.GMNOR) THEN
            MESS=' SMALL GRADIENT NORM                            '
            INDEX=4

C_______________________ Change of obj.func. is very small

         ELSE IF (DABS((FNEW-FOLD)/FOLD).LT.DMINF) THEN
            MESS=' SMALL CHANGE OF OBJ. FUN.                      '
            INDEX=5

C_______________________ Number of iterations greater than allowed

         ELSE IF (NUMITER.GE.MAXITER) THEN
            MESS=' NUMBER OF ITERATIONS IS MAXIM                  '
            INDEX=6

C_______________________ Gradient norm smaller than the one in the 
C_______________________ first iteration

         ELSE IF ((GNORM/GNORM1).LT.GMNOR1) THEN
            MESS=' GRADIENT NORM SMALLER THAN THE FIRST ONE       '
            INDEX=7

C_______________________ Process does not reach convergence (yet)

         ELSE 
            MIN_STOP=0
            INDEX=0

         ENDIF

       ELSE IF (ISUMFO.NE.0) THEN        ! BAD ITERATION

          MIN_STOP=2

C_______________________ Number of bad iterations is maximum


          IF (NITERF1.GE.NMTERF1) THEN
             MESS=' NUMBER OF BAD ITERATIONS IS MAXIMUM            '
             INDEX=1

C___________________________________Stop.Number of iterations maxim

          ELSE IF (NUMITER.GE.MAXITER) THEN
             MESS=' NUMBER OF ITERATIONS IS MAXIMUM                '
             INDEX=2

C_______________________ Process does not finish (yet)

          ELSE

             MIN_STOP=0
             INDEX=0

          ENDIF
  
       ENDIF ! ISUMFO.EQ.0

       IF (INDEX.NE.0) CALL MIN_PROCESS_INFO
     ;(FNEW     ,FOLD      ,INDEX    ,IOWPI    ,ITER1    ,MAINF     
     ;,MESS     ,NITERF1   ,NITERF2  ,NUMITER  ,XMAXIM)

       RETURN        
       END

