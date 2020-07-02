       SUBROUTINE MARQUARDT
     ;(ALF        ,COSMIN      ,DMINF       ,EPS         ,FNEW
     ;,FOLD       ,GMNOR       ,GMNOR1      ,GNORM       ,GNORM1
     ;,IDIMCOV    ,IOWIT       ,IOWPI       ,ISUMFO      ,MAINF       
     ;,MAXICOS    ,MAXITER     ,MIN_STOP    ,NBANDCOV    ,NDEVS       
     ;,NFLAGS     ,NITERF1     ,NITERF2     ,NMTERF1     ,NPAR        
     ;,NSTAT      ,NUMAX       ,NUMIN       ,NUMITER     ,NUMTOBS
     ;,NZPAR      ,OBJCON      ,OBJHED      ,OBJPAR      ,PHIMAX      
     ;,PHIMIN     ,XMARQ       ,XMAXIM      ,COVINV      ,COVPAR      
     ;,DLT_PAR    ,FOBJ_WGT    ,GRAD        ,HESS        ,HESSAUX     
     ;,IFLAGS     ,PARAUX      ,PARC        ,PARGOOD     ,PARM
     ;,PARZ        ,VJAC       ,VOBS        ,VOBSC      ,WGT_UNK
     ;,WORK        ,MEASTYP,   PHI          ,OBJGEO)

********************************************************************************
*
* PURPOSE This subroutine contains all calculations related to Marquardt's 
*         algorithm
*
* DESCRIPTION 
*
*  - Step 0: Declaration of variables
*  - PART A: Objective function has increased: BAD ITERATION       
*    - Step A.1: Echoes information of BAD ITERATION
*    - Step A.2: Recovers last good objective function and updates counters of 
*                bad iterations
*    - Step A.3: Increases Marquardt's parameter
*    - Step A.4: Recovers last good unknown and zonal parameters
*    - Step A.5: Checks STOPPING criteria (total bad iter.or number of bad 
*                iterations since the last good one)
*  - PART B: Objective function has diminished: GOOD ITER.       
*    - Step B.1: Init.counter of bad iter. since last good one
*    - Step B.2: Echoes Marquardt's process indicators
*    - Step B.3: Saves good set of unknown parameters and echoes parameters to
*                a scratch file
*    - Step B.4: Checks the goodness of the quadratic approach of the objective 
*                function and updates Marquardt's parameter accordingly
*    - Step B.5: Calculates gradient vector and its norm
*    - Step B.6: Saves first gradient norm (convergence criterion purposes)
*    - Step B.7: Checks the convergence criteria: small change in parameters, 
*                small change in obj. fun., small gradient norm, gradient norm 
*                smaller than the first one or maximum number of iter. If 
*                convergence is accepted, does not calculate a new solution of 
*                the system.
*    - Step B.8: Calculates and scales Hessian matrix. Marquardt's parameter is
*                not added yet
*  - PART C (COMMON TO A & B): Solves Marquardt's system, as many times as
*                cosin of the angle between gradient vector & parameters vector 
*                is allowed to be violated
*    - Step C.1: Solves Marquardt's system
*    - Step C.2: Checks the angle between gradient and parameter vector. If it 
*                is close o 90 deg., Marquardt's parameter is increased
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix                      
*  COVPAR                 Covariance matrix of unknown parameters
*  DLT_PAR                Array containing param's increment at every iteration.             
*  FOBJ_WGT               Array containing all objective function weights for   
*                         state variables (heads, concentrations, etc)          
*  GRAD                   Vector containing objective function gradient         
*  HESS                   Hessian matrix of objective function.                 
*  HESSAUX                Workspace for inversion of hessian matrix
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INDPAR                 Array of 0's and 1's. 0 means that the parameter is   
*                         estimated aritmethically, otherwise logarithmically.  
*  PARAUX                 Auxiliary variable to store the best computed         
*                         parameters                                            
*  PARC                   Vector containing calculated values for unknown         
*                         parameters                                            
*  PARGOOD                Vector containing last good ZONAL parameters (only for 
*                         output purposes)
*  PARM                   Vector containing measured values for unknown             
*                         parameters                                            
*  PARZ                   Array containing zonal parameters
*  VJAC                   Jacobian matrix                                       
*  VOBS                   Observation value                                     
*  VOBSC                  Value of simulated value corresponding to observation 
*  WGT_UNK                Array containing obj. fun. weights of unknown param.
*  WORK                   Workspace for residuals
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  ALF                    Reduction factor of the unknown parameters                                                                    
*  COSMIN                 XMARQ is multiplied by NUMIN if the cosinus of the    
*                         angle between the angle and the parameters increment  
*                         is less then COSMIN during MAXICOS iterations.        
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
*  IDIMCOV                Used to dimension COVINV                      
*  IOWIT                  Allows writing estimated parameters
*  IOWPI                  Allows writing Marquardt's process indicators along 
*                         the minimization process                                                        
*  ISUMFO                 If 0, obj. function diminishes (GOOD ITERATION), else
*                         increases (BAD ITERATION)                                                      
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MAXICOS                Number of successive iterations where the previous    
*                         criterion can be violatded                            
*  MAXITER                Maximum number of iterations                          
*  MIN_STOP               Control variable for stopping minimization process    
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NDEVS                  Number of observation devices
*  NFLAGS                 Maximum number of allowed flags                       
*  NITERF1                Number of failed iterations (total)
*  NITERF2                Number of failed iterations (since last good one)
*  NMTERF1                Maximum number of failed iterations                   
*  NPAR                   Total number of parameters to be estimated            
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)               
*  NUMAX                  Value to multiplcate XMARQ in apropiate iterations    
*  NUMIN                  Value to divide XMARQ in apropiate iterations         
*  NUMITER                Current iteration in inverse problem process          
*  NUMTOBS                Total number of observations                          
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  OBJCON                 Concentration contribution to objective function      
*  OBJGEO                 Geophysics contribution to objective function 
*  OBJHED                 Head contribution to objective function               
*  OBJPAR                 Parameters contribution to objective function         
*  PHIMAX                 If the above ratio is greather than PHIMAX, it is     
*                         considered a good quadratic aproximation.             
*  PHIMIN                 If the ratio between the actual change on the         
*                         objective function and its quadratic aproximation is  
*                         smaller than PHIMIN, it is considered a poor          
*                         quadratic aproximation.                               
*  XMARQ                  Initial value of Marquardts parameter (0.0)           
*  XMAXIM                 Maximum increment of parameters
*
* INTERNAL VARIABLES: SCALARS
*
*  ISOL                   Dummy counter      
*  PHI                    Coefficient depicting goodness of obj. fun. quadratic 
*                         approach                                         
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  CHECK_PAR_ANGLE        Checks the angle between gradient and param. increment 
*                         vector and updates Marquardt's parameter accordingly
*  CHECK_STOP             Checks stopping/convergence criteria   
*  COMP_PHI_QUADR         Checks the goodness of the obj. fun. quadratic approach
*                         and updates Marquardt's parameter accordingly
*  GRADIENT               Calculates gradient and its norm                               
*  HESSIANO               Calculates scaled Hessian (does not add Marq.'s param.)           
*  IO_SUB                 Flags routines
*  SOLUTION               Solves Marquardt's system
*
* HISTORY
*
*     AAR      Aug-2003   Revision and header inclusion
*     JHG      May-2020   Added geophysical observations contribution.
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NFLAGS,ISUMFO,NUMITER,MAINF,IOWPI,NITERF1,NITERF2,NUMAX
     ;         ,NPAR,MAXITER,MIN_STOP,NMTERF1,NUMIN,IDIMCOV,NBANDCOV
     ;         ,NDEVS,NSTAT,NUMTOBS,MAXICOS,IOWIT,NZPAR
     ;         ,IFLAGS(NFLAGS),MEASTYP(NUMTOBS)
                                                                 ! Real external
      REAL*8::FNEW,FOLD,OBJGEO,OBJHED,OBJCON,OBJPAR,XMARQ,XMAXIM,GNORM
     ;       ,DMINF,EPS
     ;       ,GMNOR,GMNOR1,GNORM1,ALF,PHIMAX,PHIMIN,COSMIN,PHI
     ;       ,PARC(NPAR),PARAUX(NPAR),HESS(NPAR*(NPAR+1)/2),GRAD(NPAR)
     ;       ,DLT_PAR(NPAR),WORK(2*NPAR),COVINV(IDIMCOV),FOBJ_WGT(NSTAT)
     ;       ,COVPAR(NPAR*(NPAR+1)/2),PARM(NPAR),VOBS(NUMTOBS)
     ;       ,VOBSC(NUMTOBS+NDEVS),VJAC(NUMTOBS,NPAR),WGT_UNK(NPAR)
     ;       ,HESSAUX(NPAR*(NPAR+1)/2),PARGOOD(NZPAR),PARZ(NZPAR)
                                                              ! Integer internal
      INTEGER*4 ISOL


      IF (IFLAGS(3).NE.0) CALL IO_SUB('MARQUARDT',0)

C_______________________ PART A: Objective function has increased: BAD ITERATION       

      IF (ISUMFO.NE.0) THEN

C_______________________ Step A.1: Echoes information of BAD ITERATION
      
        PRINT 2000, NUMITER,FNEW,OBJHED,OBJCON,
     ;              OBJPAR,OBJGEO,XMARQ,PHI,XMAXIM,GNORM

        IF (IFLAGS(7).NE.0) WRITE(MAINF,*) ' funcion obj. ',FNEW
        IF (IFLAGS(7).NE.0) WRITE(MAINF,2000) NUMITER,FNEW,OBJHED,
     ;             OBJCON,OBJPAR,OBJGEO,XMARQ,PHI,XMAXIM,GNORM

 2000   FORMAT(///,
     ;     '      INDICATORS OF MARQUARD METHOD; BAD ITERATION.',/,
     ;     ' ITER   FNEW     OBJHED    OBJCON    OBJPAR    OBJGEO     ',
     ;     'XMARQ    PHI(k-1) XMAXIM GNORM',/,I5,9E10.3)

        IF (IOWPI.NE.0) WRITE(71) NUMITER,ISUMFO,FNEW,OBJHED,OBJCON
     ;                           ,OBJPAR,OBJGEO,XMARQ,PHI,XMAXIM,GNORM
       
C_______________________ Step A.2: Recovers last good objective function and 
C_______________________           updates counters of bad iterations

        FNEW=FOLD
        NITERF1=NITERF1+1
        NITERF2=NITERF2+1
   
C_______________________ Step A.3: Increases Marquardt's parameter

        XMARQ=XMARQ*NUMAX+0.001

C_______________________ Step A.4: Recovers last good unknown parameters
C_______________________           and zonal parameters
           
        CALL EQUAL_ARRAY (PARC,PARAUX,NPAR)
        CALL EQUAL_ARRAY (PARZ,PARGOOD,NZPAR)

C_______________________ Step A.5: Checks STOPPING criteria (total bad iter.
C_______________________           or number of bad iterations since the last
C_______________________           good one)

        CALL CHECK_STOP       
     ;(DMINF     ,EPS        ,FNEW      ,FOLD      ,GMNOR
     ;,GMNOR1    ,GNORM      ,GNORM1    ,IOWPI     ,ISUMFO    
     ;,MAINF     ,MAXITER    ,MIN_STOP  ,NITERF1   ,NITERF2   
     ;,NMTERF1   ,NUMITER    ,XMAXIM)

        IF (IOWIT.LT.0) WRITE(70) ISUMFO,NUMITER,PARGOOD
      
      ELSE      

C_______________________ PART B: Objective function has diminished: GOOD ITER.       
      
C_______________________ Step B.1: Init.counter of bad iter. since last good one

        NITERF1=0
       
C_______________________ Step B.2: Echoes Marquardt's process indicators

        PRINT 2100, NUMITER,FNEW,OBJHED,OBJCON,
     ;             OBJPAR,XMARQ,PHI,XMAXIM,GNORM

        IF (IFLAGS(7).NE.0) WRITE(MAINF,*) ' funcion obj. ',FNEW
        IF (IFLAGS(7).NE.0) WRITE(MAINF,2100) NUMITER,FNEW,OBJHED,
     ;              OBJCON,OBJPAR,XMARQ,PHI,XMAXIM,GNORM

 2100   FORMAT(///,
     ;     '      INDICATORS OF MARQUARD METHOD; GOOD ITERATION.',/,
     ;     ' ITER   FNEW     OBJHED    OBJCON    OBJPAR    OBJGEO     ',
     ;     'XMARQ    PHI(k-1) XMAXIM GNORM',/,I5,8E10.3)
                                                                  
        IF (IOWPI.NE.0) WRITE(71) NUMITER,ISUMFO,FNEW,OBJHED,OBJCON
     ;                           ,OBJPAR,OBJGEO,XMARQ,PHI,XMAXIM,GNORM

C_______________________ Step B.3: Saves good set of unknown and zonal param.
C_______________________           and echoes parameters to scratch file
      
        CALL EQUAL_ARRAY (PARAUX,PARC,NPAR)
        CALL EQUAL_ARRAY (PARGOOD,PARZ,NZPAR)
        IF (IOWIT.NE.0) WRITE(70) ISUMFO,NUMITER,PARGOOD

C_______________________ Step B.4: Checks the goodness of the quadratic approach 
C_______________________           of the objective function and updates 
C_______________________           Marquardt's parameter accordingly

        IF (NUMITER.GT.1) CALL COMP_PHI_QUADR
     ;(ALF     ,FNEW     ,FOLD     ,NPAR     ,NUMAX
     ;,NUMIN   ,NUMITER  ,PHI      ,PHIMAX   ,PHIMIN
     ;,XMARQ   ,DLT_PAR  ,GRAD     ,HESS)

C_______________________ Step B.5: Calculates gradient vector and its norm

        CALL GRADIENT
     ;(GNORM    ,IDIMCOV  ,MAINF    ,NBANDCOV ,NDEVS    ,NFLAGS   
     ;,NPAR     ,NSTAT    ,NUMTOBS  ,COVINV   ,COVPAR   ,FOBJ_WGT 
     ;,GRAD     ,IFLAGS   ,PARC     ,PARM     ,WORK(1)  ,VJAC     
     ;,VOBS     ,VOBSC    ,WGT_UNK  ,MEASTYP)

C_______________________ Step B.6: Saves first gradient norm (convergence 
C_______________________           criterion purposes)

        IF (NUMITER.EQ.1) GNORM1=GNORM 
                    
C_______________________ Step B.7: Checks the convergence criteria: small change
C_______________________           in parameters, small change in obj. fun., 
C_______________________           small gradient norm, gradient norm smaller 
C_______________________           than the first one or maximum number of iter.
C_______________________           If convergence is accepted, does not calculate
C_______________________           a new solution of the system.

        CALL CHECK_STOP
     ;(DMINF     ,EPS        ,FNEW      ,FOLD      ,GMNOR
     ;,GMNOR1    ,GNORM      ,GNORM1    ,IOWPI     ,ISUMFO    
     ;,MAINF     ,MAXITER    ,MIN_STOP  ,NITERF1   ,NITERF2   
     ;,NMTERF1   ,NUMITER    ,XMAXIM)

        IF (MIN_STOP.EQ.1) RETURN 

C_______________________ Step B.8: Calculates and scales Hessian matrix. 
C_______________________           Marquardt's parameter is not added yet

        CALL HESSIANO
     ;(IDIMCOV  ,1          ,MAINF    ,NBANDCOV ,NFLAGS   ,NPAR     
     ;,NSTAT    ,NUMTOBS    ,COVINV   ,COVPAR   ,FOBJ_WGT ,HESS     
     ;,IFLAGS   ,VJAC       ,WGT_UNK  ,MEASTYP)

      END IF !IF ((ISUMFO.NE.O).OR.(ISUMFO.EQ.0))

C_______________________ PART C: Solves Marquardt's system, as many times as
C_______________________         cosin of the angle between gradient vector & 
C_______________________         parameters vector is allowed to be violated

      DO ISOL=1,MAXICOS

C_______________________ Step C.1: Solves Marquardt's system
           
        CALL SOLUTION 
     ;(MAINF    ,NFLAGS   ,NPAR     ,NUMITER      ,XMARQ
     ;,DLT_PAR  ,GRAD     ,HESS     ,HESSAUX      ,IFLAGS)

C_______________________ Step C.2: Checks the angle between gradient and 
C_______________________           parameter vector. If it is close o 90 deg.
C_______________________           Marquardt's parameter is increased

        CALL CHECK_PAR_ANGLE 
     ;(COSMIN  ,GNORM  ,NPAR  ,NUMIN  ,XMARQ  ,DLT_PAR  ,GRAD)

      END DO ! ISOL=1,MAXICOS

      IF (IFLAGS(3).NE.0) CALL IO_SUB('MARQUARDT',1)

      RETURN 

      END
