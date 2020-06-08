       SUBROUTINE COMP_PHI_QUADR 
     ;(ALF     ,FNEW     ,FOLD     ,NPAR     ,NUMAX
     ;,NUMIN   ,NUMITER  ,PHI      ,PHIMAX   ,PHIMIN
     ;,XMARQ   ,DLT_PAR  ,GRAD     ,HESS)

********************************************************************************
*
* PURPOSE Solves the equation system of Marquardt's method. 
*
*                     (H(k)+mu*I)Deltap(k+1)=-grad(k)
*
* DESCRIPTION 
*
* - Step 0: Declaration of variables
* - Step 1: Computes PHI (quadratic approach coefficient)
* - Step 2: Checks goodness of quadratic approach and increases/diminishes 
*           Marquardt's parameter
*
* EXTERNAL VARIABLES: ARRAYS
*
*  DLT_PAR                Array containing parameter's incr. at every iteration.             
*  GRAD                   Vector containing objective function gradient         
*  HESS                   Hessian matrix of objective function.                 
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  ALF                    Reduction coefficient of param. increment array
*  FNEW                   Objective function value in the current iteration     
*  FOLD                   Objective function computed value in last iteration   
*  NPAR                   Total number of parameters to be estimated            
*  NUMAX                  Value to multiplcate XMARQ in apropiate iterations    
*  NUMIN                  Value to divide XMARQ in apropiate iterations         
*  NUMITER                Current iteration in inverse problem process          
*  PHI                    Coefficient of goodness of the cuadratic approach
*  PHIMAX                 If the above ratio is greather than PHIMAX, it is     
*                         considered a good quadratic aproximation.             
*  PHIMIN                 If the ratio between the actual change on the         
*                         objective function and its quadratic aproximation is  
*                         smaller than PHIMIN, it is considered a poor          
*                         quadratic aproximation.                               
*  XMARQ                  Initial value of Marquardts parameter (0.0)           
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY: AMS,GGL,JCR: First coding (???)
*          AAR: Header inclusion (Aug-2003)
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

       IMPLICIT NONE
                                                              ! Integer external
       INTEGER*4 NPAR,NUMIN,NUMAX,NUMITER
                                                                 ! Real external
       REAL*8 PHIMIN,PHIMAX,FNEW,FOLD,PHI,ALF,XMARQ
     ;       ,HESS(NPAR*(NPAR+1)/2),DLT_PAR(NPAR),GRAD(NPAR)
                                                              ! Integer internal
       INTEGER*4 IDIAG,NP
                                                                 ! Real internal
       REAL*8 QUADF

C_______________________ Step 1: Computes PHI (quadratic approach coefficient)

       QUADF=0D0
       IDIAG=0
       DO NP=1,NPAR
          IDIAG=IDIAG+NP
          QUADF=QUADF-XMARQ*HESS(IDIAG)*DLT_PAR(NP)*DLT_PAR(NP)+
     ;          (1.D0-ALF/2.D0)*GRAD(NP)*DLT_PAR(NP)
       END DO
       PHI=(FNEW-FOLD)/QUADF

C_______________________ Step 2: Checks goodness of quadratic approach
C_______________________         and increases/diminishes Marq. parameter
   
 
       IF (NUMITER.GT.1)THEN

         IF (DABS(1D0-PHI).LT.PHIMIN)THEN
           XMARQ=XMARQ/NUMAX

         ELSE IF (DABS(1D0-PHI).LT.PHIMAX) THEN
           XMARQ=XMARQ/NUMIN              

         ELSE IF (DABS(1D0-PHI).GE.3.D0)THEN
           XMARQ=XMARQ*NUMIN                    !+0.001

         END IF
       END IF

       RETURN
       END

       
