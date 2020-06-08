       SUBROUTINE SOLUTION 
     ;(   MAINF    ,NFLAGS   ,NPAR     ,NUMITER      ,XMARQ      
     ;   ,PAR      ,GRAD     ,HESS     ,HESSAUX      ,IFLAGS)

********************************************************************************
*
* PURPOSE Solves the equation system of Marquardt's method. 
*
*                     (H(k)+mu*I)Deltap(k+1)=-grad(k)
*
* DESCRIPTION First, coefficient matrix and indepent term are computed, as 
*             Hessian matrix is scaled previously.
*             Second, the equation system is solved (subroutine LEQT1P)
*             Third, undesired unidentiafibility errors are identified (if this
*             is the case. Therefore, code stops automatically. 
*             Otherwise (hopefully), "real" increment of parameters is computed
*
* EXTERNAL VARIABLES: ARRAYS
*
*  GRAD                   Vector containing objective function gradient         
*  HESS                   Hessian matrix of objective function.                 
*  HESSAUX                Auxiliar matrix
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  PAR                    Parameter's increment at every iteration.             
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFLAGS                 Maximum number of allowed flags                       
*  NPAR                   Total number of parameters to be estimated            
*  NUMITER                Current iteration in inverse problem process          
*  XMARQ                  Initial value of Marquardts parameter (0.0)           
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  LEQT1P                 Solves the equation system
*
* HISTORY: AMS,GGL,JCR: First coding (???)
*          AAR: Header inclusion (Feb-2001)
*
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION HESS(NPAR*(NPAR+1)/2),HESSAUX(NPAR*(NPAR+1)/2),
     ;           GRAD(NPAR),PAR(NPAR),IFLAGS(NFLAGS)


C_____________________________Computes coefficients matrix and independent term

       I=1       
       IDIAG=0
       DO NP=1,NPAR
          IDIAG=IDIAG+NP
          PAR(NP)=-1/(2*DSQRT(HESS(IDIAG)))*GRAD(NP)
          HESSAUX(IDIAG)=1D0+XMARQ
          DO NP1=1,NP-1
             HESSAUX(I)=HESS(I)
             I=I+1  
          END DO
          I=I+1
       END DO

C----------------------------------------------Solves the system

      M=1
      IER=0
      CALL LEQT1P(HESSAUX,M,NPAR,PAR,NPAR,D1,D2,IER)

C----------------------------------------------Error message

       IF (IER.NE.0) THEN
          WRITE(MAINF,100) IER,NUMITER
 100      FORMAT(//,' ERROR',I5,
     ;       ' IN THE SOLUTION OF MARQUARDT''S  SYSTEM (SUBR. LEQT1P)'
     ;      '  AT',I3,'ITERATION',/,' PROBABLY THERE ARE SOME ',
     ;       ' IDENTIFIABILITY PROBLEMS.')
          STOP 'IDENTIFIABILITY PROBLEM'
       END IF

C_________________________Computes the "real" increment parameters value

       IDIAG=0
       DO NP=1,NPAR
          IDIAG=IDIAG+NP
          PAR(NP)=PAR(NP)/(DSQRT(HESS(IDIAG)))
       END DO

       IF(IFLAGS(9).NE.0) THEN
         WRITE(MAINF,*)' INCREMENTO DE PARAMETROS:'
         DO NP=1,NPAR
           WRITE(MAINF,*) PAR(NP)       
         END DO
       END IF

       RETURN
       END
