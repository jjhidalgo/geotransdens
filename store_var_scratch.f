        SUBROUTINE STORE_VAR_SCRATCH 
     ;(NINT     ,NTCOMP   ,NREG     ,NUMNP    ,TABSOLUT ,TIME
     ;,VAR)                                               

********************************************************************************
*
* PURPOSE  Stores last good solution at a scracteche file
*
*
* DESCRIPTION As a function of the register number, if=1 -->steady state
*                                                   ne.1 -->transient
*             If code is within steady state case, it will store VAR(NUMNP)
*             at the first register and if code is within the inner transient
*             loop, will write the current solution time (absolut) and the 
*             computed solution.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  TIME                   Observation times.                                    
*  VAR                    State variable (head, pressure or concentration)                                                      
*
* EXTERNAL VARIABLES: SCALARS
*
*  NINT                   Number of observation times                           
*  NREG                   Register number of the scratched file                                                      
*  NTCOMP                 Total number of computation steps                                                      
*  NUMNP                  Number of nodes                                       
*  TABSOLUT               Current absolut computation time                      
*
* HISTORY   AAR & GG First coding: JULY-98 (also in vacancy period)
*
********************************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      
      DIMENSION TIME(NINT),VAR(NUMNP)

C______________________________ Stores last good solution at a scratched file

      WRITE(50,REC=NREG) VAR
      
C______________________________Steady state case

      IF(NREG.EQ.1) THEN

        NTCOMP=1
        WRITE(52) TIME(1)

      ELSE

C______________________________Transient

        NTCOMP=NTCOMP+1
        WRITE(52) TABSOLUT

      END IF

      RETURN
      END

