       SUBROUTINE WRI_DERIV (INEW,INTI,NPAR,NPBV,NUMNP,VCAL,DERV)

       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION VCAL(NUMNP,NPBV),DERV(NUMNP,NPAR,2,NPBV)
    
C------------------------- Writes DERV to check derivatives

       DO IP=1,NPAR
          WRITE(69,*)'TIME STEP,IP=',INTI,IP
          DO NPB=1,NPBV
             DO J=1,NUMNP
                WRITE(69,69) VCAL(J,NPB),DERV(J,IP,INEW,NPB)
             ENDDO
          ENDDO
       ENDDO
      
 69    FORMAT(2E23.16)
       RETURN
       END
