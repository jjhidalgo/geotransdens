      SUBROUTINE ASSEMBLE_RHS_COUPLED
     &          (BCOUPLED ,BFLU     ,BTRA     ,NUMNP)

          IMPLICIT NONE

      
          INTEGER*4::I,NUMNP

          REAL*8::BCOUPLED(2*NUMNP) ,BFLU(NUMNP),BTRA(NUMNP)
      

          DO I=1,NUMNP

              BCOUPLED(2*I - 1) = BFLU(I)
              BCOUPLED(2*I) = BTRA(I)

          END DO

      END SUBROUTINE  ASSEMBLE_RHS_COUPLED