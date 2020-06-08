      SUBROUTINE UNSHUFFLE_RHS(CVEC,HVEC,NUMNP,SHUFFLED)

********************************************************************************
*     
*     PURPOSE
*     
*     Un-shuffles a vector of length 2*N into 2 vector of length N
*     
********************************************************************************

      IMPLICIT NONE

C-------------------------External

      INTEGER*4::NUMNP
      REAL*8::CVEC(NUMNP),HVEC(NUMNP),SHUFFLED(2*NUMNP)

C-------------------------Internal

      INTEGER*4::I

C-------------------------First executable statement


      DO I=1,NUMNP

          HVEC(I) = SHUFFLED(2*I-1)
          CVEC(I) = SHUFFLED(2*I)

      END DO !I=1,NN

      END SUBROUTINE UNSHUFFLE_RHS
