      SUBROUTINE CALC_FLOWDENS
     &          (BETAC    ,CAUX1    ,CREF     ,DENSREF  ,FLOW
     &          ,NODE     ,IODENS   ,NUMNP    ,NPPNP    ,PARNP)

********************************************************************************
*     Computes the density of the flow at a given node.
********************************************************************************

      IMPLICIT NONE

C-------------------- External

      INTEGER*4::IODENS   ,NODE     ,NUMNP    ,NPPNP
      REAL*8::BETAC    ,CREF     ,DENSREF  ,FLOW

      REAL*8::DENS

      REAL*8::CAUX1(NUMNP)  ,PARNP(NUMNP,NPPNP)

C-------------------- Internal
      
      REAL*8::CNODE,FLOWDENS

C-------------------- First executable statement

      IF (IODENS.EQ.1) THEN

          IF (FLOW.GT.0) THEN

              CNODE = PARNP(NODE,4)

          ELSE

              CNODE = CAUX1(NODE)

          END IF !FLOW.GT.0

          FLOWDENS = DENS(DENSREF,BETAC,CNODE,CREF)

          FLOW = FLOWDENS*FLOW

      END IF !IODENS.EQ.0

      END SUBROUTINE CALC_FLOWDENS
