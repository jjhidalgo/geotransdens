      SUBROUTINE RHS_NODE_FLUX_INV
     &          (BETAC    ,CAUDAL   ,CAUX1    ,CREF     ,DENSREF
     &          ,DERC     ,IBTCO    ,IDIMDERC ,INEW     ,IODENS
     &          ,IOLD     ,NUMNP    ,NPAR     ,THETAT)

*********************************************************************
*
* PURPOSE  Adds to inverse problem RHS the contribution of nodal
*          fluxes in nodes with mass flow boundary condition.
*
*********************************************************************

      IMPLICIT NONE

C------------------------- External
      INTEGER*4::IDIMDERC,INEW,IODENS,IOLD,NUMNP,NPAR
      REAL*8::BETAC,CREF,DENSREF,DENS,THETAT

      INTEGER*4::IBTCO(NUMNP)
      REAL*8::CAUDAL(NUMNP),CAUX1(NUMNP),DERC(NUMNP,NPAR,IDIMDERC)

C------------------------- Internal

      INTEGER*4::I,IB
      REAL*8::CAUD,DENSNOD,THETAT1

C------------------------- First executable statement

      THETAT1 = THETAT - 1D0

      DO I=1,NUMNP

          IB = IBTCO(I)
          CAUD = CAUDAL(I)

          IF (IB.EQ.2 .OR. IB.EQ.3) THEN

              IF (CAUD.GT.0) THEN

                  IF (IODENS.EQ.1) THEN

                      DENSNOD = DENS(DENSREF,BETAC,CAUX1(I),CREF)

                  ELSE

                      DENSNOD = 1D0

                  END IF !IODENS.EQ.1

                  CAUD = THETAT1*DENSNOD*CAUD

                  DERC(I,1:NPAR,INEW) = DERC(I,1:NPAR,INEW) 
     &                                 + DERC(I,1:NPAR,IOLD)*CAUD

              END IF !CAUD.GT.0

          END IF !IB.EQ.2 .OR. IB.EQ.3


      END DO !I=1,NUMNP

      END SUBROUTINE RHS_NODE_FLUX_INV