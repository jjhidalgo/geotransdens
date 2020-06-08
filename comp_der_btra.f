      SUBROUTINE COMP_DER_BTRA
     &          (CAUX1    ,DENSREF  ,DERB     ,IBTCO    ,IOCONSRC
     &          ,IODENS   ,ITPTVAR  ,NPPNP    ,NUMNP    ,PARNP
     &          ,THETAT   ,WSPECHEAT)

********************************************************************************
*
* PURPOSE
*
*  Manages the computation of the derivatives of BTRA vector w.r.t. state variables.
*
*
* DESCRIPTION
*
*  Manages the computation of the derivatives of BTRA vector w.r.t. state variables.
*  The derivative is stored nodewise.
*
*
*
********************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IOCONSRC ,IODENS   ,ITPTVAR  ,NPPNP    ,NUMNP

      REAL*8::DENSREF ,THETAT   ,WSPECHEAT

      INTEGER*4::IBTCO(NUMNP)

      REAL*8::CAUX1(NUMNP)  ,DERB(NUMNP)  ,PARNP(NUMNP,NPPNP)

C------------------------- Internal

      INTEGER*4::I
      REAL*8::ALFAX,CONC

C------------------------- First executable statement

      IF (ANY(IBTCO(1:NUMNP).EQ.5)) THEN

          DO I=1,NUMNP

              IF (IBTCO(I).EQ.5) THEN

                  ALFAX = PARNP(I,6)
                  CONC = PARNP(I,4)

C------------------------- Whe solving heat transport the heat flow
C------------------------- has to be divided by the water specific heat
C------------------------- to be consistent with the rest of the equation.

              IF (ITPTVAR.EQ.1) THEN

                  IF (IODENS.EQ.0) THEN

C------------------------- If density is constant, water density appears
C------------------------- in this term (see TRANDENS Guía rápida, 3.7).
                      ALFAX = ALFAX/(DENSREF*WSPECHEAT)

                  ELSE

                      ALFAX = ALFAX/WSPECHEAT

                  END IF !IODENS.EQ.0

             END IF !ITPTVAR.EQ.1

C--------------------------- Constant density or heat transport or
C--------------------------- solute transport with variable density and
C--------------------------- no concentration sources [alf*(w*-w)].

                  IF (IODENS.EQ.0 .OR. ITPTVAR.EQ.1 .OR.
     &           ((IODENS.EQ.1 .AND. ITPTVAR.EQ.0 .AND. IOCONSRC.EQ.0)
     &               )) THEN

                      DERB(I) = DERB(I) + ALFAX*THETAT

C--------------------------- Solute transport with variable density and
C--------------------------- concentration sources [alf*(w*-w)*(1-w)].

                  ELSE IF (IODENS.EQ.1 .AND. ITPTVAR.EQ.0
     &                .AND. IOCONSRC.EQ.1) THEN


                      DERB(I) = DERB(I)
     &                        - ALFAX*THETAT*(1D0 + CONC - 2D0*CAUX1(I))

                  END IF !IODENS.EQ.0 ...

              END IF !IBTCO(I).EQ.5

          END DO !I=1,NUMNP

      ELSE

          DERB(1:NUMNP) = 0D0

      END IF !ANY(IBTCO(1:NUMNP).EQ.5


      END SUBROUTINE COMP_DER_BTRA
