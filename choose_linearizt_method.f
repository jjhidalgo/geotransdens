      SUBROUTINE CHOOSE_LINEARIZT_METHOD
     &          (IODENS    ,ITERCHNGFL ,ITERCHNGGL ,ITERCHNGTR
     &          ,ITERCONVFL,ITERCONVTR ,ITERFL     ,ITERGL     ,ITERTR
     &          ,LINMET)

************************************************************************
*
* PURPOSE
* 
*  Manages the change of solving method for flow, transport or both.
*
*
* DESCRIPTION
*
*  Manages the computation of the  solving method according
*  current number of iterations.
*
*  The allowed combinations are:
*
*   Flow      Transport       Coupled
*  ------    -----------     ---------
*    1            1              0
*    2            2              0
*    1            2              0
*    2            1              0
*    2            2              2
*    0            0              2
*    0            1              2
*    0            2              2
*    1            0              2
*    2            0              2
*  -----------------------------------
*
*  The forbbiden combinatios are:
*
*   Flow      Transport       Coupled
*  ------    -----------     ---------
*    0            0              0
*    0            0              1
*    1            1              2
*    1            2              1
*    2            1              1
*  ----------------------------------
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IODENS       Density dependent flow.
*                 0. No
*                 1. Yes
*  ITERCHNGFL   Number of iterations after which method must change (flow).
*  ITERCHNGGL   Number of iterations after which method must change (coupled).
*  ITERCHNGTR   Number of iterations after which method must change (tpt.).
*  ITERCONVFL   Number of iterations after which local iterations must cease (flow)
*  ITERCONVTR   Number of iterations after which local iterations must cease (tpt.)
*  ITERFL       Number of iteration in current time step (flow).
*  ITERGL       Number of iteration in current time step (coupled).
*  ITERTR       Number of iteration in current time step (tpt.).
*
*
* EXTERNAL VARIABLES: ARRAYS
*
* LINMET        Array containing initial linearization method and current one.
*                 Rows:
*                   1. Flow
*                   2. Transport
*                   3. Coupled.
*
*                 Columns:
*                   1. Initial method.
*                   2. Current method.
*
*                 Values:
*                   0. Not solve (<=> Compute matrices).
*                   1. Picard
*                   2. Newton

* HISTORY
*
*   JHG         12-2003        First coding.
*
************************************************************************

      IMPLICIT NONE

      INTEGER*4::IODENS,ITERCHNGFL,ITERCHNGGL,ITERCHNGTR,ITERCONVFL,
     &           ITERCONVTR,ITERFL,ITERGL,ITERTR

      INTEGER*4::LINMET(3,2)



C------------------------- Linerarization method for flow is set.

      IF (LINMET(3,2).EQ.1) LINMET(3,2)=0

      IF (ITERFL.GT.ITERCHNGFL) THEN

C------------------------- Flow iterations can only cease if coupled
C------------------------- problem is being solved.	
          IF (ITERFL.LE.ITERCONVFL .AND. IODENS.EQ.1) THEN

              LINMET(1,2) = 0

          ELSE

              IF (LINMET(1,2).EQ.1) THEN

                  LINMET(1,2) = 2

              ELSE

                  LINMET(1,2) = 1

              END IF ! LINMET(1,2).EQ.1

          END IF ! ITERFL.LE.ITERCONVFL .AND. IODENS.EQ.1

      END IF ! ITERFL.GT.ITERCHNGFL


C------------------------- Linerarization method for trasnport is set.

      IF (ITERTR.GT.ITERCHNGTR) THEN

          IF (ITERTR.LE.ITERCONVTR .AND. IODENS.EQ.1) THEN

              LINMET(2,2) = 0

          ELSE

              IF (LINMET(2,2).EQ.1) THEN

                  LINMET(2,2) = 2

              ELSE

                  LINMET(2,2) = 1

              END IF ! LINMET(2,2).EQ.1

          END IF !ITERTR.LE.ITERCONVTR .AND. IODENS.EQ.1

      END IF ! ITERTR.GT.ITERCHNGTR


C------------------------- Linerarization method for coupled problem is set.

C------------------------- First incompatibilities are checked.
      IF (IODENS.EQ.1) THEN

          IF (LINMET(1,2).EQ.1 .AND. LINMET(2,2).EQ.1) THEN

              LINMET(3,2) = 0

          ELSE IF (LINMET(1,2).EQ.0 .OR. LINMET(2,2).EQ.0) THEN

              LINMET(3,2) = 2

C------------------------- If there is no incompatibility, the method
C------------------------- for the coupled problem is set.
          ELSE ! FLOW=1, TPT=2 or FLOW=2, TPT=1

              IF (ITERGL.GT.ITERCHNGGL) THEN

                  IF (LINMET(3,2).EQ.0) THEN

                      LINMET(3,2) = 2

                  ELSE

                      LINMET(3,2) = 0

                  END IF ! LINMET(3,2).EQ.1

              END IF ! ITERGL.GT.ITERCHNGGL

          END IF ! LINMET(1,2).EQ.1 .AND. LINMET(2,2).EQ.1

      END IF ! IODENSE.EQ.1

      END SUBROUTINE CHOOSE_LINEARIZT_METHOD

