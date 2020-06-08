      SUBROUTINE ENTDAT_DENSITY
     &          (BETAC    ,CREF     ,DENSREF  ,FILENAME ,IERROR
     &          ,INPWR    ,IODENS_INI         ,IOWAR    ,ITPTVAR
     &          ,IUPAR    ,MAINF    ,NROW     ,TEMPREF  ,VISCREF
     &          ,WSPECHEAT,WTHERMCON)
*******************************************************************************
*
*  PURPOSE
*
*    To read the density and viscosity related coeficients from 
*    the par.dat file
*  
*  FIRST CODING Luit Jan 
*
*******************************************************************************


      IMPLICIT NONE

C------------------------- External

      CHARACTER FILENAME(20)*20, LEEL*100

      INTEGER*4::IERROR     ,INPWR      ,IODENS_INI ,IOWAR
     &          ,ITPTVAR    ,IUPAR      ,MAINF      ,NROW

      REAL*8::BETAC    ,CREF     ,DENSREF  ,TEMPREF  ,VISCREF
     &       ,WTHERMCON,WSPECHEAT

C------------------------- Internal

      CHARACTER::MSG*100,LEAUX*100

C------------------------- First executable statement

C------------------------- Read density and viscosity related coeficients

      LEAUX = LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)

      READ(LEAUX,100) DENSREF,CREF,BETAC,VISCREF,TEMPREF,WTHERMCON
     &               ,WSPECHEAT


C------------------------- If solving energy transport with constant
C------------------------- density, some of the parameters values are
C------------------------- corrected.

      IF (IODENS_INI.EQ.0 .AND. ITPTVAR.EQ.1) THEN

C          IF (ABS(DENSREF - 1D0).GT.0D0) THEN
C
C              DENSREF = 1D0
C           MSG = ' ENERGY TRANSPORT WITH CONSTANT DENSITY.'
C    &            // ' REFERENCE DENSITY (DENSREF) CHANGED TO 1.0'
C           CALL ERROR
C    &            (IERROR,IOWAR,MAINF,FILENAME,MSG,NROW,0,IUPAR,0,0.0)
C
C          END IF !ABS(DENSREF - 1D0).GT.0D0

          IF (ABS(BETAC).GT.0D0) THEN

              BETAC = 0D0
              MSG = ' ENERGY TRANSPORT WITH CONSTANT DENSITY.'
     &            // ' DENSITY FACTOR (BETAC) CHANGED TO 0.0'
              CALL ERROR
     &            (IERROR,IOWAR,MAINF,FILENAME,MSG,NROW,0,IUPAR,0,0.0)

          END IF !ABS(BETAC).GT.0D0

      END IF !IODENS_INI.EQ.0 .AND. ITPTVAR.EQ.1

  100 FORMAT(7F10.3)


C------------------------- write density and viscosity related coeficients to
C------------------------- file.
 
200   FORMAT(/,
     &'REFERENCE DENSITY...................................:',F10.3,/,
     &'REFERENCE MASS FRACTION.............................:',F10.3,/,
     &'BETA FOR DENSITY....................................:',F10.3,/,
     &'REFERENCE VISCOSITY.................................:',F10.3,/,
     &'REFERENCE TEMPERATURE...............................:',F10.3,/,
     &'WATER THERMAL CONDUCTIVITY..........................:',F10.3,/,
     &'WATER SPECIFIC HEAT.................................:',F10.3,/)

      IF (INPWR.NE.0) WRITE(MAINF,200) DENSREF,CREF,BETAC,VISCREF
     &                                ,TEMPREF,WTHERMCON,WSPECHEAT
 

      END SUBROUTINE ENTDAT_DENSITY