      SUBROUTINE READ_EXTDRIFT_ZONES
     ;(IERROR   ,INPWR      ,IOWAR     ,IUGEO      ,MAINF    ,MXNZON_GS
     ;,NZON_GS  ,NEXDR_GS   ,POSZN_GS  ,EXDRZN_GS  ,FILENAME)

********************************************************************************
*
* PURPOSE Reads Card G6 (external drift for kriging purposes), concerning a
*         given group of zones.. Ext. drift can be either a locally varying 
*         mean (KTYPE_GS=2) or a external drift (up to four terms;KTYPE_GS=3).
*
* DESCRIPTION Algorithm plays as follows:
*             - Declaration of variables
*             - Writes main header at RES file (if allowed) and initialises 
*               external attributes statistics
*             - Big loop over zones of actual group. Reads card G6
*               - Step 1. Reads zone cog and value of drifts
*               - Step 2. Checks existence of the zone
*               - Step 3. Assignation and writing
*               - Step 4. Updates statistics
*             - Computation and writing of sec. attr. final statistics
*
* EXTERNAL VARIABLES: ARRAYS
*
*  POSZN_GS               Array containing coordinates of zones cog.
*  EXDRZN_GS              Array containing drifts
*  FILENAME               Array containing names for input and output
*                         data files
*
* INTERNAL VARIABLES: ARRAYS
*
*  AV                     Array containing averages of drift terms
*  EXTAUX                 Array for reading purposes
*  SS                     Array containing variance of drift terms
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data
*  INPWR                  Allows writing on MAIN FILE
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR
*  IUGEO                  GEO file unit number
*  MAINF                  Unit number of the main output file (RES.OUT)
*  MXNZON_GS              Number of zones of the most discretized group
*  NZON_GS                Number of zones defining actual group
*  NEXDR_GS               Number of external attributes to be read (1 if
*                         KTYPE_GS=2)
*
* INTERNAL VARIABLES: SCALARS
*
*  DIST2                  Distance
*  IZON                   Main loop counter dummy variable
*  IZON2                  Dummy counter variable
*  IDRIF                  Dummy counter of drift terms
*  ILEC                   Number of zone being defined (as read)
*  LEAUX                  Auxiliar string where the last read line of the
*                         current input file is stored
*  NROW                   Current record number
*  X                      X-coord for a zone-cog.
*  Y                      Y-coord for a zone-cog.
*  Z                      Z-coord for a zone-cog.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number
*  LEEL                   Returns a string value containing the current line
*                         of XXX FILE, if no coment appears.
*
* HISTORY: AAR     First coding      (Jan-2002)
*          AAR     Revision          (April-2002)
*          AAR     Revision. Inclusion of groups of zones (July-2003)
*
********************************************************************************

C______________________________________________________ Declaration of variables

      IMPLICIT NONE

      INTEGER*4 NZON_GS,MAINF,INPWR,IUGEO,IERROR,IOWAR,NEXDR_GS
     ;         ,MXNZON_GS
      REAL*8 EXDRZN_GS(MXNZON_GS,NEXDR_GS),POSZN_GS(MXNZON_GS,3)
      INTEGER*4 IZON,IZON2,NROW,IDRIF,IDUM

      CHARACTER FILENAME(18)*20,LEAUX*100,LEEL*100
      REAL*8 EXTAUX(4),XCOG,YCOG,ZCOG,X,Y,Z,DIST2,AV(4),SS(4)

C____________________ Writes main header (if allowed) and initialises statistics

      IF (INPWR.NE.0) WRITE(MAINF,2000)
 2000 FORMAT(//,17X,' EXTERNAL DRIFT AT BLOCKS ',/
     ;         ,17X,' ======== ===== == ====== ',//
     ;          ,7X,'CENTER OF GRAVITY',13X,'EXTERNAL DRIFT.',/
     ;          ,7X,'====== == =======',13X,'======== =====',/
     ;          ,9X,'X',9X,'Y',9X,'Z',12X,'DRI1.',12X,'DRI2.'
     ;          ,12X,'DRI3.',12X,'DRI4.',/
     ;          ,9X,'=',9X,'=',9X,'=',12X,'=====',12X,'====='
     ;          ,12X,'=====',12X,'=====')

      CALL ZERO_ARRAY (AV,4)
      CALL ZERO_ARRAY (SS,4)

C____________________ Reads card G8. Secondary attribute defined at all Y-blocks

      DO IZON=1,NZON_GS

                     ! Step 1. Reads zone cog and value of sec. attrib.

        LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
        READ (LEAUX,1000,ERR=9100) XCOG,YCOG,ZCOG,
     ;                             (EXTAUX(IDRIF),IDRIF=1,4)
 1000   FORMAT(8F10.0)

                     ! Step 2. Checks existence of zone


        DO IZON2=1,NZON_GS
          X=POSZN_GS(IZON2,1)
          Y=POSZN_GS(IZON2,2)
          Z=POSZN_GS(IZON2,3)
          DIST2=(X-XCOG)*(X-XCOG)+(Y-YCOG)*(Y-YCOG)+(Z-ZCOG)*(Z-ZCOG)
          IF (DIST2.LE.1.0E-6) GOTO 10
        END DO

        WRITE(MAINF,2100) IZON
 2100   FORMAT(//,' ERROR: DEFINED ZONE ',I5,' DOES NOT MATCH'
     ;            ' CALCULATED CENTERS OF GRAVITY OF ZONES'
     ;            ' BELONGING TO ACTUAL GROUP')

        WRITE(MAINF,2200) 
 2200   FORMAT(//,' LIST OF CALCULATED CENTERS OF GRAVITY',/,
     ;            ' ==== == ========== ======= == =======')

        DO IZON2=1,NZON_GS
           WRITE(MAINF,2300) IZON2,POSZN_GS(IZON2,1),POSZN_GS(IZON2,2)
     ;                            ,POSZN_GS(IZON2,3)
        END DO
 2300   FORMAT(I5,3F10.3)

        STOP ' FATAL ERROR ASSIGNING EXTERNAL DRIFT. CHECK FILE RES.OUT'

                                               ! Step 3. Assignation and writing

 10     DO IDRIF=1,NEXDR_GS
          EXDRZN_GS(IZON2,IDRIF)=EXTAUX(IDRIF)
        END DO

        IF (INPWR.NE.0) WRITE(MAINF,2400) 
     ;                  (POSZN_GS(IZON2,IDUM),IDUM=1,3)
     ;                 ,(EXDRZN_GS(IZON2,IDRIF),IDRIF=1,NEXDR_GS)
 2400   FORMAT(3E10.3,4(7X,E10.3))

C____________________________________________________ Step 4. Updates statistics

        DO IDRIF=1,NEXDR_GS
          AV(IDRIF)=AV(IDRIF)+EXDRZN_GS(IZON,IDRIF)
          SS(IDRIF)=SS(IDRIF)
     ;                 + EXDRZN_GS(IZON,IDRIF)*EXDRZN_GS(IZON,IDRIF)
        END DO
      
      END DO                                            ! Next block assignation

C_____________________________________________________ Computes final statistics

      DO IDRIF=1,4
        AV(IDRIF)=AV(IDRIF)/DFLOAT(NZON_GS)
        SS(IDRIF)=SS(IDRIF)/DFLOAT(NZON_GS)-AV(IDRIF)*AV(IDRIF)
      END DO

      IF (INPWR.NE.0) WRITE(MAINF,2500) AV(1),SS(1),AV(2),SS(2),AV(3)
     ;                                 ,SS(3),AV(4),SS(4)
 2500 FORMAT(//,17X,'EXTERNAL DRIFT STATISTICS',/
     ;         ,17X,'======== ===== ==========',//
     ;         ,26X,'AVERAGE',14X,'VARIANCE',/
     ;         ,26X,'=======',14X,'========',/
     ;         ,5X,'FIRST DRIFT:',7X,E10.3,12X,E10.3,/
     ;         ,5X,'SECOND DRIFT:',6X,E10.3,12X,E10.3,/
     ;         ,5X,'THIRD DRIFT:',7X,E10.3,12X,E10.3,/
     ;         ,5X,'FOURTH DRIFT:',6X,E10.3,12X,E10.3)

      RETURN

C______________________________________________________________ Fatal error call

 9100 CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING G6 ',NROW,1,IUGEO,2,10.1)

      RETURN
      END 

