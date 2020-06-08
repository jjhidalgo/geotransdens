      SUBROUTINE BUWEIGHT_OBS
     ;(IERROR   ,IFLAG    ,INPWR    ,IOWAR    ,IUOBS    ,MAINF
     ;,NBU      ,NROW     ,NUMTNOD  ,NUMW     ,WTOBSBU  ,FILENAME)

********************************************************************************
*
* PURPOSE
*
* This subroutines counts the number of basic unit weights read from a
* line in card E2.3.
*
* DESCRIPTION
*
* Each line of the input file can contain a maximum of 7 weights. The 7
* values read are analysed one by one (from left to right). When the
* first zero weight is reached, it is assumed that the previous weight
* is the last weight supplied by the user. If all 7 weights are
* different from zero, a new line is read. Otherwise, a flag value is
* read (IFLAG). The flag values are:
* -1 if group E2 is ended but should be read again and therefore, there are
*    more units defining current device
* -2 if group E2 is ended and should not be read again and, therefore, the last 
*    basic unit (weights) read are related to last unit defining current device
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  WTOBSBU                Weight for basic unit                                 
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  IFLAG                  End of group marker
*  INPWR                  Allows writing on MAIN FILE                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUOBS                  Unit number of OBS file                               
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NBU                    Basic unit number
*  NROW                   Current record number                                 
*  NUMTNOD                Total number of nodes used for calculating obs.       
*  NUMW                   Number of basic unit weights for current unit
*
* INTERNAL VARIABLES: SCALARS
*
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*  NW                     Number of weights counter
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR     03-2001     Revision
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER LEAUX*100,LEEL*100,FILENAME*12

      DIMENSION WTOBSBU(NUMTNOD)

      NUMW=NBU                          ! Used temporarily to remember init. NBU

      DO WHILE(IFLAG.EQ.0)

*_______________________Reads 7 weights from card E2.3

        LEAUX=LEEL(FILENAME,IUOBS,MAINF,NROW,INPWR)
        READ(LEAUX,1000,ERR=9300) WTOBSBU(NBU),
     ;WTOBSBU(NBU+1),WTOBSBU(NBU+2),WTOBSBU(NBU+3),WTOBSBU(NBU+4),
     ;WTOBSBU(NBU+5),WTOBSBU(NBU+6)
 1000   FORMAT(7F10.0)

*________________NBU is the number of the first weight read
*________________NW is the number of the last weight read

        NW=NBU+6
        DO WHILE(NBU.LE.NW)

*________________Data read is a weight => write value and go on

          IF (WTOBSBU(NBU).GT.0) THEN
            IF(INPWR.NE.0) THEN
              WRITE(MAINF,3160) NBU,WTOBSBU(NBU)
 3160         FORMAT (' BASIC UNIT  ',I5,'  WEIGHT= ',E10.5)
            ENDIF

*________________Data read is not a weight => stop and get flag value

          ELSE IF (WTOBSBU(NBU).LE.0) THEN

*________________Line contained weights. Flag has to be read from next line

            IF(NBU.NE.NW-6) THEN
              LEAUX=LEEL(FILENAME,IUOBS,MAINF,NROW,INPWR)
              READ(LEAUX,1160,ERR=9300) IFLAG
 1160         FORMAT(I5)

*________________First value read from line was a flag

            ELSE

              IFLAG=WTOBSBU(NBU)

            ENDIF

*________________Flag is OK (equal to -1 or -2)

            IF (IFLAG.EQ.-1 .OR. IFLAG.EQ.-2) THEN
              NW=NBU-1                  ! Update of counter of no. of BU weights
              IF (NW.EQ.0) CALL ERROR
     ;          (IERROR,IOWAR,MAINF,FILENAME,'ERROR: LINE EMPTY'
     ;        //' OR CONTAINS ZERO WEIGHTS',NROW,0,IUOBS,1,7.02)
              WRITE(MAINF,3165) NW
 3165         FORMAT (/,' Number of weights read: ',I5)
              GOTO 10

*________________Flag read is not OK (not equal to -1 or -2)

            ELSE
              CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;         'ERROR: NO MORE WEIGHTS. BUT NONE OF THE VALUES -1 OR -2'
     ;         //' SUPPLIED TO END CARD',NROW,0,IUOBS,1,7.02)
            ENDIF

          ENDIF                                                 ! WTOBSBU.LT.0D0
          NBU=NBU+1                        ! Number of next value to be analysed
        ENDDO                                                        ! NBU,LE.NW
      ENDDO                                                            ! IFLAG=0

   10 NUMW=NW-NUMW+1                   ! Number of BU weights read for present U
      NBU=NW+1                 ! Update NBU to number of next weights to be read

      IFLAG=IFLAG+1    ! Return flag value to initial value (if IFLAG eqauls -1)

      RETURN

 9300 CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'GENERIC FORTRAN ERROR'
     ;   //' READING CARD E2.3',NROW,0,IUOBS,1,7.02)

      RETURN

      END
