      SUBROUTINE ORDER
     ;(NITF     ,NUMTIT   ,NUMTITC  ,NOOBSIT  ,TIT)     

***********************************************************************
* PURPOSE
*
* Orders integration times.
*
* DESCRIPTION
*
* In order to avoid having to go through all integration times after
* each simulation time step, the integration times have to be ordered.
* The subroutine is called once for each device. All integration times
* for the device are ordered so that TIT(IODEVICE(ND,2)) is the first
* integration time whereas TIT(IODEVICE(ND+1,2)-1) is the last
* integration time. ND is the device number and IODEVICE(ND,2) is the
* variable used to store the number of the first integration time for
* the device.
*
* The subroutine can be used to order any variable.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  NOOBSIT                Observation number to which an integration time       
*                         belongs to                                            
*  TIT                    Integration time                                      
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  NITF                   First integration time for current device
*  NUMTIT                 Total number of integration times                     
*  NUMTITC                Last integration time for current device
*
* INTERNAL VARIABLES: SCALARS
*
*  NIT1                   Loop counter dummy variable
*  NIT2                   Loop counter dummy variable
*  NREM                   Used to remember old relation between observation 
*                         and integration time number
*  REM                    Used to remember old integration time value
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR     03-2001     Revision and header
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION TIT(NUMTIT),NOOBSIT(NUMTIT)

*_______________________Organise integration times

      DO NIT1=NUMTITC,NITF,-1
        DO NIT2=NITF,NIT1

          IF(TIT(NIT1).LT.TIT(NIT2)) THEN

C______________________________ Integration times must be exchanged

            REM=TIT(NIT2)                           ! Exchange integration times
            TIT(NIT2)=TIT(NIT1)
            TIT(NIT1)=REM
            NREM=NOOBSIT(NIT2)              ! Exchange relation with observation
            NOOBSIT(NIT2)=NOOBSIT(NIT1)
            NOOBSIT(NIT1)=NREM

          ELSE
            CONTINUE
          ENDIF
        ENDDO
      ENDDO

      RETURN

      END
