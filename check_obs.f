      SUBROUTINE CHECK_CARD11
     ;(CORRCOEF ,IERROR   ,IOINV    ,IOWAR    ,IUOBS    ,MAINF    
     ;,ND       ,NDEVS    ,NODEV    ,NROW     ,STDEV    ,IODEVICE 
     ;,FILENAME)

***********************************************************************
* PURPOSE
*
* Verifies that data supplied in card E1.1 is OK.
*
* DESCRIPTION
*
* Verifies that:
*
* 1) The device number corresponds with previously read device numbers
* 2) Data type supplied value is allowed intervals
* 3) Prints warning if user does not want to use current device for calculations
* 4) Check if spatial integration method supplied value is within allowed 
*    intervals
* 5) Check if temporal integration method supplied value is within allowed 
*    intervals
* 6) Check of type of integr. and number of integr. times. Coherence between 
*    IODEVICE(ND,4) and IODEVICE(ND,5)
* 7) Check if STDEV supplied value is negative
* 8) Check if inverse covariance matrix supplied value is within allowed 
*    intervals
* 9) Data on correlation and type of cov. matrix is coherent
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Status for calc. of obs.                    
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*                         See user's guide for further details
*
* EXTERNAL VARIABLES: SCALARS
*
*  CORRCOEF               Correlation coefficient between observations
*  IERROR                 Current number of errors on input data                
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUOBS                  Unit number of OBS file                               
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  ND                     Counter variable (loop for devices). Just for checking
*  NDEVS                  Number of devices                                     
*  NODEV                  Number of current device (input data)
*  STDEV                  Default standard deviation for current dev. observat.
*  NROW                   Current record number                                 
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR     03-2001     Revision
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER FILENAME*12

      DIMENSION IODEVICE(NDEVS+1,10)

*_______________________ Check if device number, NODEV, is in sequence

      IF (NODEV.NE.ND) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR:'
     ;  //'DEVICE NO. OUT OF SEQUENCE (CARD E1.1).',NROW,0,IUOBS,1,7.02)

*_______________________ Check if data type value is within allowed intervals

      IF (IODEVICE(ND,1).LT.1.OR.IODEVICE(ND,1).GE.7) THEN
        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: VALUE OF'
     ;          //' IODATTYP OUT OF BOUNDS (CARD E1.1).'
     ;            ,NROW,0,IUOBS,1,7.02)
      ELSE IF (IODEVICE(ND,1).GT.2.AND.IODEVICE(ND,1).LT.6) THEN
        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: VALUE OF'
     ;          //' IODATTYP NOT IMPLEMENTED YET (CARD E1.1).'
     ;            ,NROW,0,IUOBS,1,7.02)
      ELSE
        CONTINUE
      END IF

*_______________________Warning if user does not want device to be included

        IF (IODEVICE(ND,2).EQ.0) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;    'WARNING: DEVICE NOT INCLUDED IN CALCULATIONS',
     ;     NROW,0,IUOBS,0,7.02)

*_______________________ Check if spatial integration method supplied value is 
*_______________________ within allowed intervals

      IF (IODEVICE(ND,3).LT.1.OR.IODEVICE(ND,3).GT.6) THEN
        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: VALUE'
     ;          //' OF IOCALTYP OUT OF BOUNDS (CARD E1.1).',
     ;             NROW,0,IUOBS,1,7.02)
      ELSE IF (IODEVICE(ND,3).GE.4.AND.IODEVICE(ND,3).LT.6) THEN
        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: VALUE OF'
     ;          //' IOCALTYP NOT IMPLEMENTED YET (CARD E1.1).',
     ;             NROW,0,IUOBS,1,7.02)
      ELSE
        CONTINUE
      ENDIF

*_______________________ Check if temporal integration method supplied value is
*_______________________ within allowed intervals
       
      IF (IODEVICE(ND,4).LT.0.OR.IODEVICE(ND,3).GT.2) THEN
        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: VALUE'
     ;          //' OF IOINTTYP OUT OF BOUNDS (CARD E1.1).',
     ;             NROW,0,IUOBS,1,7.02)
      ELSE 
        CONTINUE
      ENDIF

*_______________________Check of type of integr. and number of integr. times
*_______________________Coherence between IODEVICE(ND,4) and IODEVICE(ND,5)

                                             ! Point-in-time or simulation times

      IF(IODEVICE(ND,4).EQ.0.OR.IODEVICE(ND,4).EQ.1) THEN        

        IF(IODEVICE(ND,5).NE.0) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;    'WARNING: INCORRECT VALUE OF NUMINT IN CARD E1.1.'
     ;    ,NROW,0,IUOBS,0,7.02)

      ELSE IF(IODEVICE(ND,4).EQ.2) THEN                         ! Fixed interval

        IF(IODEVICE(ND,5).LT.2) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;    'ERROR: TEMPORAL INTEGRATION BUT NUMINT LESS THAN 2.'
     ;    //' CARD E1.1',NROW,0,IUOBS,1,7.02)

      ELSE
        CONTINUE
      ENDIF

*_______________________Check if STDEV supplied value is negative

      IF (IOINV.GT.0 .AND. STDEV.LE.0D0) 
     ;     CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;    'ERROR: DEFAULT STANDARD DEVIATION FOR DEVICE OBERVATIONS'
     ;    //' MUST BE GREATER THAN 0D0 (CARD E1.1)',NROW,0,IUOBS,1,7.02)

*_______________________Check if inverse covariance matrix supplied value is 
*_______________________within allowed intervals

      IF (IODEVICE(ND,6).LT.1.OR.IODEVICE(ND,6).GT.4) THEN
        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: VALUE'
     ;          //' OF IOCOVTYP OUT OF BOUNDS (CARD E1.1).',
     ;             NROW,0,IUOBS,1,7.02)
      ELSE IF (IODEVICE(ND,3).GE.2.AND.IODEVICE(ND,3).LE.4) THEN
        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: VALUE OF'
     ;          //' IOCOVTYP NOT IMPLEMENTED YET (CARD E1.1).',
     ;             NROW,0,IUOBS,1,7.02)
      ELSE
        CONTINUE
      ENDIF
        
*_______________________Check if data on correlation and type of cov. matrix is
*_______________________coherent

      IF (IODEVICE(ND,6).EQ.1.AND.CORRCOEF.NE.0D0) THEN

        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;    'WARNING: CORRELATION COEFF. WILL NOT BE USED. CARD E1.1',
     ;     NROW,0,IUOBS,0,7.02)

      ELSE IF (IODEVICE(ND,6).EQ.2.AND.CORRCOEF.EQ.0D0) THEN

        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;    'WARNING: CORRELATION COEFF. IS ZERO. CARD E1.1',
     ;     NROW,0,IUOBS,0,7.02)

      ELSE
        CONTINUE
      END IF

      RETURN

      END

********************************************************************************
********************************************************************************
********************************************************************************

      SUBROUTINE CHECK_GROUP2
     ;(IERROR   ,IOCALTYP ,IOWAR    ,IUOBS    ,MAINF    ,NROW
     ;,NU       ,NUMTNOD  ,NUMTU    ,IOUTYP   ,FILENAME)

***********************************************************************
* PURPOSE
*
* Verifies data supplied in group E2.
*
* DESCRIPTION
*
* Verifies for all units that
* 1) If unit is point-defined, spatial integration cannot be requested
* 2) All supplied values are within allowed intervals
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  IOUTYP                 Unit type                                             
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  IOCALTYP               Method of spat. integr.                     
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUOBS                  Unit number of OBS file                               
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NROW                   Current record number                                 
*  NU                     Last unit defining current device
*  NUMTNOD                Total number of nodes used for calculating obs.       
*  NUMTU                  Last unit defining previous device
*
* INTERNAL VARIABLES: SCALARS
*
* N                       Dummy counter for device units
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR     03-2001     Revision
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER FILENAME*12

      DIMENSION IOUTYP(NUMTNOD)

      DO N=NUMTU+1,NU

*_______________________Checks IOUTYP (and combinations with IOCALTYP)

        IF (IOUTYP(N).EQ.1.AND.    ! Unit is a point and spat. int. is requested
     ;     (IOCALTYP.EQ.2 .OR. IOCALTYP.EQ.3)) THEN

          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: ILLEGAL COMB.'
     ;      //' OF IOUTYP AND IODEVICE(ND,3)',NROW,0,IUOBS,1,7.02)

        ELSE IF (IOUTYP(N).LT.1.OR.IOUTYP(N).GT.4) THEN

          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: VALUE'
     ;      //' OF IOUTYP. OUT OF RANGE',NROW,0,IUOBS,1,7.02)

        ELSE IF (IOUTYP(N).EQ.3.OR.IOUTYP(N).EQ.4) THEN

          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'IOUTYP OPTION'
     ;      //' NOT IMPLEMENTED YET',NROW,0,IUOBS,1,7.02)
         
        ELSE
          CONTINUE
        ENDIF

      ENDDO

      RETURN

      END

********************************************************************************
********************************************************************************
********************************************************************************

      SUBROUTINE CHECK_CARD31
     ;(IERROR   ,INPWR    ,IOCALTYP ,IOWAR    ,IUOBS    ,MAINF
     ;,NBUF     ,NROW     ,NUF      ,NUMTBU   ,NUMTNOD  ,NUMTU
     ;,BUDAT    ,IOBUTYP  ,IOCALBU  ,FILENAME) 

********************************************************************************
*
* PURPOSE
*
* Verifies that data supplied in card E3.1 is OK.
*
* DESCRIPTION
*
* Verifies that
*
* 1) User has defined device as a point and basic units with spatial extens.
* 2) For each basic unit defining unit (and so that, defining device)
*    2.1) IOBUTYP=1. BU is a point defined by its coordinates. Checks that 
*         BUDAT(4) is zero
*    2.2) IOBUTYP=2. BU is a nodal point, defined by its number (global 
*         connectivity). Checks that BUDAT(1) is not wrong and BUDAT(2:4)=0
*    2.3) IOBUTYP=3. BU is an element. Same as 2.3
*    2.4) Other cases. IOBUTYP=4 is not implemented yet: ERROR
*                      IOBUTYP>4 is out of range
* 3) If user has defined the weights related to basic units (IOCALTYP=6), then
*    a check is done over all units defining device: if user has requested and 
*    spatial integration (IOCALBU=2 or 3) and the basic unit has spatial extens.
*    COHERENCE ERROR.
*  
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  BUDAT                  Data from group E3. One column for one basic unit     
*  IOBUTYP                Basic unit type                                       
*  IOCALBU                Calculation method for basic unit                     
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOCALTYP               Type of spatial integration of current device
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR
*  IUOBS                  Unit number of OBS file                               
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NBUF                   First basic unit defining device
*  NROW                   Current record number                                 
*  NUF                    Last basic unit defining device
*  NUMTBU                 Total number of basic units
*  NUMTNOD                Total number of nodes used for calculating obs.       
*  NUMTU                  Total number of units
*
* INTERNAL VARIABLES: SCALARS
*
*  NBU                    Dummy counter for basic units
*  NU                     Dummy counter for units
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR     03-2001     Revision
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER FILENAME*12

      DIMENSION BUDAT(4,NUMTNOD),IOBUTYP(NUMTNOD),IOCALBU(NUMTNOD)

*_______________________ Writes the main header for basic units information

      IF (INPWR.NE.0) WRITE(MAINF,3200)
 3200 FORMAT(//,1X,'BASIC UNITS INFORMATION',/,
     ;          1X,'===== ===== ===========',/)

*_______________Error if user defines the device as a point and basic units with
*_______________spatial extension 

      IF (IOCALTYP.EQ.1 .AND. IOBUTYP(NBUF).GT.1) 
     ;  CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,'ERROR: IOCALTYP AND'
     ;             //' IOBUTYP DO NOT CORRESPOND',NROW,0,IUOBS,1,7.02)

*_______________________ Loop over DEVICE BASIC UNITS

      DO NBU=NBUF,NUMTBU

*_______________________Basic units: Points (maybe def. of element)

        IF (IOBUTYP(NBU).EQ.0) THEN     ! Basic unit is point (coord.) + ELEMENT
          CONTINUE
        ELSE IF (IOBUTYP(NBU).EQ.1) THEN          ! Basic unit is point (coord.)
          IF(INT(BUDAT(4,NBU)+0.5).EQ.0) THEN
            IF (INPWR.NE.0) WRITE(MAINF,3220) 
     ;        NBU,BUDAT(1,NBU),BUDAT(2,NBU),BUDAT(3,NBU)
 3220       FORMAT(1X,' BASIC UNIT ',I5,' POINT ',3F10.3)

          ELSE

            CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: WRONG VALUE '
     ;//'              OF BUDAT(4). MUST BE ZERO',NROW,0,IUOBS,1,7.02)

          ENDIF

*_______________________Basic units: Nodes

        ELSE IF(IOBUTYP(NBU).EQ.2) THEN

          IF(INT(BUDAT(1,NBU)+0.5).GT.0 .AND.
     ;           BUDAT(2,NBU).EQ.0 .AND.
     ;           BUDAT(3,NBU).EQ.0 .AND.
     ;           BUDAT(4,NBU).EQ.0) THEN
            IF (INPWR.NE.0) WRITE(MAINF,3230) NBU,INT(BUDAT(1,NBU)+0.5)
 3230       FORMAT(1X,' BASIC UNIT ',I5,' NODAL POINT ',I5)

          ELSE
            CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: WRONG '
     ;//'INFORMATION SUPPLIED IN CARD E3.1',NROW,0,IUOBS,1,7.02)
          ENDIF

*_______________________Basic units: Elements

        ELSE IF(IOBUTYP(NBU).EQ.3) THEN

          IF(INT(BUDAT(1,NBU)+0.5).GT.0 .AND.
     ;           BUDAT(2,NBU).EQ.0 .AND.
     ;           BUDAT(3,NBU).EQ.0 .AND.
     ;           BUDAT(4,NBU).EQ.0 ) THEN
            IF (INPWR.NE.0) WRITE(1,3240) NBU,INT(BUDAT(1,NBU)+0.5)
 3240       FORMAT(1X,' BASIC UNIT ',I5,' ELEMENT ',I5)
          ELSE
            CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: NEGATIVE '
     ;//'ELEMENT NO. IN CARD E3.1 (IOBUTYP(NBU)=3)',NROW,0,IUOBS,1,7.02)
          ENDIF

*_______________________Basic units: Zone (and def. of parameter type)

        ELSE IF(IOBUTYP(NBU).EQ.4) THEN

          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: OPTION OF'
     ;          //' IOBUTYP NOT IMPLEMENTED YET (CARD E3.1).'
     ;            ,NROW,0,IUOBS,1,7.02)

*_______________________Error in IOBUTYP

        ELSE
          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: VALUE'
     ; //' OF IOBUTYP(NBU) OUT OF RANGE',NROW,0,IUOBS,1,7.02)       
        ENDIF
      ENDDO
      
*_______________________Check for errors related to user defined unit weights

      IF(IOCALTYP.EQ.6) THEN
        DO NU=NUF,NUMTU
          IF((IOCALBU(NU).EQ.2 .OR. IOCALBU(NU).EQ.3) .AND.
     ;      IOBUTYP(NU).LE.2) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;      'ERROR: ILLEGAL COMBINATION OF IOCALBU (CARD E2.2) AND'
     ;      //'IOBUTYP (<=2) (CARD E3.1)',NROW,0,IUOBS,1,7.02)
        ENDDO
      ENDIF

      RETURN

      END

********************************************************************************
********************************************************************************
********************************************************************************

      SUBROUTINE CHECK_CARD41
     ;(IDIMCOV  ,IERROR   ,IOINTTYP ,IOINV    ,IOWAR    ,IUOBS    
     ;,MAINF    ,NO       ,NROW     ,NUMTOBS  ,NUMTOBSC ,TMAX     
     ;,COVINV   ,IOTINT   ,TOBS     ,FILENAME)    

********************************************************************************
*
* PURPOSE
*
* Verifies that data supplied in card E4.1 is OK.
*
* DESCRIPTION
*
* Verifies that
*
* 1) Negative standard deviation for a given obs. is supplied
* 2) IOTINT is out of range (<0 or >2)
* 3) Temporal integration requested and initial time greater then final time
* 4) Point in time observation and more than one time is defined
* 5) Two observations for the same time
* 6) Two observations for the same interval 
* 7) Do observation times exceed last simulation time?
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix
*  FILENAME               Array containing names for input and output data files
*  IOTINT                 Temporal integration (=1) or temporal averaging (=2)
*  TOBS                   Time of observation
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCOV                Used to dimension array COVINV
*  IERROR                 Current number of errors on input data
*  IOINTTYP               Device time integration option (over units)
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR
*  IUOBS                  Unit number of OBS file
*  MAINF                  Unit number of the main output file (RES.OUT)
*  NO                     Number of last observation of current device
*  NROW                   Current record number
*  NUMTOBS                Total number of observations
*  NUMTOBSC               Number of last observation for previous device
*  TMAX                   Last simulation time (TIME (NINT)
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Observations dummy counter
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR     03-2001     Revision
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER FILENAME*12

      DIMENSION COVINV(IDIMCOV),TOBS(NUMTOBS,2),IOTINT(NUMTOBS)

      DO I=NUMTOBSC+1,NO

C______________________________ Negative standard deviation for a given obs.

        IF (IOINV.GT.0 .AND. COVINV(I).LT.0.D0) 
     ;     CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;    'ERROR: NEGATIVE STANDARD DEVIATION VALUE. CARD E4.1'
     ;    ,NROW,0,IUOBS,1,7.02)

C______________________________ IOTINT is out of range (<0 or >2)

        IF (IOTINT(I).LT.0.OR.IOTINT(I).GT.2) 
     ;     CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;    'ERROR: IOTINT VALUE OUT OF RANGE. CARD E4.1'
     ;    ,NROW,0,IUOBS,1,7.02)

C______________________________ Temporal integration requested and initial time 
C______________________________ greater then final time

        IF (IOINTTYP.GT.0 .AND. TOBS(I,1).GT.TOBS(I,2))
     ;     CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR: TOBSEND'
     ;     //' SMALLER THAN OR EQUAL TO TOBS (IOINTTYP>0)'
     ;     ,NROW,0,IUOBS,1,7.02)

C______________________________ Point in time observation and more than one 
C______________________________ time is defined

        IF(IOINTTYP.EQ.0 .AND. TOBS(I,1).NE.TOBS(I,2) .AND.
     ;      TOBS(I,2).NE.0) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;      'ERROR: TOBSEND NOT EQUAL TO TOBS (IOINTTYP=0).'
     ;      ,NROW,0,IUOBS,1,7.02)

      ENDDO

      DO I=NUMTOBSC+2,NO

C______________________________ Two observations for the same time

        IF(IOINTTYP.EQ.0 .AND. TOBS(I,1).EQ.TOBS(I-1,1)) THEN
          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'TWO OBSER-'
     ;            //' VATIONS FOR SAME TIME. (CARD E4.1)'
     ;              ,NROW,0,IUOBS,0,7.02)

C______________________________ Two observations for the same interval 

        ELSE IF((IOINTTYP.EQ.1 .OR. IOINTTYP.EQ.2) .AND.
     ;      TOBS(I,1).EQ.TOBS((I-1),1) .AND.
     ;      TOBS(I,2).EQ.TOBS((I-1),2)) THEN
          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'TWO OBSER-'
     ;            //'VATIONS FOR SAME TIME INTERVAL. (CARD E4.1)'
     ;              ,NROW,0,IUOBS,0,7.02)
        ENDIF

      ENDDO


*_______________________Do observation times exceed last simulation time?

      IF(IOINTTYP.EQ.0 .AND. TOBS(NO,1).GT.TMAX) THEN
        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'OBSERVATION'
     ;          //' TIME LARGER THAN MAX. SIMULATION TIME. (CARD E4.1)'
     ;            ,NROW,0,IUOBS,0,7.02)
      ELSE IF(IOINTTYP.LT.0 .AND. TOBS(NO,2).GT.TMAX) THEN
        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'END OF'
     ;          //'OBSERVATION TIME INTERVAL LARGER THAN MAX.'
     ;          //'SIMULATION TIME. (CARD E4.1)',NROW,0,IUOBS,0,7.02)
      ENDIF

      RETURN

      END

********************************************************************************
********************************************************************************
********************************************************************************

      SUBROUTINE CHECK_ALLOC
     ;(IERROR   ,IOWAR    ,IUOBS    ,MAINF    ,NROW     ,NUMTIT   
     ;,NUMTITC  ,NUMTNOD  ,NUMTNODC ,NUMTOBS  ,NUMTOBSC ,FILENAME)

***********************************************************************
* PURPOSE
*
* Verifies that the space allocated to the variables NUMTOBS, NUMTNOD
* and NUMTIT is sufficient.
*
* DESCRIPTION
*
* While reading the input data, the variables NUMTOBSC, NUMTNODC and
* NUMTITC are used to count the number of observations, nodes and
* integration times. When all data on observations has been read, these
* three variables are compared to the allocation variables NUMTOBS,
* NUMTNOD and NUMTIT.
* If too little or too much space is allocated to one of the variables,
* it is considered an error and the program stops.
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUOBS                  Unit number of OBS file                               
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NROW                   Current record number
*  NUMTIT                 Total number of integration times
*                         Value assigned by user
*  NUMTITC                Total number of integration times
*                         Counted while reading input
*  NUMTNOD                Total number of nodes used for calculating obs.
*                         Value assigned by user
*  NUMTNODC               Total number of nodes used for calculating obs.
*                         Counted while reading input
*  NUMTOBS                Total number of observations                          
*                         Value assigned by user
*  NUMTOBSC               Total number of observations
*                         Counted while reading input
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*
* INTERNAL VARIABLES: SCALARS
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR     03-2001     Revision and header
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER FILENAME*12

      IF(NUMTOBSC.NE.NUMTOBS) THEN
          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;   'ERROR: NUMTOBSC NOT EQUAL TO NUMTOBS. CHECK DIM FILE'
     ;   ,NROW,0,IUOBS,0,7.02)
          WRITE(*,*)' NUMTOBS MUST BE EQUAL TO', NUMTOBSC
          WRITE(MAINF,*)' NUMTOBS MUST BE EQUAL TO', NUMTOBSC
          STOP ' FATAL ERROR. CHECK FILE DIM.DAT'
      END IF

      IF(NUMTNODC.NE.NUMTNOD) THEN
          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;   'ERROR: NUMTNODC NOT EQUAL TO NUMTNOD. CHECK DIM FILE'
     ;   ,NROW,0,IUOBS,0,7.02)
          WRITE(*,*)' NUMTNOD MUST BE EQUAL TO', NUMTNODC
          WRITE(MAINF,*)' NUMTNOD MUST BE EQUAL TO', NUMTNODC
          STOP ' FATAL ERROR. CHECK FILE DIM.DAT'
      END IF

      IF(NUMTITC.NE.NUMTIT) THEN
          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;   'ERROR: NUMTITC NOT EQUAL TO NUMTIT. CHECK DIM FILE'
     ;   ,NROW,0,IUOBS,0,7.02)
          WRITE(*,*)' NUMTIT MUST BE EQUAL TO', NUMTITC
          WRITE(MAINF,*)' NUMTIT MUST BE EQUAL TO', NUMTITC
          STOP ' FATAL ERROR. CHECK FILE DIM.DAT'
      END IF

      RETURN
      END
