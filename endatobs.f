      SUBROUTINE ENDATOBS
     ;(IDIMCOV  ,IERROR   ,INPWR    ,IOINV    ,IOWAR    ,IUOBS    
     ;,LMXNDL   ,MAINF    ,NDEVS    ,NROW     ,NUMEL    ,NUMNP    
     ;,NUMTIT   ,NUMTNOD  ,NUMTOBS  ,TMAX     ,AREA     ,BUDAT    
     ;,COVINV   ,DEVNAME  ,EXTNBU   ,INDEXNOD ,IOBUTYP  ,IOCALBU  
     ;,IODEVICE ,IOTINT   ,MEASTYP,IOUTYP   ,KXX ,LNNDEL,LTYPE    
     ;,NBUW     ,NOBUF    ,NOOBSIT  ,TIT      ,TOBS     ,VOBS     
     ;,WTOBSBU  ,WTOBSN   ,WTOBST   ,WTOBSU   ,X        ,Y        
     ;,Z        ,FILENAME)

********************************************************************************
*
* PURPOSE
*
* Reads observation related variables.
*
* DESCRIPTION
*
* Reads the file rootobs.dat which contains groups E1, E2, E3 and E4.
*
* Group E1 contains variables which are common for all observations
* at the device.
* Group E2 contains variables which describe the unit type and weights
* for units and basic units.
* Group E3 contains variables which describe the basic unit type and
* their location.
* Group E4 contains information on the observation values and temporal
* extension.
*
* A number of subroutines are called which make sure that the
* observation data supplied by the user is coherent. Furthermore, two
* subroutines, which determine the spatial and temporal weights used to
* determine the simulated which correspond to the observations, are
* called.
*
* Throughout the subroutine, U and BU are used in comments as
* abbreviations for unit and basic unit, respectively. IT is used for
* integration time.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  BUDAT                  Data from group E3. One column for one basic unit     
*  COVINV                 Inverse of the covariance matrix                      
*  DEVNAME                Device name                                           
*  EXTNBU                 Measure (length, area or volume) of basic unit        
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  INDEXNOD               Index relating nodes                                  
*  IOBUTYP                Basic unit type                                       
*  IOCALBU                Calculation method for basic unit                     
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Status for calc. of obs.                    
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*  IOTINT                 Temporal integration (=1) or temporal averaging (=2)  
*  IOUTYP                 Unit type                                             
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  NBUW                   Number of basic unit weights                          
*  NOBUF                  Number of first basic unit                            
*  NOOBSIT                Observation number to which an integration time       
*                         belongs to                                            
*  TIT                    Integration time                                      
*  TOBS                   Time of observation                                   
*  VOBS                   Observation value                                     
*  WTOBSBU                Weight for basic unit                                 
*  WTOBSN                 Weight for node used to calculate observation         
*  WTOBST                 Weight for integration time                           
*  WTOBSU                 Weight for unit                                       
*  X                      X-coord for a given node                              
*  Y                      Y-coord for a given node                              
*  Z                      Z-coord for a given node                              
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCOV                Used to dimension array COVINV
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUOBS                  Unit number of OBS file                               
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NDEVS                  Number of devices
*  NROW                   Current record number                                 
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NUMTIT                 Total number of integration times                     
*  NUMTNOD                Total number of nodes used for calculating obs.       
*  NUMTOBS                Total number of observations                          
*  TMAX                   Last simulation/observation time (EQUAL TO TIME(NINT))
*
* INTERNAL VARIABLES: SCALARS
*
*  CORRCOEF               Correlation coefficient between two observations
*  COV                    
*  IFLAG                  If different from 0, it writes storage partition      
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*  NBU                                                                          
*  ND                                                                           
*  NO                                                                           
*  NOB                                                                          
*  NODEV                                                                        
*  NOOBS                                                                        
*  NOUPRES                                                                      
*  NOUPREV                                                                      
*  NU                                                                           
*  NUF                                                                          
*  NUMTBU                                                                       
*  NUMTITC                                                                      
*  NUMTNODC                                                                     
*  NUMTOBSC                                                                     
*  NUMTU                                                                        
*  STDEV                                                                        
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ASS_COVINV                                                                   
*  BUWEIGHT_OBS                                                                 
*  CHECK_ALLOC                                                                  
*  CHECK_CARD11                                                                 
*  CHECK_CARD31                                                                 
*  CHECK_CARD41                                                                 
*  CHECK_GROUP2                                                                 
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*  ORDER                                                                        
*  SPAT_WEIGHT_OBS                                                              
*  TEMP_WEIGHT_OBS                                                              
*  WRIT_DEVICE                                                                  
* HISTORY
*
*     AMS        1988     First coding
*     SCR      4-1997     Revision and verification
*     AMS      1-1998     Common elimination and addition of header
*     CK      11-1999     Re-build to accomodate new obs. structure
*     AAR      7-2001     Revision and header inclusion
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)

      CHARACTER*10 DEVNAME(NDEVS)*10,LEAUX*100,LEEL*100,FILENAME(20)*20

      DIMENSION IODEVICE(NDEVS+1,10),IOUTYP(NUMTNOD),WTOBSU(NUMTNOD),
     ;IOCALBU(NUMTNOD),WTOBSBU(NUMTNOD),IOBUTYP(NUMTNOD),VOBS(NUMTOBS),
     ;TOBS(NUMTOBS,2),NOBUF(NUMTNOD+1),NBUW(NUMTNOD),NOOBSIT(NUMTIT),
     ;BUDAT(4,NUMTNOD),TIT(NUMTIT),IOTINT(NUMTOBS),COVINV(IDIMCOV)
     ;,X(NUMNP),Y(NUMNP),MEASTYP(NUMTOBS)

*_______________________Initialisations

      NUMTU=0                                               ! Total number of Us
      NUMTBU=0                                             ! Total number of BUs
      NUMTOBSC=0      ! Total number of observations (variable for verification)
      NUMTITC=0                ! Total number of ITs (variable for verification)
      NUMTNODC=0     ! Init. counter of total number of nodes (for verification)
      IFLAG=0                          ! Init. flag used to signal end of groups
      NOBUF(1)=1                                     ! First BU defining first U

*_______________________Writes header for observations

       IF (INPWR.NE.0) WRITE(MAINF,3000)
 3000  FORMAT(///,21X,'OBSERVATIONS RELATED VARIABLES ',/,
     ;            21X,'============ ======= =========')

*_______________________Begins loop for devices

      DO ND=1,NDEVS

*_______________________GROUP E1. 

*_______________________Reads card E1.1.General information on current device

        LEAUX=LEEL(FILENAME,IUOBS,MAINF,NROW,INPWR)
        READ(LEAUX,1000,ERR=9000) NODEV,DEVNAME(ND),IODEVICE(ND,1),
     ;  IODEVICE(ND,2),IODEVICE(ND,3),IODEVICE(ND,4),IODEVICE(ND,5),
     ;  STDEV,IODEVICE(ND,6),CORRCOEF,IODEVICE(ND,9)
 1000   FORMAT(I5,A10,5I5,F10.0,I5,F10.0,I5)

*_______________________Checks card E1.1

        CALL CHECK_CARD11
     ;(CORRCOEF ,IERROR   ,IOINV    ,IOWAR    ,IUOBS    ,MAINF    
     ;,ND       ,NDEVS    ,NODEV    ,NROW     ,STDEV    ,IODEVICE 
     ;,FILENAME)

*_______________________Assigns default values

        IF (IODEVICE(ND,3).EQ.0) IODEVICE(ND,3)=1         ! Point in space meas.
        IF (IODEVICE(ND,6).EQ.0) IODEVICE(ND,6)=1         ! Diagonal cov. matrix
        IF (IODEVICE(ND,9).EQ.0) IODEVICE(ND,9)=1             ! Flow/tpt prob.=1

*_______________________Writes card E1.1 for current device

        IF (INPWR.NE.0) 
     ;    WRITE(MAINF,3005) NODEV,DEVNAME(ND),
     ;IODEVICE(ND,1),IODEVICE(ND,2),IODEVICE(ND,3),IODEVICE(ND,4),
     ;IODEVICE(ND,5),STDEV,IODEVICE(ND,6),CORRCOEF,IODEVICE(ND,9)

 3005 FORMAT(////,10X,'GENERAL INFORMATION ON DEVICE NUMBER: ',I5,/,
     ;         10X,'==========================================='/,
     ;       5X,'DEVICE NAME: .............................=',A10,/,
     ;       5X,'DATA TYPE MEASURED: ......................=',I5,/,
     ;       5X,'FLAG FOR INCLUSION OF DATA: ..............=',I5,/,
     ;       5X,'METHOD OF SPATIAL INTEGRATION: ...........=',I5,/,
     ;       5X,'METHOD OF TEMPORAL INTEGRATION: ..........=',I5,/,
     ;       5X,'NUMBER OF INTEGR. TIMES: .................=',I5,/,
     ;       5X,'DEFAULT STANDARD DEVIATION: ..............=',E10.3,/,
     ;       5X,'TYPE OF INVERSE COVARIANCE MATRIX: .......=',I5,/,
     ;       5X,'COEFFICIENT FOR CORRELATION: .............=',E10.3,/,
     ;       5X,'NUMBER OF ASSOCIATED FLOW/TPT PROBLEM: ...=',I5,/)
 
*_______________________GROUP E2. Should be ommited if DEVICE===Single point. If
*_______________________this is the case, TRANSIN III methodology will be used.

        IF(IODEVICE(ND,3).NE.1) THEN              ! Device is not a single point

          NU=NUMTU                               ! Init. of counter of no. of Us
          NBU=NUMTBU+1                          ! Init. of counter of no. of BUs

          DO WHILE (IFLAG.EQ.0)                            ! End of group marker

            NU=NU+1                               ! Updates counter of no. of Us
            NBUW(NU)=0      ! Init. of counter of no. of BU weights of current U

*_______________________Reads card 2.1

            LEAUX=LEEL(FILENAME,IUOBS,MAINF,NROW,INPWR)
            READ(LEAUX,1100,ERR=9100) IOUTYP(NU)                     ! Unit type
 1100       FORMAT(I5)

            IF (IOUTYP(NU).LE.0) THEN
              IF(IOUTYP(NU).EQ.-2) IFLAG=-1              ! The group is finished
              NU=NU-1                 ! No data in line. Counter updating undone
            ELSE

              IF(IODEVICE(ND,3).EQ.6) THEN               ! User defined weights.
                                                        ! Need to read card E2.2

*_______________________Reads card 2.2. Unit weight and type of integration of
*_______________________basic units belonging to current unit.

                LEAUX=LEEL(FILENAME,IUOBS,MAINF,NROW,INPWR)
                READ(LEAUX,1120,ERR=9200) WTOBSU(NU),IOCALBU(NU)
 1120           FORMAT(F10.0,I5)

*_______________________Check whether E2 contains
*_______________________more data on the present U

                IF(WTOBSU(NU).LT.0) THEN                    ! Read next group E2
                  IF(WTOBSU(NU).EQ.-2) IFLAG=-1              ! Go on to group E3
                ELSE

                  IF(IOCALBU(NU).EQ.6) THEN  ! Users defines basic units weights

*_______________________Writes header

                    IF (INPWR.NE.0) WRITE(MAINF,1130) NU
 1130   FORMAT (//,' UNIT NUMBER ',I5,'   BASIC UNIT WEIGHT',
     ;           /,' ==== ====== ',5X,'   ===== ==== ======',/)

*_______________________Reads card E2.3 once or more. Determines for each time 
*_______________________the number of weights read

                    CALL BUWEIGHT_OBS
     ;(IERROR   ,IFLAG    ,INPWR    ,IOWAR    ,IUOBS    ,MAINF
     ;,NBU      ,NROW     ,NUMTNOD  ,NBUW(NU) ,WTOBSBU  ,FILENAME)

                  ENDIF                                       ! IOCALBU(NU).EQ.6
                ENDIF                                              ! WTOBSU <> 0
              ENDIF                                           ! IODEVICE(ND,3)=6
            ENDIF                          ! End of division depending on IOUTYP
          ENDDO                                   ! End of loop reading group E2
          IFLAG=0                           ! Return flag value to initial value

*_______________________Check group E2

        CALL CHECK_GROUP2
     ;(IERROR   ,IODEVICE(ND,3) ,IOWAR    ,IUOBS    ,MAINF    ,NROW
     ;,NU       ,NUMTNOD        ,NUMTU    ,IOUTYP   ,FILENAME)

        ENDIF                                                 ! IODEVICE(ND,3)=1

*_______________________GROUP E3. Basic units information (current device)

        NBU=NUMTBU+1                                       ! Updates BUs counter
        NU=NUMTU                 ! Init. of Us counter (last unit of prev. dev.)
        NUF=NU+1                               ! First U defining current device

*_______________________Reads card E3.1. Basic units type and information

        LEAUX=LEEL(FILENAME,IUOBS,MAINF,NROW,INPWR)
        READ(LEAUX,1190,ERR=9400)IOBUTYP(NBU),BUDAT(1,NBU),BUDAT(2,NBU),
     ;                                        BUDAT(3,NBU),BUDAT(4,NBU)
 1190   FORMAT(I5,4F10.0)

        IF(IODEVICE(ND,3).GE.2) THEN                     ! Device is not a point

          DO WHILE(IFLAG.GE.0)
            IF(IOBUTYP(NBU).GE.-1) THEN                       ! Device not ended
              IF(IOBUTYP(NBU).GT.0) THEN         ! The line contains data on BUs
                NOUPREV=NOUPRES          ! Store to which U the prev. BU belongs
                NOUPRES=NU+1             ! Store to which U the pres. BU belongs

*_______________________ Checks if basic unit types are identical , if those 
*_______________________ BU's define the same unit

                IF(NBU.NE.NUMTBU+1 .AND.                ! Comparison can be done
     ;            IOBUTYP(NBU).NE.IOBUTYP(NBU-1) .AND.  ! Different types of BUs
     ;            NOUPRES.EQ.NOUPREV)                                ! Same unit
     ;            CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'DIFF. BASIC'
     ;              //' UNIT TYPES USED FOR UNIT',NROW,0,IUOBS,1,7.02)
                NBU=NBU+1            !Update of counter of no. of BUs for device

              ELSE IF(IOBUTYP(NBU).EQ.-1) THEN     !U ended but device not ended
                NU=NU+1                          !Update of counter of no. of Us
                IOBUTYP(NU)=IOBUTYP(NBU-1)    !Assign BU type to U instead of BU

*_______________________Checks if same number of BUs supplied in E2 and E3

                IF(IOCALBU(NU).EQ.4 .AND. NBU-NUMTBU-1.NE.NBUW(NU)) 
     ;             CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR:'
     ;              //' INCONSISTENCY BETWEEN NO. OF BASIC UNITS'
     ;              //' FROM CARDS E2.3 AND E3.1',NROW,0,IUOBS,1,7.02)

                NUMTBU=NBU-1                        ! Update total number of BUs
                NOBUF(NU+1)=NUMTBU+1                  ! First BU defining next U
              ENDIF

*_______________________Reads card E3.1 again

              LEAUX=LEEL(FILENAME,IUOBS,MAINF,NROW,INPWR)
              READ(LEAUX,1200,ERR=9400)IOBUTYP(NBU),BUDAT(1,NBU),
     ;        BUDAT(2,NBU),BUDAT(3,NBU),BUDAT(4,NBU)
 1200         FORMAT(I5,4F10.0)

            ELSE IF(IOBUTYP(NBU).EQ.-2) THEN              !Data for all Us read?
              IFLAG=-1
              NBU=NBU-1                             ! Updating of counter undone

*_______________________Checks if same number of BUs are supplied in E2 and E3

              IF(IOCALBU(NU+1).EQ.4 .AND. NBU-NUMTBU.NE.NBUW(NU+1)) 
     ;           CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'ERROR:'
     ;            //' INCONSISTENCY BETWEEN NO. OF BASIC UNITS'
     ;            //' FROM CARDS E2.3 AND E3.1',NROW,0,IUOBS,1,7.02)

            ENDIF                                             ! Is device ended?
          ENDDO                                                ! End of group E3
          IFLAG=0                           ! Return flag value to initial value
        ENDIF                                            ! Device is not a point

*_______________________Updates some counters for Us and BUs

        NUMTU=NU+1                         ! Update of counter of no. of Us read
        NUMTBU=NBU                    ! Update of counter of total number of BUs
        NOBUF(NUMTU+1)=NUMTBU+1       ! First BU defining first U of next DEVICE

*_______________________Check card E3.1

        IOBUTYP(NUMTU)=IOBUTYP(NBU)          ! Assign BU type to U instead of BU

        CALL CHECK_CARD31
     ;(IERROR     ,INPWR    ,IODEVICE(ND,3) ,IOWAR    ,IUOBS    ,MAINF
     ;,NOBUF(NUF) ,NROW     ,NUF            ,NUMTBU   ,NUMTNOD  ,NUMTU
     ;,BUDAT      ,IOBUTYP  ,IOCALBU        ,FILENAME) 

*_______________________Determine nodal weights

        IODEVICE(ND,7)=NUMTNODC+1                   ! First node defining device

        CALL SPAT_WEIGHT_OBS
     ;(IODEVICE(ND,3) ,LMXNDL   ,NUF      ,NUMEL    ,NUMNP    ,NUMTBU
     ;,NUMTNOD        ,NUMTNODC ,NUMTU    ,AREA     ,BUDAT    ,EXTNBU
     ;,INDEXNOD       ,IOBUTYP  ,IOCALBU  ,IOUTYP   ,KXX      ,LNNDEL
     ;,LTYPE          ,MAINF    ,NOBUF    ,WTOBSBU  ,WTOBSN   ,WTOBSU
     ;,X              ,Y        ,Z)       

*_______________________GROUP E4. Information on the observation values and 
*_______________________temporal extension

        NO=NUMTOBSC                        ! Last observation of previous device
        IODEVICE(ND,8)=NO+1               ! First observation for current device

*_______________________If currrent device is considered for calculations,
*_______________________updates counter of number of integration times, 
*_______________________replacing IODEVINC value

        IF(IODEVICE(ND,2).NE.0) IODEVICE(ND,2)=NUMTITC+1        ! 1. IT for dev.

C------------------------- First data location for current device

        IFIRST_DATA=NUMTITC+1

        DO WHILE(IFLAG.GE.0)               ! Starts loop for device observations

          NO=NO+1                         ! Update of counter of obs. for device

*_______________________Reads card E4.1

          LEAUX=LEEL(FILENAME,IUOBS,MAINF,NROW,INPWR)
          READ(LEAUX,1300,ERR=9500) NOOBS,VAL,COV,TOBS1
     ;                             ,TOBS2,IOT
 1300     FORMAT(I5,4F10.0,I5)

*_______________________Has card ended (NOOBS pos. or not)?

          IF(NOOBS.LE.0) THEN
            IFLAG=-1            !Card ended
            NO=NO-1             !No data in line. Counter update undone
          ELSE                  !Card not ended, write and go on reading

*_______________________Once TRANSIN knows the goodness of the row, read values are
*_______________________assigned

            VOBS(NO)=VAL
            TOBS(NO,1)=TOBS1
            TOBS(NO,2)=TOBS2
            IOTINT(NO)=IOT
            IF (IODEVICE(ND,2).NE.0) THEN
               MEASTYP(NO)=IODEVICE(ND,1)
            ELSE
               MEASTYP(NO)=10
            ENDIF

*_______________________Checks if obs. st. dev. is zero. In this case, automatic
*_______________________assignment of device default standard deviation.

            IF (IOINV.GT.0) THEN
              IF(COV.EQ.0D0) THEN
                COVINV(NO)=STDEV
              ELSE
                COVINV(NO)=COV
              ENDIF
            END IF
*_______________________ Identifies observation times

            IF (IODEVICE(ND,4).EQ.0) THEN            ! Point in time observation
              NUMTITC=NUMTITC+1                            ! Total number of ITs
              TIT(NUMTITC)=TOBS(NO,1)                            ! Integral time
              NOOBSIT(NUMTITC)=NO                     ! Rel. between obs. and IT

            ELSE IF (IODEVICE(ND,4).EQ.1) THEN      ! Obs. over simulation times

              CONTINUE

            ELSE IF (IODEVICE(ND,4).EQ.2) THEN      ! Obs. over defined interval

              DO NOB=1,IODEVICE(ND,5)         ! Loop over NUMINT integrat. times
                NUMTITC=NUMTITC+1                             ! Total no. of ITs
                TIT(NUMTITC)=TOBS(NO,1)+(NOB-1)                      ! Obs. time
     ;                       *(TOBS(NO,2)-TOBS(NO,1))/(IODEVICE(ND,5)-1)
                NOOBSIT(NUMTITC)=NO                   ! Rel. between obs. and IT
              ENDDO

            ENDIF                                     ! Type of time integration
          ENDIF                                                ! Has card ended?
        ENDDO                          ! End of loop for observations for current device

*_______________________Check card E4.1. Checks all observations of current dev.

        CALL CHECK_CARD41
     ;(IDIMCOV  ,IERROR   ,IODEVICE(ND,4) ,IOINV    ,IOWAR    ,IUOBS    
     ;,MAINF    ,NO       ,NROW           ,NUMTOBS  ,NUMTOBSC ,TMAX     
     ;,COVINV   ,IOTINT   ,TOBS           ,FILENAME)    

*_______________________Organise ITs for device

        IF (IODEVICE(ND,4).EQ.0 .OR. IODEVICE(ND,4).EQ.2)
     ;  CALL ORDER
     ;(IFIRST_DATA ,NUMTIT   ,NUMTITC  ,NOOBSIT  ,TIT)     

*_______________________Determine temporal weights

        CALL TEMP_WEIGHT_OBS
     ;(IODEVICE(ND,4) ,NO      ,IODEVICE(ND,5) ,NUMTOBS  ,NUMTOBSC 
     ;,IOTINT         ,TOBS    ,WTOBST)  

*_______________________Writing information for device

*** DE MOMENTO LA DEJO, A VER QUE TAL SALE. TRAS LA VERIFICACION SERA ELIMINADA

        IF (INPWR.NE.0) CALL WRIT_DEVICE
     ;(NDEVS,ND,IODEVICE,BUDAT,NUMTNOD,WTOBSU,WTOBSBU,MAINF,NUMTU,NUF,
     ;NUMTBU,NOBUF,IOBUTYP,NO,NUMTOBSC,VOBS,TOBS,COVINV,IDIMCOV,
     ;NUMTOBS,IOTINT,WTOBSN,INDEXNOD,NUMTNODC,IOUTYP,X,Y,Z,NUMNP)

        IFLAG=0                             ! Return flag value to initial value
        NUMTOBSC=NO                    ! Update of counter of total no. of. obs.

      ENDDO                                            ! End of loop for devices

*_______________________Check that NUMTOBS,NUMTNOD and NUMTIT have
*_______________________been defined sufficiently large by comparing them to
*_______________________the values of NUMTOBSC, NUMTNODC and NUMTITC

      CALL CHECK_ALLOC
     ;(IERROR   ,IOWAR    ,IUOBS    ,MAINF    ,NROW     ,NUMTIT   
     ;,NUMTITC  ,NUMTNOD  ,NUMTNODC ,NUMTOBS  ,NUMTOBSC ,FILENAME)

*_______________________First node/obs./IT for "imaginary device"

      IODEVICE(NDEVS+1,2)=NUMTIT+1                          ! Number of first IT
      IODEVICE(NDEVS+1,7)=NUMTNOD+1                       ! Number of first node
      IODEVICE(NDEVS+1,8)=NUMTOBS+1                ! Number of first observation

*_______________________Assembly of the inverse covariance matrix

      IF (IOINV.GT.0) CALL ASS_COVINV
     ;(CORRCOEF ,IDIMCOV  ,NDEVS  ,STDEV    ,COVINV  ,IODEVICE)    

C----------------------- Stores in IODEVICE(*,10) the first "location" in the 
C----------------------- data (because IODEVICE(*,2) is modified during exec.)

       DO ND=1,NDEVS
          IODEVICE(ND,10)=IODEVICE(ND,2)
       ENDDO

*      CALL WRIT_CHR1
*     ;(COVINV,IDIMCOV,VOBS,TOBS,NUMTOBS,IODEVICE,NDEVS)

      RETURN

 9000 CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;           'GENERIC FORTRAN ERROR READING CARD E1.1'
     ;           ,NROW,0,IUOBS,1,7.03)
      RETURN
 9100 CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;           'GENERIC FORTRAN ERROR READING CARD E2.1'
     ;           ,NROW,0,IUOBS,1,7.04)
      RETURN
 9200 CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;           'GENERIC FORTRAN ERROR READING CARD E2.2'
     ;           ,NROW,0,IUOBS,1,7.05)
      RETURN
 9400 CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;           'GENERIC FORTRAN ERROR READING CARD E3.1'
     ;           ,NROW,0,IUOBS,1,7.06)
      RETURN
 9500 CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;           'GENERIC FORTRAN ERROR READING CARD E4.1'
     ;           ,NROW,0,IUOBS,1,7.07)

      RETURN
      END
