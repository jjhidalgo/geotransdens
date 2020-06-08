      SUBROUTINE READ_LOCATIONS
     ;(ICALL     ,IDIMSAM  ,IDIMVAR   ,IERROR   ,INPWR      ,IORD_PP   
     ;,IOWAR     ,IPOSITION,IUGEO     ,KTYPE    ,MAINF      ,NLECT      
     ;,NMEAS     ,NVAR     ,POSMEAS   ,TMAX     ,VMEAS      ,VSTATS     
     ;,FILENAME  ,MXNVAR_GS)    

********************************************************************************
*
* PURPOSE Read card G7 or card G8, depending on ICALL (ICALL=1 for reading 
*         sampling locations and variables values; ICALL=2 for reading pilot 
*         points location (only if IORD_PP=0). Pilot points values are assigned
*         coherently with GSLIB rountines (see Step 6)
*
* DESCRIPTION Subroutine can be summarized in several steps
*
*             - Step 0: Declaration of variables
*             - Step 1: Initialisation of statistics
*             - Step 2: Main loop over pilot points/sampling locations
*                 - Step 2.1: Two different cases.
*                             a) ICALL=1; Sampling locations. 
*                                Position and values of the different variables 
*                                are read
*                             b) ICALL=2; Pilot points. If positions are user 
*                                values, reads current card (only if IORD_PP=0)
*                                In this case, only position is read and variable 
*                                values will be coherently assigned (see Step 6).
*                 - Step 2.2: Checks correct numeration of current point
*                 - Step 2.3  Checks the number of measurements at current point
*                 - Step 2.4: For the kriging case, extensive variables do
*                             not need to be considered when prim. var. is
*                             trimmed, so that, proceeds to next location.
*                 - Step 2.5: Everything is OK. Assigns location if point is
*                             a sampling location or if is a pilot point and
*                             coordinates are user input data.
*                             Same considerations for primary variable. Also
*                             statistics are updated.
*                 - Step 2.6: Extensive variables treatment
*             - END LOOP
*             - Step 3: Calculates and writes final statistics
*             - Step 4: Stores final statistics
*             - Step 5: Checks if a variable has been trimmed at all sampling
*                       locations
*             - Step 6: Assigns variables values at pilot points such that pilot
*                       point locations are considered for the first kriging 
*                       (kriging blocks with pilot point AND sampling locations) 
*                       and such that are not considered in the second kriging 
*                       (kriging pilot points with sampling locations). Extensive
*                       variables values are chosen such that they are always 
*                       trimmed.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  POSMEAS                Array containing coordinates of pilot points and 
*                         sampling locations
*  TMAX                   Array containing maximum trimming limit for each 
*                         variable
*  VMEAS                  Array containing measurements of primary/secondary 
*                         variables and external drifts for pilot points 
*                         and sampling locations. External drifts are not read
*                         but assigned elsewhere. Variables values at pilot 
*                         points locs. are assigned coherently with GSLIB rout.
*  VSTATS                 Array containing statistics of all variables
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*
* INTERNAL VARIABLES: ARRAYS
*
*  AVERAGE                Temporarily stores averages of all variables
*  INOT                   Auxiliar array to check number of measurements
*  NSAMPLE                Temporarily stores number of sampling locations of 
*                         each variable
*  NTRIMMED               Temporarily stores number of trimmed values of each 
*                         variable                                         
*  STDEV                  Temporarily stores standard deviations of each 
*                         variable
*
* EXTERNAL VARIABLES: SCALARS
*
*  ICALL                  If 1, sampling locations definition. If 2, pilot 
*                         points definition and assignation of meas.
*  IDIMSAM                Dimension of main variables VMEAS and POSMEAS
*  IDIMVAR                Dimension of array VMEAS
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IORD_PP                If 0, reads pilot points locations      
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IPOSITION              Initial position at arrays VMEAS and POSMEAS
*                         (0 if ICALL=2, NPP_GR (number of pilot points defining 
*                         actual group if ICALL=1)
*  IUGEO                  GEO file unit number       
*  KTYPE                  Kriging type  0: Simple kriging
*                                       1: Residual kriging
*                                       2: Kriging with locally varying mean
*                                       3: Kriging with external drift (four)
*                                       4: Simple cokriging
*                                       5: Standardized ordinary cokriging
*                                       6: Traditional ordinary cokriging
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NLECT                  Number of records to be read/assigned (NMEAS if 
*                         ICALL=1; NPIPO if ICALL=2)
*  NMEAS                  Number of sampling locations
*  NVAR                   Number of variables (primary+all secondary)
*
* INTERNAL VARIABLES: SCALARS
*
*  ILECT                  Dummy counter for records
*  ISUM                   Dummy counter for checking
*  IVAR                   Dummy counter of variables                            
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*  NLOC                   Sequentiall number for reading locations
*  NROW                   Current record number                                 
*  PRIM                   Primary variable value              
*  SEC1                   First secondary variable value
*  SEC2                   Second secondary variable value
*  SEC3                   Third secondary variable value                       
*  VARMIN                 Minimum value of primary variable
*  XAUX                   Auxiliar coordinates of location                     
*  YAUX                   Auxiliar coordinates of location                      
*  ZAUX                   Auxiliar coordinates of location                      
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
* HISTORY: AAR First coding (Dec-2001)
*          AAR Revision     (Feb-2002)
*          AAR Inclusion of group of zones (July-2003)
*
********************************************************************************

C______________________________________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 INPWR,ICALL,MAINF,IORD_PP,NLECT,IUGEO,IERROR,IOWAR
     ;         ,NVAR,KTYPE,IDIMSAM,IPOSITION,NMEAS,IDIMVAR,MXNVAR_GS
                                                                 ! Real external
      REAL*8 VARMIN,TMAX(NVAR),POSMEAS(IDIMSAM,3)
     ;      ,VMEAS(IDIMSAM,IDIMVAR),VSTATS(MXNVAR_GS,4)
                                                              ! Integer internal
      INTEGER*4 IVAR,NSAMPLE(4),NTRIMMED(4),ILECT,NROW,NLOC,INOT(4),ISUM
                                                                 ! Real internal
      REAL*8 AVERAGE(4),STDEV(4),XAUX,YAUX,ZAUX,PRIM,SEC1,SEC2,SEC3
                                                                    ! Characters
      CHARACTER FILENAME(18)*20,LEAUX*100,LEEL*100

C__________________________________________ Step 1: Initialisation of statistics

      IF (ICALL.EQ.1) THEN                                  ! Sampling locatoins

        DO IVAR=1,4                                   ! IVAR=1: Primary variable
                                                ! IVAR=2,3,4 Secondary variables
          NSAMPLE(IVAR)=0                         ! Number of sampling locations
          NTRIMMED(IVAR)=0                            ! Number of trimmed values
          AVERAGE(IVAR)=0D0                                            ! Average
          STDEV(IVAR)=0D0                                   ! Standard deviation

        END DO
        VARMIN=TMAX(1)*1D5

      ELSE                                                        ! Pilot points
        VARMIN=VSTATS(1,4)
      END IF

C________________________ Step 2: Main loop over pilot points/sampling locations

      DO ILECT=1,NLECT

C__________ Step 2.1: Two different cases. 
C__________           a) Pilot points. If positions are user input values, reads
C__________              current card. In this case, only position is read and 
C__________              variable values are arbitrarily assigned.
C__________           b) Sampling locations. Current card is always read. 
C__________              Position and values of the different variables are read

        IF (ICALL.EQ.1.OR.(ICALL.EQ.2.AND.IORD_PP.EQ.0)) THEN

          LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
          READ (LEAUX,1000,ERR=9100) NLOC,XAUX,YAUX,ZAUX,PRIM,SEC1
     ;                            ,SEC2,SEC3
 1000     FORMAT(I5,8F10.0)

C_______________________ Step 2.2: Checks correct numeration of current location

          IF (NLOC.NE.ILECT) THEN                   ! Incorrect numeration. Error

            IF (ICALL.EQ.1) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ; 'INCORRECT NUMERATION: SAMPLING LOCATION. CARD G7 ',NROW,1
     ; ,IUGEO,2,10.4)

            IF (ICALL.EQ.2) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ; 'INCORRECT NUMERATION: PILOT POINTS. CARD G8 ',NROW,1
     ; ,IUGEO,2,10.4)

          END IF  ! NLOC.NE.ILECT

        END IF


        IF (ICALL.EQ.1) THEN                                ! Sampling locations

C___________________ Step 2.3 Checks the number of measurements at current point

          IF (PRIM.GE.TMAX(1)) INOT(1)=1
          IF (NVAR.GE.2.AND.SEC1.GE.TMAX(2)) INOT(2)=1
          IF (NVAR.GE.3.AND.SEC2.GE.TMAX(3)) INOT(3)=1
          IF (NVAR.GE.4.AND.SEC3.GE.TMAX(4)) INOT(4)=1

          ISUM=0
          DO IVAR=1,NVAR
            ISUM=ISUM+INOT(IVAR)
          END DO
          IF (ISUM.EQ.NVAR) THEN                         ! Location without data
            CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ; 'POINT WITHOUT MEASUREMENTS AT CARD G7 ',NROW,1
     ; ,IUGEO,0,10.4)
             GOTO 666                                            ! Next location
          END IF

C________________________ Step 2.4: For the kriging case, secondary variables do
C________________________           not need to be considered when prim. var. is
C________________________           trimmed, so that, proceeds to next location.

          IF (KTYPE.GE.3.AND.PRIM.GE.TMAX(1)) THEN
             NTRIMMED(1)=NTRIMMED(1)+1
             GOTO 666                                            ! Next location
          END IF

        END IF                            ! Only for sampling locations. Card G7

C______________________ Step 2.5: Everything is OK. Assigns location if point is
C______________________           a sampling location or if is a pilot point and
C______________________           coordinates are user input data. 
C______________________           Same considerations for primary variable. Also
C______________________           statistics are updated.

        IF (ICALL.EQ.1.OR.(ICALL.EQ.2.AND.IORD_PP.EQ.0)) THEN
          POSMEAS(IPOSITION+NLOC,1) = XAUX
          POSMEAS(IPOSITION+NLOC,2) = YAUX
          POSMEAS(IPOSITION+NLOC,3) = ZAUX
        END IF

        IF (ICALL.EQ.1) THEN                       ! Only for sampling locations
          VMEAS(IPOSITION+NLOC,1)=PRIM         ! Assigns value of primary variable

          IF (PRIM.LT.TMAX(1)) THEN               ! Primary variable not trimmed

            NSAMPLE(1)=NSAMPLE(1)+1
            AVERAGE(1)=AVERAGE(1)+PRIM
            STDEV(1)=STDEV(1)+PRIM*PRIM
            IF (PRIM.LE.VARMIN) VARMIN=PRIM   ! Ass. minimum value of prim. var.

          ELSE                                              ! Prim. var. trimmed

            NTRIMMED(1)=NTRIMMED(1)+1
          END IF                                           ! Trimmed/not trimmed

C___________ Step 2.6: Secondary variables treatment. Same way as for prim. var.

          IF (NVAR.GE.2) THEN
            VMEAS(IPOSITION+NLOC,2)=SEC1
            IF (SEC1.GE.TMAX(2)) THEN                                  ! Trimmed
              NTRIMMED(2)=NTRIMMED(2)+1 
            ELSE                                             ! Update statistics
              NSAMPLE(2)=NSAMPLE(2)+1
              AVERAGE(2)=AVERAGE(2)+SEC1
              STDEV(2)=STDEV(2)+SEC1*SEC1
            END IF                                         ! Trimmed/non trimmed

            IF (NVAR.GE.3) THEN
              VMEAS(IPOSITION+NLOC,3)=SEC2
              IF (SEC2.GE.TMAX(3)) THEN                                ! Trimmed
                NTRIMMED(3)=NTRIMMED(3)+1
              ELSE                                           ! Update statistics
                NSAMPLE(3)=NSAMPLE(3)+1
                AVERAGE(3)=AVERAGE(3)+SEC2
                STDEV(3)=STDEV(3)+SEC2*SEC2
              END IF

              IF (NVAR.GE.4) THEN
                VMEAS(IPOSITION+NLOC,4)=SEC3
                IF (SEC3.GE.TMAX(4)) THEN                              ! Trimmed
                  NTRIMMED(4)=NTRIMMED(4)+1
                ELSE                                         ! Update statistics
                  NSAMPLE(4)=NSAMPLE(4)+1
                  AVERAGE(4)=AVERAGE(4)+SEC3
                  STDEV(4)=STDEV(4)+SEC3*SEC3
                END IF

              END IF                                                 ! NVAR.GE.4
            END IF                                                   ! NVAR.GE.3
          END IF                                                     ! NVAR.GE.2
        END IF                                    ! ILECT is a sampling location
 666    CONTINUE
      END DO                                                     ! Next location

C_______________________________ Step 3: Calcultates and writes final statistics

      IF (ICALL.EQ.1) THEN

        DO IVAR=1,NVAR
          AVERAGE(IVAR)=AVERAGE(IVAR)/DMAX1(1D0,DFLOAT(NSAMPLE(IVAR)))
          STDEV(IVAR)=STDEV(IVAR)/DMAX1(1D0,DFLOAT(NSAMPLE(IVAR)))
     ;               -AVERAGE(IVAR)*AVERAGE(IVAR)
          STDEV(IVAR)=DSQRT(STDEV(IVAR))
        END DO

C_______________________________________________ Step 4: Stores final statistics

        DO IVAR=1,NVAR
          VSTATS(IVAR,1)=NSAMPLE(IVAR)
          IF (KTYPE.NE.0 .AND. IVAR.EQ.1) VSTATS(IVAR,2)=AVERAGE(IVAR)
          IF (KTYPE.EQ.0 .AND. IVAR.EQ.1 .AND.
     ;        VSTATS(1,2).NE.AVERAGE(1) .AND. IOWAR.NE.0)
     ;        WRITE(MAINF,2000) VSTATS(1,2),AVERAGE(1)
 2000         FORMAT(//,' WARNING: SUPPLIED SIMPLE KRIGING MEAN NOT'
     ;                  ' COHERENT WITH CALCULATED MEAN ON THE BASIS',/,
     ;                  ' OF MEASUREMENTS.',//,
     ;                  ' CALCULATED: ',E10.4,/,' SUPPLIED  : ',E10.4)
     
          VSTATS(IVAR,3)=STDEV(IVAR)
          IF (IVAR.EQ.1) VSTATS(IVAR,4)=VARMIN
        END DO

C_______ Step 5: Checks if a variable has been trimmed at all sampling locations

        IF (IOWAR.NE.0) THEN
          DO IVAR=1,NVAR
            IF (NTRIMMED(IVAR).EQ.NMEAS.AND.IOWAR.NE.0) 
     ;         WRITE(MAINF,2008) IVAR
          END DO
 2008     FORMAT(//,' WARNING: VARIABLE ',I2,' HAS BEEN TRIMMED AT ALL'
     ;            ' SAMPLING LOCATIONS')
        END IF
            
      ENDIF ! ICALL.EQ.1

C________ Step 6: Assigns variables values at pilot points such that pilot point
C________         locations are considered for the first kriging (kriging blocks
C________         with pilot point AND sampling locations) and such that are not
C________         considered in the second kriging (kriging pilot points with 
C________         sampling locations). Secondary variables values are chosen
C________         such that they are always trimmed.

      IF (ICALL.EQ.2) THEN
        DO ILECT=1,NLECT
          VMEAS(IPOSITION+ILECT,1)=-1D0*DFLOAT(ILECT)*1D60
          DO IVAR=2,NVAR
            VMEAS(IPOSITION+ILECT,IVAR)=1D3*TMAX(IVAR)
          END DO
        END DO
      END IF

      RETURN
C______________________________________________________________ Fatal error call

 9100  CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING CARD G7/G8 ',NROW,1,IUGEO,2,9.8)

      RETURN
      END
      
