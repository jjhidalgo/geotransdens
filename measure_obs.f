      SUBROUTINE MEASURE_OBS
     ;(IOCALTYP ,NUF      ,NUMEL    ,NUMNP    ,NUMTBU   ,NUMTNOD
     ;,NUMTU    ,AREA     ,BUDAT    ,EXTNBU   ,IOBUTYP  ,IOCALBU
     ;,IOUTYP   ,NOBUF    ,X        ,Y        ,Z)       

********************************************************************************
*
* PURPOSE
*
* Determines measures belonging all nodes which are used to describe
* the device locations.
*
* DESCRIPTION
*
* For spatial integration or averaging, it is necessary to determine
* the measures (length of a line, area of a surface and volume of a
* volume) of the part of the spatial extension of the device which
* "belongs" to each of the nodes which are used to describe device
* locations.
* By 03-2001, only the case IOUTYP=2 (Unit: Line) is implemented.
*
* Throughout the subroutine, U and BU are used in comments as
* abbreviations for unit and basic unit, respectively.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,
*                         volume for 3-D)
*  BUDAT                  Data from group E3. One column for one basic unit
*  EXTNBU                 Measure (length, area or volume) of basic unit
*  IOBUTYP                Basic unit type
*  IOCALBU                Calculation method for basic unit
*  IOUTYP                 Unit type
*  NOBUF                  Number of first basic unit
*  X                      X-coord for a given node
*  Y                      Y-coord for a given node
*  Z                      Z-coord for a given node
*
* EXTERNAL VARIABLES: SCALARS
*
*  IOCALTYP               Method of device spatial integration 
*  NUF                    First unit defining current device
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NUMTBU                 Total number of basic units                    
*  NUMTNOD                Total number of nodes used for calculating obs.       
*  NUMTU                  Total number of units
*
* INTERNAL VARIABLES: SCALARS
*
*  EXTN                   Summation dummy variable
*  EXTNU                  Summation dummy variable
*  NBU                    Dummy counter of basic units
*  NEL                    Element number (identification dummy variable)
*  NNOD1                  First node of NEL (identification dummy variable)
*  NNOD2                  Second node of NEL (identification dummy variable)
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ZERO_ARRAY                                                                   
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR     03-2001     Revision
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION EXTNBU(NUMTNOD),NOBUF(NUMTNOD+1),BUDAT(4,NUMTNOD),
     ;IOUTYP(NUMTNOD),AREA(NUMEL),IOCALBU(NUMTNOD),IOBUTYP(NUMTNOD),
     ;X(NUMNP),Y(NUMNP),Z(NUMNP)

      EXTNDEV=0                                     ! Init. of measure of device

      CALL ZERO_ARRAY (EXTNBU,NUMTBU)

      DO NU=NUF,NUMTU                          ! Loop for Us belonging to device

        EXTNU=0                               ! Init. of measure of current unit

        IF (IOUTYP(NU).EQ.1) THEN                             ! Units are points

          CONTINUE                                       ! No sense in this case

        ELSE IF(IOUTYP(NU).EQ.2) THEN                           ! Unit is a line

*_______________________ BUs: Points defining Unit (LINE)

          IF(IOBUTYP(NU).EQ.1) THEN

            DO NBU=NOBUF(NU),NOBUF(NU+1)-2 ! Loop over BU's (points of the line)

              EXTN=SQRT((BUDAT(1,NBU)-BUDAT(1,NBU+1))**2+             ! Distance
     ;                  (BUDAT(2,NBU)-BUDAT(2,NBU+1))**2+
     ;                  (BUDAT(3,NBU)-BUDAT(3,NBU+1))**2)
              EXTNBU(NBU)=EXTNBU(NBU)+EXTN/2D0              ! Updates BU measure
              EXTNBU(NBU+1)=EXTNBU(NBU+1)+EXTN/2D0          ! Updates BU measure
              EXTNU=EXTNU+EXTN                               ! Updates U measure

            ENDDO

          ENDIF

*_______________________BUs: Nodes defining Unit (LINE)

          IF(IOBUTYP(NU).EQ.2) THEN

            DO NBU=NOBUF(NU),NOBUF(NU+1)-2  ! Loop over BU's (nodes of the line)

              NNOD1=INT(BUDAT(1,NBU)+0.5)                           ! First node
              NNOD2=INT(BUDAT(1,NBU+1)+0.5)                        ! Second node
              EXTN=SQRT((X(NNOD1)-X(NNOD2))*(X(NNOD1)-X(NNOD2))+      ! Distance
     ;                  (Y(NNOD1)-Y(NNOD2))*(Y(NNOD1)-Y(NNOD2))+
     ;                  (Z(NNOD1)-Z(NNOD2))*(Z(NNOD1)-Z(NNOD2)))
              EXTNBU(NBU)=EXTNBU(NBU)+EXTN/2D0              ! Updates BU measure
              EXTNBU(NBU+1)=EXTNBU(NBU+1)+EXTN/2D0          ! Updates BU measure
              EXTNU=EXTNU+EXTN                               ! Updates U measure

            ENDDO
          ENDIF

*_______________________BUs: Elements (this option has only sense for 1D elem.
*_______________________defining a polyline)

          IF(IOBUTYP(NU).EQ.3) THEN

            DO NBU=NOBUF(NU),NOBUF(NU+1)-2  ! Loop over BU's (nodes of the line)

              NEL=INT(BUDAT(1,NBU)+0.5)          ! Identifies the element number
              EXTNBU(NBU)=AREA(NEL)                ! Direct assignment of LENGHT
              EXTNU=EXTNU+EXTNBU(NBU)                 ! Updates the measure of U

            ENDDO

          ENDIF

*_______________________BUs: Zone (only valid for point- and line-zones)

          IF(IOBUTYP(NU).EQ.4) THEN
            CONTINUE                                       ! Not implemented yet
          ENDIF

*_____________U: Surface

        ELSE IF(IOUTYP(NU).EQ.3) THEN

*_______________________BUs: Points

          IF (IOBUTYP(NU).EQ.1) THEN
            
            CONTINUE                                       ! Not implemented yet

*_______________________BUs: Nodes

            CONTINUE                                       ! Not implemented yet

*_______________________BUs: Elements (only valid for 2D elements)

          ELSE IF(IOBUTYP(NU).EQ.3) THEN

            CONTINUE                                       ! Not implemented yet

*_______________________BUs: Zone (only valid for point-, line- and area-zones)

          ELSE IF(IOBUTYP(NU).EQ.4) THEN

            CONTINUE                                       ! Not implemented yet

          ENDIF

*_____________U: Volume

        ELSE IF(IOUTYP(NU).EQ.4) THEN

*_______________________BUs: Points

          IF(IOBUTYP(NU).EQ.1) THEN

            CONTINUE                                       ! Not implemented yet

*_______________________BUs: Nodes

          ELSE IF(IOBUTYP(NU).EQ.2) THEN

            CONTINUE                                       ! Not implemented yet

*_______________________BUs: Elements (only valid for 2D elements)

          ELSE IF(IOBUTYP(NU).EQ.3) THEN

            CONTINUE                                       ! Not implemented yet

*_______________________BUs: Zones (only valid for point-, line- and area-zones)

          ELSE IF(IOBUTYP(NU).EQ.4) THEN

            CONTINUE                                       ! Not implemented yet

          ENDIF

        ENDIF

*_______________________Modification for spatial average over BUs
*_______________________Normalisation with respect to measure of U

        IF(IOCALBU(NU).EQ.3) THEN
          DO NBU=NOBUF(NU),NOBUF(NBU+1)-1
            EXTNBU(NBU)=EXTNBU(NBU)/EXTNU
          ENDDO
        ENDIF

        EXTNDEV=EXTNDEV+EXTNU                                ! Measure of device

      ENDDO

*_______________________Modification for spatial average over Us
*_______________________Normalisation with respect to measure of device

        IF(IOCALTYP.EQ.3) THEN
          DO NBU=NOBUF(NUF),NUMTBU
            EXTNBU(NBU)=EXTNBU(NBU)/EXTNDEV
          ENDDO
        ENDIF

      RETURN

      END
