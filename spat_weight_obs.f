      SUBROUTINE SPAT_WEIGHT_OBS
     ;(IOCALTYP ,LMXNDL   ,NUF      ,NUMEL    ,NUMNP    ,NUMTBU
     ;,NUMTNOD  ,NUMTNODC ,NUMTU    ,AREA     ,BUDAT    ,EXTNBU
     ;,INDEXNOD ,IOBUTYP  ,IOCALBU  ,IOUTYP   ,KXX      ,LNNDEL
     ;,LTYPE    ,MAINF    ,NOBUF    ,WTOBSBU  ,WTOBSN   ,WTOBSU
     ;,X        ,Y        ,Z)       

********************************************************************************
*
* PURPOSE
*
* Determines the spatial weights associated with all nodes used to
* define current device. 
*
* DESCRIPTION
*
* The weights are used to calculate the model-output corresponding to
* observations as a weighted sum of nodal values. This is done in
* COMP_OBS. For each device, SPATIAL_WEIGHT is called once after card
* E3.1 has been read.
* Throughout the subroutine, U and BU are used in comments as
* abbreviations for unit and basic unit, respectively.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  BUDAT                  Data from group E3. One column for one basic unit     
*  EXTNBU                 Measure (length, area or volume) of basic unit        
*  INDEXNOD               Index relating nodes                                  
*  IOBUTYP                Basic unit type                                       
*  IOCALBU                Calculation method for basic unit                     
*  IOUTYP                 Unit type                                             
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  NOBUF                  Number of first basic unit                            
*  WTOBSBU                Weight for basic unit                                 
*  WTOBSN                 Weight for node used to calculate observation         
*  WTOBSU                 Weight for unit                                       
*  X                      X-coord for a given node                              
*  Y                      Y-coord for a given node                              
*  Z                      Z-coord for a given node                              
*
* INTERNAL VARIABLES: ARRAYS
*
*  BF                     Value of basis function
*
* EXTERNAL VARIABLES: SCALARS
*
*  IOCALTYP               Method of spat. integr. (IODEVICE(ND,3))
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)
*  NUF                    Number of first unit for device
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NUMTBU                 Total number of basic units defining unit
*  NUMTNOD                Total number of nodes used for calculating obs.       
*  NUMTNODC               Counter variable. To compare with NUMTNOD
*  NUMTU                  Total number of units
*
* INTERNAL VARIABLES: SCALARS
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  BASISFUNC_OBS                                                                
*  MEASURE_OBS                                                                  
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR     03-20001    Revision
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION WTOBSN(NUMTNOD),IOCALBU(NUMTNOD),LNNDEL(NUMEL),
     ;INDEXNOD(NUMTNOD),IOUTYP(NUMTNOD),EXTNBU(NUMTNOD),WTOBSU(NUMTNOD),
     ;NOBUF(NUMTNOD+1),BUDAT(4,NUMTNOD),WTOBSBU(NUMTNOD),BF(6),
     ;IOBUTYP(NUMTNOD),KXX(LMXNDL,NUMEL),x(numnp),y(numnp)

*_______________________Device is point-defined (the usual case). TRANSIN old 
*_______________________methodology is applied

      IF (IOCALTYP.EQ.1) THEN

        XOB=BUDAT(1,NUMTBU)
        YOB=BUDAT(2,NUMTBU)
        ZOB=BUDAT(3,NUMTBU)

        CALL BASISFUNC_OBS          ! Determines basis function values for point
     ;(LMXNDL   ,MAINF    ,NEL      ,NUMEL    ,NUMNP    ,XOB
     ;,YOB      ,ZOB      ,AREA     ,BF       ,KXX      ,LTYPE
     ;,X        ,Y        ,Z)       

*_______________________ Element to which obs. point belongs to was not found

        IF (NEL.LE.0 .OR. NEL.GT.NUMEL) THEN
          WRITE(MAINF,2000) NUMTBU
          WRITE(6,2000) NUMTBU
 2000     FORMAT(//,' OBSERVATION POINT ',I5,' DOES NOT BELONG TO'
     ;              ' MODEL DOMAIN. CRITICAL STOP')
          STOP
        END IF

        NNUD=LNNDEL(NEL)     ! Num. of nodes for elem. at which point belongs to

        DO I=1,NNUD
          NUMTNODC=NUMTNODC+1                       ! Updates total no. of nodes
          WTOBSN(NUMTNODC)=BF(I)                    ! Assignment of nodal weight
          INDEXNOD(NUMTNODC)=KXX(I,NEL)                  ! Saving nodal relation
        ENDDO

*_______________________Spatial integration or spatial averaging over Us

      ELSE IF(IOCALTYP.EQ.2 .OR. IOCALTYP.EQ.3) THEN

*_______________________Determines measures belonging to BUs

        CALL MEASURE_OBS
     ;(IOCALTYP ,NUF      ,NUMEL    ,NUMNP    ,NUMTBU   ,NUMTNOD
     ;,NUMTU    ,AREA     ,BUDAT    ,EXTNBU   ,IOBUTYP  ,IOCALBU
     ;,IOUTYP   ,NOBUF    ,X        ,Y        ,Z)


*_______________________Loop for Us defining current device

        DO NU=NUF,NUMTU

*_____________U: Points

          IF (IOUTYP(NU).EQ.1) THEN

            CONTINUE                                     ! No sense in this case
            
*_____________U: Line. It can be defined by points, nodes or 1d-elements

          ELSE IF(IOUTYP(NU).EQ.2) THEN

*_______________________BUs: Points and def. of element

            IF(IOBUTYP(NU).EQ.0) THEN

C              XOB=BUDAT(1,NBU)
C              YOB=BUDAT(2,NBU)
C              ZOB=BUDAT(3,NBU)
C              NEL=INT(BUDAT(4,NBU)+0.5)

C             HOW IS THIS OPTION TO BE UNDERSTOOD?
C             Verify that the point belongs to the element
C             Write warning if this is not the case.

*_______________________BUs: Point

            ELSE IF(IOBUTYP(NU).EQ.1) THEN

*______________________  Loop over basic units of curr. unit

              DO NBU=NOBUF(NU),NOBUF(NU+1)-1   

                XOB=BUDAT(1,NBU)
                YOB=BUDAT(2,NBU)
                ZOB=BUDAT(3,NBU)

*______________________ Computes nodal weights

                CALL BASISFUNC_OBS
     ;(LMXNDL   ,MAINF    ,NEL      ,NUMEL    ,NUMNP    ,XOB
     ;,YOB      ,ZOB      ,AREA     ,BF       ,KXX      ,LTYPE
     ;,X        ,Y        ,Z)       

                NNUD=LNNDEL(NEL)  
                DO I=1,NNUD
                  NUMTNODC=NUMTNODC+1          ! Count. of total number of nodes
                  WTOBSN(NUMTNODC)=EXTNBU(NBU)*BF(I)      ! Ass. of nodal weight
                  INDEXNOD(NUMTNODC)=KXX(I,NEL)          ! Saving nodal relation
                ENDDO
              ENDDO

*_______________________BUs: Nodes

            ELSE IF(IOBUTYP(NU).EQ.2) THEN

              DO NBU=NOBUF(NU),NOBUF(NU+1)-1
                NUMTNODC=NUMTNODC+1            ! Count. of total number of nodes
                WTOBSN(NUMTNODC)=EXTNBU(NBU)              ! Ass. of nodal weight
                INDEXNOD(NUMTNODC)=INT(BUDAT(1,NBU)+0.5)  !Saving nodal relation
              ENDDO

*_______________________BUs: Elements (Only valid for 1D elements)
*_______________________(Not tested)

            ELSE IF(IOBUTYP(NU).EQ.3) THEN

              DO NBU=NOBUF(NU),NOBUF(NU+1)-1
                NEL=INT(BUDAT(1,NBU)+0.5)
                NUMTNODC=NUMTNODC+2            ! Count. of total number of nodes
                WTOBSN(NUMTNODC)=0.5*EXTNBU(NBU)          ! Ass. of nodal weight
                WTOBSN(NUMTNODC+1)=WTOBSN(NUMTNODC)       ! Ass. of nodal weight
                INDEXNOD(NUMTNODC)=KXX(1,NEL)            ! Saving nodal relation
                INDEXNOD(NUMTNODC+1)=KXX(2,NEL)          ! Saving nodal relation
              ENDDO

*_______________________BUs: Zones (Only valid for point- and line-zones)

            ELSE IF(IOBUTYP(NU).EQ.4) THEN

              CONTINUE                                     ! Not implemented yet

            ENDIF

*_____________U: Surface

          ELSE IF(IOUTYP(NU).EQ.3) THEN

              CONTINUE                                     ! Not implemented yet

*_____________U: Volumen

          ELSE IF(IOUTYP(NU).EQ.4) THEN

              CONTINUE                                     ! Not implemented yet

          ENDIF

        ENDDO

*_____________ Simple summation (of spatial average of each unit)

      ELSE IF(IOCALTYP.EQ.4) THEN

        CONTINUE                                           ! Not implemented yet

*_____________ Simple average (of spatial average of each unit)

      ELSE IF(IOCALTYP.EQ.5) THEN

        CONTINUE                                           ! Not implemented yet

*_____________User defined U weights (WTOBSU)

      ELSE IF(IOCALTYP.EQ.6) THEN

        DO NU=NUF,NUMTU

*_____________ Basic units have not spatial extension

          IF (IOCALBU(NU).EQ.1) THEN

            CONTINUE                                     ! No sense in this case

*_____________ Spatial average or integration over BU is requested

          ELSE IF(IOCALBU(NU).EQ.2 .OR. IOCALBU(NU).EQ.3) THEN

            CONTINUE                                       ! Not implemented yet

*_____________ Simple summation or simple average requested 

          ELSE IF(IOCALBU(NU).EQ.4 .OR. IOCALBU(NU).EQ.5) THEN

            CONTINUE                                       ! Not implemented yet

*_____________User defined BU weights (WTOBSBU)

          ELSE IF(IOCALBU(NU).EQ.6) THEN

*_______________________BUs: Points (and def. of element)

            IF (IOBUTYP(NU).EQ.0) THEN

              CONTINUE                                   ! No sense in this case

*_______________________BUs: Points

            ELSE IF(IOBUTYP(NU).EQ.1) THEN

              DO NBU=NOBUF(NU),NOBUF(NU+1)-1

                XOB=BUDAT(1,NBU)
                YOB=BUDAT(2,NBU)
                ZOB=BUDAT(3,NBU)

                CALL BASISFUNC_OBS
     ;(LMXNDL   ,MAINF    ,NEL      ,NUMEL    ,NUMNP    ,XOB
     ;,YOB      ,ZOB      ,AREA     ,BF       ,KXX      ,LTYPE
     ;,X        ,Y        ,Z)       

                NNUD=LNNDEL(NEL)   
                DO I=1,NNUD
                  NUMTNODC=NUMTNODC+1          ! Count. of total number of nodes
                  WTOBSN(NUMTNODC)=WTOBSU(NU)*WTOBSBU(NBU)*BF(I)   ! Ass. weight
                  INDEXNOD(NUMTNODC)=KXX(I,NEL)          ! Saving nodal relation
                ENDDO
              ENDDO

*_______________________BUs: Nodes

            ELSE IF(IOBUTYP(NU).EQ.2) THEN

              DO NBU=NOBUF(NU),NOBUF(NU+1)-1
                NUMTNODC=NUMTNODC+1            ! Count. of total number of nodes
                WTOBSN(NUMTNODC)=WTOBSU(NU)*WTOBSBU(NBU)           ! Ass. weight
                INDEXNOD(NUMTNODC)=INT(BUDAT(1,NBU)+0.5) ! Saving nodal relation
              ENDDO

*_______________________BUs: 1-D elements

            ELSE IF(IOBUTYP(NU).EQ.3) THEN

              DO NBU=NOBUF(NU),NOBUF(NU+1)-1
                NEL=INT(BUDAT(1,NBU)+0.5)                   ! Identifies element
                NUMTNODC=NUMTNODC+2            ! Count. of total number of nodes
                WTOBSN(NUMTNODC)=WTOBSU(NU)*WTOBSBU(NBU)           ! Ass. weight
                WTOBSN(NUMTNODC+1)=WTOBSN(NUMTNODC)                ! Ass. weight
                INDEXNOD(NUMTNODC)=KXX(1,NEL)            ! Saving nodal relation
                INDEXNOD(NUMTNODC+1)=KXX(2,NEL)          ! Saving nodal relation
              ENDDO

*_______________________BUs: Zones (Only valid for point-zones)

            ELSE IF(IOBUTYP(NU).EQ.4) THEN

              CONTINUE                                     ! Not implemented yet

            ENDIF                                 ! End of division (type of BU)

          ENDIF                                ! End of division (calc. over BU)

        ENDDO                                      ! Next unit of current device

      ENDIF                                ! End of division (calc. over device)

      RETURN

      END
