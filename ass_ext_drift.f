      SUBROUTINE ASS_EXT_DRIFT
     ;(IACTTYPE     ,ICALL        ,IDIMVAR_GS ,IERROR     ,IGROUP    
     ;,IOWAR        ,ISTART       ,IUGEO      ,LMXNDL     ,MAINF
     ;,MXMEASPP_GS  ,MXNZON_GS    ,MXNVAR_GS  ,NEXDR_GS
     ;,NPAREL       ,NPOINTS      ,NTYPAR     ,NUMEL      ,NUMNP
     ;,NZN_GR       ,NZPAR        ,AREA       ,COORD      ,EXDRZN_GS 
     ;,FILENAME     ,IGR_ZONE     ,INORPAR    ,KXX        ,LNNDEL
     ;,LTYPE        ,LXPAREL      ,POSMEAS_GS ,POSZN_GS   ,VMEAS_GS)
*     ,IZONMEAS_GS)

********************************************************************************
*
* PURPOSE Assigns the value of the external drifts (or locally varying mean), 
*         at a set of locations (pilot points or sampling locations).Remember 
*         that ext. drifts are defined zone to zone.
*
* DESCRIPTION Flow chart:
*
*  - Step 0: Declaration of variables
*
*  MAIN LOOP OVER LOCATIONS. 
*     - Step 1: Identifies to which group and zone point belongs to (NUMBER OF 
*               ZONE IS DEFINED ON THE BASIS ON GLOBAL NUMERATION, CONSIDERING 
*               ALL ZONES OF ACTUAL PARAMETER TYPE, AND NOT ONLY THE ZONES
*               BELONGING TO ACTUAL GROUP)
*     - Step 2: Calculates center of gravity of actual zone
*     - Step 3: Finds out number of zone in actual group (comparison of cog's)
*     - Step 4: Assigns external drifts
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  COORD                  Nodal coordinates                                     
*  EXDRZN_GS              Array containing zonal external drifts
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  IGR_ZONE               Array containing the relation zone->group (IVPAR(4))
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARZ, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IZONMEAS_GS            Array containing the measurement identifier contained in 
*                         a given zone of a given group
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  POSMEAS_GS             Array containing pilot points and sampling locations
*  POSZN_GS               Position of cog of zones belonging to actual group
*  VMEAS_GS               Array containing measurements and external drifts 
*
* INTERNAL VARIABLES: ARRAYS
*
*  BF                     Values of basis functions
*  COGELEM                Center of gravity of one element belonging to actual zone
*  COG_ZONE               Center of gravity of actual zone
*
* EXTERNAL VARIABLES: SCALARS
*
*  IACTTYPE               Type of parameters of actual group
*  ICALL                  1: meas. locations; 2: pilot points locations
*  IDIMVAR_GS             Used to dimension VMEAS_GS
*  IERROR                 Current number of errors on input data                
*  IGROUP                 Actual group
*  IOWAR                  Program echoes warnings if IOWAR.NE.0
*  ISTART                 Starting position in arrays POSMEAS_GS, VMEAS_GS
*                         ICALL=1 -> ISTART=NPP_GR; ICALL=2 -> ISTART=0
*  IUGEO                  GEO file unit number
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MXMEASPP_GS            Used to dimension POSMEAS_GS
*  MXNZON_GS              Number of zones of the most discretized group
*  MXNVAR_GS              Maximum number of variables. Last row filled
*                         in VMEAS_GS
*  NEXDR_GS               Number of external drifts in this group of zones
*  NFLAGS                 Used to dimension IFLAGS
*  NPAREL                 Number of element parameters in current problem       
*  NPOINTS                Number of points to be assigned
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZN_GR                 Number of zones defining actual group
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*
* INTERNAL VARIABLES: SCALARS
*
*  AREA_ZONE              Area of the zone to which point belongs to
*  DX2,DY2,DZ2            Dummies for distance calculations
*  I                      Dummy counter
*  IACTGROUP              Group to which point belongs to
*  IACTZONE               Zone to which point belongs to
*  IEXT                   Dummy counter of external drift terms
*  INUD                   Dummy counter of nodes
*  IPINORPAR              Pointer to array INORPAR
*  IPLXPAREL              Pointer to LXPAREL
*  IPOS                   Initial position at IGR_ZONE
*  IZON                   Dummy counter of zones
*  L                      Dummy counter of elements
*  LX                     Zone to which element belongs to
*  NEL                    Element to which point belongs to
*  NNUD                   Number of nodes of actual element
*  NODE                   Node number
*  NROW                   Current record number                                 
*  XNODE,YNODE,ZNODE      Node coordinates
*  XPOINT,YPOINT,ZPOINT   Point coordinates
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  BASISFUNC_OBS          Determines to which element a pilot points belongs to
*
* HISTORY:  AAR   First coding (Feb-2002)
*           AAR   Inclusion of groups of zones (July-2003)
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 LMXNDL,MAINF,NUMEL,NUMNP,IACTTYPE,IGROUP,NPOINTS
     ;         ,MXMEASPP_GS,IERROR,IOWAR,IUGEO,NPAREL,NTYPAR
     ;         ,NZPAR,NEXDR_GS,IDIMVAR_GS,MXNZON_GS,ISTART,ICALL
     ;         ,MXNVAR_GS,NZN_GR
     ;         ,KXX(LMXNDL,NUMEL),LTYPE(NUMEL),LXPAREL(NUMEL,NPAREL)
     ;         ,INORPAR(NTYPAR),IGR_ZONE(NZPAR)
     ;         ,LNNDEL(NUMEL)
*,IZONMEAS_GS(NZN_GR)
                                                                 ! Real external
      REAL*8 AREA(NUMEL),COORD(NUMNP,3),POSMEAS_GS(MXMEASPP_GS,3)
     ;      ,VMEAS_GS(MXMEASPP_GS,IDIMVAR_GS),EXDRZN_GS(MXNZON_GS,4)
     ;      ,POSZN_GS(MXNZON_GS,3)
                                                              ! Integer internal
      INTEGER*4 IPLXPAREL,I,IEXT,IACTZONE
     ;         ,IACTGROUP,L,LX,NNUD,INUD,NODE,IZON
                                                                 ! Real internal
      REAL*8 AREA_ZONE,XPOINT,YPOINT,ZPOINT,XNODE,YNODE,ZNODE,DX2,DY2
     ;      ,DZ2
     ;      ,COGELEM(3),COG_ZONE(3)
                                                                    ! Characters
      CHARACTER FILENAME(20)*20

C_______________________ MAIN LOOP OVER LOCATIONS

      DO I=1,NPOINTS

C_______________________ Step 1: Identifies to which group and zone the point
C_______________________         belongs to
         
         XPOINT=POSMEAS_GS(ISTART+I,1)
         YPOINT=POSMEAS_GS(ISTART+I,2)
         ZPOINT=POSMEAS_GS(ISTART+I,3)

         CALL GET_ZONE_GROUP
     ;(IACTGROUP  ,IACTTYPE   ,IACTZONE   ,IERROR   ,IGROUP   ,IOWAR
     ;,IPLXPAREL  ,ICALL      ,IUGEO      ,LMXNDL   ,MAINF    ,NPAREL   
     ;,NUMEL      ,NUMNP      ,NTYPAR     ,NZPAR    ,XPOINT   ,YPOINT   
     ;,ZPOINT     ,AREA       ,COORD      ,FILENAME ,IGR_ZONE ,INORPAR  
     ;,KXX        ,LTYPE      ,LXPAREL)

*         IF (ICALL.EQ.1) THEN    ! Measurements. Saves meas. location
*             IZONMEAS_GS(IACTZONE)=I
*         END IF

C_______________________ Step 2: Calculates center of gravity of actual zone

         AREA_ZONE=0.D0
         CALL ZERO_ARRAY(COG_ZONE,3)

         DO L=1,NUMEL

            LX=LXPAREL(L,IPLXPAREL)    
            NNUD=LNNDEL(L)

            IF (LX.EQ.IACTZONE) THEN   ! Element belongs to actual zone

C______________________ Step 2.1: Init. cog. of current element

              CALL ZERO_ARRAY(COGELEM,3)

C______________________ Step 2.2: Loop over element nodes

              DO INUD=1,NNUD

C______________________ Step 2.2.1: Identifies node and coordinates

                 NODE=KXX(INUD,L)
                 XNODE=COORD(NODE,1)
                 YNODE=COORD(NODE,2)
                 ZNODE=COORD(NODE,3)

C______________________ Step 2.2.2: Updates element cog

                 COGELEM(1)=COGELEM(1)+XNODE
                 COGELEM(2)=COGELEM(2)+YNODE
                 COGELEM(3)=COGELEM(3)+ZNODE

              END DO ! INUD=1,LNNDEL(L)

C______________________ Step 2.2.3: Calculates element cog

              COGELEM(1)=COGELEM(1)/DFLOAT(NNUD)
              COGELEM(2)=COGELEM(2)/DFLOAT(NNUD)
              COGELEM(3)=COGELEM(3)/DFLOAT(NNUD)

C______________________ Step 2.2.4: Adds contrib. of element cog. to block cog.

              COG_ZONE(1)=COG_ZONE(1)+COGELEM(1)*AREA(L)
              COG_ZONE(2)=COG_ZONE(2)+COGELEM(2)*AREA(L)
              COG_ZONE(3)=COG_ZONE(3)+COGELEM(3)*AREA(L)
              AREA_ZONE=AREA_ZONE+AREA(L)

            END IF ! LX.EQ.IACTZONE
           
         END DO ! L=1,NUMEL

C______________________ Step 2.3: Calculates zone cog

         COG_ZONE(1)=COG_ZONE(1)/AREA_ZONE
         COG_ZONE(2)=COG_ZONE(2)/AREA_ZONE
         COG_ZONE(3)=COG_ZONE(3)/AREA_ZONE

C_______________________ Step 3: Finds out actual zone (IACTZONE is global
C_______________________         numeration)

         DO IZON=1,NZN_GR
           DX2=(COG_ZONE(1)-POSZN_GS(IZON,1))**2
           DY2=(COG_ZONE(2)-POSZN_GS(IZON,2))**2
           DZ2=(COG_ZONE(3)-POSZN_GS(IZON,3))**2
           IF (DX2+DY2+DZ2 .LE. 1.E-6) GOTO 10
         END DO

C_______________________ Step 4: Assigns external drifts

 10      DO IEXT=1,NEXDR_GS
           VMEAS_GS(ISTART+I,MXNVAR_GS+IEXT)=EXDRZN_GS(IZON,IEXT)
         END DO

       END DO ! I=1,NPOINTS

*** VERIFYING EXTERNAL DRIFTS. COMPARISON WITH RESIDUAL KRIGING.
*** EXTERNAL DRIFT CAN BE DIRECTLY ASSIGNED TO COORDINATES
*** NEXT LINES SHOULD NOT BE COMMENTED. THE REST OF ROUTINE SHOULD BE

*        DO I=1,NPOINTS
*           VMEAS_GS(I+ISTART,2)=POSMEAS_GS(ISTART+I,1)
*           VMEAS_GS(I+ISTART,3)=POSMEAS_GS(ISTART+I,2)
*           VMEAS_GS(I,3)=POSMEAS_GS(I,1)*POSMEAS_GS(I,1)
*           VMEAS_GS(I,4)=POSMEAS_GS(I,1)*POSMEAS_GS(I,2)
*        END DO

      RETURN
      END

