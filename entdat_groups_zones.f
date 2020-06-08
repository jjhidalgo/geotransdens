      SUBROUTINE ENTDAT_GROUPS_ZONES
     ;(IDIMIVARIO_GS ,IDIMVAR_GS  ,IDIMWGT      ,IERROR     ,INPWR
     ;,IODIM         ,IOINV       ,IOWAR        ,ISOT       ,IUGEO
     ;,LMXNDL        ,MAINF       ,MXDISC_GS    ,MXGRPZN    ,MXMEASPP_GS
     ;,MXNPP_GS      ,MXNVAR_GS   ,MXNZON_GS    ,NGROUP_ZN  ,NPAREL
     ;,NPBMX         ,NTYPAR      ,NUMEL        ,NUMNP      ,NZPAR
     ;,AREA          ,COORD       ,COORDGR_GS   ,DIVZN_GS   ,EXDRZN_GS
     ;,FILENAME      ,INORPAR     ,IO_KG_GS     ,IOPT_GS    ,IPNT_PAR
     ;,IPOLDRIFT_GS  ,IVARIO_GS   ,IVPAR        ,KXX        ,LNNDEL
     ;,LTYPE         ,LXPAREL     ,NZONE_PAR    ,PARZ       ,POSDIS_GS
     ;,POSMEAS_GS    ,POSZN_GS    ,SEARCH_GS    ,SUPBL_GS   ,TRIM_GS
     ;,VARIO_GS      ,VMEAS_GS    ,VSTATS_GS    ,WGT_PAR    ,PARC_GS)

********************************************************************************
*
* PURPOSE Reads geostatistical data and linear combination of group of zones
*
* DESCRIPTION Summary:
*
*  - Step 0: Declaration of variables
*  - Step 1:Writes main header
*  Loop over group of zones. Writes head. for actual group
*  If group is estimated geostatistically:
*  - Step 2: Defines group geometry: cog of zones, offsets of discretization 
*            points wrt cog.and super blocks search grid (kriging)
*  - Step 3: Reads search ellipsoids data, trimming limits polynomial drift 
*            options and tolerance factor for conditional simulations
*  - Step 4: Reads external drift at zones cog. Only if kriging with locally
*            varying mean or kriging with external drift are done
*  - Step 5: Reads measurement locations. Only if unconditional simulation is
*            not done Also, checks that sampling point belongs to group. Only 
*            if ext. drift is not used. In that case, this will be checked 
*            elsewhere
*  - Step 6: If pilot point locations are fixed, they are read right now. 
*            If they are randomly drawn, will be done elsewhere. Also, checks 
*            that pilot point belongs to group. Only if ext. drift is not used. 
*            In that case, this will be checked elsewhere
*  - Step 7: Assigns external drifts to pilot point and sampling locations 
*            (zone to point assignment). Pilot points are assigned here only 
*            if the  position is fixed and inverse problem is solved
*  - Step 8: Echoes variables and drift values (point definition) + some stats
*  - Step 9: Reads corregionalization model (all variograms) of actual group 
*            of zones. If cokriging is done, coherency of the corregionalization
*            model is checked here
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  COORD                  Nodal coordinates                                     
*  COORDGR_GS             Maximum and minimum coodinates of all groups of zones
*  DIVZN_GS               Auxiliar array containing the sum of areas of elements
*                         drawing a given zone
*  EXDRZN_GS              Array containing zonal external drifts
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INORPAR                Array containing the indexes with the location        
*  IO_KG_GS               Kriging dimensions and options for geost. inv. prob. 
*                         Each row contains information of a given group of zones 
*                         Only sense if group is estimated geost. On each row:
*  IOPT_GS                General options for inverse problem. Each row contains 
*                         information of a given group of zones
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IPNT_PAR               Array contaning pointers to arrays DLT_PAR and WGT_PAR
*  IPOLDRIFT_GS           Array containing polinomial drift options
*  IVARIO_GS              Array containing integers defining all variograms
*  IVPAR                  Array containing zonal estimation options and pointers
*  IZONMEAS_GS            Array containing the measurement identifier contained in 
*                         a given zone of a given group
*  KXX                    Node numbers of every element (counterclockwise order)
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  LXPAREL                Zone number to which a given element belongs to        
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PARZ                   Array containing zonal parameters
*  POSDIS_GS              Array containing offsets of zone discretiz. points
*  POSMEAS_GS             Array containing position of sampling and pilot point 
*                         locations
*  POSZN_GS               Array containing position zones cog.
*  SEARCH_GS              Array containing geostatistical search options
*  SUPBL_GS               Array containing geostatistical search superblocks info.
*  TRIM_GS                Array containing trimming limits for all variables 
*                         and tolerance factor for conditional simulations
*  VARIO_GS               Array contaning variogram data (sill, nugget, etc)
*  VMEAS_GS               Array containing measurements values
*  VSTATS_GS              Array containing statistics of all variables
*  WGT_PAR                Array containing weights defining linear combinations
*
* INTERNAL VARIABLES: ARRAYS
*
*  CHARTYPE               Array of characters containing types of parameters
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMIVARIO_GS          Used to dimension IVARIO_GS, VARIO_GS
*  IDIMVAR_GS             NVAR_GS*NVAR_GS+NEXDR_GS
*  IDIMWGT                Used to dimension WGT_PAR, IPNT_PAR
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE
*  IODIM                  Maximum dimension of the problem                      
*  IOINV                  Inverse problem option
*  IOWAR                  Allows writing warning messages
*  ISOT                   Maximum hydraulic conductivity tensor anisotropy      
*                         degree in the problem.                                
*  IUGEO                  GEO file unit number
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MXDISC_GS              Maximum number of discretization points
*  MXGRPZN                Maximum number of groups of zones
*  MXMEASPP_GS            Maximum number of measurements+pilot points
*  MXNPP_GS               Maximum number of pilot points
*  MXNVAR_GS              Maximum number of variables in the most complex kriging 
*  MXNZON_GS              Maximum number of zones 
*  NFLAGS                 Used to dimension IFLAGS
*  NGROUP_ZN              Number of groups of zones
*  NPAREL                 Number of parameter types defined by elements
*  NPBMX                  Maximum number of flow/transport problems to be solved
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                        
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*
* INTERNAL VARIABLES: SCALARS
*
*  AUXSTRING              Auxiliar string containing the name of the param. type
*  I                      Dummy counter
*  IACTGROUP              Group to which a point belongs to
*  IACTZONE               Zone to which a point belongs to
*  IGROUP                 Dummy counter of groups
*  IO_CHECK_POS_GS        Option for checking if a meas / pp belongs to current group
*  IOGREST_GS             Estimation option of group of zones
*  IPARTYPE               Parameter type regarding group of zones
*  IPROBGR_GS             Flow/Transport defined for this group of zones
*  XPOINT                 X-coordinate of a given point
*  YPOINT                 Y-coordinate of a given point
*  ZPOINT                 Z-coordinate of a given point
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ASS_EXT_DRIFT
*  GET_ZONE_GROUP
*  READ_EXTDRIFT_ZONES
*  READ_LIN_COMB
*  READ_LOCATIONS
*  READ_SEARCH_DRIFT
*  READ_VARIOGRAMS
*  WRITE_VAR_GS
*  ZERO_ARRAY  
*  ZONE_GEOMETRY                                                                 
*
* HISTORY:  AAR  First coding (Dec-2001)
*           AAR  Revision and header (April-2002)
*           AAR  Revision and inclusion of groups of zones (July-2003)
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NGROUP_ZN,MXGRPZN,INPWR,MAINF,IODIM,LMXNDL,MXDISC_GS
     ;         ,MXNZON_GS,NTYPAR,NUMEL,NUMNP,NZPAR,NPAREL,NPBMX
     ;         ,IERROR,IOWAR,IUGEO,MXNVAR_GS,MXMEASPP_GS,IDIMVAR_GS
     ;         ,MXNPP_GS,IOINV,IDIMIVARIO_GS,IDIMWGT,ISOT
     ;         ,IOPT_GS(MXGRPZN,20),IO_KG_GS(MXGRPZN,16)
     ;         ,IVPAR(NZPAR,4),INORPAR(NTYPAR),KXX(LMXNDL,NUMEL)
     ;         ,LNNDEL(NUMEL),LXPAREL(NUMEL,NPAREL,NPBMX),LTYPE(NUMEL)
     ;         ,NZONE_PAR(NTYPAR),IPOLDRIFT_GS(9,NGROUP_ZN)
     ;         ,IPNT_PAR(IDIMWGT*NZPAR)
     ;         ,IVARIO_GS(IDIMIVARIO_GS,2,NGROUP_ZN)
*,IZONMEAS_GS(MXNZON_GS,NGROUP_ZN)
                                                                 ! Real external
      REAL*8 AREA(NUMEL),SUPBL_GS(6,NGROUP_ZN),COORD(NUMNP,3)
     ;      ,DIVZN_GS(MXNZON_GS),POSZN_GS(MXNZON_GS,3,NGROUP_ZN)
     ;      ,POSDIS_GS(MXDISC_GS,3,MXNZON_GS,NGROUP_ZN)
     ;      ,SEARCH_GS(11,NGROUP_ZN),TRIM_GS(8,NGROUP_ZN)
     ;      ,EXDRZN_GS(MXNZON_GS,4,NGROUP_ZN),COORDGR_GS(6,NGROUP_ZN)
     ;      ,POSMEAS_GS(MXMEASPP_GS,3,2,NGROUP_ZN),PARZ(NZPAR)
     ;      ,VMEAS_GS(MXMEASPP_GS,IDIMVAR_GS,2,NGROUP_ZN)
     ;      ,VSTATS_GS(MXNVAR_GS,4,NGROUP_ZN),WGT_PAR(IDIMWGT*NZPAR)
     ;      ,VARIO_GS(IDIMIVARIO_GS,8,NGROUP_ZN)
     ;      ,PARC_GS(MXNPP_GS,NGROUP_ZN)
                                                              ! Integer internal
      INTEGER*4 IGROUP,IPARTYPE,IPROBGR_GS,IOGREST_GS,I,IACTGROUP
     ;         ,IACTZONE,IPLXPAREL,NROW,IO_CHECK_POS_GS
                                                                 ! Real internal
      REAL*8 XPOINT,YPOINT,ZPOINT
                                                                    ! Characters
      CHARACTER FILENAME(20)*20,CHARTYPE(10)*3,AUXSTRING*9
     ;         ,CHARTYPE2(15)*3,LEAUX*100,LEEL*100

      DATA CHARTYPE/'TRA','STG','ARR','ART','DSL','DST','DFM','POR',
     ;              'FOD','CRD'/
      DATA CHARTYPE2/'TRA','STG','ARR','CHP','QQP','ALF','DSL','DST',
     ;               'DFM','POR','FOD','CRD','CON','AGE','DMT'/

C_______________________ Step 1:Writes main header

      IF (INPWR.NE.0) WRITE(MAINF,2000)
 2000   FORMAT(///,29X,' GEOSTATISTICAL/INTERPOLATION DATA',/
     ;            ,29X,' ============================ ====',/)

C_______________________ Loop over group of zones

      DO IGROUP=1,NGROUP_ZN

         IPARTYPE=IOPT_GS(IGROUP,1)   ! Parameter type
         IPROBGR_GS=IOPT_GS(IGROUP,12)   ! Flow/transport problem
         IOGREST_GS=IOPT_GS(IGROUP,2)    ! Estimation option
         IO_CHECK_POS_GS=IOPT_GS(IGROUP,14)    ! Checking positions
         AUXSTRING=' INTERP. '
         IF (IOGREST_GS.EQ.1) AUXSTRING=' GEOSTAT.'

C_______________________ Writes header for actual group

        IF (INPWR.NE.0) THEN
           IF (IOGREST_GS.EQ.1) 
     ;        WRITE(MAINF,2100) IGROUP,CHARTYPE(IPARTYPE),AUXSTRING
           IF (IOGREST_GS.EQ.2) 
     ;        WRITE(MAINF,2100) IGROUP,CHARTYPE2(IPARTYPE),AUXSTRING
    
 2100      FORMAT(//,30X,' GROUP OF ZONES NUMBER: ',I5,/,
     ;               30X,' ===== == ===== ======= ',//,
     ;               5X,' PARAMETER TYPE: ',A3,/,
     ;               5X,' ========= ===== ',//,
     ;               5X,' ESTIMATION TYPE:',A9,/,
     ;               5X,' ========== =====',//)
        END IF

        IF (IOGREST_GS.EQ.1) THEN                    ! Geostatistical estimation

C_______________________ Step 2: Defines group geometry: cog of zones,
C_______________________         offsets of discretization points wrt cog.
C_______________________         and super blocks search grid (kriging)

          CALL ZONE_GEOMETRY
     ;(IPARTYPE    ,IGROUP   ,INPWR      ,IODIM     ,LMXNDL    
     ;,MAINF       ,IO_KG_GS(IGROUP,6)   ,IO_KG_GS(IGROUP,7)       
     ;,IO_KG_GS(IGROUP,8)    ,MXDISC_GS  ,MXNZON_GS ,IO_KG_GS(IGROUP,9)   
     ;,IO_KG_GS(IGROUP,10)   ,IO_KG_GS(IGROUP,11)
     ;,NPAREL      ,NTYPAR   ,NUMEL      ,NUMNP     ,IOPT_GS(IGROUP,7)
     ;,NZPAR       ,AREA     ,SUPBL_GS(1,IGROUP)    ,COORD     
     ;,COORDGR_GS(1,IGROUP)  ,DIVZN_GS   ,IVPAR(1,3)
     ;,INORPAR     ,KXX      ,LNNDEL     ,LXPAREL(1,1,IPROBGR_GS)
     ;,NZONE_PAR   ,POSZN_GS(1,1,IGROUP) ,POSDIS_GS(1,1,1,IGROUP)
     ;,SUPBL_GS(4,IGROUP)) 

C_______________________ Step 3: Reads search ellipsoids data, trimming limits
C_______________________         and polynomial drift options 

          CALL READ_SEARCH_DRIFT
     ;(IO_KG_GS(IGROUP,1) ,IERROR   ,INPWR
     ;,IOWAR              ,IUGEO    ,MAINF       ,IO_KG_GS(IGROUP,13)
     ;,IOPT_GS(IGROUP,3)  ,VSTATS_GS(1,2,IGROUP) ,IPOLDRIFT_GS(1,IGROUP) 
     ;,SEARCH_GS(1,IGROUP),TRIM_GS(1,IGROUP)     ,FILENAME)

C_______________________ Step 4: Reads external drift at zones cog. Only if 
C_______________________         kriging with locally varying mean or kriging 
C_______________________         with external drift are done

          IF (IO_KG_GS(IGROUP,1).EQ.2 .OR. IO_KG_GS(IGROUP,1).EQ.3)
     ;        CALL READ_EXTDRIFT_ZONES
     ;(IERROR   ,INPWR       ,IOWAR       ,IUGEO   ,MAINF    ,MXNZON_GS
     ;,IOPT_GS(IGROUP,7)     ,IO_KG_GS(IGROUP,13)  ,POSZN_GS(1,1,IGROUP)
     ;,EXDRZN_GS(1,1,IGROUP) ,FILENAME)

C_______________________ Step 5: Reads measurement locations. Only if 
C_______________________         unconditional simulation is not done Also, 
C_______________________         checks that sampling point belongs to group. 
C_______________________         Only if ext. drift is not used. In that case,
C_______________________         this will be checked elsewhere (ASS_EXT_DRIFT
C_______________________         routine). If conditional simulation is done,
C_______________________         saves the measurement identifier.

          IF (IO_KG_GS(IGROUP,14).NE.2) THEN
   
            CALL READ_LOCATIONS
     ;(1                 ,MXMEASPP_GS        ,IDIMVAR_GS ,IERROR   
     ;,INPWR             ,IOPT_GS(IGROUP,6)  ,IOWAR   ,IOPT_GS(IGROUP,5)
     ;,IUGEO             ,IO_KG_GS(IGROUP,1) ,MAINF   ,IOPT_GS(IGROUP,4)
     ;,IOPT_GS(IGROUP,4) ,IOPT_GS(IGROUP,3)  ,POSMEAS_GS(1,1,1,IGROUP)
     ;,TRIM_GS(1,IGROUP) ,VMEAS_GS(1,1,1,IGROUP),VSTATS_GS(1,1,IGROUP)
     ;,FILENAME          ,MXNVAR_GS)           

*__________________ Checks if measurement location belongs to domain and group

            IF (IO_CHECK_POS_GS.EQ.0) THEN
               IF (IO_KG_GS(IGROUP,1).NE.2 
     ;             .AND. IO_KG_GS(IGROUP,1).NE.3) THEN

                   DO I=1,IOPT_GS(IGROUP,4)

                     XPOINT=POSMEAS_GS(IOPT_GS(IGROUP,5)+I,1,1,IGROUP)
                     YPOINT=POSMEAS_GS(IOPT_GS(IGROUP,5)+I,2,1,IGROUP)
                     ZPOINT=POSMEAS_GS(IOPT_GS(IGROUP,5)+I,3,1,IGROUP)

                     CALL GET_ZONE_GROUP
     ;(IACTGROUP  ,IPARTYPE   ,IACTZONE   ,IERROR   ,IGROUP     ,IOWAR
     ;,IPLXPAREL  ,1          ,IUGEO      ,LMXNDL   ,MAINF      ,NPAREL
     ;,NUMEL      ,NUMNP      ,NTYPAR     ,NZPAR    ,XPOINT     ,YPOINT
     ;,ZPOINT     ,AREA       ,COORD      ,FILENAME ,IVPAR(1,3) ,INORPAR
     ;,KXX        ,LTYPE      ,LXPAREL)

*                     IZONMEAS_GS(IACTZONE,IGROUP)=I

                   END DO ! I=1,IOPT_GS(IGROUP,5)

                END IF   ! IO_KG_GS...

            END IF ! IO_CHECK_POS_GS.EQ.0

          END IF ! IO_KG_GS(IGROUP,14).NE.2

C_______________________ Step 6: If pilot point locations are fixed, they are 
C_______________________         read right now. If they are randomly drawn,
C_______________________         will be done elsewhere. Also, checks that 
C_______________________         pilot point belongs to group. Only if ext.
C_______________________         drift is not used. In that case, this will be
C_______________________         checked elsewhere (ASS_EXT_DRIFT routine)

          IF (IOINV.GT.0 .AND. IOPT_GS(IGROUP,6).EQ.0) THEN 
      
             CALL READ_LOCATIONS
     ;(2                 ,MXMEASPP_GS        ,IDIMVAR_GS ,IERROR   
     ;,INPWR             ,IOPT_GS(IGROUP,6)  ,IOWAR      ,0
     ;,IUGEO             ,IO_KG_GS(IGROUP,1) ,MAINF   ,IOPT_GS(IGROUP,5)
     ;,IOPT_GS(IGROUP,4) ,IOPT_GS(IGROUP,3)  ,POSMEAS_GS(1,1,1,IGROUP)
     ;,TRIM_GS(1,IGROUP) ,VMEAS_GS(1,1,1,IGROUP),VSTATS_GS(1,1,IGROUP)
     ;,FILENAME          ,MXNVAR_GS)           

ccc ojo, modificacion 

            DO I=1,IOPT_GS(IGROUP,5)
               WRITE(779,111) POSMEAS_GS(I,1,1,IGROUP),
     ;                        POSMEAS_GS(I,2,1,IGROUP),
     ;                        POSMEAS_GS(I,3,1,IGROUP)
            END DO
111         FORMAT(3(5X,F10.3))

ccc fin del ojo

*____________ Checks if the position of the pilot point belongs to domain / group

            IF (IO_CHECK_POS_GS.EQ.0) THEN
               IF (IO_KG_GS(IGROUP,1).NE.2 
     ;             .AND. IO_KG_GS(IGROUP,1).NE.3) THEN

                   DO I=1,IOPT_GS(IGROUP,5)

                     XPOINT=POSMEAS_GS(I,1,1,IGROUP)
                     YPOINT=POSMEAS_GS(I,2,1,IGROUP)
                     ZPOINT=POSMEAS_GS(I,3,1,IGROUP)

                     CALL GET_ZONE_GROUP
     ;(IACTGROUP  ,IPARTYPE   ,IACTZONE   ,IERROR   ,IGROUP     ,IOWAR
     ;,IPLXPAREL  ,2          ,IUGEO      ,LMXNDL   ,MAINF      ,NPAREL   
     ;,NUMEL      ,NUMNP      ,NTYPAR     ,NZPAR    ,XPOINT     ,YPOINT   
     ;,ZPOINT     ,AREA       ,COORD      ,FILENAME ,IVPAR(1,3) ,INORPAR
     ;,KXX        ,LTYPE      ,LXPAREL)

                   END DO

                END IF ! IO_KG_GS...

            END IF   ! IO_CHECK_POS_GS.EQ.0

          END IF ! IOINV.GT.0 .AND. IOPT_GS(IGROUP,6).EQ.0

C_______________________ Step 7: Assigns external drifts to pilot point and 
C_______________________         sampling locations (zone to point assignment).
C_______________________         Pilot points are assigned here only if the 
C_______________________         position is fixed and inverse problem is solved

          IF (IO_KG_GS(IGROUP,1).EQ.2 .OR. IO_KG_GS(IGROUP,1).EQ.3) THEN

                                                                  ! Measurements
             IF (IO_KG_GS(IGROUP,14).NE.2) CALL ASS_EXT_DRIFT
     ;(IPARTYPE          ,1         ,IDIMVAR_GS          ,IERROR
     ;,IGROUP            ,IOWAR     ,IOPT_GS(IGROUP,5)   ,IUGEO
     ;,LMXNDL            ,MAINF     ,MXMEASPP_GS         ,MXNZON_GS
     ;,MXNVAR_GS         ,IO_KG_GS(IGROUP,13)            ,NPAREL
     ;,IOPT_GS(IGROUP,4) ,NTYPAR    ,NUMEL               ,NUMNP
     ;,IOPT_GS(IGROUP,7) ,NZPAR     ,AREA                ,COORD
     ;,EXDRZN_GS(1,1,IGROUP)        ,FILENAME
     ;,IVPAR(1,3)        ,INORPAR   ,KXX                 ,LNNDEL
     ;,LTYPE             ,LXPAREL   ,POSMEAS_GS(1,1,1,IGROUP)
     ;,POSZN_GS(1,1,IGROUP)         ,VMEAS_GS(1,1,1,IGROUP))

             IF (IOINV.GT.0 .AND. IOPT_GS(IGROUP,6).EQ.0) 
     ;          CALL ASS_EXT_DRIFT       ! Pilot points
     ;(IPARTYPE          ,2                   ,IDIMVAR_GS  ,IERROR
     ;,IGROUP            ,IOWAR               ,0           ,IUGEO
     ;,LMXNDL            ,MAINF               ,MXMEASPP_GS ,MXNZON_GS
     ;,MXNVAR_GS         ,IO_KG_GS(IGROUP,13) ,NPAREL
     ;,IOPT_GS(IGROUP,5) ,NTYPAR              ,NUMEL       ,NUMNP
     ;,IOPT_GS(IGROUP,7) ,NZPAR               ,AREA        ,COORD
     ;,EXDRZN_GS(1,1,IGROUP)                  ,FILENAME
     ;,IVPAR(1,3)        ,INORPAR             ,KXX         ,LNNDEL
     ;,LTYPE             ,LXPAREL             ,POSMEAS_GS(1,1,1,IGROUP)
     ;,POSZN_GS(1,1,IGROUP)                   ,VMEAS_GS(1,1,1,IGROUP))
*     ;,IZONMEAS_GS(1,1))

          END IF ! IO_KG_GS(IGROUP,1).EQ.2 ...

C_______________________ Step 8: Echoes variables and drift values 
C_______________________         (point definition) + some stats

          IF (INPWR.NE.0) CALL WRITE_VAR_GS
     ;(IDIMVAR_GS  ,IOINV       ,IO_KG_GS(IGROUP,14) ,IO_KG_GS(IGROUP,1)
     ;,MAINF       ,MXMEASPP_GS ,MXNVAR_GS           ,IOPT_GS(IGROUP,4)        
     ;,IOPT_GS(IGROUP,5)        ,IOPT_GS(IGROUP,3)         
     ;,POSMEAS_GS(1,1,1,IGROUP) ,VMEAS_GS(1,1,1,IGROUP)
     ;,VSTATS_GS(1,1,IGROUP))

C_______________________ Step 9: Reads corregionalization model (all variograms)
C_______________________         of actual group of zones. If cokriging is done, 
C_______________________         coherency of the corregionalization model is 
C_______________________         checked here

          CALL READ_VARIOGRAMS
     ;(IDIMIVARIO_GS         ,IERROR               ,INPWR        ,IOWAR
     ;,IUGEO                 ,IO_KG_GS(IGROUP,1)   ,MAINF
     ;,IO_KG_GS(IGROUP,12)
     ;,IOPT_GS(IGROUP,3)     ,VARIO_GS(1,2,IGROUP) ,VARIO_GS(1,4,IGROUP)
     ;,VARIO_GS(1,5,IGROUP)  ,VARIO_GS(1,6,IGROUP) ,VARIO_GS(1,7,IGROUP)
     ;,VARIO_GS(1,8,IGROUP)  ,VARIO_GS(1,1,IGROUP) ,VARIO_GS(1,3,IGROUP)
     ;,IVARIO_GS(1,2,IGROUP) ,IVARIO_GS(1,1,IGROUP),FILENAME)

C_______________________ Step 10: Reads INITIAL VALUE at pilot point locations

            IF (IOPT_GS(IGROUP,13).EQ.3) THEN
             DO I=1,IOPT_GS(IGROUP,5)
                LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
                READ(LEAUX,*) PARC_GS(I,IGROUP)
              END DO
          END IF ! IOPT_GS(IGROUP,13).EQ.3

        ELSE IF (IOGREST_GS.EQ.2) THEN                          ! Interpolation

          CALL READ_LIN_COMB
     ;(IDIMWGT            ,IERROR    ,IGROUP     ,INPWR      ,IODIM     
     ;,IOWAR              ,ISOT      ,IUGEO      ,MAINF      ,MXGRPZN 
     ;,IOPT_GS(IGROUP,11) ,NTYPAR    ,NZPAR      ,FILENAME   ,IVPAR(1,3)
     ;,INORPAR            ,IOPT_GS   ,IPNT_PAR   ,IVPAR(1,2) ,IVPAR(1,1)
     ;,NZONE_PAR          ,PARZ      ,WGT_PAR)

        END IF ! IOGREST_GS.EQ.1

      END DO ! IGROUP=1,NGROUP_ZN

      RETURN

      END
