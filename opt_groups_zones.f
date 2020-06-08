      SUBROUTINE OPT_GROUPS_ZONES
     ;(IERROR   ,INPWR    ,IOINV    ,IOWAR    ,ISOT      ,IUDIM     
     ;,IUGEO    ,MAINF    ,MXGRPZN  ,MXLINCMB ,NGROUP_ZN,NTYPAR
     ;,FILENAME ,IOPT_GS  ,IO_KG_GS)

********************************************************************************
*
* PURPOSE Reads groups of zones estimation options and dimensions related to 
*         inverse problem. Dimensions and options read here define groups of
*         zones that are estimated geostatistically or interpolated (zonal 
*         parameters defining this group are a linear combination of zonal 
*         other parameters, estimated deterministically)
*
* DESCRIPTION 
*
* - Step 0: Declaration of variables
* - Step 1: Reads FLAG variable. Not the first time. If FLAG is negative, then
*           no more groups of zones have to be read
* - Step 2: Reads geostatistical options of a group of zones. Some quick checks 
*           are done
* - Step 3: Checks geostatistical options of a group of zones. Only sense if
*           group is estimated geostatistically (not a linear combination of          
*           unknowns)
*      0) Kriging/cokriging type out of range
*      1) Number of variables used for kriging/cokriging out of range
*      2) Checks number of variables for cokriging
*      3) Warning: Kriging requested and NVAR>1 or NDMAXS>0
*      4) Number of discretization points out of range
*      5) Minimum number of samples for kriging a point/zone out of range
*      6) Maximum number of primary/secondary samples out of range
*      7) Secondary samples are not taken into account in cokriging
*      8) Maximum number of nested structures defining variogram out of range
*      9) Number of pilot points or drawing option out of range. No sense if 
*         only direct problem is solved
*      10) Number of measurements out of range. No sense if unconditional 
*          simulation is done
*      11) Number of SUPER BLOCKS out of range
*      12) Anisotropy degree out of range
*      13) Number of external drift terms
*      14) Option of conditional simulation out of range
*      15) Seed for generating conditional simulations is negative
*      16) Number of conditional simulations out of range
*      17) Number of groups of zones exceeds maximum allowed
* - Step 4: Echoes information if so desired
* - Step 5: Stores options in big arrays (useful portability)
* - Step 6: Updates number of groups of zones read 
* 
* If zonal parameters belonging to this group are a linear combination of
* unknows, reads, checks and stores the number of terms defining the linear 
* combination
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  IOPT_GS                General options for inverse problem. 
*                         Each row contains information of a given group of zones
*                         On each row:
*                         - Column 1: Identifier of parameter type
*                         - Column 2: Estimation option. (0: zones belonging to
*                                     this group are estimated deterministically;
*                                     1: geostatistically;2:linear combination)
*                         - Column 3: Number of variables used for krig./cokrig.
*                         - Column 4: Number of measurements used for k/ck
*                         - Column 5: Number of pilot points used for k/ck
*                         - Column 6: Pilot points drawing option
*                                     0: Read and fixed
*                                     1: Completely Random and variable trough
*                                        optimization process
*                                     2: Layed on a random regular mesh and 
*                                        variable (that mesh) trough opt. proc.
*                         - Column 7: Number of zones in this group
*                         - Column 8: Density of pilot points (X-direction)
*                         - Column 9: Density of pilot points (Y-direction)
*                         - Column 10: Density of pilot points (Z-direction)
*                         - Column 11: Number of terms defining linear combination
*                         - Column 12: Flow/transport problem
*                         - Column 13: Option for initial PARC
*                                      0) Equal to initial estimation at pilot 
*                                         points (=PARM)
*                                      1) The same but adding a white noise of 
*                                         given variance
*                                      2) Equal to zero
*                                      3) Read at GEO file
*                         - Column 14: Option for checking if a meas / pp location
*                                      belongs to the group of zones
*  IO_KG_GS               Kriging dimensions and options for geost. inv. prob. 
*                         Each row contains information of a given group of zones 
*                         Only sense if group is estimated geost. On each row:
*                         - Column 1: Kriging/Cokriging type
*                                     0: Simple kriging
*                                     1: Residual kriging
*                                     2: Kriging with locally varying mean
*                                     3: Kriging with external drift (up to four)
*                                     4: Simple cokriging
*                                     5: Standardized ordinary cokriging
*                                     6: Traditional ordinary cokriging
*                         - Column 2: Minimum number of samples to be considered 
*                                     in the kriging estimation. Point/Zone will 
*                                     remain unestimated if number of closest 
*                                     samples is less than this number
*                         - Column 3: Maximum number of primary variable samples 
*                                     used for kriging
*                         - Column 4: Maximum number of secondary samples used for 
*                                     kriging (all secondary variable samples are
*                                     jointly considered)
*                         - Column 5: Number of samples retain per octant. If 0, 
*                                     octant search is not performed
*                         - Column 6: Number of super-blocks in X-direction (super-
*                                     block search mesh)
*                         - Column 7: Number of super-blocks in Y-direction (super-
*                                     block search mesh)
*                         - Column 8: Number of super-blocks in Z-direction (super-
*                                     block search mesh)
*                         - Column 9: Number of discretization points of a zone in
*                                     X-direction.
*                         - Column 10: Number of discretization points of a zone in
*                                      Y-direction.
*                         - Column 11: Number of discretization points of a zone in
*                                      Z-direction.
*                                      Note: If product of last three numbers is 1, 
*                                            then point kriging is performed
*                         - Column 12: Number of nested structures defining the most
*                                      complex variogram of this group of zones 
*                                      (unique variogram if kriging is performed)
*                         - Column 13: Number of external drift terms
*                         - Column 14: Option of conditional/unconditional simulation
*                                      - 0: Conditional Estimation
*                                      - 1: Conditional simulation
*                                      - 2: Unconditional simulation
*                         - Column 15: Seed for conditional/unconditional simulations
*                         - Column 16: Number of conditional simulations
*
* INTERNAL VARIABLES: ARRAYS
*
*  CHARTYPE               Array containing name of parameter types (pilot points; 
*                         therefore, only parameter types defined by elements)
*  CHARTYPE2              Array containing name of parameter types (interpolation; 
*                         therefore, all types)
*  MXSB_GS                Array containing super block search grid definition
*  NDISC_GS               Array containing number of discr. points of zones
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOINV                  Inverse problem option                                
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  ISOT                   Maximum hydraulic conductivity tensor anisotropy      
*                         degree in the problem.                                
*  IUDIM                  DIM file unit number
*  IUGEO                  GEO file unit number
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MXGRPZN                Maximum number of group of zones (used for dim.)
*  MXLINCMB               Maximum number of unknown parameters defining a linear
*                         combination
*  NGROUP_ZN              Number of groups of zones (AS READ)
*  NLINCMB                Number of terms defining a linear combination (AS READ)
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*
* INTERNAL VARIABLES: SCALARS
*
*  FLAG                   If =>0, keep on reading groups of zones
*  I                      Dummy counter
*  ID_TP_GS               Identifier of parameter type
*  IO_CHECK_POS_GS        Option for checking if a meas / pp location
*                         belongs to the group of zones
*  IO_GREST_GS            Group Estimation option (AS READ)
*  IO_PARC_INI_GS         Option for initial value at pilot points
*  IORD_PP_GS             Pilot points drawing option
*  IOPSC_GS               Option of conditional/unconditional simulation
*  IPROBGR_GS             Flow/transport problem
*  ISEEDSC_GS             Seed for conditional simulations
*  J                      Dummy counter
*  KTYPE_GS               Kriging/simulation type
*  MX_SC_GS               Number of conditional simulations
*  MXVGM_GS               Number of nested structures defining the most complex 
*                         variogram of this group of zones
*  NZON_GS                Number of zones
*  NDMIN_GS               Minimum number of samples to consider in the kriging sys. 
*                         Point/Block will remain unestimated if number of closest 
*                         samples is below this number
*  NEXDR_GS               Number of external drift terms
*  NFM                    Auxiliar (NFORM_GS+1)
*  NMEAS_GS               Number of measurements used for k/ck
*  NMXP_GS                Maximum number of primary variable samples used for krig.
*  NMXS_GS                Maximum number of secondary samples used for kriging
*  NOCT_GS                Number of samples retain per octant. If 0, octant search 
*                         is not performed
*  NPP_GS                 Number of pilot points used for k/ck
*  NPPX_GS                Density of pilot points (X-direction)           
*  NPPY_GS                Density of pilot points (Y-direction)
*  NPPZ_GS                Density of pilot points (Z-direction)
*  NROW                   Column atread file with wrong format
*  NVAR_GS                Number of variables used for krig./cokrig.

* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*
* HISTORY
*
*     AAR      6-2003     First coding
*
********************************************************************************

      IMPLICIT NONE

C_______________________ Step 0: Declaration of variables

                                                              ! Integer external
      INTEGER*4 IERROR,INPWR,IOWAR,ISOT,IUDIM,IUGEO,MAINF,MXGRPZN
     ;         ,NGROUP_ZN,NTYPAR,MXLINCMB,IOINV
     ;         ,IOPT_GS(MXGRPZN,20)
     ;         ,IO_KG_GS(MXGRPZN,16)
                                                              ! Integer internal
      INTEGER*4 FLAG,I,IO_GREST_GS,ID_TP_GS,IORD_PP_GS
     ;         ,KTYPE_GS,MXVGM_GS,NZON_GS,NDMIN_GS
     ;         ,NEXDR_GS,NFM,NMEAS_GS,NMXP_GS,NMXS_GS
     ;         ,NOCT_GS,NPP_GS,NPPX_GS,NPPY_GS
     ;         ,NPPZ_GS,NROW,NVAR_GS,IOPSC_GS,ISEEDSC_GS
     ;         ,MX_SC_GS,NLINCMB,IPROBGR_GS,IO_PARC_INI_GS
     ;         ,IO_CHECK_POS_GS
     ;         ,MXSB_GS(3),NDISC_GS(3)
                                                                    ! Characters
      CHARACTER FILENAME(20)*20,LEAUX*100,LEEL*100,CHARTYPE(10)*3
     ;         ,CHARTYPE2(15)*3

      DATA CHARTYPE/'TRA','STG','ARR','ART','DSL','DST','DFM','POR',
     ;              'FOD','CRD'/   ! Pilot points. Only param. by elements
      DATA CHARTYPE2/'TRA','STG','ARR','CHP','QQP','ALF','DSL','DST'
     ;              ,'DFM','POR','FOD','CRD','CON','AGE','DMT'/
      
C_______________________ Step 1: Reads FLAG variable. Not first time

      NGROUP_ZN=0                          ! Total number of group of zones read 
      FLAG=0                              ! Flag to keep reading groups of zones
      MXLINCMB=0                    ! Max number of unknows defining a lin. comb.

      IF (INPWR.NE.0 .AND. NGROUP_ZN.EQ.0) WRITE(MAINF,2400)
 2400   FORMAT(///,29X,' GEOSTATISTICAL/INTERPOLATION OPTIONS',/
     ;            ,29X,' ============================ =======',/)

 10   IF (NGROUP_ZN.GT.0) THEN
         LEAUX=LEEL(FILENAME,16,MAINF,NROW,INPWR)
         READ(LEAUX,1000) FLAG
      ENDIF

      IF (FLAG.LT.0) RETURN                 ! No more groups of zones to be read

      NFM=NGROUP_ZN+1
      IF (INPWR.NE.0) WRITE(MAINF,2500) NFM
 2500   FORMAT(//,25X,' GROUP OF ZONES NUMBER: ',I5,/
     ;           ,25X,' ===== == ===== =======',/)

C_______________________ Step 2: Reads geo. options of a group of zones. Some
C_______________________         quick checks are done

                                                   ! Card G1.1. Common variables

      LEAUX=LEEL(FILENAME,16,MAINF,NROW,INPWR)
      READ(LEAUX,1000)  ID_TP_GS             ! ID. Type of parameter
     ;                 ,IO_GREST_GS          ! Estimation option
     ;                 ,NVAR_GS              ! Number of variables used
     ;                 ,NMEAS_GS             ! Number of measurements used
     ;                 ,NPP_GS               ! Number of pilot points used
     ;                 ,IORD_PP_GS           ! Option for drawing p. p.
     ;                 ,NZON_GS              ! Number of zones
     ;                 ,NPPX_GS              ! Density of pilot points X
     ;                 ,NPPY_GS              ! Density of pilot points Y
     ;                 ,NPPZ_GS              ! Density of pilot points Z
     ;                 ,IPROBGR_GS           ! Flow/transport problem
     ;                 ,IO_PARC_INI_GS       ! Opt. for init. value of PARC
     ;                 ,IO_CHECK_POS_GS      ! Opt. for checking meas and pp loc
 1000 FORMAT(16I5)

* 1) ID of parameter type or of group of zones out of range

      IF (ID_TP_GS.LE.0.OR.ID_TP_GS.GT.NTYPAR) THEN
         WRITE(MAINF,2000) NGROUP_ZN+1
         WRITE(*,2000) NGROUP_ZN+1
 2000    FORMAT(//,' PARAMETER TYPE OUT OF RANGE READING GROUP OF'
     ;             ' ZONES:',I5,' FORCED STOP')
         STOP
      END IF

* 2) Estimation option out of range

      IF (IO_GREST_GS.LT.0. OR. IO_GREST_GS.GT.2) THEN
         WRITE(MAINF,2100) NGROUP_ZN+1
         WRITE(*,2100) NGROUP_ZN+1
 2100    FORMAT(//,' ESTIMATION OPTION OF GROUP OF ZONES:',I5,
     ;             ' CANNOT <0 OR >2. FORCED STOP')
         STOP
      END IF

      IF (IO_GREST_GS.EQ.1) THEN    ! Group is estimated geostatistically

                          ! Card G.2. Kriging options and variogram dimensions

         LEAUX=LEEL(FILENAME,16,MAINF,NROW,INPWR)
         READ(LEAUX,1100) KTYPE_GS             ! Kriging/simulation type
     ;                   ,NDMIN_GS             ! Min number of samples
     ;                   ,NMXP_GS              ! Max of primary variable
     ;                   ,NMXS_GS              ! Max od secondary varS.
     ;                   ,NOCT_GS              ! Octant search
     ;                   ,(MXSB_GS(I),I=1,3)   ! Num of superblocks XYZ
     ;                   ,(NDISC_GS(I),I=1,3)  ! Num of discr. points XYZ
     ;                   ,MXVGM_GS             ! Num of nested structures
     ;                   ,NEXDR_GS             ! Num of external drifts
     ;                   ,IOPSC_GS             ! Option for cond. simul.
     ;                   ,ISEEDSC_GS           ! Seed for cond. simul.
     ;                   ,MX_SC_GS             ! Number of cond. simul.

 1100    FORMAT(14I4,I10,I4)

C_______________________ Step 3: Checks geo. options of a group of zones. 

* 0) Checks kriging type

         IF (KTYPE_GS.LT.0.OR.KTYPE_GS.GT.6) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,'KRIGING TYPE OUT OF RANGE.'
     ;,NROW,1,IUGEO,2,9.0)  

*** PROVISIONAL ***

         IF (KTYPE_GS.GT.3) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,'COKRIGING NOT PROGRAMMED YET'
     ;,NROW,1,IUGEO,2,9.0)  

* 1) Number of variables for kriging out of range
       
         IF (NVAR_GS.LT.1 .OR. NVAR_GS.GT.4) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'NUMBER OF VARIABLES USED FOR KRIGING OUT OF RANGE (<1 OR >4)'
     ;,NROW,1,IUGEO,2,9.1)

* 2) Checks number of variables for cokriging

         IF (KTYPE_GS.GE.4 .AND. NVAR_GS.EQ.1) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'COKRIGING REQUESTED AND NVAR_GS=1. SHOULD BE >1'
     ;,NROW,1,IUGEO,2,9.2)

* 3) Warning: Kriging requested and NVAR>1 or NDMAXS>0

         IF (KTYPE_GS.LE.3 .AND. 
     ;      (NVAR_GS.GT.1.OR.NMXS_GS.NE.0)) THEN 
            CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' KRIGING REQUESTED AND NVAR OR NDMXS_GS' 
     ;//' ARE OUT OF RANGE. RESETED TO 1 AND 0 RESPECTIVELY.'
     ;,NROW,1,IUGEO,0,9.3)

            NVAR_GS=1
            NMXS_GS=0
         END IF

* 4) Number of discretization points out of range

         IF (NDISC_GS(1)*NDISC_GS(2)*NDISC_GS(3).LE.0) 
     ;      CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'CONCEPTUAL ERROR. NUMBER OF BLOCK DISCRET. POINTS IS ZERO'
     ;  ,NROW,2,IUGEO,2,9.4)

* 5) Minimum number of samples for kriging a block out of range
  
         IF (NDMIN_GS.LE.0) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'CONCEPTUAL ERROR. MIN. NUMBER OF SAMPLES FOR INTERP. IS ZERO'
     ;,NROW,2,IUGEO,2,9.5)

* 6) Maximum number of primary/secondary samples out of range

         IF (NMXP_GS.LT.NDMIN_GS .OR.
     ;   (NVAR_GS.GT.1.AND.NMXS_GS+NMXP_GS.LT.NDMIN_GS)) 
     ;     CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'CONCEPTUAL ERROR. ' 
     ;//'MAXIMUM NUMBER OF SAMPLES FOR INTERP. IS OUT OF RANGE',
     ;   NROW,2,IUGEO,2,9.6)

* 7) Secondary samples are not taken into account in cokriging

         IF (KTYPE_GS.GE.4 .AND. KTYPE_GS.LE.6 .AND. 
     ;    ((NMXS_GS/2).LE.NVAR_GS.OR.NMXS_GS.EQ.0)) THEN
           IF (IOWAR.NE.0) WRITE(MAINF,2200)
 2200  FORMAT(//,' WARNING: WITH TRADITIONAL ORDINARY COKRIGING THE SUM'
     ;      ,/,' OF THE WEIGHTS APPLIED TO EACH SECONDARY DATA IS ZERO.' 
     ;      ,/,' WITH NMXS_GS SET LOW AND NVAR_GS LARGE THE SECONDARY'
     ;      ,/,' DATA WILL NOT CONTRIBUTE TO THE ESTIMATE',//)
         END IF

* 8) Maximum number of nested structures defining variogram out of range

         IF (MXVGM_GS.LT.1) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'MAXIMUM NUMBER OF VARIOGRAM NESTED STRUCTURES OUT OF RANGE.', 
     ;  NROW,2,IUGEO,2,9.8)

* 9) Number of pilot points or drawing option out of range. No sense if
*    only direct problem is being solved

         IF (IOINV.GT.0) THEN

           IF (NPP_GS.LT.1) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'NUMBER OF PILOT POINTS OUT OF RANGE.', NROW,2,IUGEO,2,9.9)

           IF (IORD_PP_GS.LT.0 .OR. IORD_PP_GS.GT.2) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'OPTION FOR DRAWING PILOT POINTS OUT OF RANGE.', 
     ;  NROW,2,IUGEO,2,9.9)

           IF (IORD_PP_GS.EQ.2) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'OPTION FOR DRAWING PILOT POINTS UNIFORMLY NOT READY.', 
     ;  NROW,2,IUGEO,2,9.9)
         END IF ! IOINV.GT.0

* 10) Number of measurements out of range. No sense if unconditional 
*     simulation is done

         IF (IOPSC_GS.LE.1 .AND. NMEAS_GS.LE.0) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'NUMBER OF MEASUREMENTS OUT OF RANGE.', 
     ;  NROW,2,IUGEO,2,9.10)

* 11) Number of SUPER BLOCKS out of range

         IF (MXSB_GS(1).EQ.0 .OR. MXSB_GS(2).EQ.0 .OR. MXSB_GS(3).EQ.0) 
     ;      CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'CONCEPTUAL ERROR. NUMBER OF SUPER BLOCKS IS ZERO'
     ;  ,NROW,2,IUGEO,2,9.11)

* 12) Anisotropy degree out of range
         
         IF (ISOT.GT.1) CALL ERROR          
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' CURRENT VERSION OF GEOEST. INV. PROB. IS LIMITED TO ISOT=1.',
     ;  NROW,2,IUDIM,0,9.12)

* 13) Number of external drift terms

         IF ((KTYPE_GS.EQ.2 .OR. KTYPE_GS.EQ.3) .AND. NEXDR_GS.LE.0) 
     ;     CALL ERROR          
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' NUMBER OF EXTERNAL DRIFT TERMS OUT OF RANGE',
     ;  NROW,2,IUGEO,2,9.13)

         IF (KTYPE_GS.EQ.3  .AND. NEXDR_GS.GT.3) CALL ERROR          
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' NUMBER OF EXTERNAL DRIFT TERMS OUT OF RANGE.'
     ;//' ACTUALLY LIMITED TO THREE TERMS. SORRY',
     ;  NROW,2,IUGEO,2,9.13)

         IF ((KTYPE_GS.NE.2. AND. KTYPE_GS.NE.3).AND.NEXDR_GS.NE.0) THEN
            CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' KRIGING WITH EXT. ATTRIB. NOT REQUESTED AND NUMBER OF EXT.'
     ;//' DRIFT TERMS IS NOT ZERO. IT IS SET TO ZERO',
     ;  NROW,2,IUGEO,0,9.13)
            NEXDR_GS=0
         END IF

         IF (KTYPE_GS.EQ.2.AND.NEXDR_GS.NE.1) THEN
             CALL ERROR          
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' KRIGING WITH VARYING MEAN REQUIRED AND NUMBER OF EXT.'
     ;//' ATTRIBUTE TERMS IS NOT ONE. IT MUST BE ONE',
     ;  NROW,2,IUGEO,0,9.13)
           NEXDR_GS=1
         END IF

* 14) Option of conditional simulation our of range

         IF (IOPSC_GS.LT.0 .OR. IOPSC_GS.GT.2) CALL ERROR          
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' OPTION OF SIMULATION OUT OF RANGE (<0 OR >2)',
     ;  NROW,2,IUGEO,2,9.14)

*** PROVISIONAL ***

         IF (IOPSC_GS.EQ.2) CALL ERROR          
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' OPTION OF UNCONDITIONAL SIMULATION NOT ENCODED YET',
     ;  NROW,2,IUGEO,2,9.14)

* 15) Seed for conditional simulations is zero or negative

         IF (IOPSC_GS.GT.0.AND.ISEEDSC_GS.LE.0) CALL ERROR 
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' SEED FOR CONDITIONAL SIMULATIONS OUT OF RANGE (<=0)',
     ;  NROW,2,IUGEO,2,9.15)

* 16) Number of conditional simulations
           
         IF (IOPSC_GS.GT.0.AND.MX_SC_GS.LE.0) CALL ERROR 
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' NUMBER OF CONDITIONAL SIMULATIONS OUT OF RANGE',
     ;  NROW,2,IUGEO,2,9.16)

         IF (IOPSC_GS.EQ.0) MX_SC_GS=0

* 17) Option for initial value of PARC out of range

         IF (IO_PARC_INI_GS.LT.0 .OR. IO_PARC_INI_GS.GT.3) CALL ERROR 
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' OPTION FOR INITIAL VALUES AT PILOT POINTS OUT OF RANGE',
     ;  NROW,2,IUGEO,2,9.16)

         IF (IO_PARC_INI_GS.NE.0 .AND. IORD_PP_GS.NE.0) CALL ERROR 
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' OPTION FOR INITIAL VALUES AT PILOT POINTS AND'
     ;//' OPTION FOR RANDOM DRAWING OF PILOT POINTS NOT COHERENT',
     ;  NROW,2,IUGEO,2,9.16)


      ELSE IF (IO_GREST_GS.EQ.2) THEN 

        LEAUX=LEEL(FILENAME,16,MAINF,NROW,INPWR)
        READ(LEAUX,1000) NLINCMB             ! Number of terms defining comb.

       ! Checks that number of terms is >1. If it is one, estimation option
       ! should be 0

        IF (NLINCMB.LT.0) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' NUMBER OF TERMS DEFINING A LINEAR COMBINATION OUT OF ORDER'//
     ;  ' (<=0)',
     ;  NROW,2,IUDIM,0,9.12)

        IF (NLINCMB.EQ.1) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  ' NUMBER OF TERMS DEFINING A LINEAR COMBINATION OUT OF ORDER'//
     ;  '. SHOULD BE >1. CHANGE ESTIMATION OPTION OF THIS GROUP ',
     ;  NROW,2,IUDIM,0,9.12)

       ! Updates max. number of parameters in a linear comb. (used for dimens.)

         IF (NLINCMB.GT.MXLINCMB) MXLINCMB=NLINCMB

      END IF ! IO_GREST_GS.EQ.0
          
* 17) Checks if number of group of zones exceeds maximum allowed

      IF (NGROUP_ZN.GT.MXGRPZN-1) THEN
          WRITE(MAINF,2300)
          WRITE(6,2400)
 2300     FORMAT(//,' NUMBER OF GROUPS OF ZONES EXCEEDED'
     ;           //'. CONTACT TECHNICAL SERVICE. SORRY')
          STOP
      END IF

C_______________________ Step 4: Echoes information if so desired

      IF (INPWR.NE.0) THEN
        IF (IO_GREST_GS.EQ.0) THEN  ! PARAMETER ESTIMATED DETERMINIST.
          WRITE(MAINF,2601) CHARTYPE(ID_TP_GS),IO_GREST_GS
 2601     FORMAT(//,10X,' GENERAL OPTIONS',/,10X,' ======= =======',//,
     ; 5X,'TYPE OF PARAMETER................................... =',4X
     ;    ,A3,/,
     ; 5X,'ESTIMATION OPTION................................... =',I5)

        ELSE IF (IO_GREST_GS.EQ.1) THEN  ! PARAMETER ESTIMATED GEOEST.
          WRITE(MAINF,2600) CHARTYPE(ID_TP_GS),IO_GREST_GS,NVAR_GS
     ;                     ,NMEAS_GS,NPP_GS,IORD_PP_GS,NZON_GS,NPPX_GS
     ;                     ,NPPY_GS,NPPZ_GS,IPROBGR_GS,IO_PARC_INI_GS

 2600     FORMAT(//,10X,' GENERAL OPTIONS',/,10X,' ======= =======',//,
     ; 5X,'TYPE OF PARAMETER................................... =',4X
     ;    ,A3,/,
     ; 5X,'ESTIMATION OPTION................................... =',I5,/,
     ; 5X,'NUMBER OF VARIABLES (PRIM+ALL SEC.)................. =',I5,/,
     ; 5X,'TOTAL NUMBER OF SAMPLING LOCATIONS.................. =',I5,/,
     ; 5X,'NUMBER OF PILOT POINTS.............................. =',I5,/,
     ; 5X,'OPTION FOR RANDOM DRAWING OF PILOT POINTS........... =',I5,/,
     ; 5X,'NUMBER OF ZONES..................................... =',I5,/,
     ; 5X,'DENSITY OF PILOT POINTS X-DIRECTION................. =',I5,/,
     ; 5X,'DENSITY OF PILOT POINTS Y-DIRECTION................. =',I5,/,
     ; 5X,'DENSITY OF PILOT POINTS Z-DIRECTION................. =',I5,/,
     ; 5X,'FLOW-TRANSPORT PROBLEM.............................. =',I5,/,
     , 5X,'OPTION FOR INITIAL VALUE AT PILOT POINTS............ =',I5)



          WRITE(MAINF,2700) KTYPE_GS,NDMIN_GS,NMXP_GS
     ;                     ,NMXS_GS,NOCT_GS
     ;                     ,(MXSB_GS(I),I=1,3)
     ;                     ,(NDISC_GS(I),I=1,3),MXVGM_GS
     ;                     ,NEXDR_GS,IOPSC_GS,ISEEDSC_GS
     ;                     ,MX_SC_GS

 2700     FORMAT(//,10X,' KRIGING OPTIONS',/,10X,' ======= =======',//,
     ; 5X,'KRIGING TYPE........................................ =',I5,/,
     ; 5X,'MINIMUM NUMBER OF SAMPLES FOR KRIGING A BLOCK....... =',I5,/,
     ; 5X,'MAX. NUM. OF (PRIM.+PI.POINTS) FOR KRIGING A BLOCK.. =',I5,/,
     ; 5X,'MAX. NUM. OF SEC. SAMPLES FOR KRIGING A BLOCK....... =',I5,/,
     ; 5X,'NUMBER OF SAMPLES RETAINED PER SEARCH OCTANT........ =',I5,/,
     ; 5X,'NUMBER OF SUPER BLOCKS (COARSE SEARCH GRID).(X-AXIS) =',I5,/,
     ; 5X,'NUMBER OF SUPER BLOCKS (COARSE SEARCH GRID).(Y-AXIS) =',I5,/,
     ; 5X,'NUMBER OF SUPER BLOCKS (COARSE SEARCH GRID).(Z-AXIS) =',I5,/,
     ; 5X,'NUMBER OF BLOCK DISCRETIZATION POINTS (X AXIS)...... =',I5,/,
     ; 5X,'NUMBER OF BLOCK DISCRETIZATION POINTS (Y AXIS)...... =',I5,/,
     ; 5X,'NUMBER OF BLOCK DISCRETIZATION POINTS (Z AXIS)...... =',I5,/,
     ; 5X,'MAX. NUM. OF VARIOGRAM NESTED STRUCTURES............ =',I5,/,
     ; 5X,'NUMBER OF EXTERNAL ATTRIBUTE TERMS.................. =',I5,/,
     ; 5X,'OPTION OF CONDITIONAL SIMULATION.................... =',I5,/,
     ; 5X,'SEED FOR CONDITIONAL SIMULATIONS.................... =',I10,/
     ; 5X,'NUMBER OF CONDITIONAL SIMULATIONS................... =',I5)

        ELSE IF (IO_GREST_GS.EQ.2) THEN   ! Linear combination of unknowns
 
          WRITE (MAINF,2800) CHARTYPE2(ID_TP_GS),NLINCMB
 2800     FORMAT(
     ; 5X,'TYPE OF PARAMETER................................... =',4X
     ;    ,A3,/,5X,'NUMBER OF TERMS DEFINING LINEAR'
     ;              ' COMBINATION......... =',I5)

        END IF ! Type of estimation
      END IF ! INPWR.NE.0

C_______________________ Step 5: Stores read options in geostatistical arrays

      IOPT_GS(NFM,1)=ID_TP_GS
      IOPT_GS(NFM,2)=IO_GREST_GS

      IF (IO_GREST_GS.EQ.0) THEN       ! Deterministic. Continues
        CONTINUE
      ELSE IF (IO_GREST_GS.EQ.1) THEN  ! Only if pilot points are used

        IOPT_GS(NFM,3)=NVAR_GS
        IOPT_GS(NFM,4)=NMEAS_GS
        IOPT_GS(NFM,5)=NPP_GS
        IOPT_GS(NFM,6)=IORD_PP_GS
        IOPT_GS(NFM,7)=NZON_GS
        IOPT_GS(NFM,8)=NPPX_GS
        IOPT_GS(NFM,9)=NPPY_GS
        IOPT_GS(NFM,10)=NPPZ_GS
        IOPT_GS(NFM,12)=IPROBGR_GS
        IOPT_GS(NFM,13)=IO_PARC_INI_GS
        IOPT_GS(NFM,14)=IO_CHECK_POS_GS

        IO_KG_GS(NFM,1)=KTYPE_GS
        IO_KG_GS(NFM,2)=NDMIN_GS
        IO_KG_GS(NFM,3)=NMXP_GS
        IO_KG_GS(NFM,4)=NMXS_GS
        IO_KG_GS(NFM,5)=NOCT_GS
        IO_KG_GS(NFM,6)=MXSB_GS(1)
        IO_KG_GS(NFM,7)=MXSB_GS(2)
        IO_KG_GS(NFM,8)=MXSB_GS(3)
        IO_KG_GS(NFM,9)=NDISC_GS(1)
        IO_KG_GS(NFM,10)=NDISC_GS(2)
        IO_KG_GS(NFM,11)=NDISC_GS(3)
        IO_KG_GS(NFM,12)=MXVGM_GS
        IO_KG_GS(NFM,13)=NEXDR_GS
        IO_KG_GS(NFM,14)=IOPSC_GS
        IO_KG_GS(NFM,15)=ISEEDSC_GS
        IO_KG_GS(NFM,16)=MX_SC_GS


      ELSE ! Parameter is a linear combination of zonal paramaters

         IOPT_GS(NFM,11)=NLINCMB

      END IF ! IOPT_GS(NFM,2).EQ.1

C_______________________ Step 6: Updates number of group of zones read 

      NGROUP_ZN=NGROUP_ZN+1
      GOTO 10                                                ! Reads IFLAG again

      RETURN
      END

