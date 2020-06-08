      SUBROUTINE GEO_CALC
     ;(IDIMCROSS_GS   ,IDIMIVARIO_GS ,IDIMVAR_GS
     ;,IDIMWGT        ,IDIMWORK      ,IDIMZONPP_GS    ,IERROR          
     ;,IOINV          ,LMXNDL          ,MAINF
     ;,MXCLOSE_GS     ,MXDISC_GS     ,MXGRPZN         ,MXKRIG_GS       
     ;,MXMEASPP_GS    ,MXNPP_GS      ,MXNPRIM_GS      ,MXNVAR_GS       
     ;,MXNZON_GS      ,MXROT_GS      ,MXSAM_GS        ,MXSB_GS         
     ;,MXZONPP_GS     ,NFLAGS        ,NGROUP_ZN       ,NPAR            
     ;,NPARDET        ,NPAREL        ,NTYPAR          ,NUMEL           
     ;,NUMNP          ,NUMITER       ,NZPAR           ,NWRITE          
     ;,AREA           ,CLOSESAM_GS   ,COORD           ,COORDGR_GS     
     ;,COVPAR         ,CROSSCOV_GS   ,DATASC_GS       ,ESTKRIG_GS      
     ;,EXDRZN_GS      ,FILENAME      ,ICHECK_GS       ,ICROSSCOV_GS    
     ;,IFLAGS         ,INDPAR        ,INORPAR         ,IO_KG_GS        
     ;,IOLG_PAR       ,IOPT_GS       ,IOWRITE         ,IPNT_PAR        
     ;,IPOLDRIFT_GS   ,ISOZ          ,ISUPBL_GS       ,IVARIO_GS       
     ;,IVPAR          ,IZN_NPP_GS    ,IZN_PP_GS       ,KRIGAUX_GS      
     ;,KRISOL_GS      ,KRISYS_GS     ,KXX             ,LDIM_GS         
     ;,NUMSB_GS       ,NZONE_PAR     ,PAR_WGT         ,PARC            
     ;,PARM           ,PARZ          ,POSDIS_GS       ,POSDISAUX_GS    
     ;,POSMEAS_GS     ,POSZN_GS      ,ROTMAT_GS       ,SEARCH_GS       
     ;,SUPBL_GS       ,TRIM_GS       ,VARIO_GS        ,VMEAS_GS        
     ;,VSTATS_GS      ,WGT_PAR       ,WGT_UNK         ,WORK            
     ;,ZNWGT_GS       ,IACTSIMU      ,IDIMDATASC_GS   ,PARC_GS
     ;,IFLAG_SIMUL    ,COVPAR_GR)

********************************************************************************
*
* PURPOSE
*
* This routine performs all initializations and calculations related to
* geostatistics.
*
* DESCRIPTION
*
* - Step 0: Declaration of variables
* - Step 1: Recovers covar. matrix from its inverse. May be it will change due 
*           to random pilot points 
* - Step 2: For all groups of zones, draws randomly pilot points positions and 
*           assigns external drifts at those positions. Also, "measurements" 
*           are assigned
* - Step 3: Initializes geostatistical arrays
* - Step 4: Conditional estimation/simulation. OUTPUTS
*           1) Pilot points covariance matrix (COVPAR)
*           2) Conditional expectation value at pilot point locations. (PARM)
*           3) Kriging weights defining value of zonal parameters at groups 
*              of zones defined geostatistically (WGT_PAR)
*           4) Initial guess of zonal parameters at groups of zones defined 
*              geostatistically (conditional est. / sim.). (PARZ)
* - Step 5: Corrects array PARZ if estimation is performed logarithmically
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CLOSESAM_GS            Auxiliar array (GSLIB requirement) to store closest
*                         sampling points to an estimation point 
*  COORD                  Nodal coordinates                                     
*  COORDGR_GS             Maximum and minimum coodinates of all groups of zones
*  COVPAR                 A priori covariance matrix of inverse problem unknowns
*  COVPAR_GR              Block of COVPAR related to a given group of zones
*  CROSSCOV_GS            Used to store cross covariance matrix of pilot points
*                         Will be used in the calculation of covariance matrix
*  DATASC_GS              Used to store data for conditional simulation (updated
*                         sequentially)
*  ESTKRIG_GS             Array used to store estimated values through kriging
*  EXDRZN_GS              Array that stores values of the external drifts zonal
*                         values
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  ICHECK_GS              Array used to check if all pilot points have been used 
*  ICROSSCOV_GS           Used to store variable type of samples parameterizing a 
*                         zone. Will be used in the calc. of cov. matrix (COKRIGING)
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INDPAR                 Array of 0's and 1's. 0 means that the parameter is   
*                         estimated aritmethically, otherwise logarithmically.  
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARZ, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
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
*  IOLG_PAR               Array containing all logarithmic options of           
*                         estimated parameters                                  
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
*                         - Column 11; Number of terms defining linear combination
*                         - Column 12. Flow/transport problem
*  IOWRITE                Array containing all output options                   
*  IPNT_PAR               Array of pointer of components of DLT_PAR used in the 
*                         parameterization of a given zonal value
*  IPOLDRIFT_GS           Array containing polynomial drift options (residual kriging)
*                         for all groups of zones
*  ISOZ                   Anisotropy of every transmissivity zone               
*  ISUPBL_GS              Array used (GSLIB requirements) for super block search
*  IVARIO_GS              Array containing integers defining all variograms
*  IVPAR                  Vector containing estimation index for all            
*                         parameters  
*  IZN_NPP_GS             Used to store the number of pilot points used to param.
*                         a zonal value
*  IZN_PP_GS              Used to store which pilot points are used to param. a zonal 
*                         value
*  KRIGAUX_GS             GSLIB requirements. Workspace array
*  KRISOL_GS              GSLIB requirements. Workspace array
*  KRISYS_GS              GSLIB requirements. Workspace array
*  KXX                    Node numbers of every element (counterclockwise order)
*  LDIM_GS                Array containing dimension of transmisivity zones
*  NUMSB_GS               GSLIB req. for superblock search
*  NZONE_PAR              Array containing the number of zones of all parameters
*  PAR_WGT                Array containing objective function weights for       
*                         all estimated parameters                              
*  PARC                   Vector containing calculated values for all parameters
*  PARM                   Vector containing measured/expected values for all 
*                         parameters                                            
*  PARZ                   Vector containing calculated values for all parameters                                            
*  POSDIS_GS              Array containing offsets of the discretization points (wrt. 
*                         center of gravity of all zones)
*  POSDISAUX_GS           GSLIB requirements. Workspace array
*  POSMEAS_GS             Array containing all measurement locations
*  POSZN_GS               Array containing positions of center of gravity of zones
*  ROTMAT_GS              GSLIB requirements. Workspace array for rotation matrices
*  SEARCH_GS              GSLIB requirements. Workspace array for locations search
*  SUPBL_GS               Array containing position of center of gravity and size 
*                         of the first super block
*  TRIM_GS                Array containing trimming limmits
*  VARIO_GS               Array contaning variogram data (sill, nugget, etc)
*  VMEAS_GS               Array containing measurements values
*  VSTATS_GS              Array containing statistics of all variables
*  WGT_PAR                Array containing weights defining linear combinations
*  WGT_UNK                Array containing parameters objectivefunctions weights 
*                         of problems unknowns
*  WORK                   Workspace array
*  ZNWGT_GS               Array used to store temporarily the kriging weights of the
*                         pilot points used in the parameterization of zonal values
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCROSS_GS           Used to dimension CROSSCOV_GS
*  IDIMIVARIO_GS          Used to dimension IVARIO_GS, VARIO_GS
*  IDIMVAR_GS             NVAR_GS*NVAR_GS+NEXDR_GS
*  IDIMWGT                Used to dimension WGT_PAR, IPNT_PAR
*  IDIMZONPP_GS           Maximum among maximum number of pilot points and 
*                         maximum number of zones
*  IERROR                 Current number of errors on input data                
*  IOINV                  Inverse problem option
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MXCLOSE_GS             Maximum number of closest samples allowed
*  MXDISC_GS              Maximum number of discretization points of a given zone
*  MXGRPZN                Maximum number of groups of zones
*  MXKRIG_GS              Used to dimension kriging systems arrays
*  MXMEASPP_GS            Maximum among maximum number of meas. loc. and max. number
*                         of pilot points
*  MXNPP_GS               Maximum number of pilot points
*  MXNPRIM_GS             Maximum number of primary var+pilot point locations used in
*                         the parameterization of a zonal value
*  MXNVAR_GS              Maximum number of variables in the most complex kriging 
*  MXNZON_GS              Maximum number of zones 
*  MXROT_GS               Used to dimension ROTMAT_GS
*  MXSAM_GS               Used to dimension KRIGAUX_GS           
*  MXSB_GS                Maximum number of superblocks 
*  MXZONPP_GS             Maximum among maximum number of zones and max. number of
*                         pilot points
*  NFLAGS                 Used to dimension IFLAGS
*  NGROUP_ZN              Number of groups of zones
*  NPAR                   Number of parameters to be estimated
*  NPARDET                Number of parameters to be estimated deterministically
*  NPAREL                 Number of element parameters in current problem       
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                        
*  NUMITER                Actual inverse problem iteration
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  NWRITE                 Number of output options (used for dimensioning)      
*
* INTERNAL VARIABLES: SCALARS
*
*  D1                     LINV1P requirements
*  D2                     LINV1P requirements
*  IER                    LINV1P requirements
*  IGR                    Actual group of zones
*  IPIPO                  Dummy counter of pilot points
*  IVAR                   Dummy counter of variables
*  IZPAR                  Dummy counter of zonal parameters
*  IGROUP                 Group to which a zone belongs to
*  IFIRSTPOS              First useful position. Initialization of arrays
*  ISVARIABLE             Boolean flagging variable positions of the pilot points 
*                         at any of the groups
*  ISTOTPIPO              Dummy counter of pilot points
*  NPP_GR                 Number of pilot points of a given group
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ASS_EXT_DRIFT
*  COND_EST_SIM
*  CORRECT_PARZ
*  EQUAL_ARRAY
*  IO_SUB
*  LINV1P
*  RANDOM_PIPO
*  ZERO_ARRAY
*  ZERO_ARRAY_I
*
* HISTORY: AAR. First coding (Sept-2002)
*
********************************************************************************

C___________________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 MAINF,NPAR,NUMITER,NGROUP_ZN,IOINV,MXGRPZN,LMXNDL,NPAREL
     ;         ,MXMEASPP_GS,NFLAGS,NTYPAR,NUMEL,NUMNP,NZPAR,IDIMVAR_GS
     ;         ,IERROR,NWRITE,MXNZON_GS,MXNVAR_GS,IDIMWGT,NPARDET,ISOZ
     ;         ,IDIMCROSS_GS,IDIMIVARIO_GS,IDIMZONPP_GS
     ;         ,MXCLOSE_GS,MXDISC_GS,MXKRIG_GS,MXNPP_GS,LDIM_GS
     ;         ,MXNPRIM_GS,MXROT_GS,MXSAM_GS,MXSB_GS,MXZONPP_GS,IDIMWORK
     ;         ,IACTSIMU,IDIMDATASC_GS
     ;         ,IOPT_GS(MXGRPZN,20),IO_KG_GS(NGROUP_ZN,16)
     ;         ,IVPAR(NZPAR,4),INORPAR(NTYPAR),KXX(LMXNDL,NUMEL)
     ;         ,LTYPE(NUMEL),LXPAREL(NUMEL,NPAREL),IOWRITE(NWRITE)
     ;         ,LNNDEL(NUMEL),IPNT_PAR(NZPAR*IDIMWGT),NZONE_PAR(NTYPAR)
     ;         ,ICHECK_GS(MXNPP_GS),ICROSSCOV_GS(MXNPRIM_GS,MXZONPP_GS)
     ;         ,INDPAR(NPAR),IOLG_PAR(NTYPAR),IPOLDRIFT_GS(9,NGROUP_ZN)
     ;         ,ISUPBL_GS(MXSB_GS*8,3),IZN_NPP_GS(MXNZON_GS)
     ;         ,IVARIO_GS(IDIMIVARIO_GS,2,NGROUP_ZN),IFLAGS(NFLAGS)
     ;         ,IZN_PP_GS(IDIMZONPP_GS),NUMSB_GS(MXSB_GS)
     ;         ,IFLAG_SIMUL(MXZONPP_GS)
                                                                 ! Real external
      REAL*8 WORK(IDIMWORK),COVPAR(NPAR*(NPAR+1)/2),AREA(NUMEL)
     ;      ,COORD(NUMNP,3),COORDGR_GS(6,NGROUP_ZN),PARM(NPAR)
     ;      ,POSMEAS_GS(MXMEASPP_GS,3,2,NGROUP_ZN),PARC(NPAR)
     ;      ,EXDRZN_GS(MXNZON_GS,4,NGROUP_ZN),WGT_PAR(NZPAR*IDIMWGT)
     ;      ,POSZN_GS(MXNZON_GS,3,NGROUP_ZN),TRIM_GS(8,NGROUP_ZN)
     ;      ,VMEAS_GS(MXMEASPP_GS,IDIMVAR_GS,2,NGROUP_ZN)
     ;      ,CLOSESAM_GS(MXCLOSE_GS,2),DATASC_GS(IDIMDATASC_GS,18)
     ;      ,CROSSCOV_GS(IDIMCROSS_GS,MXZONPP_GS,2)
     ;      ,KRIGAUX_GS(MXSAM_GS,8),KRISOL_GS(MXKRIG_GS,3)
     ;      ,KRISYS_GS(MXKRIG_GS*MXKRIG_GS),PARZ(NZPAR)
     ;      ,POSDIS_GS(MXDISC_GS,3,MXNZON_GS,NGROUP_ZN)
     ;      ,POSDISAUX_GS(MXDISC_GS,3),ROTMAT_GS(MXROT_GS,3,3)
     ;      ,PAR_WGT(NTYPAR),WGT_UNK(NPAR),SEARCH_GS(11,NGROUP_ZN)
     ;      ,SUPBL_GS(6,NGROUP_ZN),VSTATS_GS(MXNVAR_GS,4,NGROUP_ZN)
     ;      ,VARIO_GS(IDIMIVARIO_GS,8,NGROUP_ZN)
     ;      ,ZNWGT_GS(MXNPRIM_GS*MXNZON_GS),ESTKRIG_GS(MXZONPP_GS)
     ;      ,PARC_GS(MXNPP_GS,NGROUP_ZN)
     ;      ,COVPAR_GR(MXNPP_GS*(MXNPP_GS+1)/2)

                                                              ! Integer internal
      INTEGER*4 IER,IGR,IPIPO,IVAR,IZPAR,IGROUP,IFIRSTPOS,NPP_GR
     ;         ,ITOTPIPO,PIPOVARIABLE
                                                                 ! Real internal
      REAL*8 D1,D2
      CHARACTER FILENAME(20)*20

      IF (IFLAGS(3).NE.0) CALL IO_SUB('GEO_CALC',0)

C___________________________ Step 1: Recovers covar. matrix from its inverse. 
C___________________________         May be it will change due to random pilot 
C___________________________         points. Only for 2nd or more iteration 
C___________________________        (if pilot points vary)or at first iteration 
C___________________________        for 2nd or more simulations

*** OJO, DEBERIA REINVERTIRSE POR BLOQUES TAMBIEN... DE MOMENTO LO DEJO ASI

      IF ( (NUMITER.GT.1 .AND. PIPOVARIABLE.NE.0) .OR. 
     ;     (IACTSIMU.GT.1 .AND. NUMITER.EQ.1) ) THEN
         CALL EQUAL_ARRAY (WORK,COVPAR,NPAR*(NPAR+1)/2)
         CALL LINV1P (WORK,NPAR,COVPAR,D1,D2,IER)
         IF (IER.EQ.129) THEN
            WRITE(MAINF,2200)
 2200       FORMAT(//,' ERROR INVERTING A PRIORI COVARIANCE MATRIX OF'
     ;             ' DETERMINISTICALLY ESTIMATED PARAMETERS.',/,' IT IS'
     ;             ' NOT POSITIVE DEFINITE. FORCED STOP, SORRY',/)
            STOP ' CRITICAL STOP. CHECK FILE RES.OUT'
         END IF
      END IF

*** FIN DEL OJO

C___________________________ Step 2: For all groups of zones, draws randomly 
C___________________________         pilot points positions and assigns
C___________________________         external drifts at those positions. 
C___________________________         Also, "measurements" are assigned

      DO IGR=1,NGROUP_ZN

         IF (IOINV.GT.0. AND. IOPT_GS(IGR,6).GT.0) THEN

            CALL RANDOM_PIPO
     ;(IOPT_GS(IGR,1)    ,IGR         ,IOPT_GS(IGR,6)    ,LMXNDL    
     ;,MAINF             ,MXMEASPP_GS ,NFLAGS            ,NPAREL        
     ;,IOPT_GS(IGR,5)    ,NTYPAR      ,NUMEL             ,NUMNP                
     ;,NZPAR             ,AREA        ,COORD       ,COORDGR_GS(1,IGR)
     ;,IFLAGS            ,IVPAR(1,3)  ,INORPAR           ,KXX          
     ;,LTYPE             ,LXPAREL     ,POSMEAS_GS(1,1,1,IGR))


            CALL ASS_EXT_DRIFT
     ;(IOPT_GS(1,IGR)    ,2          ,IDIMVAR_GS         ,IERROR
     ;,IGR               ,IOWRITE(2) ,0                  ,16
     ;,LMXNDL            ,MAINF      ,MXMEASPP_GS        ,MXNZON_GS
     ;,MXNVAR_GS         ,IO_KG_GS(IGR,13)
     ;,NPAREL            ,IOPT_GS(IGR,5)                 ,NTYPAR
     ;,NUMEL             ,NUMNP      ,IOPT_GS(IGR,7)     ,NZPAR
     ;,AREA              ,COORD      ,EXDRZN_GS(1,1,IGR) ,FILENAME
     ;,IVPAR(1,3)        ,INORPAR    ,KXX
     ;,LNNDEL            ,LTYPE      ,LXPAREL    ,POSMEAS_GS(1,1,1,IGR)
     ;,POSZN_GS(1,1,IGR) ,VMEAS_GS(1,1,1,IGR))            
*     ;,IZONMEAS_GS(1,1))

            DO IPIPO=1,IOPT_GS(IGR,5)
               VMEAS_GS(IPIPO,1,1,IGR)=-1D0*DFLOAT(IPIPO)*1D60
               DO IVAR=2,IOPT_GS(IGR,3)
                  VMEAS_GS(IPIPO,IVAR,1,IGR)=1D3*TRIM_GS(IVAR,IGR)
               END DO
            END DO

         END IF ! IOPT_GS(IGR,6).GT.0

      END DO ! IGR=1,NGROUP_ZN

C___________________________ Step 3: Initializes geostatistical arrays

      IF (IOINV.GT.0) THEN

         DO IZPAR=1,NZPAR

            IGROUP=IVPAR(IZPAR,3)
            IF (IGROUP .LE. 0) CYCLE
            IF (IGROUP.LE.NGROUP_ZN .AND. IOPT_GS(IGROUP,1).EQ.1) THEN
                   
               IFIRSTPOS=(IZPAR-1)*IDIMWGT+1
               IF (IOPT_GS(IGROUP,6).NE.0) THEN   ! Variable pilot points
                 CALL ZERO_ARRAY_I (IPNT_PAR(IFIRSTPOS),IDIMWGT)
                 CALL ZERO_ARRAY (WGT_PAR(IFIRSTPOS),IDIMWGT)
               END IF

            END IF
         
         END DO

         ITOTPIPO=0
         DO IGR=1,NGROUP_ZN

           IF (IOPT_GS(IGR,6).NE.0) THEN   ! Variable pilot points

              NPP_GR=IOPT_GS(IGR,5)
              IF (IGR.GT.1) ITOTPIPO=ITOTPIPO+NPP_GR

              IFIRSTPOS=NPARDET*(NPARDET+1)/2+ITOTPIPO+1
              CALL ZERO_ARRAY (COVPAR(IFIRSTPOS),NPP_GR*(NPP_GR+1)/2)

              IFIRSTPOS=NPARDET+ITOTPIPO+1
              CALL ZERO_ARRAY (PARM(IFIRSTPOS),NPP_GR)
              CALL ZERO_ARRAY (PARC(IFIRSTPOS),NPP_GR)

           END IF

         END DO



      END IF

C___________________________ Step 4: Conditional estimation/simulation. OUTPUTS
C___________________________         1) Pilot points covariance matrix (COVPAR)
C___________________________         2) Conditional expectation value at pilot 
C___________________________            point locations. (PARM)
C___________________________         3) Kriging weights defining value of zonal
C___________________________            parameters at groups of zones defined 
C___________________________            geostatistically (WGT_PAR)
C___________________________         4) Initial guess of zonal parameters at 
C___________________________            groups of zones defined geostatistically
C___________________________            (conditional est. / sim.). (PARZ)

      CALL COND_EST_SIM
     ;(IDIMCROSS_GS    ,IDIMIVARIO_GS    ,IDIMVAR_GS
     ;,IDIMWGT         ,IDIMWORK         ,IDIMZONPP_GS   ,IOINV          
     ;,MAINF           ,MXCLOSE_GS       ,MXDISC_GS      ,MXGRPZN        
     ;,MXKRIG_GS       ,MXMEASPP_GS      ,MXNPP_GS       ,MXNPRIM_GS     
     ;,MXNVAR_GS       ,MXNZON_GS        ,MXROT_GS       ,MXSAM_GS       
     ;,MXSB_GS         ,MXZONPP_GS       ,NFLAGS         ,NGROUP_ZN      
     ;,NPAR            ,NPARDET          ,NTYPAR         ,NUMITER        
     ;,NZPAR           ,NZONE_PAR(1)     ,PIPOVARIABLE   ,CLOSESAM_GS
     ;,COVPAR          ,ESTKRIG_GS       ,EXDRZN_GS      
     ;,ICHECK_GS       ,ICROSSCOV_GS     ,IFLAGS         ,IVPAR(1,3)  
     ;,INDPAR          ,INORPAR          ,IO_KG_GS       ,IOLG_PAR       
     ;,IOPT_GS         ,IVPAR(1,2)       ,IPNT_PAR       ,IVPAR(1,1)     
     ;,IPOLDRIFT_GS    ,ISOZ             ,ISUPBL_GS      ,IVARIO_GS      
     ;,IZN_NPP_GS      ,IZN_PP_GS        ,KRIGAUX_GS     ,KRISOL_GS      
     ;,KRISYS_GS       ,LDIM_GS          ,NUMSB_GS       ,NZONE_PAR      
     ;,PAR_WGT         ,PARC             ,PARM           ,PARZ           
     ;,POSDIS_GS       ,POSDISAUX_GS     ,POSMEAS_GS     ,POSZN_GS       
     ;,ROTMAT_GS       ,SEARCH_GS        ,SUPBL_GS       ,TRIM_GS        
     ;,VARIO_GS        ,VMEAS_GS         ,VSTATS_GS      ,WGT_PAR        
     ;,WGT_UNK         ,WORK             ,ZNWGT_GS
     ;,IACTSIMU        ,DATASC_GS        ,IDIMDATASC_GS  ,PARC_GS
     ;,IFLAG_SIMUL     ,COVPAR_GR)


C___________________________ Step 5: Corrects array PARZ if estimation is 
C___________________________         performed logarithmically

      IF (PIPOVARIABLE.EQ.1 .OR. NUMITER.EQ.1) THEN
         CALL CORRECT_PARZ
     ;(MXGRPZN  ,NTYPAR    ,NUMITER   ,NZPAR   ,NZONE_PAR(1)  ,INORPAR
     ;,IOPT_GS  ,ISOZ      ,IVPAR     ,LDIM_GS ,PARZ)

C_______________________ Step 6: Echoes complete array PARZ (GSLIB format)
C_______________________         and PARC

         IF (IFLAGS(23).LT.0) THEN
            IF (IACTSIMU.EQ.1 .AND. NUMITER.EQ.1)
     ;         OPEN(UNIT=999,FILE='INITIAL_FIELDS.DAT',STATUS='UNKNOWN')
            WRITE(666,2701) IACTSIMU
            WRITE(999,2702) IACTSIMU
 2701       FORMAT(//,' VALORES INICIALES DE PARZ EN LA SIMULACION:',I5
     ;             ,/,' ======= ========= == ==== == == ==========',/)
 2702       FORMAT(//,' VALORES INICIALES DE PARZ EN LA SIMULACION:',I5
     ;             ,/,' ======= ========= == ==== == == ==========',/,
     ;                ' 1',/,' PARZ',/)
            WRITE(999,2601) (PARZ(IZPAR),IZPAR=1,NZPAR)
            WRITE(666,2602) (PARZ(IZPAR),IZPAR=1,NZPAR)
 2601       FORMAT(E10.4)
 2602       FORMAT(7E10.4)

         END IF ! IFLAGS...
      END IF

      IF (NUMITER.EQ.1) THEN 
         WRITE(MAINF,998) IACTSIMU
         WRITE(6,999) IACTSIMU
      END IF

 998  FORMAT(//,' RESULTS OF CONDITIONAL SIMULATION / ESTIMATION: '
     ;    ,I5,/,' ======= == =========== ========== = ===========',//)
 999  FORMAT(//,' END OF GEOSTATISTICAL CALCULATIONS. UPDATING'
     ;          ' FIELDS.',/,' SIMULATION / ESTIMATION NUMBER',I5,//)
      IF (IFLAGS(3).NE.0) CALL IO_SUB('GEO_CALC',0)

      RETURN
      END

