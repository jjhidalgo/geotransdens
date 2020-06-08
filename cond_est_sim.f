      SUBROUTINE COND_EST_SIM
     ;(IDIMCROSS_GS    ,IDIMIVARIO_GS    ,IDIMVAR_GS
     ;,IDIMWGT         ,IDIMWORK         ,IDIMZONPP_GS   ,IOINV
     ;,MAINF           ,MXCLOSE_GS       ,MXDISC_GS      ,MXGRPZN
     ;,MXKRIG_GS       ,MXMEASPP_GS      ,MXNPP_GS       ,MXNPRIM_GS
     ;,MXNVAR_GS       ,MXNZON_GS        ,MXROT_GS       ,MXSAM_GS
     ;,MXSB_GS         ,MXZONPP_GS       ,NFLAGS         ,NGROUP_ZN
     ;,NPAR            ,NPARDET          ,NTYPAR         ,NUMITER
     ;,NZPAR           ,NZTRA            ,PIPOVARIABLE   ,CLOSESAM_GS
     ;,COVPAR          ,ESTKRIG_GS       ,EXDRZN_GS      
     ;,ICHECK_GS       ,ICROSSCOV_GS     ,IFLAGS         ,IGR_ZONE
     ;,INDPAR          ,INORPAR          ,IO_KG_GS       ,IOLG_PAR
     ;,IOPT_GS         ,IPNT_END         ,IPNT_PAR       ,IPNT_START
     ;,IPOLDRIFT_GS    ,ISOZ             ,ISUPBL_GS      ,IVARIO_GS
     ;,IZN_NPP_GS      ,IZN_PP_GS        ,KRIGAUX_GS     ,KRISOL_GS
     ;,KRISYS_GS       ,LDIM_GS          ,NUMSB_GS       ,NZONE_PAR
     ;,PAR_WGT         ,PARC             ,PARM           ,PARZ
     ;,POSDIS_GS       ,POSDISAUX_GS     ,POSMEAS_GS     ,POSZN_GS
     ;,ROTMAT_GS       ,SEARCH_GS        ,SUPBL_GS       ,TRIM_GS
     ;,VARIO_GS        ,VMEAS_GS         ,VSTATS_GS      ,WGT_PAR
     ;,WGT_UNK         ,WORK             ,ZNWGT_GS
     ;,IACTSIMUL       ,DATASC_GS        ,IDIMDATASC_GS  ,PARC_GS
     ;,IFLAG_SIMUL     ,COVPAR_GR)

C____________________________ Declaration of variables and DEBUG file

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NFLAGS,NUMITER,MXGRPZN,NGROUP_ZN,MXMEASPP_GS,IDIMVAR_GS
     ;         ,MXZONPP_GS,IDIMCROSS_GS,MXNPP_GS,MXNVAR_GS
     ;         ,MAINF,MXDISC_GS,MXKRIG_GS,IDIMIVARIO_GS,MXNPRIM_GS
     ;         ,MXNZON_GS,MXCLOSE_GS,IDIMZONPP_GS,MXSB_GS,MXROT_GS
     ;         ,MXSAM_GS,NPARDET,NPAR,NTYPAR,IOINV,IDIMWGT,NZPAR,NZTRA
     ;         ,PIPOVARIABLE,IDIMWORK,IACTSIMUL,IDIMDATASC_GS
     ;         ,IFLAGS(NFLAGS),IOPT_GS(MXGRPZN,20),IO_KG_GS(MXGRPZN,16)
     ;         ,ICROSSCOV_GS(MXNPRIM_GS,MXZONPP_GS)
     ;         ,IZN_NPP_GS(MXNZON_GS),IPOLDRIFT_GS(9,NGROUP_ZN)
     ;         ,IVARIO_GS(IDIMIVARIO_GS,2,NGROUP_ZN),INDPAR(NPAR)
     ;         ,ISUPBL_GS(MXSB_GS*8,3),NUMSB_GS(MXSB_GS)
     ;         ,IOLG_PAR(NTYPAR),ICHECK_GS(MXNPP_GS),INORPAR(NTYPAR)
     ;         ,NZONE_PAR(NZPAR),IGR_ZONE(NZPAR),IPNT_START(NZPAR)
     ;         ,IPNT_END(NZPAR),ISOZ(NZTRA),LDIM_GS(NZTRA)
     ;         ,IZN_PP_GS(IDIMZONPP_GS),IPNT_PAR(NZPAR*IDIMWGT)
     ;         ,IFLAG_SIMUL(MXZONPP_GS)
                                                                 ! Real external
      REAL*8 POSMEAS_GS(MXMEASPP_GS,3,2,NGROUP_ZN)
     ;      ,TRIM_GS(8,NGROUP_ZN)
     ;      ,VMEAS_GS(MXMEASPP_GS,IDIMVAR_GS,2,NGROUP_ZN)
     ;      ,ESTKRIG_GS(MXZONPP_GS),DATASC_GS(IDIMDATASC_GS,18)
     ;      ,VSTATS_GS(MXNVAR_GS,4,NGROUP_ZN),SEARCH_GS(11,NGROUP_ZN)
     ;      ,SUPBL_GS(6,NGROUP_ZN),KRISYS_GS(MXKRIG_GS*MXKRIG_GS)
     ;      ,VARIO_GS(IDIMIVARIO_GS,8,NGROUP_ZN),KRISOL_GS(MXKRIG_GS,3)
     ;      ,ZNWGT_GS(MXNPRIM_GS*MXNZON_GS),CLOSESAM_GS(MXCLOSE_GS,2)
     ;      ,ROTMAT_GS(MXROT_GS,3,3),KRIGAUX_GS(MXSAM_GS,8)
     ;      ,COVPAR(NPAR*(NPAR+1)/2),WORK(IDIMWORK)
     ;      ,PAR_WGT(NTYPAR),PARC(NPAR),PARM(NPAR),WGT_UNK(NPAR)
     ;      ,POSDIS_GS(MXDISC_GS,3,MXNZON_GS,NGROUP_ZN)
     ;      ,EXDRZN_GS(MXNZON_GS,4,NGROUP_ZN),POSDISAUX_GS(MXDISC_GS,3)
     ;      ,POSZN_GS(MXNZON_GS,3,NGROUP_ZN),PARC_GS(MXNPP_GS,NGROUP_ZN)
     ;      ,WGT_PAR(NZPAR*IDIMWGT),PARZ(NZPAR)
     ;      ,COVPAR_GR(MXNPP_GS*(MXNPP_GS+1)/2)

      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: CROSSCOV_GS

                                                              ! Integer external
      INTEGER*4 NVAR_GR,NMEAS_GR,NPP_GR,NZN_GR,IO_RDPP_GR,IUSEPP
     ;         ,KTYPE_GR,NDMIN_GR,NMXP_GR,NMXS_GR,NOCT_GR,MXSBX_GR,IER
     ;         ,MXSBY_GR,MXSBZ_GR,MXNST_GR,NEXDR_GR,IGR,NDISC_GR
     ;         ,MXROT_GR,MAXSAM_GR,MXSB_GR,IDEBUG,MAX_1,MAX_2,MXKRIG_GR
     ;         ,NTOTALPP,NZEROS,IDIMCOV_GR,IOPSC_GR,IPIPO,IPOS
     ;         ,IACTTYPE,NONUSEPP,IZON,IDPIPO,IDIM
     ;         ,IPINORPAR,IPNZPAR,ISTART,IPOSDLT_PAR,ISUM,IADD,ICOMPO
     ;         ,NRESTRI
     ;         ,IO_PARC_INI_GR
                                                                 ! Real internal
      REAL*8 R_DUMMY,D1,D2,GAUPROB,AUXPARZ, dMaxCov
                                                                 ! Integer internal
      INTEGER*4 JPOS,ISEED,ILOOP,ILOOP2
     ;         ,IPOSCOVA_GR,IPOSCOVA


      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: daCovParSec

      IF (IFLAGS(21).GT.0.OR.IFLAGS(22).GT.0.OR.IFLAGS(23).GT.0) THEN ! On demand

         IF (IACTSIMUL.EQ.1 .AND. NUMITER.EQ.1) THEN      ! Beginning of the process
            OPEN(UNIT=666,FILE='CE_CS_INFO.OUT',STATUS='UNKNOWN')
         END IF

         IF (NUMITER.EQ.1) THEN
            WRITE(666,1000) IACTSIMUL
 1000       FORMAT(//,' GEOSTATISTICS INFORMATION. CONDITIONAL' 
     ;             ' ESTIMATION / SIMULATION NUMBER:',I5,/
     ;             ' ============= ============ =========== =========='
     ;             ' = ========== =======',/)
         END IF
         WRITE(666,1001) NUMITER
 1001    FORMAT(//,' ITERATION NUMBER: ',I5,/, ' ========= =======')
      END IF

      ALLOCATE(CROSSCOV_GS(IDIMCROSS_GS,MXZONPP_GS,2))
      CROSSCOV_GS = 0.0D0

C____________________________ LOOP over groups of zones, considering only those
C____________________________ estimated geostatistically

      PIPOVARIABLE=0       ! Flag marking variable position of pilot points
      NTOTALPP=0           ! Total number of pilot points
      DO IGR=1,NGROUP_ZN
         IF (IOPT_GS(IGR,2).EQ.1) THEN   ! Estimated geostatistically

         IF (IFLAGS(21)+IFLAGS(22)+IFLAGS(23).GT.0) 
     ;      WRITE(666,1002) IGR
 1002       FORMAT(//,' INFORMATION RELATED TO GROUP OF ZONES :',I5,/,
     ;                ' =========== ======= == ===== == ===== =')
               
C____________________________ Step 1: Identification of useful variables of current group

             NVAR_GR=IOPT_GS(IGR,3)   ! Numer of variables (primary + all second.)
             NMEAS_GR=IOPT_GS(IGR,4)  ! Total number of measurements
             NPP_GR=IOPT_GS(IGR,5)    ! Number of pilot points          
             NZN_GR=IOPT_GS(IGR,7)    ! Number of zones

             IO_RDPP_GR=IOPT_GS(IGR,6)    ! Option for variable location of pilot p.
             IF (IO_RDPP_GR.NE.0) PIPOVARIABLE=1  ! Updates flag if loc. is variable

             IO_PARC_INI_GR=IOPT_GS(IGR,13)  ! Option for initial value of pilot p.

             KTYPE_GR=IO_KG_GS(IGR,1)        ! Kriging type
             NDMIN_GR=IO_KG_GS(IGR,2)    ! Min. num. of samples for kriging a block
             NMXP_GR=IO_KG_GS(IGR,3)     ! Max. num. of primary variable samples 
             NMXS_GR=IO_KG_GS(IGR,4)     ! Max. num. of sec. variable samples 
             NOCT_GR=IO_KG_GS(IGR,5)     ! Octant search

             MXSBX_GR=IO_KG_GS(IGR,6)    ! Num. of superblocks (X-axis)
             MXSBY_GR=IO_KG_GS(IGR,7)    ! Num. of superblocks (Y-axis)
             MXSBZ_GR=IO_KG_GS(IGR,8)    ! Num. of superblocks (Z-axis)
             MXSB_GR=                    ! Number of superblocks
     ;           IO_KG_GS(IGR,6)*IO_KG_GS(IGR,7)*IO_KG_GS(IGR,8)

             MXNST_GR=IO_KG_GS(IGR,12)   ! Num. of variog. nested structures
             NEXDR_GR=IO_KG_GS(IGR,13)   ! Num. of external drift terms

             NDISC_GR=                   ! Num. of discretization points per block
     ;           IO_KG_GS(IGR,9)*IO_KG_GS(IGR,10)*IO_KG_GS(IGR,11)

             MXROT_GR=MXNST_GR*NVAR_GR*NVAR_GR+1  ! Maximum number of rotation mat.
             MAXSAM_GR=NMXP_GR+NMXS_GR   ! Maximum number of samples to search

             IF (IGR.GT.1) NTOTALPP=NTOTALPP+IOPT_GS(IGR-1,5)  ! Number of pilot p.
             IOPSC_GR=IO_KG_GS(IGR,14)    ! Option for conditional simulation
             IACTTYPE=IOPT_GS(IGR,1)      ! Param. type of zones of this group

             IPINORPAR=1                  ! Pointer to INORPAR of that type
             IPNZPAR=1                    ! Pointer to NZONE_PAR of that type
             IF (IACTTYPE.EQ.2) THEN      ! Storage coeff.
                IPINORPAR=7
                IPNZPAR=IACTTYPE
             ELSE IF (IACTTYPE.EQ.3 .OR. IACTTYPE.EQ.4) THEN  ! Areal recharge
                IPINORPAR=8
                IPNZPAR=3
             ELSE IF (IACTTYPE.GT.4) THEN ! Others
                IPINORPAR=IACTTYPE+7
                IPNZPAR=IACTTYPE+2
             END IF
             IPOS=INORPAR(IPINORPAR)          ! Initial position at IGR_ZONE

C____________________________ Step 2: CE / CS of pilot points to measurements. 
C____________________________         OUTPUTS : - Prior information of pilot points
C____________________________                   - Posible Initial value of pilot points. 
C____________________________                     It can be read from GEO file if the 
C____________________________                     "RESTART" pmethodology is used
C____________________________                   - Inverse of cova matrix of pilot points

             IF ( (IOINV.GT.0. AND. NPP_GR.GT.0) ! Inverse problem at this group
     ;      .AND. (IO_RDPP_GR.NE.0               ! Variable position of pilot p.
     ;       .OR.  NUMITER.EQ.1) ) THEN         ! At first iteration, ALWAYS

C____________________________ Step 2.A) CE of pilot points to measurements. ALWAYS 
C____________________________           done, as the weights and cross-cova. matrix 
C____________________________           of pilot points to measurements are required 
C____________________________           for the calculation of the covariance matrix
C____________________________           of parameters. OUTPUTS: weights and cross-cova 
C____________________________           matrix (needed later) and posible initial value
C____________________________           and prior information of pilot points

                IDEBUG=0
                IF (IFLAGS(21).EQ.1) THEN
                   IDEBUG=1
                   WRITE(666,1100)
                END IF
 1100           FORMAT(//,' CONDITIONAL ESTIMATION OF PILOT POINTS'
     ;                    ' TO MEASUREMENTS',/,
     ;                    ' =========== ========== == ===== ======' 
     ;                    ' == ============',/)

                                 ! Pilot point and sampling locations and meas. 
                                 ! values are saved (arrays VMEAS_GS,
                                 ! POSMEAS_GS are changed in GSLIB routines).
                                 ! Main outputs of current step are also init.

                CALL EQUAL_ARRAY (POSMEAS_GS(1,1,2,IGR)
     ;                 ,POSMEAS_GS(1,1,1,IGR),3*MXMEASPP_GS)

                CALL EQUAL_ARRAY (VMEAS_GS(1,1,2,IGR)
     ;                 ,VMEAS_GS(1,1,1,IGR),MXMEASPP_GS*IDIMVAR_GS)

                CALL ZERO_ARRAY (ESTKRIG_GS,MXZONPP_GS)
                CALL ZERO_ARRAY_I (ICROSSCOV_GS,MXNPRIM_GS*MXZONPP_GS)
                CALL ZERO_ARRAY (CROSSCOV_GS,IDIMCROSS_GS*MXNPP_GS*2)

                                 ! Calculates maximum number of kriging / cokriging equations

                MAX_1=NVAR_GR*NMEAS_GR+NVAR_GR   
                MAX_2=NMEAS_GR+10+NEXDR_GR
                MXKRIG_GR=MAX0(MAX_1,MAX_2) ! Max kriging dimensions

                                 ! Performs the conditional estimation

                IF (KTYPE_GR.LE.4) THEN  ! Kriging

                   CALL KRIGING
     ;(IDEBUG                ,IDIMCROSS_GS        ,1                 
     ;,1                     ,1                   ,KTYPE_GR
     ;,MAINF                 ,MXKRIG_GS           ,MXROT_GR             
     ;,MXSB_GR               ,MXDISC_GS           ,NMEAS_GR             
     ;,1                     ,IVARIO_GS(1,1,IGR)  ,NPP_GR               
     ;,NEXDR_GR              ,NMXP_GR             ,NDMIN_GR             
     ;,NOCT_GR               ,NRESTRI             ,MXSBX_GR             
     ;,MXSBY_GR              ,MXSBZ_GR            ,VARIO_GS(1,1,IGR)    
     ;,SEARCH_GS(1,IGR)      ,SEARCH_GS(9,IGR)    ,SEARCH_GS(10,IGR)    
     ;,SEARCH_GS(11,IGR)     ,SEARCH_GS(3,IGR)    ,SEARCH_GS(4,IGR)     
     ;,VSTATS_GS(1,2,IGR)    ,SUPBL_GS(1,IGR)     ,SUPBL_GS(4,IGR)      
     ;,SUPBL_GS(2,IGR)       ,SUPBL_GS(5,IGR)     ,SUPBL_GS(3,IGR)      
     ;,SUPBL_GS(6,IGR)       ,VARIO_GS(1,4,IGR)   ,VARIO_GS(1,5,IGR)    
     ;,VARIO_GS(1,6,IGR)     ,VARIO_GS(1,7,IGR)   ,VARIO_GS(1,8,IGR)    
     ;,CLOSESAM_GS(1,1)      ,KRISOL_GS(1,1)      ,KRISOL_GS(1,2)       
     ;,KRISYS_GS             ,CROSSCOV_GS(1,1,1)   
     ;,VMEAS_GS(NPP_GR+1,1,2,IGR)              
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+1,2,IGR)
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+2,2,IGR)      
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+3,2,IGR)
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+4,2,IGR)         
     ;,IPOLDRIFT_GS(1,IGR)   ,IVARIO_GS(1,2,IGR)  ,ISUPBL_GS(1,1)        
     ;,ISUPBL_GS(1,2)        ,ISUPBL_GS(1,3)      ,NUMSB_GS              
     ;,ICROSSCOV_GS(1,1)     ,0.0D0               ,VARIO_GS(1,2,IGR)     
     ;,ROTMAT_GS             ,VMEAS_GS(1,2,2,IGR) ,VMEAS_GS(1,3,2,IGR)   
     ;,VMEAS_GS(1,4,2,IGR)   ,VMEAS_GS(1,5,2,IGR) ,VARIO_GS(1,3,IGR)     
     ;,CLOSESAM_GS(1,2)      ,TRIM_GS(1,IGR)      ,KRIGAUX_GS(1,4)       
     ;,KRIGAUX_GS(1,5)       ,KRISOL_GS(1,3)      ,CROSSCOV_GS(1,1,2)    
     ;,POSMEAS_GS(NPP_GR+1,1,2,IGR)               ,POSMEAS_GS(1,1,2,IGR) 
     ;,KRIGAUX_GS(1,1)       ,0.0D0                
     ;,POSMEAS_GS(NPP_GR+1,2,2,IGR)               ,POSMEAS_GS(1,2,2,IGR)
     ;,KRIGAUX_GS(1,2)       ,0.0D0               
     ;,POSMEAS_GS(NPP_GR+1,3,2,IGR)               ,POSMEAS_GS(1,3,2,IGR)
     ;,KRIGAUX_GS(1,3)       ,0.0D0               ,ESTKRIG_GS
     ;,POSMEAS_GS(NPP_GR+1,1,1,IGR) ,POSMEAS_GS(NPP_GR+1,2,1,IGR) 
     ;,POSMEAS_GS(NPP_GR+1,3,1,IGR) ,MXNPRIM_GS   ,MXMEASPP_GS
     ;,R_DUMMY)

                ELSE ! Cokriging
                   CONTINUE   ! Not yet
                END IF ! Kriging / Cokriging

C____________________________ Step 2.B) Calculation of CE covariance matrix. ALWAYS done. 
C____________________________           If algorithm is in CE mode, this will  be the 
C____________________________           covariance matrix of parameters. In CS mode, 
C____________________________           it will be corrected later. At this step, only the block
C____________________________           of CE cova. matrix related to current group of zones
C____________________________           is calculated

                IDIMCOV_GR=NPP_GR*(NPP_GR+1)/2
                CALL ZERO_ARRAY(WORK,NPAR*(NPAR+1)/2)
                CALL ZERO_ARRAY(COVPAR_GR,MXNPP_GS*(MXNPP_GS+1)/2)

                IDEBUG = 0
                CALL COVA_KRIG_MATRIX
     ;(IDIMCOV_GR         ,IDIMCROSS_GS ,1           ,IDEBUG
     ;,MXROT_GR           ,MXDISC_GS    ,1           ,IVARIO_GS(1,1,IGR) 
     ;,NPP_GR             ,NRESTRI      ,VARIO_GS(1,1,IGR)
     ;,0                  ,WORK         ,COVPAR_GR
     ;,CROSSCOV_GS(1,1,1) ,IVARIO_GS(1,2,IGR)        
     ;,WORK(NPP_GR*(NPP_GR+1)/2+1)      ,0.0D0       ,VARIO_GS(1,2,IGR)         
     ;,ROTMAT_GS          ,VARIO_GS(1,3,IGR)         ,CROSSCOV_GS(1,1,2) 
     ;,POSMEAS_GS(1,1,1,IGR)            ,POSMEAS_GS(1,2,1,IGR)            
     ;,POSMEAS_GS(1,3,1,IGR))

                IF (IOPSC_GR.NE.0) THEN

C____________________________ Step 2.C) CS of pilot points to measurements on demand.
C____________________________           OUTPUTS : - Prior information of pilot points
C____________________________                     - Posible initial value of pilot points
C____________________________                     - Correction of the kriging matrix = COVPAR

                   IDEBUG=0
                   IF (IFLAGS(22).EQ.1) THEN
                      IDEBUG=1
                      WRITE(666,1101)
                   END IF
 1101              FORMAT(//,' CONDITIONAL SIMULATION OF PILOT POINTS'
     ;                       ' TO MEASUREMENTS',/,
     ;                       ' =========== ========== == ===== ======' 
     ;                       ' == ============',/)

                                 ! Pilot point and sampling locations and meas. 
                                 ! values are saved (arrays VMEAS_GS,
                                 ! POSMEAS_GS are changed in GSLIB routines).
                                 ! Main outputs of current step are also init.

                   CALL EQUAL_ARRAY (POSMEAS_GS(1,1,2,IGR)
     ;                     ,POSMEAS_GS(1,1,1,IGR),3*MXMEASPP_GS)
 
                   CALL EQUAL_ARRAY (VMEAS_GS(1,1,2,IGR)
     ;                     ,VMEAS_GS(1,1,1,IGR),MXMEASPP_GS*IDIMVAR_GS)

                   CALL ZERO_ARRAY (ESTKRIG_GS,MXZONPP_GS)
                   CALL ZERO_ARRAY_I(ICROSSCOV_GS,MXNPRIM_GS*MXZONPP_GS)
                   CALL ZERO_ARRAY (CROSSCOV_GS,IDIMCROSS_GS*MXNPP_GS*2)

                                 ! Calculates maximum number of kriging / cokriging equations
                   MAX_1=NVAR_GR*(NMEAS_GR+NPP_GR)+NVAR_GR   
                   MAX_2=NMEAS_GR+NPP_GR+10+NEXDR_GR
                   MXKRIG_GR=MAX0(MAX_1,MAX_2) ! Max kriging dimensions

                                 ! Performs the conditional simulation

                   CALL CONDITIONAL_SIMULATION
     ;(NPP_GR        ,NMEAS_GR     ,MAINF
     ;,MXNPRIM_GS   ,IDIMCROSS_GS
     ;,IO_KG_GS(IGR,15)            ,IACTSIMUL
     ;,IDEBUG        ,KTYPE_GR    ,MXKRIG_GS    ,MXROT_GR    ,MXSBX_GR
     ;,MXSBY_GR      ,MXSBZ_GR    ,MXDISC_GS    ,MXSB_GR     ,NDISC_GR
     ;,IDIMIVARIO_GS ,NEXDR_GR    ,NMXP_GR      ,NDMIN_GR    ,NOCT_GR
     ;,IDIMVAR_GS    ,MXCLOSE_GS  ,MXMEASPP_GS  ,MXZONPP_GS  ,MAXSAM_GR
     ;,NUMSB_GS      ,ICROSSCOV_GS,IFLAG_SIMUL
     ;,IPOLDRIFT_GS(1,IGR)  ,IVARIO_GS(1,1,IGR) ,ISUPBL_GS   ,ESTKRIG_GS   
     ;,POSMEAS_GS(1,1,1,IGR),POSMEAS_GS(1,2,1,IGR),POSMEAS_GS(1,3,1,IGR)
     ;,VMEAS_GS(NPP_GR+1,1,2,IGR) ,TRIM_GS(1,IGR)      
     ;,POSMEAS_GS(NPP_GR+1,1,2,IGR) ,POSMEAS_GS(NPP_GR+1,2,2,IGR)
     ;,POSMEAS_GS(NPP_GR+1,3,2,IGR) ,DATASC_GS(1,9)
     ;,DATASC_GS(1,1)   ,DATASC_GS(1,4)
     ;,VMEAS_GS(1,MXNVAR_GS+1,2,IGR)
     ;,VMEAS_GS(1,MXNVAR_GS+2,2,IGR)
     ;,VMEAS_GS(1,MXNVAR_GS+3,2,IGR)
     ;,VMEAS_GS(1,MXNVAR_GS+4,2,IGR)
     ;,VARIO_GS(1,1,IGR)    ,SEARCH_GS(1,IGR)    
     ;,VSTATS_GS(1,1,IGR)  
     ;,SUPBL_GS(1,IGR)
     ;,KRISOL_GS     ,CLOSESAM_GS ,KRISYS_GS    ,CROSSCOV_GS ,0.0
     ;,ROTMAT_GS     ,KRIGAUX_GS  ,POSDISAUX_GS ,NMEAS_GR+NPP_GR
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+1,2,IGR)
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+2,2,IGR)
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+3,2,IGR)
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+4,2,IGR)
     ;,DATASC_GS(1,11),DATASC_GS(1,16))

C____________________________  Correction of initial cova matrix of parameters. 
C____________________________  If CS is used, multiplies *2 the CE cova matrix.

                    DO ILOOP=1,NPP_GR*(NPP_GR+1)/2
                      COVPAR_GR(ILOOP)=COVPAR_GR(ILOOP)*2.0D0
                    END DO

                END IF ! IOPSC_GR.NE.0

C____________________________ Step 2.D) Inversion and storage of the cova. matrix
C____________________________           of current group (COVPAR_GR). It is a block 
C____________________________           of the global block diagonal cova. matrix 
C____________________________           of parameters(a priori) 

C____________________________           Inversion. Standardizes first to get rid 
C____________________________           of some big numbers

               ALLOCATE(daCovParSec(NPP_GR*(NPP_GR+1)/2))
               daCovParSec = COVPAR_GR        ! Keeps a back-up, just in case LINV1P fails

               CALL ZERO_ARRAY(WORK,IDIMWORK)
               CALL EQUAL_ARRAY (WORK,COVPAR_GR,NPP_GR*(NPP_GR+1)/2)

               dMaxCov = MAXVAL(COVPAR_GR)
               WORK = WORK / dMaxCov          ! Scales covariance matrix

               CALL LINV1P (WORK,NPP_GR,COVPAR_GR,D1,D2,IER)

      
               IF (IER.EQ.129) THEN
                  WRITE(*,2000) IGR
                  WRITE(MAINF,2000) IGR
 2000             FORMAT(//,' ERROR INVERTING A PRIORI COVARIANCE'
     ;                   ' MATRIX OF PARAMETERS OF GROUP OF ZONES ',I5,/
     ;                   ,' IT IS NOT POSITIVE DEFINITE.',/
     ;                   ,' I WILL USE ONLY ITS DIAGONAL. SORRY',/)


                  COVPAR_GR = 0.0D0
                  DO iLoop = 1, npp_gr
                     jPos = iLoop * (iLoop + 1) / 2
                     COVPAR_GR(jPos) = 1.0D0 / daCovParSec(jPos) 
                  END DO

               ELSE        ! Did not fail. Full covariance matrix is re-scaled

                  COVPAR_GR = COVPAR_GR / dMaxCov
                                 
               END IF 
      
               DEALLOCATE(daCovParSec)

C____________________________           Storage at COVPAR

               NZEROS=NPARDET+NTOTALPP

               DO ILOOP=1,NPP_GR          ! Loop over rows
                  DO ILOOP2=1,ILOOP       ! Loop over columns (symmetric)
                     IPOSCOVA_GR=ILOOP*(ILOOP+1)/2+ILOOP2-ILOOP  ! Cova. matrix of current group

                     IF (ILOOP.EQ.1 .AND. ILOOP2.EQ.1)                
     ;                  IPOSCOVA=NZEROS*(NZEROS+1)/2+NZEROS+1

                        COVPAR(IPOSCOVA)=COVPAR_GR(IPOSCOVA_GR)

                     IF (ILOOP2.EQ.ILOOP) THEN
                        IPOSCOVA=IPOSCOVA+NZEROS+1
                     ELSE
                        IPOSCOVA=IPOSCOVA+1
                     END IF

                  END DO   ! ILOOP2=1,ILOOP
               END DO      ! ILOOP=1,NPP_GR


C____________________________ Step 2.E) Saves local variables onto global

C____________________________   - Step 2.E.1) Weights of objective function of parameters
C____________________________                 and log-estimation options values

                IF (NUMITER.EQ.1) THEN
                   DO IPIPO=1,NPP_GR
                      WGT_UNK(NPARDET+NTOTALPP+IPIPO)=PAR_WGT(IACTTYPE)
                      INDPAR(NPARDET+NTOTALPP+IPIPO)=IOLG_PAR(IACTTYPE)
                   END DO ! IPIPO=1,NPP_GR
                END IF ! NUMITER.EQ.1

C____________________________   - Step 2.E.2) Prior information of pilot points

*** OJO2, PARA LA VERIFICACION DE DERIVADAS
*                ESTKRIG_GS(1)=ESTKRIG_GS(1)+1.D-6

                CALL EQUAL_ARRAY
     ;             (PARM(NPARDET+NTOTALPP+1),ESTKRIG_GS,NPP_GR)

C____________________________   - Step 2.E.3) Initial value of pilot points
C____________________________         - a) Previously calculated
C____________________________         - b) Previously calculated + white noise
C____________________________         - c) Set to zero
C____________________________         - d) Read at GEO file

                IF (IO_PARC_INI_GR.EQ.0) THEN   ! Set to calculated values

                   CALL EQUAL_ARRAY 
     ;             (PARC(NPARDET+NTOTALPP+1),ESTKRIG_GS,NPP_GR)

                ELSE IF (IO_PARC_INI_GR.EQ.1) THEN  ! Calc. values + White Noise

                   ISEED=IO_KG_GS(IGR,15)
                   DO IPIPO=1,NPP_GR
                      CALL NORMAL(ISEED,1.0D0,GAUPROB,0.0D0)
                      PARC(NPARDET+NTOTALPP+IPIPO)=
     ;                   ESTKRIG_GS(IPIPO)+GAUPROB
                   END DO

                ELSE IF (IO_PARC_INI_GR.EQ.2) THEN  ! Set to zero

                   CALL ZERO_ARRAY 
     ;             (PARC(NPARDET+NTOTALPP+1),NPP_GR)

                ELSE IF (IO_PARC_INI_GR.EQ.3) THEN  ! Read at GEO file
    
                   CALL EQUAL_ARRAY
     ;             (PARC(NPARDET+NTOTALPP+1),PARC_GS(1,IGR),NPP_GR)

                END IF ! IO_PARC_INI_GR.EQ....

                                       ! Saves pilot point values as meas. First part of the VMEAS vector (1:NPP_GR)

                CALL EQUAL_ARRAY
     ;            (VMEAS_GS(1,1,1,IGR),PARC(NPARDET+NTOTALPP+1),NPP_GR)

             END IF  ! (IOINV.GT.0. AND. NPP_GR.GT.0) . AND. ...

C____________________________ Step 3: CE / CS of zones to measurements. Only if they are not read at PAR file
C____________________________         OUTPUT : - Initial value of zonal parameters

             IF (IO_PARC_INI_GR.NE.3) THEN
                IF ( (IOINV.LT.0 .AND. NUMITER.EQ.1) .OR.   ! Direct problem OR
     ;             ( (IOINV.GT.0 .AND. NPP_GR.GT.0) .AND.   ! Inverse problem
     ;               (NUMITER.EQ.1 .OR. IO_RDPP_GR.NE.0) ) ) THEN  

                                 ! Pilot point and sampling locations and meas. 
                                 ! values are saved (arrays VMEAS_GS,
                                 ! POSMEAS_GS are changed in GSLIB routines).
                                 ! Main outputs of current step are also init.

                   CALL EQUAL_ARRAY (POSMEAS_GS(1,1,2,IGR)
     ;                     ,POSMEAS_GS(1,1,1,IGR),3*MXMEASPP_GS)

                   CALL EQUAL_ARRAY(VMEAS_GS(1,1,2,IGR)
     ;                     ,VMEAS_GS(1,1,1,IGR),MXMEASPP_GS*IDIMVAR_GS)

                   CALL ZERO_ARRAY (ESTKRIG_GS,MXZONPP_GS)
                   CALL ZERO_ARRAY (ZNWGT_GS,MXNPRIM_GS*MXNZON_GS)
                   CALL ZERO_ARRAY_I (IZN_PP_GS,IDIMZONPP_GS)
                   CALL ZERO_ARRAY_I (IZN_NPP_GS,MXNZON_GS)
                   CALL ZERO_ARRAY_I (ICHECK_GS,MXNPP_GS)

                   IF (IOPSC_GR.EQ.0) THEN   ! CONDITIONAL ESTIMATION

                      IDEBUG=0
                      IF (IFLAGS(21).EQ.2) THEN
                         IDEBUG=1
                         WRITE(666,1102) 
                      END IF
 1102                 FORMAT(//,' CONDITIONAL ESTIMATION OF ZONES'
     ;                          ' TO MEASUREMENTS',/,
     ;                          ' =========== ========== == =====' 
     ;                          ' == ============',/)

                                 ! Calculates maximum number of kriging / cokriging equations

                      MAX_1=NVAR_GR*NMEAS_GR+NVAR_GR
                      MAX_2=NMEAS_GR+10+NEXDR_GR
                      MXKRIG_GR=MAX0(MAX_1,MAX_2)

                      IF (KTYPE_GR.LE.4) THEN   ! Kriging

                         CALL KRIGING
     ;(IDEBUG                ,IDIMCROSS_GS        ,0
     ;,0                     ,0                   ,KTYPE_GR
     ;,MAINF                 ,MXKRIG_GS           ,MXROT_GR       
     ;,MXSB_GR               ,MXDISC_GS           ,NMEAS_GR
     ;,NDISC_GR              ,IVARIO_GS(1,1,IGR)  ,NZN_GR         
     ;,NEXDR_GR              ,NMXP_GR             ,NDMIN_GR      
     ;,NOCT_GR               ,NRESTRI             ,MXSBX_GR      
     ;,MXSBY_GR              ,MXSBZ_GR            ,VARIO_GS(1,1,IGR)         
     ;,SEARCH_GS(1,IGR)      ,SEARCH_GS(9,IGR)    ,SEARCH_GS(10,IGR)    
     ;,SEARCH_GS(11,IGR)     ,SEARCH_GS(3,IGR)    ,SEARCH_GS(4,IGR)     
     ;,VSTATS_GS(1,2,IGR)    ,SUPBL_GS(1,IGR)     ,SUPBL_GS(4,IGR)      
     ;,SUPBL_GS(2,IGR)       ,SUPBL_GS(5,IGR)     ,SUPBL_GS(3,IGR)      
     ;,SUPBL_GS(6,IGR)       ,VARIO_GS(1,4,IGR)   ,VARIO_GS(1,5,IGR)    
     ;,VARIO_GS(1,6,IGR)     ,VARIO_GS(1,7,IGR)   ,VARIO_GS(1,8,IGR)    
     ;,CLOSESAM_GS(1,1)      ,KRISOL_GS(1,1)      ,KRISOL_GS(1,2)       
     ;,KRISYS_GS             ,CROSSCOV_GS(1,1,1)  
     ;,VMEAS_GS(NPP_GR+1,1,2,IGR)         
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+1,2,IGR)              
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+2,2,IGR)      
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+3,2,IGR)
     ;,VMEAS_GS(NPP_GR+1,MXNVAR_GS+4,2,IGR)         
     ;,IPOLDRIFT_GS(1,IGR)   ,IVARIO_GS(1,2,IGR)  ,ISUPBL_GS(1,1)        
     ;,ISUPBL_GS(1,2)        ,ISUPBL_GS(1,3)      ,NUMSB_GS              
     ;,ICROSSCOV_GS(1,1)     ,POSDIS_GS(1,1,1,IGR),VARIO_GS(1,2,IGR)     
     ;,ROTMAT_GS             ,EXDRZN_GS(1,1,IGR)  ,EXDRZN_GS(1,2,IGR)   
     ;,EXDRZN_GS(1,3,IGR)    ,EXDRZN_GS(1,4,IGR)  ,VARIO_GS(1,3,IGR)     
     ;,CLOSESAM_GS(1,2)      ,TRIM_GS(1,IGR)      ,KRIGAUX_GS(1,4)       
     ;,KRIGAUX_GS(1,5)       ,KRISOL_GS(1,3)      ,CROSSCOV_GS(1,1,2)
     ;,POSMEAS_GS(NPP_GR+1,1,2,IGR)               ,POSZN_GS(1,1,IGR)   
     ;,KRIGAUX_GS(1,1)       ,POSDISAUX_GS(1,1)     
     ;,POSMEAS_GS(NPP_GR+1,2,2,IGR)               ,POSZN_GS(1,2,IGR)         
     ;,KRIGAUX_GS(1,2)       ,POSDISAUX_GS(1,2)  
     ;,POSMEAS_GS(NPP_GR+1,3,2,IGR)               ,POSZN_GS(1,3,IGR)     
     ;,KRIGAUX_GS(1,3)       ,POSDISAUX_GS(1,3)
     ;,ESTKRIG_GS            ,POSMEAS_GS(NPP_GR+1,1,1,IGR)
     ;,POSMEAS_GS(NPP_GR+1,2,1,IGR) ,POSMEAS_GS(NPP_GR+1,3,1,IGR)
     ;,MXNPRIM_GS            ,MXMEASPP_GS         ,R_DUMMY)

                      ELSE       ! Cokriging
                         CONTINUE
                      END IF     ! Kriging / Cokriging

                   ELSE          ! CONDITIONAL / UNCONDITIONAL SIMULATION

                      IDEBUG=0
                      IF (IFLAGS(22).EQ.2) THEN
                         IDEBUG=1
                         WRITE(666,1103)
                      END IF
 1103                 FORMAT(//,' CONDITIONAL SIMULATION OF ZONES'
     ;                          ' TO MEASUREMENTS',/,
     ;                          ' =========== ========== == =====' 
     ;                          ' == ============',/)

                      IF (IOPSC_GR.EQ.1) THEN   ! CONDITIONAL SIMULATION

                                 ! Calculates maximum number of kriging / cokriging equations

                         MAX_1=NVAR_GR*(NMEAS_GR+NZN_GR)+NVAR_GR  ! Max kriging dim
*                         MAX_2=NMEAS_GR+NMEAS_GR+NZN_GR+10+NEXDR_GR
                         MAX_2=NMEAS_GR+NZN_GR+10+NEXDR_GR
                         MXKRIG_GR=MAX0(MAX_1,MAX_2)

                         CALL CONDITIONAL_SIMULATION
     ;(NZN_GR      ,NMEAS_GR    ,MAINF
     ;,MXNPRIM_GS  ,IDIMCROSS_GS
     ;,IO_KG_GS(IGR,15)           ,IACTSIMUL
     ;,IDEBUG        ,KTYPE_GR    ,MXKRIG_GS    ,MXROT_GR    ,MXSBX_GR
     ;,MXSBY_GR      ,MXSBZ_GR    ,MXDISC_GS    ,MXSB_GR     ,NDISC_GR
     ;,IDIMIVARIO_GS ,NEXDR_GR    ,NMXP_GR      ,NDMIN_GR    ,NOCT_GR
     ;,IDIMVAR_GS    ,MXCLOSE_GS  ,MXMEASPP_GS  ,MXZONPP_GS  ,MAXSAM_GR
     ;,NUMSB_GS      ,ICROSSCOV_GS,IFLAG_SIMUL
     ;,IPOLDRIFT_GS(1,IGR)  ,IVARIO_GS(1,1,IGR) ,ISUPBL_GS   ,ESTKRIG_GS   
     ;,POSZN_GS(1,1,IGR)    ,POSZN_GS(1,2,IGR)  ,POSZN_GS(1,3,IGR)
     ;,VMEAS_GS(NPP_GR+1,1,2,IGR)        ,TRIM_GS(1,IGR)      
     ;,POSMEAS_GS(NPP_GR+1,1,2,IGR)      ,POSMEAS_GS(NPP_GR+1,2,2,IGR)
     ;,POSMEAS_GS(NPP_GR+1,3,2,IGR)
     ;,DATASC_GS(1,9)     
     ;,DATASC_GS(1,1)   ,DATASC_GS(1,4)
     ;,EXDRZN_GS(1,1,IGR),EXDRZN_GS(1,2,IGR),EXDRZN_GS(1,3,IGR)
     ;,EXDRZN_GS(1,4,IGR)
     ;,VARIO_GS(1,1,IGR)    ,SEARCH_GS(1,IGR)    
     ;,VSTATS_GS(1,1,IGR)  
     ;,SUPBL_GS(1,IGR)
     ;,KRISOL_GS     ,CLOSESAM_GS ,KRISYS_GS    ,CROSSCOV_GS ,POSDIS_GS
     ;,ROTMAT_GS     ,KRIGAUX_GS  ,POSDISAUX_GS ,NMEAS_GR+NZN_GR
     ;,VMEAS_GS(NPP_GR+1,2,2,IGR) ,VMEAS_GS(NPP_GR+1,3,2,IGR)
     ;,VMEAS_GS(NPP_GR+1,4,2,IGR) ,VMEAS_GS(NPP_GR+1,5,2,IGR)       
     ;,DATASC_GS(1,11),DATASC_GS(1,16))

                      ELSE                      ! UNCONDITIONAL SIMULATION

                      END IF                    ! CONDITIONAL / UNCONDITIONAL SIM.

                   END IF  ! CONDITIONAL ESTIMATION / SIMULATION


C____________________________ Step 3.A: Estimated / simulated values are saved onto PARZ

                   ICOMPO=0
                   DO IZON=1,NZONE_PAR(IPNZPAR)
                      IF (IGR_ZONE(IPOS+IZON).EQ.IGR) THEN  ! Zone belongs to group
                         ICOMPO=ICOMPO+1
                         PARZ(IPOS+IZON)=ESTKRIG_GS(ICOMPO)
                         IF (IACTTYPE.EQ.1) THEN    ! Special case: Isotropic T
                            IF (ISOZ(IZON).EQ.1) THEN
                               DO IDIM=1,LDIM_GS(IZON)
                                  PARZ(INORPAR(IDIM)+IZON) =
     ;                               ESTKRIG_GS(ICOMPO)
                               END DO ! IDIM=1,LDIM(NZ)
                            END IF ! ISOZ(IZON).EQ.1
                         END IF ! IACTTYPE.EQ.1
                      END IF ! IGR_ZONE(IPOS+IZON).EQ.IGR

                   END DO ! IZON=1,NZONE_PAR(IPNZPAR)

                END IF ! (IOINV.LT.0 .AND. NUMITER.EQ.1) .OR. ...

             ELSE   ! IO_PARC_INI_GR.EQ.3 : PARZ/PARC are given

                IF (NUMITER.EQ.1) THEN
                  ICOMPO=0
                  DO IZON=1,NZONE_PAR(IPNZPAR)
                      IF (IGR_ZONE(IPOS+IZON).EQ.IGR) THEN  ! Zone belongs to group
                         ICOMPO=ICOMPO+1
                         IF (IOLG_PAR(IACTTYPE).NE.0) THEN  ! Log-estimation
                            AUXPARZ=DLOG10(PARZ(IPOS+IZON))
                         ELSE
                            AUXPARZ=PARZ(IPOS+IZON)
                         END IF         

                         PARZ(IPOS+IZON)=AUXPARZ
                         IF (IACTTYPE.EQ.1) THEN    ! Special case: Isotropic T
                            IF (ISOZ(IZON).EQ.1) THEN
                               DO IDIM=1,LDIM_GS(IZON)
                                  PARZ(INORPAR(IDIM)+IZON) =AUXPARZ
                               END DO ! IDIM=1,LDIM(NZ)
                            END IF ! ISOZ(IZON).EQ.1
                         END IF ! IACTTYPE.EQ.1
                      END IF ! IGR_ZONE(IPOS+IZON).EQ.IGR

                   END DO ! IZON=1,NZONE_PAR(IPNZPAR)

                 END IF ! NUMITER.EQ.1

             END IF ! IO_PARC_INI_GR.EQ.3 : PARC, PARZ are given

****   END OF STEP 3. CE/CS OF ZONES TO MEASUREMENTS

C____________________________ Step 4: CE of zones to pilot points and measurements. 
C____________________________         OUTPUT : - Weights of zonal param. to pilot points

             IF ( (IOINV.GT.0. AND. NPP_GR.GT.0) ! Inverse problem at this group
     ;      .AND. (IO_RDPP_GR.NE.0               ! Variable position of pilot p.
     ;       .OR.  NUMITER.EQ.1) ) THEN         ! At first iteration, ALWAYS

                                 ! Pilot point and sampling locations and meas. 
                                 ! values are saved (arrays VMEAS_GS,
                                 ! POSMEAS_GS are changed in GSLIB routines).
                                 ! Main outputs of current step are also init.

                CALL EQUAL_ARRAY (POSMEAS_GS(1,1,2,IGR)
     ;                     ,POSMEAS_GS(1,1,1,IGR),3*MXMEASPP_GS)

                CALL EQUAL_ARRAY(VMEAS_GS(1,1,2,IGR),VMEAS_GS(1,1,1,IGR)
     ;                     ,MXMEASPP_GS*IDIMVAR_GS)

                CALL ZERO_ARRAY (ESTKRIG_GS,MXZONPP_GS)
                CALL ZERO_ARRAY (ZNWGT_GS,MXNPRIM_GS*MXNZON_GS)
                CALL ZERO_ARRAY_I (IZN_PP_GS,IDIMZONPP_GS)
                CALL ZERO_ARRAY_I (IZN_NPP_GS,MXNZON_GS)
                CALL ZERO_ARRAY_I (ICHECK_GS,MXNPP_GS)
                
                IF (KTYPE_GR.LE.4) THEN   ! Kriging

                   IDEBUG=0
                   IF (IFLAGS(21).EQ.3) THEN
                      IDEBUG=1
                      WRITE(666,1104)
                   END IF
 1104              FORMAT(//,' CONDITIONAL ESTIMATION OF ZONES'
     ;                       ' TO PILOT POINTS AND MEASUREMENTS',/,
     ;                       ' =========== ========== == =====' 
     ;                       ' == ===== ====== === ============',/)

                                 ! Calculates maximum number of kriging / cokriging equations

                   MAX_1=NVAR_GR*(NPP_GR+NMEAS_GR)+NVAR_GR   ! Max kriging dim
                   MAX_2=NPP_GR+NMEAS_GR+10+NEXDR_GR
                   MXKRIG_GR=MAX0(MAX_1,MAX_2)

                   CALL KRIGING
     ;(IDEBUG                ,IDIMCROSS_GS        ,0
     ;,1                     ,1                   ,KTYPE_GR      
*     ;,MAINF                 ,MXKRIG_GR           ,MXROT_GR       
     ;,MAINF                 ,MXKRIG_GS           ,MXROT_GR       
     ;,MXSB_GR               ,MXDISC_GS           ,NPP_GR+NMEAS_GR
     ;,NDISC_GR              ,IVARIO_GS(1,1,IGR)  ,NZN_GR         
     ;,NEXDR_GR              ,NMXP_GR             ,NDMIN_GR      
     ;,NOCT_GR               ,NRESTRI             ,MXSBX_GR      
     ;,MXSBY_GR              ,MXSBZ_GR            ,VARIO_GS(1,1,IGR)         
     ;,SEARCH_GS(1,IGR)      ,SEARCH_GS(9,IGR)    ,SEARCH_GS(10,IGR)    
     ;,SEARCH_GS(11,IGR)     ,SEARCH_GS(3,IGR)    ,SEARCH_GS(4,IGR)     
     ;,VSTATS_GS(1,2,IGR)    ,SUPBL_GS(1,IGR)     ,SUPBL_GS(4,IGR)      
     ;,SUPBL_GS(2,IGR)       ,SUPBL_GS(5,IGR)     ,SUPBL_GS(3,IGR)      
     ;,SUPBL_GS(6,IGR)       ,VARIO_GS(1,4,IGR)   ,VARIO_GS(1,5,IGR)    
     ;,VARIO_GS(1,6,IGR)     ,VARIO_GS(1,7,IGR)   ,VARIO_GS(1,8,IGR)    
     ;,CLOSESAM_GS(1,1)      ,KRISOL_GS(1,1)      ,KRISOL_GS(1,2)       
     ;,KRISYS_GS             ,CROSSCOV_GS(1,1,1)  
     ;,VMEAS_GS(1,1,2,IGR)         
     ;,VMEAS_GS(1,MXNVAR_GS+1,2,IGR)              
     ;,VMEAS_GS(1,MXNVAR_GS+2,2,IGR)      
     ;,VMEAS_GS(1,MXNVAR_GS+3,2,IGR)
     ;,VMEAS_GS(1,MXNVAR_GS+4,2,IGR)         
     ;,IPOLDRIFT_GS(1,IGR)   ,IVARIO_GS(1,2,IGR)  ,ISUPBL_GS(1,1)        
     ;,ISUPBL_GS(1,2)        ,ISUPBL_GS(1,3)      ,NUMSB_GS              
     ;,ICROSSCOV_GS(1,1)     ,POSDIS_GS(1,1,1,IGR),VARIO_GS(1,2,IGR)     
     ;,ROTMAT_GS             ,EXDRZN_GS(1,1,IGR)  ,EXDRZN_GS(1,2,IGR)   
     ;,EXDRZN_GS(1,3,IGR)    ,EXDRZN_GS(1,4,IGR)  ,VARIO_GS(1,3,IGR)     
     ;,CLOSESAM_GS(1,2)      ,TRIM_GS(1,IGR)      ,KRIGAUX_GS(1,4)       
     ;,KRIGAUX_GS(1,5)       ,KRISOL_GS(1,3)      ,CROSSCOV_GS(1,1,2)
     ;,POSMEAS_GS(1,1,2,IGR) ,POSZN_GS(1,1,IGR)   
     ;,KRIGAUX_GS(1,1)       ,POSDISAUX_GS(1,1)     
     ;,POSMEAS_GS(1,2,2,IGR) ,POSZN_GS(1,2,IGR)         
     ;,KRIGAUX_GS(1,2)       ,POSDISAUX_GS(1,2)  
     ;,POSMEAS_GS(1,3,2,IGR) ,POSZN_GS(1,3,IGR)     
     ;,KRIGAUX_GS(1,3)       ,POSDISAUX_GS(1,3)
     ;,ESTKRIG_GS            ,POSMEAS_GS(1,1,1,IGR)
     ;,POSMEAS_GS(1,2,1,IGR) ,POSMEAS_GS(1,3,1,IGR)
     ;,MXNPRIM_GS            ,MXMEASPP_GS         ,R_DUMMY)

                ELSE  ! Cokriging
                   CONTINUE
                END IF ! Kriging / Cokriging

C____________________________ Step 4.A: Saves weights and related pilot points

                NONUSEPP=0
                DO IZON=1,NZN_GR
                   IUSEPP=0
                   DO ILOOP=1,NMXP_GR
                      IF (ICROSSCOV_GS(ILOOP,IZON).GT.0 .AND. 
     ;                    ICROSSCOV_GS(ILOOP,IZON).LE.NPP_GR) THEN
                         IUSEPP=IUSEPP+1
                         IZN_NPP_GS(IZON) = IZN_NPP_GS(IZON) + 1
                      END IF ! ICROSSCOV_GS(IPOS,IZON).GT.0
                   END DO ! ILOOP=1,NMEAS_GR + NPP_GR
                   IF (IUSEPP.EQ.0) NONUSEPP = NONUSEPP + 1
                END DO ! IZON=1,NZN_GR

! 100            IPOS=INORPAR(IPINORPAR)          ! Initial position at IGR_ZONE
                ISTART=IPOS*IDIMWGT              ! Initial position at IPNT_PAR
                IPOSDLT_PAR=NPARDET+NTOTALPP

C_______________________ Step 4.A.1: Loop over zones belonging to this group

                ISUM=0
                DO IZON=1,NZONE_PAR(IPNZPAR)
                   IF (IGR_ZONE(IPOS+IZON).EQ.IGR) THEN  ! Zone belongs to group

C_______________________ Step 4.A.1.1: Assigns first and last useful components
C_______________________               at IPNT_PAR. 

                      ISUM=ISUM+1
                      IPNT_START(IPOS+IZON)=ISTART+(IZON-1)*IDIMWGT+1
                      IPNT_END(IPOS+IZON)=
     ;                    IPNT_START(IPOS+IZON)+IZN_NPP_GS(ISUM)-1

C_______________________ Step 4.A.1.2: Loop over pilot points used on the 
C_______________________               parameterization, assigning the value of
C_______________________               IPNT_PAR & WGT_PAR

                      IADD=0
                      DO ICOMPO=1,NMXP_GR
                         IDPIPO=ICROSSCOV_GS(ICOMPO,ISUM)

                         IF (IDPIPO.GT.0. AND. IDPIPO.LE.NPP_GR) THEN

                            IADD=IADD+1
                            IPNT_PAR(IPNT_START(IPOS+IZON)-1+IADD)=
     ;                         IPOSDLT_PAR+IDPIPO
                            WGT_PAR(IPNT_START(IPOS+IZON)-1+IADD)=
     ;                         CROSSCOV_GS(IDPIPO,ISUM,2)

                         END IF ! IZN_PP_GS(ICOMPO).NE.0           

                      END DO ! ICOMPO=1,IZN_NPP_GS(IZON)

                   END IF ! IGR_ZONE(IPOS+IZON).EQ.IGR

                END DO ! IZON=1,NZONE_PAR(IPNZPAR)

C____________________________ Step 5: Checks coherency of the methodology

C____________________________   - Step 5.A: A zonal parameter is not parameterized using
C____________________________               pilot points. It will not vary during optimiz.
C____________________________               process

                IF (NONUSEPP.NE.0) THEN
                   WRITE(MAINF,1400) IGR,NONUSEPP
                   WRITE(MAINF,1500)
                   WRITE(6,1400) IGR,NONUSEPP
 1400              FORMAT(//,' WARNING DRAWING INITIAL GUESS OF ZONES'
     ;                      ,' BELONGING TO GROUP: ',I5,/,I5
     ;                      ,' ZONE(S) DOES NOT USE PILOT POINTS.'
     ;                       ' CHECK FILE RES.OUT',/)
 1500              FORMAT(/,5X' ZONE',14X,'X',14X,'Y',14X,'Z',/
     ;                     ,5X' ====',14X,'=',14X,'=',14X,'=')

                   DO IZON=1,NZN_GR
                     IF (IZN_NPP_GS(IZON).EQ.0) THEN
                        WRITE(MAINF,1600) IZON
     ;                      ,(POSZN_GS(IZON,ILOOP,IGR),ILOOP=1,3)
                     END IF ! IZN_NPP_GS(IZON).EQ.0
 1600                FORMAT(5X,I5,3(5X,F10.3))
                   END DO ! IZON=1,NZON_GR
                END IF ! NONUSEPP.NE.0

C____________________________   - Step 5.B: All pilot points of this group are used. It can 
C____________________________               cause singularity of hessian matrix if prior 
C____________________________               information is not used

                DO IZON=1,NZN_GR
                   DO ILOOP=1,NMXP_GR
                      IDPIPO=ICROSSCOV_GS(ILOOP,IZON)
                      IF (IDPIPO.GT.0 .AND. IDPIPO.LE.NPP_GR .AND. 
     ;                   ICHECK_GS(IDPIPO).EQ.0) ICHECK_GS(IDPIPO)=1
                   END DO ! ILOOP=1,NMXP_GR
                END DO ! IZON=1,NZN_GR

                NONUSEPP=0

                DO ILOOP=1,NPP_GR
                   IF (ICHECK_GS(ILOOP).EQ.0) NONUSEPP=NONUSEPP+1
                END DO ! ILOOP=1,NPP_GR

                IF (NONUSEPP.NE.0) THEN
                   WRITE(MAINF,1700) IGR,NONUSEPP
                   WRITE(MAINF,1800)
                   WRITE(*,1700) IGR,NONUSEPP
                   WRITE(*,1800)
 1700              FORMAT(//,' WARNING SOLVING KRIGING OF ZONES'
     ;                       '  BELONGING'
     ;                      ,' TO GROUP: ',I5,/,I5,' PILOT POINT(S)'
     ;                      ,' IS NOT USED.'
     ;                    ,/,' TRY TO INCREASE PRIMARY VARIABLE SEARCH'
     ;                      ,' RADIUS OR REMOVE IT. IF YOU DO NOT'
     ;                       ' CONSIDER'
     ;                       ' PRIOR INFORMATION, IT WILL LEAD TO A'
     ;                       ' SINGULARITY OF JACOBIAN MATRIX',/)
 1800              FORMAT(/,5X'POINT',14X,'X',14X,'Y',14X,'Z',/
     ;                     ,5X'=====',14X,'=',14X,'=',14X,'=')
 
                   DO ILOOP=1,NPP_GR
                      IF (ICHECK_GS(ILOOP).EQ.0) THEN
                         WRITE(MAINF,1600) ILOOP,
     ;                      (POSMEAS_GS(ILOOP,ILOOP2,1,IGR),ILOOP2=1,3)
                         WRITE(*,1600) ILOOP,
     ;                      (POSMEAS_GS(ILOOP,ILOOP2,1,IGR),ILOOP2=1,3)
                      END IF ! ICHECK_GS(ILOOP).EQ.0
                   END DO ! ILOOP=1,NPP_GR
               
                END IF ! NONUSEPP.NE.0

             END IF ! (IOINV.GT.0. AND. NPP_GR.GT.0) ...

          END IF  ! Estimated geostatistically

      END DO  ! IGR=1,NGROUP_ZN

C____________________________ Step 6) Echoes inverse of the a priori
C____________________________         covariance matrix of parameters
C____________________________         (determ. + pilot points)

      IF (IOINV.GT.0) THEN

         IF (NUMITER.EQ.1 .AND.
     ;      (IFLAGS(23).EQ.2 .OR. IFLAGS(23).EQ.3) ) THEN
            WRITE(666,2100) NUMITER
 2100       FORMAT(//,22X,'INVERTED A PRIORI COVARIANCE MATRIX.'
     ;                    ' ITERATION:',I5,/,22X,'========'
     ;                    ' = ====== ========== ======= ==========',/)
            WRITE(666,3000) COVPAR
         END IF ! IFLAGS(22).NE.0

      END IF ! IOINV.GT.0

 3000 FORMAT(7E10.4)

      DEALLOCATE(CROSSCOV_GS)
      RETURN
      END
     













