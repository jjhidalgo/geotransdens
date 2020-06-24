       PROGRAM TRANSDENS
****************************************************************************
*
* PURPOSE
*     Makes the partition of two big arrays, RV and IV
*
* DESCRIPTION
*
*     Calls LECDIM to read the dimensions and then makes the partitions.
*     Most variables starting with ID are used only to position arrays
*     inside the working space of RV and IV and are not included in this
*     header.
*
* INTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output
*                         data files
*  FOBJ_WGT               Array containing all objective function weights for
*                         state variables (heads, concentrations, etc)
*  IFLAGS                 Array with different writing options. Used mainly for
*                         debugging.
*  INORPAR                Array containing the indexes with the location
*                         of the different paramaters in arrays PARZ, PARM,
*                         IVPAR, NFTPAR, STPAR and FNTPAR
*  IOLG_PAR               Array containing all logarithmic options of
*                         estimated parameters
*  IOWRITE                Array containing all output options
*  IPAR_DIR               Array containing all integer direct problem
*                         parameters
*  IPAR_INV               Array containing all integer inverse problem
*                         parameters
*  IV                     Integer array used to reserve most of the necessary
*                         space for all the integer variables of the problem.
*                         Only some small arrays are not included in it:
*                         IFLAGS, INORPAR, IOLG_PAR, IOWRITE ,IPAR_DIR,
*                         IPAR_INV and NZONE_PAR
*  KV                     Character array used to reserve most of the necessary
*                         space for all the string variables of the problem.
*  NZONE_PAR              Array containing the number of zones of all
*                         parameters
*  PAR_DIR                Array containing all real direct problem
*                         parameters
*  PAR_INV                Array containing all real inverse problem
*                         parameters
*  PAR_WGT                Array containing objective function weights for
*                         all estimated parameters
*  RV                     Real array used to reserve most of the necessary
*                         space for all the real variables of the problem.
*                         Only some small arrays are not included in it:
*                         PAR_DIR, PAR_INV, PAR_WGT and FOBJ_WGT
*
* INTERNAL VARIABLES: SCALARS
*
*  IDIMAFLU               Used to dimension array AFLU
*  IDIMBB                 Used to dimension array BIBI
*  IDIMDFLU               Used to dimension array DFLU
*  IDIMDQ                 Used to dimension array DAT_VD (second dimension)
*  IDIMDTRA               Used to dimension array DTRA (second dimension)
*  IDIMFNT                Used for dimensioning array FNT, it coincides with
*                         NFNT if this latter is not zero
*  IDIMHESS               Used to dimension array HESS. It is equal to
*                         NPAR*(NPAR+1)/2
*  IDIMQ                  Used to dimension array QXYZ
*  IERROR                 Current number of errors on input data
*  IOBALC                 Same as IOBALH for solute mass balance
*  IOBALGC                Same as IOBALH for solute mass balance
*  IOBALGH                If non zero, the code computeds a global mass balance
*  IOBALH                 Flow mass balance computation option
*  IOCNSF                 Scheme for storage term in flow problem
*  IOCNST                 Scheme for mass storage term in transport problem
*  IOCRITRAP              Option of treatement on the direct problem convergence
*                         criteria
*  IODIM                  Maximum dimension of any element included
*                         in the problem
*  IOEQT                  Type of problem to be solved
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1
*  IOFLSAT                Indicates the possibility that one part of the domain
*                         reaches unsaturated state.
*  IOFMLF                 Flow Formulation number
*  IOFMLT                 Transport formulation number
*  IOFOBJ                 If 0, the code solves completely the direct problem
*                         at each iteration of the inverse problem. If sets
*                         1 it computes completely the simulation only for
*                         succesfull iterations of the inverse problem.
*  IOINV                  Inverse problem option
*  IOLD                   Counts the location of the different integer arrays
*                         in the space defined by array IV
*  IOOBS                  If it equals 1, all nodal points are observation
*                         points. Otherwise they are a set of specified points
*                         (usually 0)
*  IOPART                 If 1, the list of problem partitions is written
*  IOPINITC               Option for the extrapolation of concentrations
*                         at the next time step in the Newton process.
*  IOPINITH               Option for the extrapolation of heads or pressures
*                         at the next time step in the Newton process.
*  IOPRHED                Indicates whether the flow state variable state is
*                         preasure (set to 1) or head (set to 0)
*  IORTS                  Transport regime
*  IOSTC                  If it equals 1, standard deviation of measured
*                         concentration errors is taken constant (usually 0)
*  IOSTH                  If it equals 1, standard deviation of measured head
*                         errors is taken constant (usually 0)
*  IOSUCHUM               Indicates if measures are given in terms of pressure
*                         or piezometric heads (set to 0), or in terms of
*                         saturation degree
*  IOTRLI                 Idem to IOFLLI, in the case of transport.
*  IOTRS                  Flow regime
*  IPROCESS               Controls if it4s the first time that subroutine
*                         LECDIM is called.
*  ISOT                   Maximum hydraulic conductivity tensor anisotropy
*                         degree in the problem.
*  ITT                    Variable needed for routine INIT_CPUTIME.
*                         It is not used
*  KOLD                   Counts the location of the different integer arrays
*                         in the space defined by array KV
*  LMXNDL                 Maximum number of nodes per element
*  MAINF                  Unit number of the main output file (RES.OUT)
*  MAXNEOP                Used to reserve some space. It is equal to
*                         MAX (NUMEL,NUMNP,NUOBS,NPAR)
*  NBAND                  Band with (maximumdifferecne between the numbers of
*                         two nodes
*                         belonguing to the same element
*  NBAND1                 Used to dimension. It is equal to NBAND+1
*  NBAND2                 Used to dimension. It is equal to 2*NBAND+1
*  NFNL                   Total number of non-linear functions required
*                         in the problem
*  NFNT                   Number of time functions used for describing time
*                         dependence of all transient parameters
*  NINT                   Number of observation times
*  NOLD                   Counts the location of the different integer arrays
*                         in the space defined by array RV
*  NPAR                   Total number of parameters to be estimated
*  NPAREL                 Number of element parameters in current problem
*  NPARF                  Number of transient parameters to be estimated
*  NPARFPRG               Number of uncertain generic parameter zones involved
*                         in the non-linear flow inverse problem
*  NPARNP                 Number of nodal parameters in current problem
*  NPARPRG                Total number of uncertain generic parameter zones
*                         involved in the inverse problem
*  NTDMT                  If diferent form zero, number of terms in matrix
*                         diffusion(if zero, no diffusion)
*  NUMEL                  Number of elements
*  NUMNP                  Number of nodes
*  NUOBS                  Number of observation points
*  NZALF                  Number of leakance zones
*  NZARR                  Number of areal recharge zones
*  NZCHP                  Number of prescribed head zones
*  NZDMT                  Number of matrix diffusion zones
*  NZFOD                  Number of zones of first order decay
*  NZFOF                  Number of zones of formation factor      
*  NZPAR                  Total number of zones for all nodal and element
*                         parameters including each transmissivity tensor
*                         component.
*  NZQQP                  Number of prescribed flow zones
*  NZSTG                  Number of storage Coefficient zones
*  NZTRA                  Number of transmissivity zones
*  TIMECPU                Cpu time in hour, minutes and seconds
*  ZERO                   Real number 0.0 in double precision
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  CPUTIME                Computes the cpu time consumed
*  INIT_CPUTIME           Initializes the counter of cpu time. It is useless
*                         in most machines except VMS
*  INIT_INTEG             Defines some integer variables
*  LECDIM                 Reads dimensions, options and scalar parameters of
*                         the current problem
*  PRINCIPAL              It is actually the routine that runs the code
*  WRI_PART               Writes the partitions of arrays RV and IV
*
* HISTORY
*
*     AMS       3-1997     First coding
*     AMS       1-1998     Revision
*****************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*14 FILENAME(20)*20
       INCLUDE 'COMMON.FOR'
       INCLUDE 'MAIN_COMM.FOR'
       DATA ZERO/0.D0/

C--------------- Total dimension for program TRANSIN


C--------------- Parameters to dimension arrays used to group scalar
C--------------- variables

       PARAMETER (NTYPAR=30, NPARALG=40, NWRITE=20, NFLAGS=50,
     ;            NSTAT=10 , NPPNP=12  , NPPEL=17 , NOPTS=50,
     ;            NTYPEL=6 , MAXPG=6   , MXGRPZN=1000)

       DIMENSION NZONE_PAR(NTYPAR),IOLG_PAR(NTYPAR,2),IPAR_INV(NPARALG),
     ;   IPAR_DIR(NPARALG),IOWRITE(NWRITE),IFLAGS(NFLAGS),
     ;   PAR_WGT(NTYPAR),PAR_INV(NPARALG),PAR_DIR(NPARALG),
     ;   FOBJ_WGT(NSTAT),INORPAR(NTYPAR),IOPTS(NOPTS),
     ;   LINMET(3,2)
     ;   ,IOPT_GS(MXGRPZN,20),IO_KG_GS(MXGRPZN,16)

C------------------------- FIRST EXECUTABLE STATEMENT.

       WRITE (*,10)
   10  FORMAT(/,'GEOTRANSDENS v0.1 (beta version)'
     &       ,/,'==============================',/)


C--------------- Initializes cputime counter to zero

       CALL INIT_CPUTIME (ITT)

C--------------- Initializes some variables

       IPROCESS=0

C--------------- Reads dimensions and problem options

       CALL LECDIM
     ;(IDIMFNT  ,IERROR   ,IOPTS(18),IOPTS(22),IOPTS(21),IOPTS(20)
     ;,IOPTS(19),IOPTS(17),IOCNSF   ,IOCNST   ,IOCRITRAP
     ;,IODIM    ,IOEQT    ,IOFLLI   ,IOFLSAT  ,IOFMLF
     ;,IOFMLT   ,IOINV    ,IOPART   ,IOPINITC
     ;,IOPINITH           ,IOPTS(30),IORTS    ,IOPTS(28),IOPTS(29)
     ;,IOSUCHUM ,IOTRLI   ,IOTRS    ,IOVAR    ,IPROCESS
     ;,ISOT     ,LMXNDL   ,MAINF
     ;,NBAND    ,NBANDCOV ,NBLCVP
     ;,NDEVS
     ;,NFLAGS   ,NFNL     ,NFNT     ,NINT     ,NMAXF
     &,NMAXT
     ;,NOPTS    ,NPAR     ,NPARALG  ,NPARF    ,NPARFPRG
     ;,NPARPRG  ,NPBFL    ,NPBTP    ,NSTAT    ,NTDMT
     ;,NTYPAR   ,NUMEL    ,NUMNP    ,NUMTIT   ,NUMTNOD  ,NUMTOBS
     ;,NWRITE   ,FOBJ_WGT ,IFLAGS   ,IOLG_PAR ,IOPTS
     ;,IOWRITE  ,IPAR_DIR ,IPAR_INV ,NZONE_PAR,PAR_DIR
     ;,PAR_INV  ,PAR_WGT  ,FILENAME
     ;!NUEVOS
     ;,IODENS_INI      ,IODIRECT    ,IOSPARSE  ,ITPTVAR ,LINMET
     &,IOCONSRC
     ;,IOINV_GS ,MXGRPZN   ,MXLINCMB ,NGROUP_ZN ,IO_KG_GS ,IOPT_GS)


C------------------------- Provisionally sets the COMMON variables to use
C------------------------- in the part of the code that uses the COMMON.FOR

       CALL PROVISIONAL
     ; (NPARALG  ,NSTAT    ,NTYPAR   ,NWRITE
     ; ,FOBJ_WGT ,IOLG_PAR ,IOWRITE  ,IPAR_DIR
     ; ,IPAR_INV ,NZONE_PAR,PAR_DIR  ,PAR_INV  ,PAR_WGT)

       IDNFNT=IDIMFNT          !PROVISIONAL.................
       IDIMTRAC=ISOT           !PROVISIONAL..............

C--------------- Sets some auxiliar variables used to partition

       NZDMT = NZONE_PAR(16)
       NZTRA = NZONE_PAR(1)
       NZSTG = NZONE_PAR(2)
       NZARR = NZONE_PAR(3)
       NZCHP = NZONE_PAR(4)
       NZQQP = NZONE_PAR(5)
       NZALF = NZONE_PAR(6)
       NZFOD = NZONE_PAR(11)
       NZPRG = NZONE_PAR(14)
       IOSMFL=IOPTS(28)
       IOSMTP=IOPTS(29)

       IOBALDC=IOPTS(22)
       IOBALDH=IOPTS(21)
       IOBALC=IOPTS(18)
       IOBALH=IOPTS(17)
       IOBALGC=IOPTS(20)
       IOBALGH=IOPTS(19)

C--------------- Initialization of integer variables used to dimension


       CALL INIT_INTEG
     ;(IDIMAFLU ,IDIMBB   ,IDIMCOV  ,IDIMCROSS_GS ,IDIMDFLU   ,IDIMDQ
     ;,IDIMDTRA ,IDIMFNT  ,IDIMHESS ,IDIMQ    ,IDIMVAR_GS ,IOCNSF
     ;,IOCNST   ,IODENS_INI,IODIM   ,IOEQT    ,IOFLLI     ,IOFLSAT
     ;,IOINV    ,IOPART   ,IORTS    ,IOTRLI   ,IOTRS
     ;,IPARTRA  ,ISOT     ,LMXNDL   ,MAINF
     ;,MAXNEOP
     ;
     ;,NBAND    ,NBAND1   ,NBAND2   ,NBANDCOV   ,NDEVS
     ;,NFNT     ,NINT
     ;,NPAR     ,NPAREL   ,NPARNP   ,NTYPAR     ,NUMEL
     ;,NUMNP    ,NUMTIT   ,NUMTNOD  ,NUMTOBS  ,NZALF
     ;,NZARR    ,NZCHP    ,NZONE_PAR(13) ,NZONE_PAR(12),NZONE_PAR(9)
     ;,NZDMT    ,NZONE_PAR(7),NZFOD    ,NZPAR ,NZONE_PAR(10)
     ;,NZPRG    ,NZQQP    ,NZSTG    ,NZTRA    ,INORPAR
     ;,IDIMWORK ,IOSMTP   ,IOSMFL   ,NPBTP    ,NPBFL
     ;!nuevos
     ;,IDIMCFLU ,IDIMATRA
     ;,IOSPARSE ,IAFLUDSC_ROWS, IAFLUDSC_COLS,IATRADSC_ROWS
     ;,IATRADSC_COLS  ,IA_COUPLED_DSC_ROWS ,IA_COUPLED_DSC_COLS
     ;,ICAN_CN  ,IALW_CN  ,IPAR_DIR ,NPARALG, ITYPAFLU, ITYPBFLU
     ;,ITYPCFLU ,ITYPDFLU ,ITYPATRA ,ITYPDTRA
     ;,ITYPBTRA ,ITYPFLUDSC,ITYPTRADSC,ITYPCOUPLDSC, ITYPDERIV
     ;,IDIMDENS ,IOPTS(31),IDIMGRAVEL,NZONE_PAR(18),NZONE_PAR(19)
     &,IDIMDERH,IDIMDERC
     ; ,IDIMIVARIO_GS  ,IDIMWGT        ,IDIMZONPP_GS   ,IOINV_GS
     ; ,MXCLOSE_GS     ,MXDISC_GS      ,MXGRPZN        ,MXKRIG_GS
     ; ,MXLINCMB       ,MXMEASPP_GS    ,MXNPP_GS       ,MXNPRIM_GS
     ; ,MXNST_GS       ,MXNVAR_GS      ,MXNZON_GS      ,MXROT_GS
     ; ,MXSAM_GS       ,MXSB_GS        ,MXSC_GS        ,MXVGM_GS
     ; ,MXZONPP_GS     ,NGROUP_ZN      ,IOPT_GS        ,IO_KG_GS
     ; ,IDIMDATASC_GS)

C--------------- Initializes index counter of array RV. It accounts for all
C--------------- real arrays

       NOLD=1

C--------------- Real arrays related to model parameters

       IDPARZ=NOLD                                      !PARZ
       IDPARM=IDPARZ+NZPAR                              !PARM
       NOLD=IDPARM+NZPAR
       IDSTPAR=NOLD                                     !STPAR
       IF (IOINV.GT.0) NOLD=IDSTPAR+NZPAR

       IDCFPAREL=NOLD                                   !CFPAREL
       IDCFPARNP=IDCFPAREL+NPAREL*NUMEL                 !CFPARNP
       NOLD=  IDCFPARNP+NPARNP*NUMNP

C--------------- Element or nodal values of flow and transport parameters

       IDPARNP=NOLD                                     !IDPARNP
       IDPAREL=IDPARNP+NUMNP*NPPNP*
     ;         MAX(1,IOPTS(29)*NPBTP,IOPTS(28)*NPBFL)   !IDPAREL
       NOLD=IDPAREL+NUMEL*NPPEL*MAX(1,IOPTS(29)*NPBTP,IOPTS(28)*NPBFL)

C--------------- A Priori Covariance matrix of parameters

       IDCOVPAR=NOLD                                      !COVPAR
       IF (IOINV.GT.0) NOLD=IDCOVPAR+NPAR*(NPAR+1)/2

       IDSOURCE=NOLD                                      !SOURCE
       IF (IOSMFL.NE.0 .OR. IOSMTP.NE.0) NOLD=IDSOURCE+2*NUMNP

C--------------- Partition of flow equation arrays

       IF (IOEQT.NE.2) THEN

C--------------- Auxiliar variable to reserve space in some flow arrays
C--------------- (AFLU, etc)

          NPBFLDIM=MAX(1,IOSMFL*NPBFL)

C--------------- Gravity flow


          IF (IOFLLI.EQ.1 .OR. IODENS_INI.EQ.1) THEN
             IDGRAV=NOLD+NUMNP                  !GRAVity
             IDGRAVEL=IDGRAV+3                  !GRAVEL
             NOLD=IDGRAVEL+3*NUMEL
          ELSE
             IDGRAV=NOLD
             IDGRAVEL=NOLD
          ENDIF

C--------------- Discretization arrays

          IDAFLU=NOLD                                   !AFLU
          IDBFLU=IDAFLU+NUMEL*IDIMAFLU* NPBFLDIM         !BFLU
          IDHCALIT = IDBFLU+NUMNP*NPBFLDIM
          IDHCALAN=IDHCALIT+NUMNP*NPBFLDIM                     !HCALAN
          IDDFLU=IDHCALAN+NUMNP*NPBFLDIM                       !DFLU
          IDALFA=IDDFLU+NUMEL*IDIMDFLU*NPBFLDIM                ! ALFA
          IDHAUX1=IDALFA+NUMNP*NPBFLDIM              !HAUX1
          IDHAUX2=IDHAUX1+NUMNP*NPBFLDIM                     !HAUX2
          NOLD=IDHAUX2+NUMNP*NPBFLDIM


C-------------- density dependency matrixes
          IF (IODENS_INI.EQ.1) THEN
             IDCFLU=NOLD
             NOLD=IDCFLU + NUMEL*IDIMCFLU
          ELSE
             IDCFLU=NOLD
          ENDIF

C---------------Flow System Matrix
          IDAFLUDSC=NOLD
          NOLD = NOLD + IAFLUDSC_ROWS * IAFLUDSC_COLS * NPBFLDIM

c-------------- factorized flow matrix for watsolv
          IDAFLUDSCF = NOLD
          IF(IOSPARSE.EQ.1 .AND. IPAR_DIR(22).EQ.1)
     ;       NOLD = NOLD + IPAR_DIR(24)*IAFLUDSC_COLS


C---------------solution vector
          IDSOLUTION = NOLD
          NOLD = NOLD + (ICAN_CN+1)*NUMNP

C--------------- Inverse problem related arrays

          IDDERH=NOLD             !DERH
          IF (IOINV.GT.0) THEN

C--------------- Computes the number of calibrated parameters.
C--------------- When in variable density problems NPARF=NPAR.
C--------------- However, the value of NPARF is not changed yet because
C--------------- it is used to check input data (i.e., PAR file).
C--------------- The NPARF=NPAR change is made at the end of ENTDAT.

              NPARAUX = 0
              IF (IODENS_INI.GT.0) THEN

                  NPARAUX = NPAR

              ELSE

                  NPARAUX = NPARF

              END IF

              NOLD=IDDERH+IDIMDERH*NPARAUX*NUMNP*NPBFLDIM

          END IF !IOINV.GT.0
C--------------- Flow mass balance

          IDBM_ND_FL = NOLD                            ! BM_ND_FL
          IF (IOBALDH.NE.0) NOLD = NOLD + NUMNP*8*2

          IDBM_ZN_FL=NOLD                              ! BM_ZN_FL

          IF (IOBALH.NE.0 .OR. IOBALGH.NE.0) THEN

              NOLD=NOLD+( NZONE_PAR(2)*(IODENS_INI+1)
     &                   +NZONE_PAR(3)
     &                   +NZONE_PAR(4)
     &                   +NZONE_PAR(5)
     &                   +NZONE_PAR(6)
     &                   +NZONE_PAR(13)*IODENS_INI)*2*NPBFL

        END IF !IOBALH.NE.0 .OR. IOBALGH.NE.0

       ELSE                                           ! IOEQT .EQ. 2

          IDGRAV=NOLD
          IDGRAVEL=NOLD
          IDAFLU=NOLD
          IDAFLUDSC=NOLD
          IDBFLU=NOLD
          IDHCALIT = NOLD
          IDDFLU=NOLD
          IDALFA=NOLD
          IDHAUX1=NOLD
          IDHAUX2=NOLD
          IDDERH=NOLD
          IDBM_ND_FL=NOLD
          IDBM_ZN_FL=NOLD
          IDLEAK_CONT = NOLD

          !density dependency arrays
          IDCFLU= NOLD

       END IF                                 !END ANSWER IOEQT .NE. 2


C--------------- Inverse problem general arrays

       IDGRAD=NOLD
       IF (IOINV.GT.0) THEN
          IDHESSAUX=IDGRAD+NPAR                     ! HESSAUX
          IDHESS=IDHESSAUX+IDIMHESS                 ! HESS
          IDDLT_PAR=IDHESS+IDIMHESS                 ! DLT_PAR
          IDPARAUX=IDDLT_PAR+NPAR
          IDCOVINV=IDPARAUX+NPAR                    ! COVINV
          IDVJAC=IDCOVINV+IDIMCOV                   ! VJAC
          IDPARC=IDVJAC+NPAR*NUMTOBS                ! PARC
          IDPARM=IDPARC+NPAR                        ! PARM
          IDWGT_PAR=IDPARM+NPAR                     ! WGT_PAR
          IDDERIV=IDWGT_PAR+IDIMWGT*NZPAR           ! DERIV
          IDWGT_UNK=IDDERIV+NPAR                    ! WGT_UNK
          IDPARGOOD=IDWGT_UNK+NPAR                  ! PARGOOD
          NOLD=IDPARGOOD+NZPAR
       ELSE
          IDGRAD=NOLD
          IDHESSAUX=NOLD
          IDHESS=NOLD
          IDDLT_PAR=NOLD
          IDPARAUX=NOLD
          IDCOVINV=NOLD
          IDVJAC=NOLD
          IDPARC=NOLD
          IDPARM=NOLD
          IDWGT_PAR=NOLD
          IDDERIV=NOLD
          IDWGT_UNK=NOLD
          IDPARGOOD=NOLD
       END IF

C--------------- Coordinates and finite element integral arrays

       IDCOORD=NOLD                            !COORD
       IDBIBI=IDCOORD+3*NUMNP                  !BIBI
       IDVOLNOD=IDBIBI+IDIMBB*NUMEL              !VOLNOD
       IDAREA=IDVOLNOD+NUMNP                     !AREA
       IF (IOINV.NE.0) THEN
          NOLD=IDAREA+NUMEL
       ELSE
          NOLD=IDAREA+NUMEL*LMXNDL
       END IF

       IDGRDFF=NOLD                            !GRDFF
       NOLD=IDGRDFF+IODIM*LMXNDL*NUMEL

C--------------- Auxiliar arrays

       IDWORK=NOLD                             !WORK
       NOLD=IDWORK+IDIMWORK

C--------------- Real time arrays

       IF (IOTRS+IORTS.NE.0) THEN
          IDFNT=NOLD                           !FNT
          IDTIME=IDFNT+IDIMFNT*(NINT+1)           !TIME
          NOLD=IDTIME+NINT
       ELSE
          IDFNT=NOLD
          IDTIME=NOLD
       END IF

       IDCAUDAL=NOLD                          !CAUDAL
       NOLD=IDCAUDAL+NUMNP*NPBFLDIM

       IDCONCFLOW = NOLD                      !CONCENTRATION SOURCES
       IF (IOCONSRC.EQ.1) THEN
           NOLD = NOLD + NUMNP
       END IF !IOCONSRC.EQ.1

C--------------- Darcy velocity and related variables
       IDVD=NOLD                              !VD
       IDQXYZ=IDVD+IODIM*NUMEL*NPBFLDIM       !QXYZ
       IDXNORVD=IDQXYZ+IDIMQ*NUMEL*NPBFLDIM   !XNORVD
       NOLD=IDXNORVD+NUMEL*NPBFLDIM

      IF (IODENS_INI.EQ.1 .OR. IOTRLI.EQ.1 .OR. IOINV.EQ.3) THEN !provisional
          IDDQDFLU = NOLD
          NOLD = NOLD + NUMEL*LMXNDL*LMXNDL   !DQDFLU

          IDDQDTRA = NOLD
          NOLD = NOLD + NUMEL*LMXNDL*LMXNDL   !DQDTRA

          IDDVDH = NOLD
          NOLD = NOLD + NUMEL*IODIM*LMXNDL    !DVDH

          IDDVDC = NOLD
          IF (IODENS_INI.EQ.1) NOLD = NOLD + NUMEL*IODIM*LMXNDL  !DVDC

      ELSE

          IDDQDFLU = NOLD
          IDDQDTRA = NOLD
          IDDVDH = NOLD
          IDDVDC = NOLD

      END IF

C---------------Consistent velocity arrays
       IF (IODENS_INI.EQ.1) THEN
          IDPOINTWEIGHT=NOLD
          IDGRADLOC = IDPOINTWEIGHT + NTYPEL * MAXPG
          IDBUOYANCY = IDGRADLOC  + IODIM*LMXNDL*MAXPG
          IDDBUOYANCY = IDBUOYANCY + IODIM*LMXNDL*NUMEL
          IDGP_COORD = IDDBUOYANCY +IODIM*LMXNDL*LMXNDL*NUMEL
          NOLD=IDGP_COORD+8*6*3
       ELSE
          IDPOINTWEIGHT=NOLD
          IDGRADLOC = NOLD
          IDBUOYANCY  = NOLD
          IDDBUOYANCY = NOLD
          IDGP_COORD= NOLD
       ENDIF

C--------------- Real transport arrays

       IF (IOEQT.NE.1) THEN

C--------------- Auxiliar variable to reserve space in some transport arrays
C--------------- (ATRA, etc)

          NPBTPDIM=MAX(1,IOSMTP*NPBTP)

C--------------- Water content of a given element

          IDWATVOL=NOLD
          IF (IOPTS(31).LE.1) THEN        ! WATVOL by elements
             NOLD=IDWATVOL+NUMEL*3*NPBTPDIM     !In time k, k+th, k+1-th
          ELSE                            ! WATVOL by nodes
             NOLD=IDWATVOL+LMXNDL*NUMEL*3*NPBTPDIM
          ENDIF

C--------------- Derivatives of water content w. r. t. concentratrion.

          IDDWDH = NOLD

          IF (IOTRLI.NE.0 .OR. IODENS_INI.EQ.1) THEN

              IF (IOPTS(31).LE.1) THEN        ! WATVOL by elements
                  NOLD = IDDWDH + NUMEL*NPBTPDIM
              ELSE                            ! WATVOL by nodes
                  NOLD = IDDWDH + 2*LMXNDL*NUMEL*NPBTPDIM
              END IF

          END IF

C--------------- Aquifer thickness

          IDACTH=NOLD                               !ACTH
          NOLD=NOLD+NUMEL

C--------------- Discretization arrays

          IDATRA=NOLD                                !ATRA
          IDBTRA=IDATRA+NUMEL*IDIMATRA*NPBTPDIM
          IDCCALIT=IDBTRA+MAXNEOP*NPBTPDIM           !CCALIT
          NOLD=IDCCALIT+NUMNP*NPBTPDIM
          IDCCALAN=NOLD                              !CCALAN
          IF (IORTS.NE.0) NOLD=NOLD+NUMNP*NPBTPDIM
          IDDTRA=NOLD                                !DTRA
          IF (IORTS.NE.0) THEN
             NOLD=NOLD+NUMEL*IDIMDTRA*NPBTPDIM
          ENDIF

          IDCAUX1=NOLD                               !CAUX1
          IDCAUX2=IDCAUX1+NUMNP*NPBTPDIM             !CAUX2
          NOLD=IDCAUX2+NUMNP*NPBTPDIM

C-----------ATRADSC

          IDATRADSC = NOLD
          NOLD = IDATRADSC + IATRADSC_ROWS*IATRADSC_COLS*NPBTPDIM

c-------------- factorized flow matrix for watsolv
          IDATRADSCF = NOLD
          IF(IOSPARSE.EQ.1 .AND. IPAR_DIR(22).EQ.1)
     ;       NOLD = NOLD + IATRADSC_COLS * IPAR_DIR(24)

          IDDERC=NOLD                               !DERC
          IF (IOINV.GT.0) THEN

              NOLD = IDDERC+IDIMDERC*NPAR*NUMNP*NPBTPDIM
              IDDVDP=NOLD                            !DVDP

             IF (NPARF.NE.0) NOLD=IDDVDP+IODIM*NPARF*NUMEL

          ELSE

             IDDVDP=NOLD

          ENDIF

C--------------- Transport mass balance

          IDBM_ND_TT=NOLD                                 ! BM_ND_TT

          IF (IOBALDC.NE.0) NOLD=NOLD+NUMNP*12*2

          IDBM_ZN_TT=NOLD                                 ! BM_ZN_TT
          IF (IOBALC.NE.0 .OR. IOBALGC.NE.0)
     ;       NOLD=NOLD+( NZONE_PAR(10)+NZONE_PAR(11) +          !NZPOR+NZFOD
     ;                   NZONE_PAR(17)+NZONE_PAR(13)*6 +        !NZZOR+NZCOE*6
     &                   NZONE_PAR(18)+NZONE_PAR(16))*2  *NPBTP !NZCLK+NZMATDIF

       ELSE                                               ! IOEQT .EQ. 1
          IDWATVOL=NOLD
          IDACTH=NOLD
          IDCJAC=NOLD
          IDDERC=NOLD
          IDCAUX1=NOLD
          IDCAUX2=NOLD
          IDDVDP=NOLD
          IDATRA=NOLD
          IDATRADSC=NOLD
          IDATRADSCF=NOLD
          IDBTRA=NOLD
          IDCCALIT=NOLD
          IDDTRA=NOLD
          IDCCALAN=NOLD
          IDBM_ZN_TT=NOLD                                 ! BM_ZN_TT
          IDBM_ND_TT=NOLD                                 ! BM_ND_TT
          IDDWDH = NOLD
       END IF                                  ! END ANSWER IOEQT .NE. 1

C--------------- Observation arrays

       IF (IOINV.GT.0.OR.IOWRITE(5).NE.0.OR.IOWRITE(6).NE.0
     ;               .OR.IOWRITE(11).NE.0.OR.IOWRITE(12).NE.0) THEN

C--------------- Observation values, calculated values and derivatives

          IDVOBS=NOLD                               !VOBS
          IDVOBSC=IDVOBS+NUMTOBS                    !VOBSC
          IDDVOBS=IDVOBSC+NUMTOBS+NDEVS             !DVOBS
          NOLD=IDDVOBS
          IF (IOINV.GT.0) NOLD=IDDVOBS+3*NPAR

C--------------- Observation times and integration times

          IDTOBS=NOLD                               !TOBS
          IDTIT=IDTOBS+2*NUMTOBS                    !TIT
          NOLD=IDTIT+NUMTIT

C--------------- Basic units

          IDBUDAT=NOLD                               !BUDAT
          IDEXTNBU=IDBUDAT+4*NUMTNOD                 !EXTNBU
          NOLD=IDEXTNBU+NUMTNOD

C--------------- Weights

          IDWTOBSBU=NOLD                             !WTOBSBU
          IDWTOBSU=IDWTOBSBU+NUMTNOD                 !WTOBSU
          IDWTOBSN=IDWTOBSU+NUMTNOD                  !WTOBSN
          IDWTOBST=IDWTOBSN+NUMTNOD                  !WTOBST
          NOLD=IDWTOBST+NUMTOBS

       ELSE

          IDVOBSC=NOLD                              !VOBSC
          IDVOBS=NOLD                               !VOBS
          IDDVOBS=NOLD                              !DVOBS
          IDTOBS=NOLD                               !TOBS
          IDTIT=NOLD                                !TIT
          IDBUDAT=NOLD                              !BUDAT
          IDEXTNBU=NOLD                             !EXTNBU
          IDWTOBSBU=NOLD                            !WTOBSBU
          IDWTOBSU=NOLD                             !WTOBSU
          IDWTOBSN=NOLD                             !WTOBSN
          IDWTOBST=NOLD                             !WTOBST

       END IF



C--------------- Non linear arrays



C--------------- FLOW

       IF (IOFLLI.NE.0 .OR. IODENS_INI.EQ.1) THEN
          IDHBASE=NOLD                                       !HBASE
          IDHPREV1=IDHBASE+NUMEL                             !HPREV1
          IDHPREV2=IDHPREV1+NUMNP                            !HPREV2
          IDDFLUDFLU=IDHPREV2+NUMNP                          !DFLUDFLU
          IDDBFLUDFLU=IDDFLUDFLU+LMXNDL*LMXNDL*NUMEL*NPBFLDIM !DBFLUDFLU
          NOLD=IDDBFLUDFLU +NUMNP*NPBFLDIM

          IDDBFLUDTRA=NOLD                             !DBFLUDTRA

          IF (IODENS_INI.EQ.1 .AND. ICAN_CN.EQ.1) THEN

              NOLD = IDDBFLUDTRA+NUMNP
              IDDERVISC = NOLD
              NOLD = NOLD + NUMEL

          ELSE

              IDDERVISC = NOLD

          END IF

C---------------  Stores nodal variables at k time



          IDDPARELDH=NOLD                     ! DPARELDH
          IDDNODALRH=IDDPARELDH+NUMEL*NPPEL
     ;              *MAX(1,IOPTS(29)*NPBTP,IOPTS(28)*NPBFL)
          NOLD=IDDNODALRH+4*NUMNP
       ELSE                                   ! IOFLLI=0
          IDHBASE=NOLD
          IDHPREV1=NOLD
          IDDFLUDFLU=NOLD
          IDHPREV2=NOLD
          IDDNODALRH=NOLD                     ! DNODALRH
          IDDPARELDH = NOLD
          IDDBFLUDFLU=NOLD                !DBFLUDFLU
          IDDBFLUDTRA=NOLD                             !DBFLUDTRA
          IDDERVISC = NOLD
       ENDIF                                  ! Ends IF IOFLLI.NE.0

C--------------- Transport
       IF (IOTRLI.NE.0 .OR. IODENS_INI.EQ.1 .OR. IOINV.EQ.3) THEN

        IDHINI = NOLD
          IDDAT_VD = NOLD + NUMNP
          NOLD = IDDAT_VD + IODIM*IDIMDQ*NUMEL*NPBFLDIM
       ELSE
          IDHINI = NOLD
          IDDAT_VD = NOLD
       ENDIF

       IF (IOTRLI.NE.0 .OR. IODENS_INI.EQ.1)THEN
          IDCBASE=NOLD                                 !CBASE
          IDCPREV1= IDCBASE+NUMEL                      !CPREV1
          IDCPREV2= IDCPREV1+NUMNP                     !CPREV2
          IDDTRADTRA=IDCPREV2+NUMNP                    !DTRADTRA
          NOLD=IDDTRADTRA+LMXNDL*LMXNDL*NUMEL*NPBTPDIM



      !The derivatives of parameters to concentration
          IDDPARELDC = NOLD

          NOLD = IDDPARELDC+NUMEL*NPPEL
     ;              *MAX(1,IOPTS(29)*NPBTP,IOPTS(28)*NPBFL)


       ELSE

          IDCBASE=NOLD
          IDCPREV1=NOLD
          IDCPREV2=NOLD
          IDDTRADTRA=NOLD

          IDDPARELDC=NOLD
       ENDIF

C---------------Coupled flow and transport
       IDA_COUPLED_DSC=NOLD
       IF (IODENS_INI.EQ.1 .AND. ICAN_CN.GT.0) THEN  !if the coupled system is
                                              !going to be solved

           NOLD = IDA_COUPLED_DSC +  IA_COUPLED_DSC_ROWS
     ;            *IA_COUPLED_DSC_COLS

          !factorized matrix for watsolv
          IDA_COUPLED_DSCF=NOLD
          IF(IOSPARSE.EQ.1 .AND. IPAR_DIR(22).EQ.1) THEN
               NOLD = IDA_COUPLED_DSCF + IA_COUPLED_DSC_COLS
     ;                *IPAR_DIR(24)*2
          ENDIF

          !the derivatives matrices
          IDDFLUDTRA = NOLD
          IDDTRADFLU = IDDFLUDTRA+LMXNDL*LMXNDL*NUMEL
          NOLD = IDDTRADFLU + LMXNDL*LMXNDL*NUMEL

          IDBCOUPLED = NOLD


          NOLD = IDBCOUPLED + 2*NUMNP

       ELSE    !Newton's method for coupled system not going to be used
          IDDFLUDTRA = NOLD
          IF (IODENS_INI.EQ.1) NOLD = NOLD +LMXNDL*NUMEL !provisional. only necesary
          IDDTRADFLU = NOLD                          !if tpte uses newton
          IDA_COUPLED_DSCF=NOLD
          IDBCOUPLED = NOLD
       ENDIF


C--------------- Common non linear variables (flow or transport)

      IF (IOFLLI.NE.0. OR .IOTRLI.NE.0)THEN
        IDDTMXDS=NOLD                                       !DTMXDS
        IDPARACD=IDDTMXDS+NINT                              !PARACD
        NOLD=IDPARACD+3*NFNL
      ELSE
        IDDTMXDS=NOLD
        IDPARACD=NOLD
      ENDIF

C------------density is always dimensioned.

      IDDENSITY = NOLD
      NOLD = IDDENSITY + MAX(NUMEL,NUMNP) !Por exceso

C------------density dependency arrays

      IF (IODENS_INI.EQ.1) THEN
         IDVISCOSITY = NOLD
         IDDELTAITER  = IDVISCOSITY + NUMEL
         NOLD = IDDELTAITER + NUMNP
      ELSE
         IDVISCOSITY = NOLD
         IDDELTAITER = NOLD
      ENDIF


C--------------- Geoestatistical inverse problem variables

       IF (IOINV_GS.NE.0) THEN

                             ! Variogram definition arrays and rotation matrices

         IDVARIO_GS=NOLD
         IDROTMAT_GS=IDVARIO_GS+IDIMIVARIO_GS*8*NGROUP_ZN
         NOLD=IDROTMAT_GS+3*3*MXROT_GS

        ! Sampling and pilot point locations and variables values. Zonal ext. att.

         IDPOSMEAS_GS=NOLD                                          !POSMEAS_GS
         IDVMEAS_GS=IDPOSMEAS_GS+MXMEASPP_GS*3*NGROUP_ZN*2          !VMEAS
         IDEXDRZN_GS=IDVMEAS_GS+MXMEASPP_GS*IDIMVAR_GS*NGROUP_ZN*2  !SECTZN_GS

                                     ! Zonal positions and discretization points

         IDPOSDIS_GS=IDEXDRZN_GS+MXNZON_GS*4*NGROUP_ZN              !POSDIS_GS
         IDPOSDISAUX_GS=IDPOSDIS_GS+MXDISC_GS*3*MXNZON_GS*NGROUP_ZN !POSDAUX_GS
         IDPOSZN_GS=IDPOSDISAUX_GS+MXDISC_GS*3                      !POSZN_GS
         IDDIVZN_GS=IDPOSZN_GS+MXNZON_GS*3*NGROUP_ZN                !DIVZN_GS
         NOLD=IDDIVZN_GS+MXNZON_GS

                                   ! Var. Statistics, search and trimming limits

         IDVSTATS_GS=NOLD                                           !VSTATS_GS
         IDTRIM_GS=IDVSTATS_GS+MXNVAR_GS*4*NGROUP_ZN                !TRIM_GS
         IDSEARCH_GS=IDTRIM_GS+8*NGROUP_ZN              !SEARCH_GS
         NOLD=IDSEARCH_GS+12*NGROUP_ZN

                                           ! Estimation, weights and covariances

         IDESTKRIG_GS=NOLD                                          !ESTKRIG_GS
         IDZNWGT_GS=IDESTKRIG_GS+MXZONPP_GS                         !ZNWGT_GS
         !IDCROSSCOV_GS=IDZNWGT_GS+MXNPRIM_GS*MXNZON_GS              !CROSSCOV_GS
         !IDESTPARZ_GS=IDCROSSCOV_GS+IDIMCROSS_GS*MXZONPP_GS*2       !ESTPARZ_GS 
         IDESTPARZ_GS = IDZNWGT_GS+MXNPRIM_GS*MXNZON_GS              !PARZ_GS
         NOLD=IDESTPARZ_GS+MXNZON_GS*MXSC_GS*NGROUP_ZN

                                           ! Conditional simulations

         IDDATASC_GS=NOLD                                             !DATA_SC
         NOLD=IDDATASC_GS+IDIMDATASC_GS*18

                                           ! GSLIB requirements

         IDSUPBL_GS=NOLD                                            !SUPBL_GS
         IDKRISYS_GS=IDSUPBL_GS+6*NGROUP_ZN                         !KRISYS_GS
         IDKRISOL_GS=IDKRISYS_GS+MXKRIG_GS*MXKRIG_GS                !KRISOL_GS
         IDCLOSESAM_GS=IDKRISOL_GS+MXKRIG_GS*3                      !CLOSESAM_GS
         IDKRIGAUX_GS=IDCLOSESAM_GS+MXCLOSE_GS*2                    !KRIGAUX_GS
         NOLD=IDKRIGAUX_GS+MXSAM_GS*8

                                           ! Maximum/minimum group coordinates
         IDCOORDGR_GS=NOLD
         NOLD=IDCOORDGR_GS+6*NGROUP_ZN

                                           ! PARC at all groups

         IDPARC_GS=NOLD
         NOLD=IDPARC_GS+MXNPP_GS*NGROUP_ZN

                                           ! Part of COVPAR related to a given group of zones

         IDCOVPAR_GR=NOLD
         NOLD=IDCOVPAR_GR+MXNPP_GS*(MXNPP_GS+1)/2

       ELSE

         IDVARIO_GS=NOLD
         IDROTMAT_GS=NOLD
         IDPOSMEAS_GS=NOLD
         IDVMEAS_GS=NOLD
         IDEXDRZN_GS=NOLD
         IDPOSDIS_GS=NOLD
         IDPOSDISAUX_GS=NOLD
         IDPOSZN_GS=NOLD
         IDDIVZN_GS=NOLD
         IDVSTATS_GS=NOLD
         IDTRIM_GS=NOLD
         IDSEARCH_GS=NOLD
         IDESTKRIG_GS=NOLD
         IDZNWGT_GS=NOLD
         IDCROSSCOV_GS=NOLD
         IDESTPARZ_GS=NOLD
         IDDATASC_GS=NOLD
         IDSUPBL_GS=NOLD
         IDKRISYS_GS=NOLD
         IDKRISOL_GS=NOLD
         IDCLOSESAM_GS=NOLD
         IDKRIGAUX_GS=NOLD
         IDCOORDGR_GS=NOLD
         IDPARC_GS=NOLD
         IDCOVPAR_GR=NOLD

       END IF

       IDDTPREVINV = NOLD
       NOLD = NOLD + NINT

C--------------- Statistical analysis arrays

       IF (IOVAR.NE.0) THEN
          IDDEVICESTAT=NOLD
          IDEIGENVEC=IDDEVICESTAT+7*NDEVS
          IDRESID=IDEIGENVEC+MAX(NPAR*NPAR,2*NPAR)
          IDRESIDPAR = IDRESID + NUMTOBS
          NOLD=IDRESID+2*NPAR
       ELSE
          IDDEVICESTAT=NOLD
          IDEIGENVEC=NOLD
          IDRESID=NOLD
          IDRESIDPAR = NOLD
       END IF


C--------------- String arrays

       KOLD=1
       IDDEVNAME=KOLD                 !DEVNAME
       KOLD=IDDEVNAME+10*NDEVS

C--------------- Integer arrays

       IOLD=1

       IF (IOEQT.NE.2) THEN

C--------------- Flow boundary conditions

          IDIBCOD=IOLD                              !IBCOD

C--------------- Transmisivity

          IDISOZ=IDIBCOD+NUMNP*NPBFL                       !ISOZ
          IOLD=IDISOZ+NZTRA
       ELSE
          IDIBCOD=IOLD
          IDISOZ=IOLD
       END IF

C--------------- Parameter arrays

       IDNFTPAR=IOLD                              !NFTPAR
       IDNFNLPAR=IDNFTPAR+NZPAR                   !NFNLPAR
       IDIVPAR=IDNFNLPAR+NZPAR                    !IVPAR
       IOLD=IDIVPAR+4*NZPAR

       IF (IOINV.GT.0) THEN
         IDIDMBLCVP=IOLD                           !IDMBLCVP
         IDIPOS=IDIDMBLCVP+NPAR                    !IPOS
         IDIPNT_PAR=IDIPOS+NPAR                    !IPNT_PAR
         IOLD=IDIPNT_PAR+IDIMWGT*NZPAR
       ELSE
         IDIDMBLCVP=IOLD
         IDIPOS=IOLD
         IDIPNT_PAR=IOLD
       END IF

       IDLXPAREL=IOLD                             !LXPAREL
       IDIXPARNP=IDLXPAREL+NPAREL*NUMEL*MAX(NPBFL,NPBTP)    !IXPARNP
       IOLD=IDIXPARNP+NPARNP*NUMNP*MAX(NPBFL,NPBTP)

C--------------- Discretization arrays

       IDKXX=IOLD                              !KXX
       IDLTYPE=IOLD+LMXNDL*NUMEL                 !LTYPE
       IDLNNDEL=IDLTYPE+NUMEL                       !LNNDEL
       IDLDIM=IDLNNDEL+NUMEL                     !LDIM
       IOLD=IDLDIM+NUMEL

C---------------  Sparse storage adjacency arrays


       IF (IOSPARSE.EQ.1) THEN          !if we use sparse storage

          IF (MAX(IAFLUDSC_ROWS,IATRADSC_ROWS).GT.0) THEN

             IDIAD_S=IOLD
             IDIADD_S=IDIAD_S + NUMNP*IPAR_DIR(19)
             IDIADN_S=IDIADD_S + NUMNP
             IOLD = IDIADN_S + NUMNP

             IF (IPAR_DIR(22).EQ.1) THEN  !if we use matrix preconditioning

                IDIAFD_S = IOLD
                IDIAFDD_S = IDIAFD_S + IPAR_DIR(24)*NUMNP
                IDIAFDN_S = IDIAFDD_S + NUMNP
                IOLD = IDIAFDN_S + NUMNP

             ELSE
                IDIAFD_S = IOLD
                IDIAFDD_S = IOLD
                IDIAFDN_S = IOLD
             ENDIF
          ELSE
                !if we solve either flow or transport separately
             IDIAD_S=IOLD       !adjacency arrays for either flow or tpt
             IDIADD_S=IOLD
             IDIADN_S=IOLD

             IDIAFD_S=IOLD      !adjacency array for factorized matrix of flow or tpt
             IDIAFDD_S=IOLD
             IDIAFDN_S=IOLD
          ENDIF


          IF (IA_COUPLED_DSC_ROWS.NE.0)THEN !if we solve coupled system

             IDIAD_D=IOLD
             IDIADD_D=IDIAD_D + 4*NUMNP*IPAR_DIR(19)
             IDIADN_D=IDIADD_D + 2*NUMNP
             IOLD = IDIADN_D + 2*NUMNP

             IF (IPAR_DIR(22).EQ.1) THEN !if we use matrix preconditioning

                IDIAFD_D = IOLD
                IDIAFDD_D = IDIAFD_D + 4*NUMNP*IPAR_DIR(24)
                IDIAFDN_D = IDIAFDD_D + 2*NUMNP
                IOLD = IDIAFDN_D + 2*NUMNP
             ELSE
                IDIAFD_D = IOLD
                IDIAFDD_D = IOLD
                IDIAFDN_D = IOLD
             ENDIF

          ELSE    !if we solve the coupled system
            IDIAD_D=IOLD       !adjacency arrays for coupled flow & tpt
            IDIADD_D=IOLD
            IDIADN_D=IOLD

            IDIAFD_D=IOLD      !adjacency arrays for factorized matrix
            IDIAFDD_D=IOLD     !of coupled flow &tpt
            IDIAFDN_D=IOLD
          ENDIF

       ELSE
            IDIAD_S=IOLD       !adjacency arrays for either flow or tpt
            IDIADD_S=IOLD
            IDIADN_S=IOLD

            IDIAFD_S=IOLD      !adjacency array for factorized matrix of flow or tpt
            IDIAFDD_S=IOLD
            IDIAFDN_S=IOLD

            IDIAD_D=IOLD       !adjacency arrays for coupled flow & tpt
            IDIADD_D=IOLD
            IDIADN_D=IOLD

            IDIAFD_D=IOLD      !adjacency arrays for factorized matrix
            IDIAFDD_D=IOLD     !of coupled flow &tpt
            IDIAFDN_D=IOLD

       ENDIF       !if we use sparse storage



C--------------- Transport arrays

       IF (IOEQT.NE.1) THEN

C--------------- Transport boundary conditions

          IDIBTCO=IOLD                              !IBTCO
          IOLD=IDIBTCO+NUMNP*NPBTP
       ELSE
          IDIBTCO=IOLD
       END IF

       IDINDPAR=IOLD                               !INDPAR
       IF (IOINV.GT.0) THEN
          IOLD=IOLD+NPAR
       END IF

C--------------- Integer arrays related to measurements and to dimensions
C--------------- NUMTNOD, NUMTOBS and NUMTIT

       IF (IOINV.GT.0.OR.IOWRITE(5).NE.0.OR.IOWRITE(6).NE.0
     ;               .OR.IOWRITE(11).NE.0.OR.IOWRITE(12).NE.0) THEN

          IDIODEVICE=IOLD                             !IODEVICE
          IOLD=IDIODEVICE+10*(NDEVS+1)

C--------------- Observations related variables dimensioned NUMTNOD

C--------------- Basic units, units, nodes and location of elements

          IDINDEXNOD=IOLD                             !INDEXNOD
          IDIOBUTYP=IDINDEXNOD+NUMTNOD                !IOBUTYP
          IDIOCALBU=IDIOBUTYP+NUMTNOD                 !IOCALBU
          IDIOUTYP=IDIOCALBU+NUMTNOD                  !IOUTYP
          IDNBUW=IDIOUTYP+NUMTNOD                     !NBUW
          IDNOBUF=IDNBUW+NUMTNOD                      !NOBUF
          IOLD=IDNOBUF+NUMTNOD+1

C--------------- Observations related variables dimensioned NUMTOBS or
C--------------- NUMTIT. Integration times

          IDNOOBSIT=IOLD                              !NOOBSIT
          IDIOTINT=IDNOOBSIT+NUMTIT                   !IOTINT
          IDMEASTYP=IDIOTINT+NUMTOBS                  !IDMEASTYP
          IOLD=IDMEASTYP+NUMTOBS


       ELSE
          IDIODEVICE=IOLD
          IDINDEXNOD=IOLD
          IDIOBUTYP=IOLD
          IDIOCALBU=IOLD
          IDIOUTYP=IOLD
          IDNBUW=IOLD
          IDNOBUF=IOLD
          IDNOOBSIT=IOLD
          IDIOTINT=IOLD
          IDMEASTYP = IOLD
       END IF

C--------------- Time arrays

       IDKINT=IOLD                                  !KINT
       IF (IOTRS+IORTS.NE.0) IOLD=IOLD+NINT
       IDISOLEQ=IOLD                                  !ISOLEQ
       IF (IOTRS+IORTS.NE.0) IOLD=IOLD+NINT*4

C--------------- Auxiliary array

       IDINTAUX=IOLD                                 !INTAUX
       IOLD=IOLD+MAXNEOP

C--------------- Non linear integer arrays

       IF (IOFLLI.NE.0. OR .IOTRLI.NE.0)THEN
          IDNFNLTIP=IOLD                                       !NFNLTIP
          IDNFNLPRG=IDNFNLTIP+NFNL                              !NFNLPRG
          IOLD=IDNFNLPRG+8*NFNL
       ELSE
          IDNFNLTIP=IOLD
          IDNFNLPRG=IOLD
       ENDIF

C----------------Consistent velocity arrays
       IF (IODENS_INI.EQ.1) THEN
           IDLCOORD = IOLD
           IDIPARTNER = IDLCOORD + 6*3*6
           IOLD = IDIPARTNER +6*3*6
       ELSE
           IDLCOORD=IOLD
           IDIPARTNER=IOLD
       ENDIF


C--------------- Geoestatistical inverse problem integer arrays

       IF (IOINV_GS.NE.0) THEN

                                   ! Integer array defining variograms structure
         IDIVARIO_GS=IOLD                                   !IVARIO_GS
         IOLD=IDIVARIO_GS+IDIMIVARIO_GS*2*NGROUP_ZN

                                                              ! Polynomial trend

         IDIPOLDRIFT_GS=IOLD                                !IPOLDRIFT_GS
         IOLD=IDIPOLDRIFT_GS+9*NGROUP_ZN


                                                    ! Block interpolation arrays

         IDIZN_PP_GS=IOLD                                   !IZN_PP_GS
         IDIZN_NPP_GS=IDIZN_PP_GS+IDIMZONPP_GS              !IZN_NPP_GS
         IDICHECK_GS=IDIZN_NPP_GS+MXNZON_GS                 !ICHECK_GS
         IDIZONMEAS_GS=IDICHECK_GS+MXNPP_GS                 !IZONMEAS_GS
         IOLD=IDIZONMEAS_GS+MXNZON_GS*NGROUP_ZN

                                                 ! Covariance calculation arrays

         IDICROSSCOV_GS=IOLD                                !ICROSSCOV_GS
         IOLD=IDICROSSCOV_GS+MXNPRIM_GS*MXZONPP_GS

                                                 ! GSLIB requirements

         IDIVARCK_GS=IOLD                                   !IVARCK_GS
         IDISUPBL_GS=IDIVARCK_GS+MXKRIG_GS                  !ISUPBL_GS
         IDNUMSB_GS=IDISUPBL_GS+MXSB_GS*8*3                 !NUMSB_GS
         IOLD=IDNUMSB_GS+MXSB_GS
                                        ! Dimension of transmissivity elements
         IDLDIM_GS=IOLD
         IOLD=IDLDIM_GS+NZONE_PAR(1)

         IDIFLAG_SIMUL=IOLD
         IOLD=IDIFLAG_SIMUL+MXZONPP_GS

       ELSE

         IDIVARIO_GS=IOLD
         IDIPOLDRIFT_GS=IOLD
         IDIZN_PP_GS=IOLD
         IDIZN_NPP_GS=IOLD
         IDICHECK_GS=IOLD
         IDICROSSCOV_GS=IOLD
         IDIVARCK_GS=IOLD
         IDISUPBL_GS=IOLD
         IDNUMSB_GS=IOLD
         IDLDIM_GS=IOLD
         IDIFLAG_SIMUL=IOLD
         IDIZONMEAS_GS=IOLD

       END IF

C--------------- Statistical analysis arrays

       IDOBSCLASS=IOLD
       IF (IOVAR.NE.0) THEN
         IDITYPEPAR=IDOBSCLASS+NSTAT*(NDEVS+1)
         IDIZPAR=IDITYPEPAR+2*NPAR
         IOLD = IDIZPAR + NPAR
       ELSE
         IDITYPEPAR=IOLD
         IDIZPAR=IOLD
       END IF


       LASTII=IOLD
       LASTIR=NOLD

C--------------- Writes variables location in arrays RV, IV and KV

       IF (IOPART.NE.0) CALL WRI_PART
     ;(IDPARM     ,IDSTPAR  ,IDCFPAREL
     ;,IDCFPARNP ,IDPARNP    ,IDPAREL  ,IDCOVPAR
     ;,IDSOURCE  ,IDGRAV   ,IDGRAVEL
     ;,IDAFLU    ,IDBFLU     ,IDHCALAN
     ;,IDDFLU    ,IDALFA     ,IDHAUX1  ,IDHAUX2
     ;,IDCFLU    ,IDAFLUDSC,IDAFLUDSCF
     ;,IDDERH    ,IDBM_ND_FL ,IDBM_ZN_FL
     ;!------------------------------------------ Inverse problem general arrays
     ;,IDGRAD    ,IDHESSAUX  ,IDHESS   ,IDDLT_PAR
     ;,IDPARAUX   ,IDCOVINV  ,IDVJAC   ,IDPARC
     ;,IDWGT_PAR  ,IDDERIV  ,IDWGT_UNK ,IDPARGOOD
     ;!---------------------------Coordinates and finite element integral arrays
     ;
     ;,IDCOORD   ,IDBIBI     ,IDVOLNOD  ,IDAREA
     ;,IDGRDFF
     ;!-------------------------------------Auxiliar arrays and Real time arrays
     ;,IDWORK    ,IDFNT      ,IDTIME
     ;,IDCAUDAL  ,IDCONCFLOW
     ;!------------------------------------ Darcy velocity and related variables
     ;,IDVD      ,IDQXYZ     ,IDXNORVD
     ;,IDDQDFLU  ,IDDQDTRA   ,IDDVDH    ,IDDVDC
     ;!-----------------------------------------------Consistent velocity arrays
     ;,IDPOINTWEIGHT,IDGRADLOC
     ;,IDGP_COORD
     ;!---------------------------------------------------------transport arrays
     ;,IDWATVOL  ,IDDWDH     ,IDACTH
     ;,IDDERC    ,IDCAUX1    ,IDCAUX2    ,IDDVDP
     ;,IDATRA    ,IDATRADSC  ,IDATRADSCF ,IDBTRA
     ;,IDDTRA    ,IDCCALAN   ,IDBM_ZN_TT ,IDBM_ND_TT
     ;
     ;!------------------------------------------------------ Observation arrays
     ;,IDVOBSC   ,IDVOBS     ,IDDVOBS  ,IDTOBS
     ;,IDTIT     ,IDBUDAT    ,IDEXTNBU ,IDWTOBSBU
     ;,IDWTOBSU  ,IDWTOBSN   ,IDWTOBST
     ;!---------------------------------------------------------- FLOW nonlinear
     ;,IDHBASE   ,IDHPREV1   ,IDDFLUDFLU,IDDBFLUDFLU
     ;,IDHPREV2  ,IDDNODALRH ,IDDPARELDH
     ;,IDHINI    ,IDDAT_VD
     &,IDCBASE   ,IDCPREV1   ,IDCPREV2
     ;,IDCCALIT  ,IDDTRADTRA
     ;!-----------------------------------------------Coupled flow and transport
     ;,IDA_COUPLED_DSC,IDA_COUPLED_DSCF,IDDFLUDTRA   ,IDDTRADFLU
     ;,IDDPARELDC     ,IDBCOUPLED   ,IDDBFLUDTRA
     ;,IDDERVISC
     ;
     ;!------------------------- Common non linear variables (flow or transport
     ;,IDHCALIT  ,IDDTMXDS     ,IDPARACD
     ;
     ;!------------------------------------------------density dependency arrays
     ;,IDDENSITY ,IDVISCOSITY, IDDELTAITER
     ;
     ;!--------------------------------------------- Statistical analysis arrays
     ;,IDDEVICESTAT,IDEIGENVEC,IDRESID ,IDRESIDPAR
     ;
     ;!----------------------------------------------------------- String arrays
     ;,IDDEVNAME
     ;
     ;!================================================================
     ;!--------------- Integer arrays
     ;,IDIBCOD  ,IDISOZ    ,IDNFTPAR    ,IDNFNLPAR
     ;,IDIVPAR  ,IDIDMBLCVP
     ;,IDLXPAREL   ,IDIXPARNP
     ;
     ;!--------------------------------------------------- Discretization arrays
     ;,IDKXX     ,IDLTYPE   ,IDLNNDEL    ,IDLDIM
     ;,IDIBTCO   ,IDINDPAR
     ;!----------------------------------------------------Sparse storage arrays
     ;,IDIAD_S   ,IDIADD_S  ,IDIADN_S    ,IDIAD_D
     ;,IDIADD_D  ,IDIADN_D  ,IDIAFD_S    ,IDIAFDD_S
     ;,IDIAFDN_S ,IDIAFD_D  ,IDIAFDD_D   ,IDIAFDN_D
     ;
     ;!---------------- Integer arrays related to measurements and to dimensions
     ;!---------------- NUMTNOD, NUMTOBS and NUMTIT
     ;,IDIODEVICE,IDINDEXNOD ,IDIOBUTYP  ,IDIOCALBU
     ;,IDIOUTYP  ,IDNBUW     ,IDNOBUF    ,IDNOOBSIT
     ;,IDIOTINT
     ;
     ;!--------------- Time arrays & Auxiliary array & Non linear integer arrays
     ;,IDKINT    ,IDISOLEQ  ,IDINTAUX, IDNFNLTIP
     ;,IDNFNLPRG ,IDLCOORD,IDIPARTNER
     ;
     ;!-------------------------------------------- Statistical analysis arrays
     ;
     ;,IDOBSCLASS,IDITYPEPAR,IDIZPAR
     ;
     ;!------------------------------------------------------------------others
     ;,LASTII      ,LASTIR       ,MAINF)
      
ccc modificacion 

       OPEN(UNIT = 778, FILE = 'transin_output_optimum.out', 
     ;      STATUS = 'UNKNOWN')
       OPEN(UNIT = 779, FILE = 'pilot_point_locs_transin.out', 
     ;      STATUS = 'UNKNOWN')

ccc fin de la modificacion

C--------------- Checks partition size

       IF (NOLD.GT.IRMAX .OR. IOLD.GT.IIMAX .OR. KOLD.GT.IKMAX) THEN
          WRITE(6,1000) IRMAX,NOLD,IIMAX,IOLD,IKMAX,KOLD
          WRITE(MAINF,1000) IRMAX,NOLD,IIMAX,IOLD,IKMAX,KOLD

 1000     FORMAT(25X,' PARTITION ERROR  ',/,25X,'-------------------',
     ;     /,10X,'IRMAX IS',I9,' IT MUST BE GREATHER OR EQUAL THAN',I9,
     ;     /,10X,'IIMAX IS',I9,' IT MUST BE GREATHER OR EQUAL THAN',I9,
     ;     /,10X,'IKMAX IS',I9,' IT MUST BE GREATHER OR EQUAL THAN',I9,
     ;     /,' EDIT THE MAIN PROGRAM AND CHANGE THE STATEMENT ',
     ;     /,10X,'PARAMETER (IRMAX=...,IIMAX=...,IKMAX=...)',
     ;     /,10X,'THEN COMPILE AND LINK THE PROGRAM')

          STOP ' NOT ENOUGH SPACE'
       ENDIF



C--------------- Initializes to zero all variables

C--------------- Real

       DO I=1,IRMAX
          RV(I)=0D0 !ZERO
       ENDDO

C--------------- Integer

       DO I=1,IIMAX
          IV(I)=0
       ENDDO

C--------------- Character

       DO I=1,IKMAX
          KV(I)='          '
       ENDDO

C--------------- Calls main subprogram

      CALL PRINCIPAL
     ;(!--------------------------------------------------dimensioning variables
     ;IAFLUDSC_COLS        ,IAFLUDSC_ROWS          ,IATRADSC_COLS
     ;,IATRADSC_ROWS       ,IA_COUPLED_DSC_COLS    ,IA_COUPLED_DSC_ROWS
     ;,IDIMAFLU            ,IDIMDFLU               ,IDIMCFLU
     ;,IDIMATRA
     &,IDIMDTRA            ,IDIMDQ    ,IDIMBB
     ;,IDIMCOV             ,IDIMFNT                ,IDIMGRAVEL,IDIMHESS
     ;,IDIMQ               ,IDIMWORK  ,IDIMDERC    ,IDIMDERH
     ;,NBAND               ,NBAND1                 ,NBANDCOV
     ;,NDEVS               ,NFLAGS                 ,NUMTOBS,NBLCVP
     ;
     ;!------------------------------------------some new dimensioning variables
     ;,IPAR_DIR(19),IPAR_DIR(24),MAX(NPBTP,NPBFL)  ,MAXPG ,NTYPEL
     ;,NZPRG    ,NUMTIT  ,NUMTNOD
     ;,NZTRA
     ;,IDIMDENS  ,NMAXF   ,NMAXT
     ;
     ;!-------------------------------------------------------------io-variables
     ;,IOCNSF   ,IOCNST   ,IOCONSRC     ,IOCRITRAP    ,IODENS_INI,IODIM
     ;,IOEQT    ,IOFLLI   ,IOFLSAT      ,IOFMLF       ,IOFMLT
     ;,IOINV    ,IOPRHED  ,IOPINITC     ,IOPINITH     ,IORTS     ,IOSMFL
     &,IOSMTP   ,IOSPARSE ,IODIRECT
     ;,IOTRLI   ,IOTRS        ,IOVAR        ,ITPTVAR
     ;
     ;!-------------------------------------------------------------------arrays
     ;,FOBJ_WGT  ,IFLAGS    ,INORPAR   ,IOLG_PAR  ,IOPTS     ,IOWRITE
     ;,IPAR_DIR  ,IPAR_INV  ,LINMET
     ;,NZONE_PAR ,PAR_DIR   ,PAR_INV   ,PAR_WGT
     ;
     ;!--------------------------------------------------------general variables
     ;,FILENAME ,IERROR   ,ISOT      ,LMXNDL   ,MAINF    ,NFNL
     ;,NFNT     ,NINT     ,NOPTS     ,NPAR     ,NPARALG  ,NPAREL
     ;,NPARF    ,NPARFPRG ,NPARNP    ,NPBFL    ,NPBTP
     ;,NPPEL    ,NPPNP    ,NSTAT     ,NTDMT    ,NTYPAR   ,NUMEL
     ;,NUMNP    ,NWRITE   ,NZPAR
     ;!-------------------------------------------geostatistical problem scalars
     ;,MXROT_GS      ,MXMEASPP_GS
     ;,MXNPP_GS    ,IDIMVAR_GS  ,MXNZON_GS     ,MXDISC_GS
     ;,MXNVAR_GS   ,MXNPRIM_GS  ,IDIMCROSS_GS  ,MXSAM_GS
     ;,MXKRIG_GS   ,MXCLOSE_GS  ,IDIMIVARIO_GS ,IDIMZONPP_GS
     ;,MXSB_GS     ,IDIMWGT     ,IOINV_GS
     ;,MXGRPZN     ,NGROUP_ZN     ,IOPT_GS
     ;,IO_KG_GS    ,MXSC_GS     ,IDIMDATASC_GS
     ;!-------------------------------------------------------------Matrix types
     ;,ITYPAFLU  ,ITYPCFLU  ,ITYPDFLU
     ;,ITYPATRA  ,ITYPDTRA  ,ITYPFLUDSC
     &,ITYPTRADSC,ITYPCOUPLDSC
     ;!--------------------------------------------------------------real arrays
     ;,RV(IDPARC)    ,RV(IDPARM)     ,RV(IDSTPAR)  ,RV(IDCFPAREL)
     ;,RV(IDCFPARNP) ,RV(IDPARNP)    ,RV(IDPAREL)  ,RV(IDCOVPAR)
     ;,RV(IDSOURCE)  ,RV(IDGRAV)   ,RV(IDGRAVEL)
     ;,RV(IDAFLU)    ,RV(IDBFLU)     ,RV(IDHCALAN)
     ;,RV(IDDFLU)    ,RV(IDALFA)     ,RV(IDHAUX1)  ,RV(IDHAUX2)
     ;,RV(IDCFLU)    ,RV(IDAFLUDSC),RV(IDAFLUDSCF)
     ;,RV(IDDERH)    ,RV(IDBM_ND_FL) ,RV(IDBM_ZN_FL)
     ;,RV(IDSOLUTION)
     ;!------------------------------------------ Inverse problem general arrays
     ;,RV(IDGRAD)    ,RV(IDHESSAUX)  ,RV(IDHESS)   ,RV(IDDLT_PAR)
     ; ,RV(IDPARAUX)   ,RV(IDCOVINV) ,RV(IDVJAC)
     ;!---------------------------Coordinates and finite element integral arrays
     ;
     ;,RV(IDCOORD)   ,RV(IDBIBI)     ,RV(IDVOLNOD)  ,RV(IDAREA)
     ;,RV(IDGRDFF)
     ;!-------------------------------------Auxiliar arrays and Real time arrays
     ;,RV(IDDAT_VD)  ,RV(IDWORK)     ,RV(IDFNT)     ,RV(IDTIME)
     ;,RV(IDCAUDAL)  ,RV(IDCONCFLOW)
     ;!------------------------------------ Darcy velocity and related variables
     ;,RV(IDVD)      ,RV(IDQXYZ)     ,RV(IDXNORVD)
     ;,RV(IDDQDFLU)  ,RV(IDDQDTRA)   ,RV(IDDVDH)  ,RV(IDDVDC)
     ;!-----------------------------------------------Consistent velocity arrays
     ;,RV(IDPOINTWEIGHT),RV(IDGRADLOC), RV(IDBUOYANCY),RV(IDDBUOYANCY)
     ;,RV(IDGP_COORD)
     ;!---------------------------------------------------------transport arrays
     ;,RV(IDWATVOL)  ,RV(IDACTH)
     ;,RV(IDDERC)    ,RV(IDCAUX1)    ,RV(IDCAUX2)    ,RV(IDDVDP)
     ;,RV(IDATRA)    ,RV(IDATRADSC)  ,RV(IDATRADSCF) ,RV(IDBTRA)
     ;,RV(IDCCALAN)  ,RV(IDDTRA)     ,RV(IDBM_ZN_TT)
     ;,RV(IDBM_ND_TT),RV(IDDWDH)
     ;
     ;!------------------------------------------------------ Observation arrays
     ;,RV(IDVOBS)    ,RV(IDVOBSC)    ,RV(IDDVOBS)  ,RV(IDTOBS)
     ;,RV(IDTIT)     ,RV(IDBUDAT)    ,RV(IDEXTNBU) ,RV(IDWTOBSBU)
     ;,RV(IDWTOBSU)  ,RV(IDWTOBSN)   ,RV(IDWTOBST)
     ;!---------------------------------------------------------- FLOW nonlinear
     ;,RV(IDHBASE)   ,RV(IDHINI)     ,RV(IDHPREV1)
     ;,RV(IDDFLUDFLU),RV(IDDBFLUDFLU)
     ;,RV(IDHPREV2)  ,RV(IDDNODALRH) ,RV(IDDPARELDH)
     ;
     ;!------------------------------------------------------Transport nonlinear
     ;,RV(IDCBASE)   ,RV(IDCPREV1)    ,RV(IDCPREV2)
     ;,RV(IDCCALIT)  ,RV(IDDTRADTRA)
     ;
     ;!-----------------------------------------------Coupled flow and transport
     ;,RV(IDA_COUPLED_DSC),RV(IDA_COUPLED_DSCF),RV(IDDFLUDTRA)
     ;,RV(IDDTRADFLU),RV(IDDPARELDC), RV(IDBCOUPLED)
     ;,RV(IDDBFLUDTRA),RV(IDDERVISC)
     ;
     ;!------------------------- Common non linear variables (flow or transport)
     ;,RV(IDHCALIT)  ,RV(IDDTMXDS)     ,RV(IDPARACD)
     ;
     ;!------------------------------------------------density dependency arrays
     ;,RV(IDDENSITY) ,RV(IDVISCOSITY),RV(IDDELTAITER)
     ;
     ;!------------------------------- Geostatistical inverse problem variables
     ;,RV(IDPARZ)         ,RV(IDWGT_PAR) ,RV(IDDERIV)  ,IV(IDIPOS)
     ;,IV(IDIPNT_PAR)     ,RV(IDVARIO_GS)          ,RV(IDROTMAT_GS)
     ;,RV(IDPOSMEAS_GS)   ,RV(IDVMEAS_GS)          ,RV(IDPOSDIS_GS)
     ;,RV(IDPOSDISAUX_GS) ,RV(IDPOSZN_GS)          ,RV(IDDIVZN_GS)
     ;,RV(IDVSTATS_GS)    ,RV(IDTRIM_GS)           ,RV(IDSEARCH_GS)
     ;,RV(IDESTKRIG_GS)   ,RV(IDZNWGT_GS)          ,RV(IDCROSSCOV_GS)
     ;,RV(IDDATASC_GS)         ,RV(IDSUPBL_GS)
     ;,RV(IDKRISYS_GS)    ,RV(IDKRISOL_GS)         ,RV(IDCLOSESAM_GS)
     ;,RV(IDKRIGAUX_GS)   ,RV(IDEXDRZN_GS)         ,RV(IDCOORDGR_GS)
     ;,IV(IDIVARIO_GS)    ,IV(IDIPOLDRIFT_GS)
     ;,IV(IDIZN_PP_GS)    ,IV(IDIZN_NPP_GS)        ,IV(IDICROSSCOV_GS)
     ;,IV(IDISUPBL_GS)         ,IV(IDNUMSB_GS)
     ;,MXZONPP_GS         ,IV(IDICHECK_GS)         ,IV(IDLDIM_GS)
     ;,RV(IDWGT_UNK)      ,RV(IDPARGOOD)           ,RV(IDPARC_GS)
     ;,IV(IDIFLAG_SIMUL)  ,RV(IDCOVPAR_GR)
     ;                    ,RV(IDDTPREVINV)
     ;
     ;!--------------------------------------------- Statistical analysis arrays
     ;,RV(IDDEVICESTAT),RV(IDEIGENVEC),RV(IDRESID),RV(IDRESIDPAR)
     ;
     ;!----------------------------------------------------------- String arrays
     ;,KV(IDDEVNAME)
     ;
     ;!================================================================
     ;!--------------- Integer arrays
     ;,IV(IDIBCOD)  ,IV(IDISOZ)    ,IV(IDNFTPAR)    ,IV(IDNFNLPAR)
     ;,IV(IDIVPAR)  ,IV(IDIDMBLCVP),IV(IDLXPAREL)   ,IV(IDIXPARNP)
     ;
     ;!--------------------------------------------------- Discretization arrays
     ;,IV(IDKXX)     ,IV(IDLTYPE)   ,IV(IDLNNDEL)    ,IV(IDLDIM)
     ;,IV(IDIBTCO)   ,IV(IDINDPAR)
     ;!----------------------------------------------------Sparse storage arrays
     ;,IV(IDIAD_S)   ,IV(IDIADD_S)  ,IV(IDIADN_S)    ,IV(IDIAD_D)
     ;,IV(IDIADD_D)  ,IV(IDIADN_D)  ,IV(IDIAFD_S)    ,IV(IDIAFDD_S)
     ;,IV(IDIAFDN_S)  ,IV(IDIAFD_D)  ,IV(IDIAFDD_D)   ,IV(IDIAFDN_D)
     ;
     ;!---------------- Integer arrays related to measurements and to dimensions
     ;!---------------- NUMTNOD, NUMTOBS and NUMTIT
     ;,IV(IDIODEVICE),IV(IDINDEXNOD) ,IV(IDIOBUTYP)  ,IV(IDIOCALBU)
     ;,IV(IDIOUTYP)  ,IV(IDNBUW)     ,IV(IDNOBUF)    ,IV(IDNOOBSIT)
     ;,IV(IDIOTINT),IV(IDMEASTYP)
     ;
     ;!--------------- Time arrays & Auxiliary array & Non linear integer arrays
     ;,IV(IDKINT)    ,IV(IDISOLEQ)  , IV(IDNFNLTIP)
     ;,IV(IDNFNLPRG) ,IV(IDIPARTNER),IV(IDLCOORD)
     ;
     ;
     ;!-------------------------------------------- Statistical analysis arrays
     ;
     ;,IV(IDOBSCLASS),IV(IDITYPEPAR))

C-------------------------------------------------- Computes the total cputime

       CALL CPUTIME (MAINF)

C-------------------------------------------------- Writes code version

       WRITE (MAINF,10)

       END
