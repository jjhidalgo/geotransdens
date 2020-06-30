      SUBROUTINE SIM_JAC
     ;(A_COUPL_DSC,A_COUPL_DSCF,ACTH     ,AFLU     ,AFLUDSC  
     ;,AFLUDSCF  ,ALFA     
     &,AREA      ,ATRA      ,ATRADSC     ,ATRADSCF ,BCOUPLED
     &,BETAC     ,BFLU      ,BIBI                  ,BM_ND_FL
     &,BM_ND_TT  ,BM_ZN_FL  ,BM_ZN_TT    ,BTRA
     &,BUOYANCY  ,DBUOYANCY  
     &,CAUDAL    ,CAUX1     ,CAUX2       ,CCALAN   ,CCALIT
     &,CFLU      ,CFPAREL   ,CFPARNP     ,CONCFLOW ,COORD    ,CPREV1
     &,CPREV2    ,CREF      ,DAT_VD      ,DBFLUDFLU,DBFLUDTRA 
     &,DELTAITER   ,DENSITY  ,DENSREF
     &,DER_VISC 
     &,DERC      ,DERH      ,DFLU        ,DFLUDFLU ,DFLUDTRA ,DNODALRH 
     &,DPARELDC  ,DPARELDH  ,DQDFLU      ,DQDTRA 
     &,DTMXDS    ,DTPREVINV ,DTRA        ,DTRADFLU ,DTRADTRA ,DVDC    
     &,DVDH      ,DVOBS     ,DWDH        ,FILENAME
     &,FNT       ,GP_COORD  ,GRAVEL      ,GRDFF    ,HAUX1    ,HAUX2   
     &,HBASE     ,HCALAN    ,HCALIT      ,HPREV1   ,HPREV2
     &,IAD_S     ,IADD_S      ,IADN_S
     &,IAD_D     ,IADD_D      ,IADN_D
     &,IAFD_S    ,IAFDD_S     ,IAFDN_S
     &,IAFD_D    ,IAFDD_D     ,IAFDN_D
     &,IA_COUPLED_DSC_COLS ,IA_COUPLED_DSC_ROWS
     &,IAFLUDSC_COLS,IAFLUDSC_ROWS,IATRADSC_COLS,IATRADSC_ROWS,IBCOD
     &                      ,IBTCO       ,IDIMAFLU ,IDIMATRA,IDIMBB
     &,IDIMCFLU  ,IDIMDENS  ,IDIMDERC    ,IDIMDERH ,IDIMDFLU,IDIMDQ
     &,IDIMDTRA  ,IDIMFNT   ,IDIMGRAVEL,IDIMQ       ,IDIMWORK
     &,IFLAGS    ,INDEXNOD  ,INORPAR     ,IOBMCMP
     &,IOCONSRC    ,IODENS_INI
     &,IODEVICE  ,IODIM
     &,IODIRECT  ,IOEQT     ,IOFLLI      ,IOFLSAT  ,IOFMLF  ,IOFMLT
     &           ,IOINV     ,IOPINVDT    ,IOPTS
     &,IORTS     ,IOTRLI    ,IOTRS       
     &,IOWRITE   ,IPAR_DIR  ,IPARTNER    ,ISOLEQ   ,ISOT    ,ISOZ
     &,ISPARSE   ,ITERGLMX  ,ITPTVAR 
     &,ITYPAFLU  ,ITYPCFLU    ,ITYPDFLU
     &,ITYPATRA  ,ITYPDTRA  
     &,ITYPTRADSC,ITYPACOUPLDSC          ,ITYPFLUDSC
     &,IVPAR
     &,IXPARNP   ,KINT      ,KXX         ,LDIM
     &,LINMET    ,LMXNDL    ,LNNDEL      ,GRADLOC  ,LTYPE   
     &,LXPAREL   ,MAXITER     ,MAXNB    ,MAXNBF
     &,MAXNN     ,MAXPG
     &,MIN_STOP  ,NBAND       ,NBAND1   
     &,NDEVS                ,NFLAGS      ,NFNL
     &,NFNLPAR   ,NFNLPRG   ,NFNLTIP     ,NFTPAR   ,NINT    ,NMAXF
     &,NMAXT     ,NOOBSIT
     &,NOPTS     ,NPAR      ,NPARALG     ,NPAREL   ,NPARF   ,NPARNP
     &,NPBFL     ,NPBMX     ,NPBTP       ,NPPEL    ,NPPNP   ,NROW
     &,NTDMT     ,NTYPAR    ,NTYPEL      ,NUMEL    ,NUMITER ,NUMNP
     &,NUMTIT    ,NUMTNOD   ,NUMTOBS     ,NWRITE   ,NZONE_PAR
     &,NZPAR     ,NZPRG     ,NZTRA       ,PAR_DIR  ,PARACD  ,PARZ
     &,PAREL     ,PARNP     ,POINTWEIGHT,PRGC      ,QXYZ    ,SOLUTION
     &,SOURCE    ,TIME      ,TIT         ,TOBS     ,VAR_REF ,VD
     &,VISCOSITY ,VISCREF   ,VJAC        ,VOBSC    ,WATVOL  ,WORK
     &,WTOBSN    ,WTOBST    ,XNORVD,DVDP,IOLG_PAR,IOCTRA,HINI,WSPECHEAT
     &,WTHERMCON
     ;,IDIMWGT   ,WGT_PAR  ,IPNT_PAR,IPOS     ,DERIV
     &,NDEVGEO   ,PARNAME)

********************************************************************************
*
* PURPOSE  To manage the resolution of flow and transport through time steps,
*          compute mass balance, compute RHS of inverse problem.
*
* DESCRIPTION 
*
* VARIABLES
*
*  DELTACGL      Increment of mass fraction between two global iterations 
*              Used for checking convergence in the picard scheme
*  DELTAHGL      Increment of mass fraction between two global iterations 
*              Used for checking convergence in the picard scheme
*  DELTAITER   Array used to store updates of flow or transport to compare 
*              the solutions of two global iterations
*  IDELHGL     Node for which the head update between 2 global iterations was the 
*              largest
*  IDELCGL     Node for which the mass fraction update between 2 global iterations 
*              was the largest
*  INEW        Index to locate the RHS of the derivatives of head
*                 at current time with respect to parameters inside
*                 array DERH
*  INEWT       Location in array DERC
*  IOBMCMP     Mass balance option
*                 0. Do not compute mass balance.
*                 1. Compute mass balance.
*  IOLD        Index to locate the derivatives of head with 
*                 respect to parameters at the previous time step
*  IOLDT       Third dimension of array DERC used to store           
*                 derivatives of C w.r.t. parameters at previous time   
*                 step 
*
*
*  IODENS       Density dependent flow.
*                 0. No
*                 1. Yes
*  IODIV        Divergence control.
*                 0. No divergence.
*                 1. Solving process has diverged.
*                    Simulation stops.
*                 The divergence of the process takes place
*                  when the solution of the simultion has 
*                  not converged and it is not possible to
*                  reduce the time increment.
*  IOEQT        Type of problem to be solved.
*                 0. Check only.
*                 1. Flow
*                 2. Transport
*                 3. Flow and transport
*  IOFLLI       Linear flow problem option
*                 0. Linear
*                 1. Non linear.
*  IOINV        Inverse problem options:
*                 0. Only simulation.
*                 1. Estimation of flow parameters.
*                 2. Estimation of transport parameters.
*                 3. Estimation of flow and transport parameters.
*  IOTRLI       Linear transport problem option
*                 0. Linear
*                 1. Non linear.
*  IOTRS        Flow regime.
*  ISOLEQ       Array containing the type of head/concentration
*                 solution desired for the user at each obs. time
*                 0. Transient
*                 1. Steady state.
*                 2. Read initial conditions
*                 3. Null initial conditions
*                 4. Not solve.
*  ISOLFL       Status of flow problem in current time step.
*                 0. Not solved.
*                 1. Steady flow solved.
*                 2. Transient flow solved.
*  ISOLTR       Status of transport problem in current time step.
*                 0. Not solved.
*                 1. Steady transport solved.
*                 2. Transient transport solved.
*  ITTERTOTGL   Cumulative iteration counter for coupled problem
*  ITTERTOTFL   Cumulative iteration counter for flow problem
*  ITTERTOTTR   Cumulative iteration counter for transport.
*  LINMET       Array containing initial linearization method and current one.
*                 Rows:
*                   1. Flow
*                   2. Transport
*                   3. Coupled.
*
*                 Columns:
*                   1. Initial method.
*                   2. Current method.
*
*                 Values:
*                   0. Not solve (<=> Compute matrices).
*                   1. Picard
*                   2. Newton
*  MIN_STOP     Control variable for stopping minimization process.
*                 0. Objective function diminishes, so code performs
*                    a new Marquardt's iteration.
*                 1. Code accepts convergence and stops.
*                 2. Code performs the last simulation and stops.
*  NPBFL        Number of flow problems.
*  NTCOMP       Number of computation times carried out
*  NTRNRF       Number of reduction time steps (flow) / Total number of N-R process time reductions                                                   
*  NTRNRT       Number of reduction time steps (transport)  
*  NUMITER      Number of iterations of inverse problem.
*  TINC         Current time increment
*  TOLD         Previous time increment that lead to convergence.
********************************************************************************
      IMPLICIT NONE

      INCLUDE 'MAIN_COMM.FOR'
      INCLUDE 'COMMON_DMT.FOR'
      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER*4::IA_COUPLED_DSC_COLS,IA_COUPLED_DSC_ROWS,IAFLUDSC_COLS
     &          ,IAFLUDSC_ROWS      ,IATRADSC_COLS      ,IATRADSC_ROWS
     &          ,IDELCGL  ,IDELHGL
     &          ,IDIMAFLU ,IDIMATRA ,IDIMBB   ,IDIMCFLU ,IDIMDENS
     &          ,IDIMDERC ,IDIMDERH ,IDIMDFLU ,IDIMDQ   ,IDIMDTRA
     &          ,IDIMFNT  ,IDIMGRAVEL
     &          ,IDIMQ    ,IDIMWORK ,IENTRY   ,IERROR   ,II
     &          ,INARR    ,INCOE    ,INCON    ,INDENDDT 
     &          ,INDSSTR  ,INEW     ,INEWT    ,INTI     ,IOBMCMP
     &          ,IOCONVF  ,IOCONVGL ,IOCONVT
     &          ,IOCONSRC ,IDIMWGT  ,IPNT_PAR ,IPOS
     &          ,IODENS   ,IODENS_INI         ,IODIM    ,IODIRECT
     &          ,IODIV    ,IOEQT    ,IOFLSAT  ,IOFLLI   ,IOFMLF
     &          ,IOFMLT   ,IOINV    ,IOITERFLEND
     &          ,IOITERGLEND        ,IOITERTREND        ,IOLD
     &          ,IOLDT    ,IOPINVDT ,IOPL     ,IORTS
     &          ,IOTRLI   ,IOTRS    ,IOWNRR,IOCTRA   
     &          ,IPBFL    ,IPBTP
     &          ,IPROB    ,IREDTIMC ,IREDTIMH ,IRESCMAX
     &          ,ISOLFL   ,ISOLTR   ,ISOT     ,ISPARSE  ,ISUMFO
     &          ,ITERFL   ,ITERGL   ,ITERGLMX ,ITERM
     &          ,ITERTOTFL,ITERTOTGL,ITERTOTTR,ITERTR
     &          ,ITPTVAR  ,ITYPACOUPLDSC      ,ITYPAFLU ,ITYPATRA
     &          ,ITYPCFLU ,ITYPDFLU ,ITYPDTRA ,IOCALMAT
     &          ,ITYPFLUDSC         ,ITYPTRADSC
     &          ,IUCAL    ,LMXNDL   ,MAXITER  ,MAXNB
     &          ,MAXNBF   ,MAXNN    ,MAXPG    ,MIN_STOP
     &          ,MXNRTV   ,NBAND    ,NBAND1   ,NCONVI
     &          ,NCONVITR ,NDEVS
     &          ,NFLAGS   ,NFNL     ,INTIPREV ,IZ       ,IPOSIC
     &          ,NINT     ,NOOBSIT  ,NOPTS    ,NPAR     ,NPARALG
     &          ,NPAREL   ,NPARF    ,NPARNP   ,NPBFL    ,NPBMX
     &          ,NPBTP    ,NPPEL    ,NPPNP    ,NROW     ,NTCOMP
     &          ,NTDMT    ,NTRNRF   ,NTRNRT   ,NTRNRV   ,NTYPAR
     &          ,NTYPEL   ,NUCNVAR  ,NUMEL    ,NUMITER  ,NUMNP
     &          ,NUMTIT   ,NUMTNOD  ,NUMTOBS  ,NWRITE   ,NZPAR
     &          ,NZPRG    ,NZTRA    ,NMAXF    ,NMAXT

      LOGICAL::FILEEX

      INTEGER*4::IAD_D(2*MAXNB,2*NUMNP)
     ;          ,IAD_S(MAXNB,NUMNP)
     &          ,IADD_D(2*NUMNP)              ,IADD_S(NUMNP)
     &          ,IADN_D(2*NUMNP)              ,IADN_S(NUMNP)
     &          ,IAFD_D(2*MAXNBF,NUMNP*2)
     ;          ,IAFD_S(MAXNBF,NUMNP)
     &          ,IAFDD_D(NUMNP*2)             ,IAFDD_S(NUMNP)
     &          ,IAFDN_D(NUMNP*2)             ,IAFDN_S(NUMNP)
     &          ,IBCOD(NUMNP,NPBFL)           
     &          
     ;          ,IBTCO(NUMNP,NPBTP)
     &          ,IFLAGS(NFLAGS)               ,INDEXNOD(NUMTNOD)
     &          ,INORPAR(NTYPAR)              ,IODEVICE(NDEVS+1,10)
     &          ,IOLG_PAR(NTYPAR,2)           ,IOPTS(NOPTS)
     &          ,IOWRITE(NWRITE)              ,IPAR_DIR(NPARALG)
     &          ,IPARTNER(6,3,6)              ,ISOLEQ(NINT,4)
     &          ,ISOZ(NZTRA)                  ,IVPAR(NZPAR)
     &          ,IXPARNP(NUMNP,NPARNP,NPBMX)  ,KINT(NINT)
     &          ,KXX(LMXNDL,NUMEL)
     &          ,LDIM(NUMEL)                  ,LINMET(3,2)
     &          ,LNNDEL(NUMEL)                ,LTYPE(NUMEL)               
     &          ,LXPAREL(NUMEL,NPAREL,NPBMX)  ,NFNLPAR(NZPAR)
     &          ,NFNLPRG(8,NFNL)              ,NFNLTIP(NFNL)                
     &          ,NFTPAR(NZPAR)                ,NZONE_PAR(NTYPAR)

      REAL*8::BETAC    ,CREF     ,DELTAC_SJ,DELTAH_SJ
     &       ,DELTACOLD,DELTAHOLD  
     &       ,DELTAC_TR  ,DELTAH_FL   ,DENSREF  
     &       ,DRELCMX,DRELHMX
     &       ,DTAVGIO  ,DTIMEF   ,DTIMET
     &       ,DTINITIAL,DTMAXIO  ,DTPREVINV,RESCMAX  ,RESCMAXOLD
     &       ,RESHMAX  ,RESHMAXOLD
     &       ,RESHMAXGL,RESCMAXGL,RESCMAXGLOLD,RESHMAXGLOLD
     &       ,TABSOLUT ,THETAF   ,THETAT   ,TICAL
     &       ,TICALAN  ,TIMEINT  ,TINC     ,TINCINI  ,TINCLAST
     &       ,TINTERVOBS         ,TOLD     ,VAR_REF  ,VISCREF
     &       ,WSPECHEAT,WTHERMCON,WGT_PAR  ,DERIV


      REAL*8::A_COUPL_DSC(IA_COUPLED_DSC_ROWS,IA_COUPLED_DSC_COLS)
     &       ,A_COUPL_DSCF(MAXNB,2*NUMNP)
     &       ,ACTH(NUMEL)       ,AFLU(NUMEL,IDIMAFLU,NPBFL)
     &       ,AFLUDSC(IAFLUDSC_ROWS,IAFLUDSC_COLS,NPBFL)
     &       ,AFLUDSCF(MAXNBF,NUMNP)
     &       ,ALFA(NUMNP,NPBFL)                       ,AREA(NUMEL)
     &       ,ATRA(NUMEL,IDIMATRA,NPBTP)
     &       ,ATRADSC(IATRADSC_ROWS,IATRADSC_COLS,NPBTP)
     &       ,ATRADSCF(MAXNBF,NUMNP)    ,BCOUPLED(2*NUMNP)
     &       ,BFLU(NUMNP,NPBFL)         ,BIBI(IDIMBB,NUMEL)
     &       ,BTRA(NUMNP,NPBTP)
     &       ,CAUDAL(NUMNP,NPBFL)       ,CAUX1(NUMNP,NPBTP)
     &       ,CAUX2(NUMNP,NPBTP)        
     &       ,CCALAN(NUMNP,NPBTP)       ,CCALIT(NUMNP,NPBTP)
     &       ,CFLU(NUMEL,IDIMCFLU)      ,CFPAREL(NUMEL,NPAREL)
     &       ,CFPARNP(NUMNP,NPARNP)     ,CONCFLOW(NUMNP)
     &       ,COORD(NUMNP,3)            
     &       ,CPREV1(NUMNP,NPBTP)       ,CPREV2(NUMNP,NPBTP)
     &       ,DAT_VD(IODIM,IDIMDQ,NUMEL,NPBFL)
     &       ,DBFLUDFLU(NUMNP,NPBFL)    ,DBFLUDTRA(NUMNP)
     &       ,DENSITY(IDIMDENS)         ,DELTAITER(NUMNP)
     &       ,DER_VISC(NUMEL)           ,DERC(NUMNP,NPAR,IDIMDERC,NPBTP)
     &       ,DERH(NUMNP,NPARF,IDIMDERH,NPBFL)
     &       ,DFLU(NUMEL,IDIMDFLU,NPBFL)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL,NPBFL)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL)    ,DNODALRH(NUMNP,4)
     &       ,DPARELDC(NPPEL,NUMEL)      ,DPARELDH(NPPEL,NUMEL)
     &       ,DQDFLU(NUMEL,LMXNDL*LMXNDL),DQDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DTMXDS(NINT)              ,DTRA(NUMEL,IDIMDTRA,NPBTP)
     &       ,DTRADFLU(NUMEL,LMXNDL*LMXNDL)    
     &       ,DTRADTRA(NUMEL,LMXNDL*LMXNDL,NPBTP)
     &       ,DVDC(LMXNDL,IODIM,NUMEL)  ,DVDH(LMXNDL,IODIM,NUMEL)
     &       ,DVOBS(NPAR,3)
     &       ,FNT(IDIMFNT,NINT)         ,GP_COORD (6,8,3)
     &       ,GRAVEL(IDIMGRAVEL,MIN(3*IDIMGRAVEL,3))
     &       ,GRDFF(IODIM,LMXNDL,NUMEL) ,HAUX1(NUMNP,NPBFL)
     &       ,HAUX2(NUMNP,NPBFL)        ,HBASE(NUMEL)
     &       ,HCALAN(NUMNP,NPBFL)
     &       ,HCALIT(NUMNP,NPBFL)       ,HINI(NUMNP)
     &       ,HPREV1(NUMNP,NPBFL)       ,HPREV2(NUMNP,NPBFL)
     &       ,GRADLOC (IODIM,LMXNDL,MAXPG),PAR_DIR(NPARALG)
     &       ,PARACD(3,NFNL)            ,PARZ(NZPAR)
     &       ,PAREL(NUMEL,NPPEL,NPBMX)  ,PARNP(NUMNP,NPPNP,NPBMX)
     &       ,POINTWEIGHT(MAXPG,NTYPEL) ,PRGC(NZPRG)
     &       ,QXYZ(IDIMQ,NUMEL,NPBFL)   ,SOLUTION!(*)
     &       ,SOURCE(NUMNP,2)             ,TIME(NINT)
     &       ,TIT(NUMTIT)               ,TOBS(NUMTOBS,2)
     &       ,VD(IODIM,NUMEL,NPBFL)     ,VISCOSITY(NUMEL)
     &       ,VJAC(NUMTOBS,NPAR)        ,VOBSC(NUMTOBS+NDEVS)
     &       ,WATVOL(MAX(1,(IOPTS(31)-1)*LMXNDL),NUMEL,3,NPBTP)
     &       ,WORK(NUMNP,NBAND1+2,NPBMX)     ,WTOBSN(NUMTNOD)
     &       ,WTOBST(NUMTOBS)           ,XNORVD(NUMEL,NPBFL)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL)
     &       ,DBUOYANCY(IODIM,LMXNDL*LMXNDL,NUMEL)
     &       ,DWDH(MAX(1,(IOPTS(31)-1)*2*LMXNDL),NUMEL,NPBTP)
     &       ,DVDP(IODIM,NPARF,NUMEL)
     &       ,BM_ZN_FL(NMAXF,2,NPBFL) ,BM_ZN_TT(NMAXT,2,NPBTP)
     &       ,BM_ND_FL(NUMNP,8,2)     ,BM_ND_TT(NUMNP,12,2)

      CHARACTER*20::FILENAME
      CHARACTER::PARNAME(NPAR)*4

      INTEGER*4 IREDTIMGL,NUMDIVCGL,NUMDIVHGL,NCONVIFL
     &     ,ISYMETRIC,ND,NFL_SIM,NTP_SIM,I_REC,INCLK,IORECATRA,IOFIRST
     &     ,NDEVGEO

C------------------------- First executable statement

C------------------------- If requiered by user, the entrance in the
C------------------------- subroutine is written in the main unit.

      IF (IFLAGS(3).EQ.1) CALL IO_SUB('SIM_JAC',0)

       MAINF=25       ! Esto es necesario para evitar conflictos con pasar como 
                      ! argumento MAINF y en el COMMON _DMT Se puede arreglar mas 
                      ! finamente pero requiere mas toqueteos. Como MAINF no cambia de 
                      ! valor tampoco tiene mayor importancia

      !inicialize some variables
      IERROR = 0
      II =0
      INARR= 4   !index for areal recharge in array variables (LXPAREL etc)       
      INCOE = 10 !index for external concentration by element in lxparel etc
      INCON = 7  !index for external concentration by nodes in lxparnp etc
      INCLK = 10 !index for conc. leakage coeficient by nodes in lxparnp etc
      IUCAL = 15 !index for INI file unit

      NFL_SIM = MAX(1,IOPTS(28)*NPBFL)  ! # of simultaneous flow pb 
      NTP_SIM = MAX(1,IOPTS(29)*NPBTP)  ! # of simultaneous tpt pb
      MXNRTV = 0
      NTRNRV = 0

      ITERFL = 0
      ITERTR = 0

      IF (IODIRECT.NE.0) THEN
         MAXNB=1  ! para que no pete balance
         MAXNBF=1  ! para que no pete JAC_C
      ENDIF

C------------------------- Initialization of time step control variables

      DTIMEF = 0D0
      DTIMET = 0D0
      TINCLAST = 0D0


      IOFIRST = 0

      ISYMETRIC = 0

      IDELHGL = 0
      IDELCGL = 0


C------------------------- If inverse problem second or further 
C------------------------- iteration is being solved or if mass balance
C------------------------- is to be computed along with inverse problem,
C------------------------- initial conditions must be recovered.

      IF (NUMITER.GT.1 .OR. (IOINV.GT.0 .AND. IOBMCMP.EQ.1) ) THEN

          REWIND(15) ! INI File.
      
C------------------------- Zeros in call are to avoid writing again
C------------------------- previous information in file RES

          CALL ENDATINICOND
     & (IERROR   ,0        ,IORTS    ,IOTRS    ,0         ,15
     & ,MAINF    ,NFL_SIM  ,NTP_SIM  ,NUMNP    ,CCALIT      ,FILENAME
     & ,HCALIT)


      END IF

C------------------------- Stores in HCALAN and CCALAN the first values

      IF (IOEQT.NE.2) CALL EQUAL_ARRAY
     ;                    (HCALAN,HCALIT,NUMNP*MAX(1,IOPTS(28)*NPBFL))
      IF (IOEQT.NE.1) CALL EQUAL_ARRAY
     ;                    (CCALAN,CCALIT,NUMNP*MAX(1,IOPTS(29)*NPBTP))

C---------------------------------------------
C---------------------------------------------
C----- Begins global initialization process --
C---------------------------------------------
C---------------------------------------------


C------------------------- Initializes time steps counter and arrays 
C------------------------- related to non-linear problem. 

      NTCOMP=0        !NUMBER OF COMPUTATIONS TIMES.

C------------------------- If there is any non-linear problem...
      IF ((IOFLLI+IOTRLI).NE.0) THEN

C------------------------- ...initialices the total time reduction number
C------------------------- Updates files containing a simulation solution 
C------------------------- which will be used for initializing during the 
C------------------------- current simulation

          NTRNRF = 0
          NTRNRT = 0
          NTRNRV = 0
          MXNRTV = MAX(IPAR_DIR(1),IPAR_DIR(2))

          CALL UPDATE_FILES_SCRATCH
     &        (IOFLLI   ,IOINV    ,IPAR_DIR(5),IPAR_DIR(4),IOTRLI,ISUMFO
     &        ,MAXITER  ,NUMITER)

      END IF !((IOFLLI+IOTRLI).NE.0)

C------------------------- Output options
C-------------------------IF IOPL >0 something is plotted.
C------------------------- IOPL = IOPLH + IOPLC + IOCMH +IOCMC

      IOPL = IOWRITE(5) + IOWRITE(6) + IOWRITE(11) + IOWRITE(12)

C------------------------- If inverse problem is solved,
C------------------------- initializes arrays related to
C------------------------- generic observations.
      
      IF ((IOINV.LT.0 .AND. IOPL.GE.1) .OR. 
     &    (IOINV.GT.0 .AND. MIN_STOP.NE.1)) THEN

          CALL ZERO_ARRAY(VOBSC,NUMTOBS) 

      END IF !(IOINV.LT.0 .AND. IOPL.GE.1) .OR. ...

      IF (IOINV.GT.0 .AND. MIN_STOP.NE.1) THEN

          CALL ZERO_ARRAY(DVOBS,3*NPAR)      !simulated value.
          CALL ZERO_ARRAY(VJAC,NPAR*NUMTOBS) !Observations' derivatives.

      END IF !IOINV.GT.0 .AND. MIN_STOP.NE.1


C------------------------- Initializes auxiliar variables
C------------------------- related to time increment and
C------------------------- DERH and DERC array (the index
C------------------------- that controls where values at 
C------------------------- current time are stored)

      TINC=0D0        !Current time increment.
      TOLD=0D0        !Previous time increment that lead to convergence.
      IOLD=1          !Index to manage DERH 
      INEW=2          !Index to manage DERH
      IOLDT=1         !Index to manage DERC
      INEWT=2         !Index to manage DERC


C------------------------- Initializes divergence indicator.
       IODIV=0

C------------------------- Initializes cumulative iteration
C--------------------------counters.
      ITERTOTFL = 0
      ITERTOTTR = 0
      ITERTOTGL = 0


C------------------------- Initial water volumen if density dependent flow
C------------------------- or transport is solved (redundant, since density
C------------------------- dependent flow implies to solve transport).
      IF (IODENS_INI.EQ.1. OR. IOEQT.GT.1) THEN
          CALL WATVOL_INI
     &        (INORPAR  ,IOPTS    ,IOPTS(31),LMXNDL   ,LXPAREL
     &        ,NPAREL   ,NOPTS    ,NPBMX    ,NPBTP    ,NTYPAR
     &        ,NUMEL    ,NZPAR    ,LNNDEL   ,ACTH     ,CFPAREL
     &        ,PARZ     ,WATVOL)

      END IF

C-------------------------Initializes HINI for transport inverse problem.

      IF (IOINV.GE.3) THEN

          HINI(:) = HCALIT(:,1)

      END IF !IOINV.GE.3

      !VISCOSITY(:) = VISCREF

C-------------------------------------------
C-------------------------------------------
C----- Ends global initialization process --
C-------------------------------------------
C-------------------------------------------


C-------------------------------------------
C-------------------------------------------
C----- Begins time loop --------------------
C-------------------------------------------
C-------------------------------------------


      DO INTI = MAX(0,IOPTS(50)) ,NINT-1  

          PRINT*,""
          PRINT 10, INTI,NINT-1
          PRINT*,""
   10     FORMAT ('TIME STEP ',I5,' OF ',I5)

C------------------------- If divergence happens, exits time loop.

          IF (IODIV.GT.0) EXIT

          INDENDDT = 0
          IENTRY=0
          IOCONVGL = 0
          INDSSTR=MIN(INTI,1)  ! 0 in first SS at TIME(1) and 1 in the rest
C______________________________prevent that transient flow is calculated when 
C______________________________this is not necesary 

C------------------------- TABSOLUT is initialized to TIME(1) every step.
C------------------------- It seems to be useful to prevent errors when INTI=0,
C------------------------- since when INTI>0, its value is changen in UPDATE_TIME.

          TABSOLUT=TIME(1)


          IF (INDSSTR.EQ.0) THEN
             THETAF = 1D0
             THETAT = 1D0
          ELSE
             THETAF = PAR_DIR(29)
             THETAT = PAR_DIR(30)
          END IF !INDSSTR.EQ.0

C-------------------------Initializes maximum residuals

          RESHMAXGL = HUGE(1D0)
          RESCMAXGL = HUGE(1D0)


C----------------------- Loop over the calculation times between two 
C----------------------- observation times         

          DO WHILE (INDENDDT.EQ.0 .AND. (1-INDSSTR)*IOCONVGL.EQ.0
     &             .AND. IODIV.EQ.0)

                      !in the transient state, 1 - indsstr = 0 and the second
                      !condition is always satisfied. In the steady state, the
                      !second condition is satisfied when the model has not converged
              
C------------------------- Initializes "local" iteration counters.

              ITERGL = 0
              IOITERGLEND = 0

C------------------------- Indexes related to solved problems
C------------------------- are reseted for new time step since
C------------------------- neither flow nor transport have been
C------------------------- solved yet.

              ISOLFL = 0
              ISOLTR = 0

C------------------------- Assesses which problems are to be solved.
C------------------------- This assumes that current interval solution
C------------------------- (steady state, initial condition, etc.) is
C------------------------- adressed within each subroutine (FLOW or
C------------------------- TRANSPORT).

              IOCONVF = 0
              IOCONVT = 0

              IF (IOEQT.EQ.1) THEN !Only flow
                  IOCONVT = 1 
              ELSE IF (IOEQT.EQ.2) THEN !Only transport
                  IOCONVF = 1
              END IF

C------------------------- Initializes linearization method

              LINMET(1:3,2) = LINMET(1:3,1)



C------------------------- If we solve with coupled newton and inti = 0
C------------------------- and we have initial conditions for either flow
C------------------------- or transport but not for the other, then we should
C------------------------- solve for the missing state variable, not by solving
C------------------------- the coupled system but by solving a single system.

              IF (INTI .EQ.0 .AND. LINMET(3,2).EQ.2) THEN

                  IF(IOCONVF.EQ.1 .AND. IOCONVT.EQ.0) THEN

                      LINMET(3,2) = 1
                      IF (LINMET(2,2).EQ.0) LINMET(2,2) = 2

                  ELSE IF(IOCONVF.EQ.0 .AND. IOCONVT.EQ.1) THEN

                      LINMET(3,2) = 1
                      IF (LINMET(1,2).EQ.0) LINMET(1,2) = 2

                  END IF

              END IF


C------------------------- Iterates through FLOW and TRANSPORT
C------------------------- until both have converged or diverged
C------------------------- ("converge" means "get values for state variable",
C------------------------- no matter wether they come from the solution of
C------------------------- an equation, read as initial conditions or taken from
C------------------------- steady state).

              IOCONVGL = 0
              IODIV = 0

              NUMDIVCGL = 0
              NUMDIVHGL = 0

              IREDTIMH = 0
              IREDTIMC = 0
              IREDTIMGL = 0

              IF (INTI.GE.1) THEN
C------------------------- Updates time.
C------------------------- Udates parameters associated with the time increment.
C------------------------- It takes in to account if a reduction or increment of the
C------------------------- step has to be done.
C------------------------- Computes the value of IODIV, since last resort to achieve
C------------------------- convergence is to reduce time increment.
C------------------------- Includes matrix diffusion terms (UPD_DT_MT subroutine).
                  IOWNRR = 0  !provisional : no se utiliza iownrr


                  CALL UPDATE_TIME
     &                (IODENS     ,IOFLLI     ,IOTRLI     ,IPAR_DIR(3)
     &                ,IREDTIMH   ,IREDTIMC   ,NUCNVAR    ,INDENDDT
     &                ,IENTRY     ,INTI       ,NINT       ,NUMITER
     &                ,NPARALG    ,IOPINVDT   ,IOINV      ,KINT
     &                ,TINC       ,TOLD       ,DTMAXIO    ,PAR_DIR(3)
     &                ,PAR_DIR(2) ,DTAVGIO    ,TINTERVOBS ,TICALAN
     &                ,DTINITIAL  ,TICAL      ,DTIMEF     ,DTIMET
     &                ,TABSOLUT   ,TIME       ,DTPREVINV  ,DTMXDS
     &                ,PAR_DIR    ,TINCINI    ,TINCLAST   ,IOEQT
     &                ,ISOLEQ     ,MXNRTV     ,NTRNRV     ,MAINF)

                  IREDTIMGL = 0

C------------------------- Changes in matrix diffusion dimensionless param.
C------------------------- in response to a change in time increment

         IF (NTDMT .NE.0. AND. IOEQT. NE. 1. AND. 
     ;     ( (DABS(TINC-TOLD).GE.1.D-6*TINC).OR.IENTRY.EQ.1) 
     ;                                          .AND. INTI.GE.1 ) THEN
            IF (ISOLEQ(INTI,2).LE.1) CALL UPD_DT_DMT (TINC)
         ENDIF

C------------------------- Asseses what method will be used to solve the problem.
C------------------------- This has to be done every iteration because one might want to
C------------------------- switch.

                  CALL CHOOSE_LINEARIZT_METHOD
     &            (IODENS_INI ,IPAR_DIR(8),IPAR_DIR(14),IPAR_DIR(11)
     &            ,IPAR_DIR(9),IPAR_DIR(12),ITERFL     ,ITERGL
     &            ,ITERTR    ,LINMET)

              END IF !(INTI.GE.1)

              DELTAH_SJ= 1E20
              DELTAC_SJ= 1E20
              DELTAHOLD = 1E20
              DELTACOLD = 1E20

              DO WHILE (IOCONVGL.EQ.0.AND.IREDTIMGL.EQ.0
     &                 .AND.IOITERGLEND.EQ.0)
          

C------------------------- Flow problems are solved if needed.

                  IF (IOCONVF.EQ.0) THEN

                      DO IPROB=1,MAX(1,IOPTS(28)*NPBFL)

C------------------------- IPBFL is the number of the current flow problem

                          IF (IOPTS(28).EQ.1) THEN
                              IPBFL=IPROB
                          ELSE
                              IPBFL=MAX(1,ISOLEQ(MAX(1,INTI),3))
                          ENDIF
                      
C------------------------- Only the first porblem can be a coupled one.

                          IF (IPROB.GT.1) THEN
                              IODENS = 0
                          ELSE
                              IODENS = IODENS_INI
                          END IF


C------------------------- Solves/computes matrices of flow equation.

                          CALL FLOW_EQN
     &(IAFLUDSC_COLS   ,IAFLUDSC_ROWS   ,IDIMAFLU
     &,IDIMBB          ,IDIMCFLU        ,IDIMDENS     ,IDIMDFLU
     &,IDIMFNT         ,IDIMGRAVEL      ,IDIMWORK     
     &,IOCONVF         ,IODIRECT
     &,IENTRY          ,IERROR          ,INCLK        ,INCON   ,INDENDDT
     &,INTI            ,IOCONSRC        ,IODENS
     &,IODIM           ,IOFLLI          ,IOFLSAT      ,IOFMLF   ,IOFMLT
     &                 ,IOINV           ,IOITERFLEND  ,ITPTVAR  ,IOTRLI
     &,IOTRS           ,IOPTS(31)       ,IOWRITE(2)
     &,IREDTIMH        ,ISOLEQ          ,ISOLFL       ,ISOT         
     &,ISPARSE         ,IPAR_DIR(8)     ,IPAR_DIR(14) ,IPAR_DIR(11)
     &,IPAR_DIR(9)     ,IPAR_DIR(12)    ,ITERFL       ,ITERGL
     &,ITERTOTFL       ,ITERTR          
     ;,ITYPAFLU        ,ITYPFLUDSC      ,ITYPCFLU     ,ITYPDFLU 
     ;,LMXNDL          ,MAINF
     &,MAXNB           ,MAXNBF          ,NBAND       
     &,NBAND1                           ,NCONVI       ,NCONVIFL
     &                 ,NFLAGS          ,NFNL         ,NINT               
     &,NPAREL          ,NPARALG         ,NPARNP       ,NPBTP
     &,NPPEL           ,NPPNP           ,NTRNRF       ,NTYPAR
     &,NUMEL           ,NUMITER         ,NUMNP        ,NWRITE
     &,NZTRA           ,NZPAR           ,NZPRG
     &!REAL SCALARS
     &,BETAC          ,CREF             ,DENSREF        ,PAR_DIR(11)
     &,DTIMEF         ,DTIMET
     &,TABSOLUT       ,THETAF           ,THETAT   
     &,TICAL          ,TICALAN
     &,TINC           ,TINCINI          ,TINCLAST     ,TINTERVOBS
     &,TOLD           ,VISCREF          ,VAR_REF      ,PAR_DIR(8)
     &,FILENAME    
     &!INTEGER ARRAYS
     &,IBCOD(1,IPBFL)                                ,IBTCO   ,IFLAGS
     &,INORPAR        ,IPAR_DIR         ,ISOZ        ,IXPARNP(1,1,IPBFL) 
     &,KXX            ,LDIM             ,LINMET       ,LNNDEL
     &,LTYPE          ,LXPAREL(1,1,IPBFL),NFNLPAR      ,NFNLTIP
     &,NFNLPRG        ,NFTPAR           ,NROW         ,NZONE_PAR
     &,IAD_S,IADD_S,IADN_S,IOWRITE
     ;,IAFD_S         ,IAFDD_S          ,IAFDN_S
     &!REAL ARRAYS
     &,AFLU(1,1,IPROB)  ,AFLUDSC(1,1,IPROB),AFLUDSCF
     &,ALFA(1,IPROB)  
     &,AREA           ,BFLU(1,IPROB)    ,BIBI
     &,CAUX1(1,1)     ,CAUX2(1,1)       ,COORD     
     &,CCALAN(1,1)    ,CCALIT(1,1)      ,CFLU         ,CFPAREL
     &,CFPARNP        ,CPREV1           ,CPREV2       ,DENSITY
     &,DER_VISC       ,DNODALRH         ,DFLU(1,1,IPROB)
     &,DFLUDFLU       ,DFLUDTRA         ,DBFLUDFLU(1,IPROB),DBFLUDTRA
     &,DPARELDH       ,DPARELDC
     &,FNT            ,GRAVEL           ,HAUX1(1,IPROB)
     &,HAUX2(1,IPROB) ,HBASE            ,HCALAN(1,IPROB)
     &,HCALIT(1,IPROB),HPREV1(1,IPROB)
     &,HPREV2(1,IPROB),PARZ
     &,PARACD         ,PAR_DIR    ,PAREL(1,1,IPROB)  ,PARNP(1,1,IPROB)
     &,PRGC
     &,SOLUTION       ,TIME             ,VISCOSITY    ,WATVOL
     &,DELTAITER      ,DELTAH_FL        ,IDELHGL      ,RESHMAX
     &,RESHMAXOLD     ,WORK
     &,BUOYANCY       ,DBUOYANCY        ,IPARTNER
     &,GRADLOC        ,GP_COORD         ,POINTWEIGHT  ,MAXPG ,NTYPEL
     &,grdff          ,DWDH(1,1,1)      ,CONCFLOW     ,ATRA  ,DTRA
     &,IDIMATRA       ,IDIMDTRA         ,ITYPATRA     ,ITYPDTRA,CAUDAL
     &,IORTS)

 
                      END DO !IPROB


                  END IF  !If (IOCONVF = 0)             

               IF (ISOLFL.NE.0) THEN

                  IPBFL=MAX(1,ISOLEQ(MAX(1,INTI),3)) ! # of flow problem

C------------------------- Velocities and fluxes are computed if 
C------------------------- transport is to be solved or if they 
C------------------------- are required for output purposes.

                  IF (IOCONVF.EQ.1 .AND.
     &               (IOCONVT.EQ.0 .OR. IOWRITE(16).NE.0)) THEN

                      CALL COMVEL
     &                    (AREA     ,BUOYANCY ,COORD    ,DBUOYANCY
     &                    ,DPARELDC ,DPARELDH ,DVDC     ,DVDH
     &                    ,GP_COORD ,GRADLOC  ,GRAVEL   ,GRDFF
     &                    ,HAUX1    ,IDIMQ    ,IODENS   ,IODIM
     &                    ,IOFLLI   ,ISOZ     ,KXX      ,LDIM
     &                    ,IFLAGS   ,LINMET   ,LMXNDL   ,LNNDEL
     ;                    ,LTYPE    ,LXPAREL(1,1,IPBFL) ,MAXPG
     &                    ,NPAREL   ,NPPEL
     &                    ,NTYPEL   ,NUMEL    ,NUMNP    ,NZTRA
     &                    ,PAREL    ,POINTWEIGHT        ,NFLAGS
     ;                    ,QXYZ     ,VD       ,XNORVD)


                  END IF !IOCONVF.EQ.1.AND. ...

                  IF (IOCONVF.EQ.1 .AND. IOCONVT.EQ.0) THEN
                      if(iflags(21).eq.1) write(7,*) 'INTI ',inti
                      CALL COMFLOW
     &                (AFLU    ,AREA     ,BETAC    ,BUOYANCY
     ;                ,CAUDAL
     &                ,CAUX1    ,CAUX2    ,CFLU     ,COORD    ,CREF
     &                ,DBUOYANCY,DENSITY  ,DENSREF  ,DFLU
     ;                ,DFLUDFLU
     &                ,DFLUDTRA ,DPARELDH
     &                ,GP_COORD ,GRADLOC  ,GRAVEL   ,GRDFF
     ;                ,HAUX1
     &                ,HAUX2    ,IBCOD(1,IPBFL)     ,IDIMAFLU
     ;                ,IDIMDFLU           ,IFLAGS
     &                ,IODENS   ,IODIM    ,IOFLLI   ,ISOLFL   ,ISOZ
     &                ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &                ,LXPAREL(1,INARR,IPBFL)       ,LXPAREL(1,1,IPBFL)
     ;                ,MAXPG    ,NFLAGS
     &                ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL   ,NUMEL
     &                ,NUMNP    ,NZTRA    ,PAREL
     ;                ,PARNP
     &                ,POINTWEIGHT        ,THETAT  ,CONCFLOW
     &                ,IOCONSRC,IBTCO(1,1))


C------------------------- Derivatives of fluxes w. r. t. state
C------------------------- variable are computed if transport is
C------------------------- linearized by N-R method.

                      IF (LINMET(2,2).EQ.2 .OR. LINMET(3,2).EQ.2 
     &                .OR. IOINV.GE.3) THEN

                         CALL COMP_DER_FLOW
     &                       (AFLU     ,BETAC    ,CAUDAL   ,CFLU
     &                       ,CREF     ,DENSREF  ,DFLU     ,DFLUDFLU
     &                       ,DFLUDTRA ,DNODALRH ,DQDFLU   ,DQDTRA
     &                       ,PAR_DIR(31)        ,PAR_DIR(32)
     &                       ,HAUX1    ,IBCOD    ,IDIMAFLU ,IDIMDFLU
     &                       ,0        ,IODENS   ,ISOLFL   ,KXX
     &                       ,LINMET   ,LMXNDL   ,LNNDEL   ,NPPNP
     &                       ,NUMEL    ,NUMNP    ,PARNP    ,TINC)

                      END IF !LINMET(2,2).EQ.2 .OR. ...

                  ENDIF  ! IOCONVF.EQ.1 .AND. IOCONVT.EQ.0
               ENDIF ! ISOLFL.NE.0

                  IF (IOCONVF.EQ.1 .AND. IOCONVT.EQ.0) THEN

C------------------------- Transport problems are solved if needed.

                      IF (INTI.EQ.1) THEN
                          INTIPREV=1
                      ELSE
                          INTIPREV=INTI-1
                      ENDIF

C------------------------- Sets variable IOCALMAT to 1 if ATRA and DTRA 
C------------------------- arrays need to be computed.
C------------------------- ATRA and DTRA arrays have to be computed when
C------------------------- either: its the first time that transport is
C------------------------- solved, there has been a change of tpt problem,
C------------------------- or flow has just been calculated

                      IF ( (IENTRY.EQ.1 .AND. INTI.EQ.1) .OR. 
     ;                   ( ISOLEQ(MAX(1,INTI),4).NE.
     ;                      ISOLEQ(MAX(1,INTIPREV),4) )  .OR.
     ;                                        ISOLFL.GT.0 ) THEN

                          IOCALMAT=1

C--------------------- This is done for matrix diffussion

                          do iz=1,NZONE_PAR(16)
                              iposic=iv(IN_IV_DMT + 4*iz-3)
                              rv(iposic+11)=0.d0
                          enddo

                          ICH_IT_MAR=1! We have just started a new inverse problem iteration

                      ENDIF

                      IF (IOPTS(29).NE.0) IPBFL=1 ! Several simultaneous of flow and one tpt,
                                                  ! the tpt is solved with the FIRST flow prob.
                      DO IPROB=1, MAX(1,IOPTS(29)*NPBTP)

                          IF (IOPTS(29).EQ.1) THEN
                              IPBTP=IPROB
                          ELSE
                              IPBTP=MAX(1,ISOLEQ(MAX(1,INTI),4))
                          ENDIF

C------------------------- Only the first problem can be a coupled one.

                          IF (IPROB.GT.1) THEN
                              IODENS = 0
                          ELSE
                              IODENS = IODENS_INI
                          END IF


                          CALL TRANSPORT
     &          (AREA       ,ATRA(1,1,IPROB)
     &          ,ATRADSC(1,1,IPROB)   ,ATRADSCF
     &          ,BETAC
     &          ,BIBI       ,BTRA(1,IPROB),CAUDAL    
     &          ,CAUX1(1,IPROB)     
     &          ,CAUX2(1,IPROB)       ,CCALAN(1,IPROB)  
     &          ,CCALIT(1,IPROB)      ,CFLU
     &          ,CFPAREL  ,CFPARNP    ,CPREV1    ,CPREV2
     &          ,CREF     ,PAR_DIR(35),DAT_VD
     &          ,PAR_DIR(12)
     &          ,DELTAC_TR   ,DELTAITER  ,DENSITY
     &          ,DENSREF  ,DPARELDC   ,DPARELDH  ,DQDFLU
     &          ,DQDTRA   ,DRELCMX    ,PAR_DIR(34),DTIMET
     &          ,DTRA(1,1,IPROB),DTRADFLU   ,DTRADTRA,DVDH ,DVDC,DWDH
     &          ,PAR_DIR(31),PAR_DIR(32),FILENAME  ,FNT
     &          ,GRDFF
     &          ,IAD_S    ,IADD_S     ,IADN_S  ,IAFD_S    ,IAFDD_S
     &          ,IAFDN_S  
     &          ,IATRADSC_COLS
     &          ,IATRADSC_ROWS,IBTCO(1,IPBTP)   ,IDELCGL
     &          ,IDIMATRA
     &          ,IDIMBB   ,IDIMCFLU   ,IDIMDENS,IDIMDFLU,IDIMDQ   
     &          ,IDIMDTRA
     &          ,IDIMFNT    ,IDIMQ   ,IDIMWORK
     &          ,IFLAGS   ,INARR  ,INCLK    ,INCON
     &          ,INDENDDT ,1          ,INDSSTR ,INORPAR
     &          ,IOWRITE(1),INTI      ,IOCONSRC
     &          ,IOCONVT  ,IODENS     ,IODIM   ,IODIRECT
     &          ,IOFMLT   ,IOINV      ,IOITERTREND,IPAR_DIR(5),IOPTS(30)
     &          ,IORTS    ,IOTRLI      ,IOTRS  ,IOPTS(31) ,iowrite(2)
     &          ,iowrite(17),IPAR_DIR   ,IPROB   ,IREDTIMC  ,IRESCMAX
     &          ,ISOLEQ   ,ISOLTR  ,ISPARSE   ,IPAR_DIR(8)
     &          ,IPAR_DIR(14)         ,IPAR_DIR(11)   ,IPAR_DIR(9)
     &          ,IPAR_DIR(12)         ,ITERFL  ,ITERGL    ,ITERM
     &          ,ITERTOTTR  ,ITERTR  ,ITPTVAR   ,ITYPATRA
     &          ,ITYPTRADSC ,ITYPCFLU ,ITYPDTRA,IUCAL ,IOCALMAT
     &          ,IXPARNP(1,1,IPBTP)   ,KXX        ,LDIM    ,LINMET    
     &          ,LMXNDL     ,LNNDEL   ,LTYPE      ,LXPAREL(1,1,IPBTP)    
     &          ,MAINF      ,MAXNB
     &          ,MAXNBF
     &          ,NBAND    ,NBAND1     ,NCONVI  ,NCONVITR
     &          ,NFLAGS   ,NFNL       ,NFNLPAR ,NFNLPRG   ,NFNLTIP
     &          ,NFTPAR   ,NINT       ,NPARALG ,NPAREL    ,NPARNP
     &          ,NPPEL    ,NPPNP      ,NROW      ,NTDMT
     &          ,NTYPAR   ,NUMEL      ,NUMITER ,NUMNP     ,NZONE_PAR
     &          ,NZPAR    ,NZPRG      ,PAR_DIR ,PARACD    ,PARZ
     &          ,PAREL(1,1,IPROB)     ,PARNP(1,1,IPROB)   ,PRGC
     &          ,QXYZ     ,RESCMAX    ,RESCMAXOLD,PAR_DIR(7),SOLUTION
     &          ,SOURCE(1,2),TABSOLUT ,PAR_DIR(30),TICAL  ,TICALAN
     &          ,TIME     ,TINC       ,TINCINI ,TINCLAST  ,TINTERVOBS
     &          ,TOLD     ,VD      ,WATVOL(1,1,1,IPROB),WORK(1,1,IPROB)
     &          ,WSPECHEAT,WTHERMCON  ,XNORVD   ,PAR_DIR(9),ACTH,NPBTP)


C------------------------- If NO density dependent flow is being solved

                          IF (IODENS.EQ.0) THEN

C------------------------- Global convergence is calculated in case the first
C------------------------- flow or transport problem is non-linear (rho=cte.)
      
                              IOCONVGL = IOCONVF*IOCONVT
                              IREDTIMGL = MAX(IREDTIMH,IREDTIMC)
                              IF (NTRNRV.GT.MXNRTV) IODIV = 1

                          ELSE
C------------------------- If density dependent flow IS BEING solved

                              CALL COUPLED_FLOW_TRANSPORT
     &          (A_COUPL_DSC   ,A_COUPL_DSCF  ,AFLU          ,AREA
     &          ,ATRA          ,BCOUPLED      ,BETAC         ,BFLU
     &          ,BTRA          ,CAUDAL        ,CAUX1         ,CAUX2
     &          ,CCALAN        ,CCALIT        ,CFLU          ,CREF
     &          ,DBFLUDFLU     ,DBFLUDTRA     ,DELTAC_SJ     ,DELTAC_TR
     &          ,DELTACOLD     ,DELTAH_FL     ,DELTAH_SJ     ,DELTAHOLD
     &          ,DELTAITER     ,DENSREF       ,DFLU
     &          ,DFLUDFLU      ,DFLUDTRA      ,DQDFLU        ,DQDTRA
     &          ,DRELCMX       ,DRELHMX       ,DTRA          ,DTRADFLU
     &          ,DTRADTRA      ,HAUX1         ,HAUX2         ,HCALAN
     &          ,HCALIT        ,IA_COUPLED_DSC_COLS
     &          ,IA_COUPLED_DSC_ROWS          ,IAD_D         ,IAD_S
     &          ,IADD_D        ,IADN_D        ,IADN_S        ,IAFD_D
     &          ,IAFDD_D       ,IAFDN_D       ,IBCOD         ,IBTCO
     &          ,IDELCGL       ,IDELHGL       ,IDIMAFLU      ,IDIMATRA
     &          ,IDIMCFLU      ,IDIMDFLU      ,IDIMDTRA      ,IDIMWORK
     &          ,IFLAGS        ,INDENDDT      ,INDSSTR       ,INTI
     &          ,IOCONSRC      ,IOCONVF       ,IOCONVGL      ,IOCONVT
     &          ,IODENS        ,IODIRECT      ,IOFLLI
     &          ,IOITERGLEND   ,IOPTS         ,IORTS         ,IOTRLI
     &          ,IOTRS         ,IOWRITE       ,IPAR_DIR      ,IREDTIMC
     &          ,IREDTIMGL     ,IREDTIMH      ,ISOLFL        ,ISOLTR
     &          ,ISPARSE       ,ITERGL        ,ITERGLMX      ,ITERM
     &          ,ITERTOTGL     ,ITPTVAR       ,ITYPACOUPLDSC ,ITYPCFLU
     &          ,ITYPDFLU      ,KXX           ,LINMET        ,LMXNDL
     &          ,LNNDEL        ,MAINF         ,MAXNB         ,MAXNBF
     &          ,MAXNN         ,NBAND         ,NBAND1        ,NCONVIFL
     &          ,NCONVITR      ,NFLAGS        ,NOPTS         ,NPARALG
     &          ,NPBFL         ,NPBMX         ,NPBTP         ,NPPEL
     &          ,NPPNP         ,NTYPAR        ,NUMDIVCGL     ,NUMDIVHGL
     &          ,NUMEL         ,NUMNP         ,NWRITE        ,NZONE_PAR
     &          ,PAR_DIR       ,PAREL         ,PARNP         ,RESCMAX
     &          ,RESCMAXGL     ,RESCMAXGLOLD  ,RESCMAXOLD    ,RESHMAX
     &          ,RESHMAXGL     ,RESHMAXGLOLD  ,RESHMAXOLD    ,SOLUTION
     &          ,SOURCE        ,TINC          ,WATVOL        ,WORK
     &          ,WSPECHEAT)

                          ENDIF !(IODENS.EQ.0)

C------------------------- When there are simultaneous transport problems,
C------------------------- if the first transport problem is non-linear
C------------------------- and has not converged, the rest of simultaneous
C------------------------- problems cannot be solved, then the transport 
C------------------------- problems loop is interrupted. There

                          IF (NTP_SIM.GT.1) THEN
                              IF(IREDTIMGL.EQ.1 .OR. IOCONVGL.EQ.0)THEN

                                  WRITE(*,60)
   60                             FORMAT
     &                            (/,'First non-linear transport '
     &                              ,'problem did not converge.'
     &                              ,/,'Exiting simultaneous transport '
     &                              ,'problems loop.',/)
                                  EXIT 

                              ELSE

                                  WRITE(*,61)
   61                             FORMAT
     &                            (/,'The first non-linear transport '
     &                              ,'has converged.',/
     &                              ,'Simultaneous transport loop'
     &                              ,' continues',/)

                              END IF !IREDTIMGL.EQ.1 .OR. IOCONVGL.EQ.0
                          END IF !NTP_SIM.GT.1

                      END DO !IPROB=1, MAX(1,IOPTS(29)*NPBTP)

                  END IF  !IOCONVF.EQ.1 .AND. IOCONVT.EQ.0

C------------------------- Global convergence is calculated in case the first
C------------------------- flow or transport problem is non-linear (rho=cte.)
c------------------------- no está muy claro, parece que es para inti=0.

                      IOCONVGL = IOCONVF*IOCONVT
                      IREDTIMGL = MAX(IREDTIMH,IREDTIMC)
                      IF (IREDTIMGL.GT.0) INDENDDT = 0
                      !!IOITERGLEND = 1
                      IF (NTRNRV.GT.MXNRTV .OR.IOITERGLEND.EQ.1)IODIV=1
                      NCONVI = MIN(NCONVIFL,NCONVITR)

              END DO !WHILE IOCONVGL.EQ.0.AND.IODIV.EQ.0


C------------------------- Once both flow and transport have converged...

              IF (IOCONVGL.EQ.1) THEN

C------------------------- Updates density, viscosity and water content.

                  IF (IODENS_INI.EQ.1) THEN

                      CALL CALC_DENS
     &                    (BETAC    ,CAUX1    ,CREF    ,DENSREF
     &                    ,DENSITY  ,IDIMDENS ,1       ,KXX
     &                    ,LMXNDL   ,LNNDEL   ,NUMEL    ,NUMNP)

                      CALL CALC_VISC
     &                    (CAUX1    ,DER_VISC ,ITPTVAR  ,KXX     ,LMXNDL
     &                    ,LNNDEL   ,NUMEL    ,NUMNP    ,VAR_REF
     &                    ,VISCOSITY)

                  END IF !IODENS_INI.EQ.1

                  IF (IOPTS(31).GT.0) THEN

                     CALL CALC_WATVOL
     &                   (HCALAN   ,HCALIT   ,INTI     ,IOPTS(31)
     &                   ,KXX      ,LMXNDL   ,LNNDEL   ,NPBTP
     &                   ,NPPEL    ,NUMEL    ,NUMNP    ,PAREL
     &                   ,THETAT   ,WATVOL)
                  END IF !IOPTS(31).GT.0

C------------------------- Write head or concentration field to temporary file

                  CALL WRITE_TEMP_STATE_VARS
     &                (CCALIT    ,HCALIT    ,INDENDDT  ,INDSSTR
     &                ,INTI      ,IOEQT     ,IOWRITE(8),IOWRITE(7)
     &                ,ISOLFL    ,ISOLTR    ,NINT      
     &                ,MAX(1,IOPTS(28)*NPBFL)
     &                ,MAX(1,IOPTS(29)*NPBTP)          ,NUMNP)

C------------------------- current time interval is computed
C_________________________¿¿since it might have change due to time reduction???
                  IF (INTI.NE.0) THEN
                      TIMEINT=TIME(INTI+1)-TIME(INTI)
                  END IF
*** EL INTI DE LA PRIMERA LINEA ES INDSSTR (SI INTI=0, ESTADO ESTACIONARIO,
*** POR TANTO, INDSSTR=0)

C-------------------------Write velocities.

          IF (IOBMCMP.EQ.1 .AND. IOWRITE(16).NE.0 .AND.
     &        INDENDDT.EQ.1 .AND.
     &        INTI.NE.0 .AND. IOEQT.NE.0 .AND.
     &        (MOD(INTI,IOWRITE(16)).EQ.0 .OR. INTI.EQ.1 .OR.
     &         INTI.EQ.NINT)) THEN

              CALL WRITE_VELOCITY_VMSH
     &                (FILENAME ,INTI     ,IODIM    ,IOFIRST  ,KXX
     &                ,LMXNDL   ,LNNDEL   ,NINT     ,NUMEL    ,NUMNP
     &                ,TIME     ,VD       ,COORD(1,1)
     &                ,COORD(1,2))

              IOFIRST = 1

          END IF !INTI.NE.0 .AND. IOEQT.NE.1...

C---------------------------------------------
C---------------------------------------------
C----- Begins mass balance -------------------
C---------------------------------------------
C---------------------------------------------

C------------------------ Mass balance is computed if required.

                  IF (IOBMCMP.EQ.1) THEN

C______________________________ Time interval between observation times

                      IF (INTI.NE.0) TIMEINT=TIME(INTI+1)-TIME(INTI)

C------------------------ States if ATRA contains recharge
C------------------------ (only when Picard method is used for transport)
C------------------------ (computed again, it may have changed inside TRANSPORT)
                      IF (LINMET(2,2).EQ.1
     &                   .OR. LINMET(3,2).EQ.1 )THEN
                          IORECATRA = 1
                      ELSE
                          IORECATRA = 0
                      END IF !LINMET(2,2).EQ.1 ...

                      CALL BALANCE
     &              (ACTH     ,AFLU     ,AREA     ,ATRA     ,BETAC
     &              ,BM_ND_FL ,BM_ND_TT ,BM_ZN_FL ,BM_ZN_TT ,BUOYANCY
     &              ,CAUDAL   ,CAUX1    ,CAUX2    ,CCALAN   ,CCALIT
     &              ,CFLU     ,COORD    ,CREF     ,DBUOYANCY,TINC
     &              ,DENSITY  ,DENSREF  ,DFLU     ,DFLUDFLU ,DFLUDTRA
     &              ,DPARELDH ,DTRA     ,GP_COORD ,GRADLOC  ,GRAVEL
     &              ,GRDFF    ,HCALAN   ,HCALIT   ,HAUX1    ,HAUX2
     &              ,I_REC    ,IAD_S    ,IADN_S   ,IBCOD    ,IBTCO
     &              ,IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU ,IDIMDTRA
     &              ,IENTRY   ,IFLAGS   ,INDENDDT ,INTI     ,IOPTS(18)
     &              ,IOPTS(22),IOPTS(21),IOPTS(20),IOPTS(19),IOPTS(17)
     &              ,IOCONSRC ,IODENS   ,IODIM    ,IOEQT    ,IOFLLI
     &              ,IOPTS(28),IOPTS(29),IORECATRA,IOPTS(31),ISOLEQ
     &              ,ISOLFL   ,ISOLTR   ,ISOZ     ,ITPTVAR  ,ITYPAFLU
     &              ,ITYPATRA ,ITYPCFLU ,ITYPDFLU ,ITYPDTRA ,IXPARNP
     &              ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &              ,LXPAREL  ,MAXNB    ,MAXPG    ,NFLAGS   ,NINT
     &              ,NMAXF    ,NMAXT    ,NPARALG  ,NPAREL   ,NPARNP
     &              ,NPBFL    ,NPBMX    ,NPBTP    ,NPPEL    ,NPPNP
     &              ,NTYPAR   ,NTYPEL   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &              ,NZTRA    ,PAR_DIR  ,PAREL    ,PARNP    ,POINTWEIGHT
     &              ,TABSOLUT ,TIME     ,WSPECHEAT,WATVOL)

                  END IF ! IOBMCMP=1


C---------------------------------------------
C---------------------------------------------
C----- Ends mass balance -------------------
C---------------------------------------------
C---------------------------------------------


C------------------------- If solving, checking, or solving flow, then
C------------------------- stores last solution and current time (absolut)

                  IF(IOEQT.NE.2.AND.IOTRS.NE.0) THEN
                      CALL STORE_VAR_SCRATCH
     &              (NINT   ,NTCOMP ,NTCOMP+1  ,NUMNP  ,TABSOLUT
     &              ,TIME   ,HCALIT)
                  END IF

                  IF(IOEQT.NE.1.AND.IORTS.NE.0) THEN
                      CALL STORE_VAR_SCRATCH
     &            (NINT   ,NTCOMP ,NTCOMP+1    ,NUMNP  ,TABSOLUT
     &            ,TIME   ,CCALIT)
                  END IF

C---------------------------------------------
C---------------------------------------------
C----- Begins Inverse problem RHS computing --
C---------------------------------------------
C---------------------------------------------


C------------------------- RHS of inverse problem computed if inverse problem
C------------------------- is to be solved and minimization must not stop.

                  IF (IOINV.GT.0 .AND. MIN_STOP.NE.1) THEN

C------------------------- If flow parameters are to be estimated
C------------------------- and flow has been solved, RHS for inverse 
C------------------------- flow problem is computed for all ¿simultaneous? 
C------------------------- problems.
                      IF (IOINV.NE.2 .AND. ISOLFL.GE.1) THEN

                          II=2-ISOLFL

                          DO IPROB=1,MAX(1,IOPTS(28)*NPBFL)


                              IF (IODENS.EQ.1 .AND. IPROB.EQ.1) THEN
                                  ISYMETRIC = 0
                              ELSE
                                  ISYMETRIC=1
                              END IF !IODENS.EQ.1 .AND. IPROB.EQ.1

                              IF (IFLAGS(25).GT.0) THEN

                                  WRITE (747,*)
     &                               'ITERACION ',NUMITER, ' INTI ',INTI
                              END IF !IFLAGS(25).GT.0


                               CALL JAC_H
     &          (ACTH     ,AFLU     ,AFLUDSC  ,AFLUDSCF ,ALFA
     &         ,AREA      ,BETAC    ,BIBI               ,BUOYANCY
     &         ,CAUX1     ,CAUX2    ,CCALAN   ,CCALIT   ,CFLU
     &         ,CFPAREL   ,CFPARNP  ,COORD    ,CREF     ,DBFLUDFLU
     &         ,DBFLUDTRA ,DENSITY  ,DENSREF  ,DERC     ,DERH
     &         ,DFLU      ,DFLUDFLU ,DFLUDTRA ,DTIMEF   ,PAR_DIR(31)
     &         ,PAR_DIR(32),FNT      ,GP_COORD ,GRADLOC  ,HAUX1
     &         ,HAUX2     ,HBASE    ,HCALAN   ,HCALIT   ,IAD_S
     &         ,IADD_S    ,IADN_S   ,IAFLUDSC_COLS      ,IAFLUDSC_ROWS
     &         ,IAFD_S    ,IAFDD_S  ,IAFDN_S  ,IBCOD
     &                    ,1        ,IDIMAFLU ,IDIMBB   ,IDIMCFLU
     &         ,IDIMDENS  ,IDIMDERC ,IDIMDERH ,IDIMDFLU ,IDIMFNT
     &         ,IDIMWORK  ,IODIRECT ,IFLAGS   ,6-II     ,4-II
     &         ,2-II      ,1-II     ,INEW     ,INEWT    ,INORPAR
     &         ,4-II      ,INTI     ,0        ,IODENS   ,IODIM
     &         ,IOFLLI    ,IOFMLF   ,IOFMLT             ,IOINV
     &         ,IOLD      ,IOLDT    ,IOTRLI   ,IOPTS(31),IPAR_DIR
     &         ,ISOZ      ,ISPARSE  ,ISYMETRIC,ITYPAFLU ,ITYPCFLU
     &         ,ITYPDFLU  ,IVPAR    ,IXPARNP  ,KXX      ,LDIM
     &         ,LMXNDL    ,LNNDEL   ,LTYPE    ,LXPAREL  ,MAINF
     &         ,MAXNB     ,MAXNBF   ,MAXPG    ,NBAND1
     &                    ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &         ,NFNLTIP   ,NFTPAR   ,NINT     ,NPARALG  ,NPAREL
     &         ,NPARF     ,NPARNP   ,NPAR     ,NPPEL    ,NPPNP
     &         ,NTYPAR    ,NTYPEL   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &         ,NZPAR     ,NZTRA    ,PARACD   ,PARZ     ,PAREL
     &         ,PARNP     ,POINTWEIGHT        ,SOLUTION ,THETAF
     &         ,THETAT    ,TINC     ,TINTERVOBS         ,WORK
     ;         ,IDIMWGT   ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV   ,NPAR)

                          END DO !IPROB=1,MAX(1,IOPTS(28)*NPBFL)

C------------------------- INEW and IOLD are interchanged only if there is no
C------------------------- density dependent flow.

                          IF (IODENS.EQ.0) THEN

                              INEW = 3 - INEW
                              IOLD = 3 - IOLD

                          END IF !IODENS.EQ.0


C------------------------- Else if flow parameters are to be estimated
C------------------------- and second or further time interval is solved,...

                      ELSE IF (IOINV.NE.2 .AND. INTI.GE.1) THEN

C------------------------- ...if initial conditions have been read in this 
C------------------------- interval derivatives must be set to zero (to use
C------------------------- in the next interval)

                          IF (ISOLEQ(INTI,1).EQ.2 .OR.
     &                        ISOLEQ(INTI,1).EQ.3) THEN
                             DERH(:,:,:,1:MAX(1,IOPTS(28)*NPBFL) )=0.D0
                          END IF

                      ENDIF ! IOINV.NE.2 .AND. ISOLFL.GE.1


C------------------------- Transport inverse problem? .AND. Transport has been
C------------------------- solved in the current time?

                      IF (IOINV.NE.1 .AND. ISOLTR.GE.1) THEN

C------------------------- Loop over simultaneous problems

                          DO IPROB=1,MAX(1,IOPTS(29)*NPBTP)


C------------------------- IPBTP is the number of the current transport problem

                            IF (IOPTS(29).EQ.1) THEN
                               IPBTP=IPROB
                            ELSE
                               IPBTP=MAX(1,ISOLEQ(MAX(1,INTI),4))
                            ENDIF
                            IPBFL=ISOLEQ(MAX(1,INTI),3)

                            IF (IFLAGS(25).GT.0) THEN  
                                          
                                WRITE (746,*)
     &                               'ITERACION ',NUMITER, ' INTI ',INTI

                            END IF !IFLAGS(25).GT.0

                              CALL JAC_C
     ;          (ACTH     ,AFLU     ,AREA     ,ATRA(1,1,IPROB)
     ;          ,ATRADSC(1,1,IPROB)
     ;          ,ATRADSCF ,BETAC    ,BIBI     ,BUOYANCY ,CAUDAL
     ;          ,CAUX1(1,IPROB)    ,CAUX2(1,IPROB)    ,CCALAN(1,IPROB)
     ;          ,CCALIT(1,IPROB)   ,CFPAREL
     ;          ,CFPARNP  ,COORD    ,CREF     ,DAT_VD   ,DENSITY
     ;          ,DENSREF  ,DERC(1,1,1,IPROB)
     ;          ,DERC(1,1,1,MIN(IOPTS(30)*IPROB+1,NPBTP))
     ;          ,DERH     ,DFLU     ,DTIMEF   ,DTRA(1,1,IPROB)     
     ;          ,DTRADTRA
     ;          ,DVDP     ,PAR_DIR(31)        ,PAR_DIR(32)
     ;          ,FNT      ,GP_COORD ,GRADLOC  ,GRDFF    ,HAUX1
     ;          ,HAUX2    ,HBASE    ,HCALAN   ,HCALIT   ,HINI
     ;          ,IAD_S    ,IADD_S   ,IADN_S   ,IATRADSC_COLS
     ;          ,IATRADSC_ROWS      ,IAFD_S   ,IAFDD_S  ,IAFDN_S
     ;          ,IATRADSC_COLS      ,IATRADSC_ROWS      ,IBCOD
     ;          ,IBTCO(1,IPBTP)     ,1        ,IDIMAFLU ,IDIMATRA
     ;          ,IDIMBB
     ;          ,IDIMDERC ,IDIMDERH ,IDIMDFLU ,IDIMDQ   ,IDIMDTRA
     ;          ,IDIMFNT  ,IDIMQ    ,IDIMWORK ,IODIRECT ,IFLAGS
     ;          ,INEW     ,INEWT    ,INORPAR  ,INTI     ,0
     ;          ,IOCTRA   ,IODENS   ,IODIM    ,IOFLLI   ,IOFMLF
     ;          ,IOFMLT   ,IOINV    ,IOLD     ,IOLDT    ,IOLG_PAR
     ;          ,IOPTS(30),IOTRLI   ,IOPTS(31),IPAR_DIR ,IPROB
     ;          ,ISOLFL   ,ISOLTR   ,ISOZ     ,ITERM    ,ITPTVAR
     ;          ,ITYPAFLU ,ITYPATRA ,ITYPDFLU ,ITYPDTRA ,IVPAR
     ;          ,IXPARNP(1,1,IPBFL) ,IXPARNP(1,1,IPBTP) ,KXX      ,LDIM
     ;          ,LMXNDL   ,LNNDEL
     ;          ,LTYPE    ,LXPAREL(1,1,IPBTP)  ,LXPAREL(1,1,IPBFL) 
     ;          ,MAINF
     ;          ,MAXNB    ,MAXNBF   ,MAXPG    ,NBAND1   ,NFLAGS
     ;          ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR
     ;          ,NINT     ,NPAR     ,NPARALG  ,NPAREL   ,NPARF
     ;          ,NPARNP   ,NPBTP    ,NPPEL    ,NPPNP    ,NTYPAR
     ;          ,NTYPEL   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     ;          ,NZTRA    ,PARACD   ,PARZ     ,PAREL(1,1,IPROB)
     ;          ,PARNP(1,1,IPROB)
     ;          ,POINTWEIGHT        ,QXYZ     ,SOLUTION ,THETAF
     ;          ,THETAT   ,TINC     ,TINTERVOBS         ,VD
     ;          ,WATVOL(1,1,1,IPROB)   ,WORK(1,1,IPROB)
     ;          ,WSPECHEAT,WTHERMCON
     ;          ,WORK(1,NBAND1+1,IPROB),WORK(1,NBAND1+2 ,IPROB)
     ;          ,IDIMWGT   ,WGT_PAR  ,IPNT_PAR ,IPOS    ,DERIV
     &          ,IOCONSRC  ,DTRADFLU ,CFLU     ,DFLUDFLU,DFLUDTRA
     &          ,DNODALRH  ,DQDFLU   ,DQDTRA   ,LINMET  ,DPARELDH)


                          ENDDO   ! IPROB=1,MAX(1,IOPTS(29)*NPBTP)

C------------------------- INEWT and IOLDT are interchanged only if there is no
C------------------------- density dependent flow, else it is done after JAC_COUPLED.

                          IF (IODENS.EQ.0) THEN

                              INEWT = 3 - INEWT
                              IOLDT = 3 - IOLDT

                          END IF !IODENS.EQ.0


C------------------------- Coupled flow and transport inverse problem

                          IF (IODENS.EQ.1) THEN

                              CALL JAC_COUPLED
     &                 (A_COUPL_DSC        ,A_COUPL_DSCF       ,BCOUPLED
     &                 ,DERC     ,DERH     ,IA_COUPLED_DSC_COLS
     &                 ,IA_COUPLED_DSC_ROWS,IAD_D    ,IADD_D   ,IADN_D
     &                 ,IAFD_D  ,IAFDD_D   ,IAFDN_D  ,IDIMDERC ,IDIMDERH
     &                 ,IDIMWORK,IFLAGS    ,INEW     ,INEWT    ,INTI
     &                 ,IODIRECT,IPAR_DIR  ,ITERM    ,MAINF    ,MAXNB
     &                 ,MAXNBF  ,NBAND1    ,NFLAGS   ,NPAR     ,NPARALG
     &                 ,NUMNP   ,PAR_DIR   ,SOLUTION ,WORK     ,PARNAME
     &                 ,FILENAME)

C------------------------- INEW and IOLD and  INEWT and IOLDT are interchanged.

                              INEW = 3 - INEW
                              IOLD = 3 - IOLD

                              INEWT = 3 - INEWT
                              IOLDT = 3 - IOLDT

                          END IF !IODENS.EQ.1 .AND.IPROB.EQ.1

                     ELSE IF (IOINV.NE.1 .AND. INTI.GE.1) THEN

C------------------------- If initial conditions have been read in this 
C------------------------- interval derivatives must be set to zero (to use
C------------------------- in the next interval)

                        IF (ISOLEQ(INTI,2).EQ.2 .OR. 
     ;                      ISOLEQ(INTI,2).EQ.3)
     ;                        CALL ZERO_ARRAY(DERC,NUMNP*NPAR*2*
     ;                                         MAX(1,IOPTS(29)*NPBTP))
                     ENDIF   ! IOINV.NE.1 .AND. ISOLTR.GE.1

                  ENDIF    ! IOINV.GT.0 .AND. MIN_STOP.NE.1    Derivatives

C---------------------------------------------
C---------------------------------------------
C----- Ends Inverse problem RHS computing ----
C---------------------------------------------
C---------------------------------------------


C------------------------- Writes state variable and derivatives to check these

                  IF (IFLAGS(4).EQ.1 .AND. IOINV.GT.0
     &                               .AND. INDENDDT.EQ.1)THEN
                      IF (IOINV.EQ.1) THEN
                          CALL WRI_DERIV
     &                  (IOLD,INTI,NPARF,NPBFL,NUMNP,HCALIT,DERH)
                      ELSE IF (IOINV.EQ.2) THEN
                          CALL WRI_DERIV
     &                  (IOLDT,INTI,NPAR,NPBTP,NUMNP,CCALIT,DERC)
                      ELSE
                          CALL WRI_DERIV
     &                  (IOLD,INTI,NPARF,NPBFL,NUMNP,HCALIT,DERH)
                          CALL WRI_DERIV
     &                  (IOLDT,INTI,NPAR,NPBTP,NUMNP,CCALIT,DERC)
                      ENDIF
                  ENDIF

C------------------------- Calculating "observations"

                  IF ((IOINV.LT.0 .AND. IOPL.GE.1)
     &                .OR. (IOINV.GT.0 .AND. MIN_STOP.NE.1) )THEN
              
                        IPBFL=MAX(1,ISOLEQ(MAX(1,INTI),3)) !FLOW PROBLEM
                        IPBTP=MAX(1,ISOLEQ(MAX(1,INTI),4)) !TRANSPORT PROBLEM

                      CALL COMP_OBS
     ;(81         ,INEW     ,IOINV    ,IOLD     ,IOWRITE(6)
     ;,IOWRITE(5) ,NDEVS    ,NPAR     ,NPARF    ,NPBFL    ,NPBTP    
     ;,NUMNP      ,NUMTIT   ,NUMTNOD  ,NUMTOBS  ,TABSOLUT ,TINC     
     ;,CCALIT     ,CCALAN   ,DERC     ,DERH     ,DVOBS    ,HCALIT     
     ;,HCALAN     ,INDEXNOD ,IODEVICE ,NOOBSIT  ,TIT      ,TOBS     
     ;,VJAC       ,VOBSC    ,WTOBSN   ,WTOBST   ,INEWT    ,IOLDT
     ;,IOPTS(28)  ,IOPTS(29),IPBFL    ,IPBTP    ,TIME     ,NINT
     ;,IDIMDERH   ,IDIMDERC ,IFLAGS   ,NFLAGS)

                      IF (NDEVGEO.GT.0) THEN

C------------------------- Writes mesh for ERT if not already done.
                         IF (.NOT. FILEEX('HYDROMESH.vtk')) THEN

                            OPEN(123, File='HYDROMESH.vtk'
     &                              , Status='UNKNOWN')

                           CALL Write_VTK_Unstructured_Grid
     &        (COORD   ,0       ,KXX      ,LMXNDL   ,LNNDEL   ,LTYPE
     &        ,NUMEL   ,NUMNP   ,'HYDROMESH'        ,123)

                           CLOSE(123)
                           END IF  !No HYDDROMESH.vtk
  
                           CALL  COMP_OBS_GEO 
     ;     (CCALIT   ,ACTH     ,CFPAREL  ,DERC     ,DVOBS    ,IDIMDERC
     ;     ,IDIMWGT  ,IPNT_PAR ,IOLG_PAR ,IFLAGS   ,INEWT    ,INORPAR
     ;     ,IODEVICE ,IOINV    ,IOPTS(29),IPBTP    ,IVPAR    ,LXPAREL  
     ;     ,NDEVGEO  ,NDEVS    ,NFLAGS   ,NOOBSIT  ,NPAR     ,NPAREL
     ;     ,NPBMX    ,NPBTP    ,NTYPAR   ,NUMEL    ,NUMNP    ,NUMTIT
     ;     ,NUMTOBS  ,NZONE_PAR(19)      ,NZPAR    ,PARZ     ,TABSOLUT
     ;     ,TIT      ,VJAC     ,VOBSC)
                      END IF !NDEVGEO.GT.0
                      
                  END IF !IOINV.LT.0 .AND. IOPL.GE.1 ...

C-------------------------Updates WATVOL for the new time step
C-------------------------Only in INT>0 and transport is going to be solved
C-------------------------(INTI=0 initial WATVOl is needed)

          IF (INTI.GT.0 .AND. IOEQT.GT.1) THEN
             CALL WATVOL_UPDATE
     &           (HCALIT ,HCALAN ,IOPTS(31),KXX
     &           ,LMXNDL ,LNNDEL ,NPBTP  ,NPPEL    ,NUMEL
     &           ,NUMNP  ,PAREL  ,WATVOL)

          END IF !INTI.GT.0 .AND. IOEQT.GT.1

C------------------------ Assign arrays for next time step

                  IF (IOFLLI.NE.0 .AND. INTI.NE.0) THEN

                      CALL EQUAL_ARRAY
     &                    (HPREV2,HPREV1,NUMNP*MAX(1,IOPTS(28)*NPBFL))

                      CALL EQUAL_ARRAY
     &                    (HPREV1,HCALAN,NUMNP*MAX(1,IOPTS(28)*NPBFL))

                  END IF !IOFLLI.NE.0 .AND. INTI.NE.0

                  IF (IOTRLI.NE.0 .AND. INTI.NE.0) THEN

                      CALL EQUAL_ARRAY
     &                    (CPREV2,CPREV1,NUMNP*MAX(1,IOPTS(29)*NPBTP))

                   CALL EQUAL_ARRAY
     &                 (CPREV1,CCALAN,NUMNP*MAX(1,IOPTS(29)*NPBTP))

                  END IF !IOTRLI.NE.0 .AND. INTI.NE.0

                  IF (IOEQT.NE.2) THEN

                      CALL EQUAL_ARRAY
     &                    (HCALAN,HCALIT,NUMNP*MAX(1,IOPTS(28)*NPBFL))

                      CALL EQUAL_ARRAY
     &                    (HAUX1,HCALIT,NUMNP*MAX(1,IOPTS(28)*NPBFL))

                      CALL ZERO_ARRAY
     &                    (HAUX2,NUMNP*MAX(1,IOPTS(28)*NPBFL))

                  END IF !IOEQT.NE.2

                  IF (IOEQT.NE.1)  THEN

                      CALL EQUAL_ARRAY
     &                    (CCALAN,CCALIT,NUMNP*MAX(1,IOPTS(29)*NPBTP))

                      CALL EQUAL_ARRAY
     &                    (CAUX1,CCALIT,NUMNP*MAX(1,IOPTS(29)*NPBTP))

                      CALL ZERO_ARRAY
     &                    (CAUX2,NUMNP*MAX(1,IOPTS(29)*NPBTP))

                  END IF !IOEQT.NE.1

              END IF !(IOCONVGL.EQ.1)   
          ENDDO  !WHILE (INTI.NE.0 .AND. INDENDDT.EQ.0)

      END DO ! INTI=0,NINT-1
C-------------------------------------------
C-------------------------------------------
C----- Ends time loop ----------------------
C-------------------------------------------
C-------------------------------------------

      
C------------------------- Returns index of first integration times to
C------------------------- original values. Only if they were updated
C------------------------- at COMP_OBS

      IF ((IOINV.LT.0 .AND. IOPL.GE.1) .OR. 
     &    (IOINV.GT.0 .AND. MIN_STOP.NE.1)) THEN

C----------------------- Recovers initial value of IODEVICE(*,2)

         DO ND=1,NDEVS
           IODEVICE(ND,2)=IODEVICE(ND,10)
         ENDDO

      END IF


      IF (IODIV.EQ.1) STOP 'DIVERGENCE'


      IF(IFLAGS(3).EQ.1) CALL IO_SUB('SIM_JAC',1)
      RETURN

      END SUBROUTINE SIM_JAC
