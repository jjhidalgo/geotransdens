      SUBROUTINE TRANSPORT
     &          (AREA       ,ATRA    ,ATRADSC   ,ATRADSCF
     &          ,BETAC
     &          ,BIBI     ,BTRA       ,CAUDAL  ,CAUX1
     &          ,CAUX2    ,CCALAN     ,CCALIT  ,CFLU
     &          ,CFPAREL  ,CFPARNP    ,CPREV1    ,CPREV2
     &          ,CREF     ,DABSMTR    ,DAT_VD
     &          ,DCITMX
     &          ,DELTACGL ,DELTAITER  ,DENSITY
     &          ,DENSREF  ,DPARELDC,DPARELDH  ,DQDFLU
     &          ,DQDTRA   ,DRELCMX    ,DRELMXTR,DTIMET
     &          ,DTRA     ,DTRADFLU   ,DTRADTRA,DVDH      ,DVDC,DWDH
     &          ,EPSFLU   ,EPSTRA     ,FILENAME,FNT
     &          ,GRDFF
     &          ,IAD_S    ,IADD_S     ,IADN_S  ,IAFD_S    ,IAFDD_S
     &          ,IAFDN_S  ,IATRADSC_COLS
     &          ,IATRADSC_ROWS,IBTCO   ,IDELCGL   ,IDIMATRA
     &          ,IDIMBB   ,IDIMCFLU,IDIMDENS   ,IDIMDFLU,IDIMDQ,IDIMDTRA
     &          ,IDIMFNT  ,IDIMQ   ,IDIMWORK
     &          ,IFLAGS   ,INARR      ,INCLK    ,INCON
     &          ,INDENDDT ,INDFLTR    ,INDSSTR ,INORPAR
     &          ,INPWR    ,INTI       ,IOCONSRC
     &          ,IOCONVT  ,IODENS     ,IODIM   ,IODIRECT
     &          ,IOFMLT   ,IOINV      ,IOITERTREND,IOPINITC  ,IORDCH
     &          ,IORTS    ,IOTRLI     ,IOTRS   ,IOVRWC    ,IOWAR
     &          ,IOWNR    ,IPAR_DIR   ,IPROB   ,IREDTIMC  ,IRESCMAX
     &          ,ISOLEQ   ,ISOLTR  ,ISPARSE   ,ITERCHNGFL
     &          ,ITERCHNGGL           ,ITERCHNGTR         ,ITERCONVFL
     &          ,ITERCONVTR           ,ITERFL  ,ITERGL    ,ITERM
     &          ,ITERTOTTR ,ITERTR  ,ITPTVAR   ,ITYPATRA
     &          ,ITYPTRADSC,ITYPCFLU  ,ITYPDTRA
     &          ,IUCAL    ,IOCALMAT
     &          ,IXPARNP  ,KXX        ,LDIM    ,LINMET    ,LMXNDL
     &          ,LNNDEL   ,LTYPE      ,LXPAREL ,MAINF     ,MAXNB
     &          ,MAXNBF   ,NBAND   ,NBAND1    ,NCONVI
     &          ,NCONVITR
     &          ,NFLAGS   ,NFNL       ,NFNLPAR ,NFNLPRG   ,NFNLTIP
     &          ,NFTPAR   ,NINT       ,NPARALG ,NPAREL    ,NPARNP
     &          ,NPPEL    ,NPPNP      ,NROW    ,NTDMT
     &          ,NTYPAR   ,NUMEL      ,NUMITER ,NUMNP     ,NZONE_PAR
     &          ,NZPAR    ,NZPRG      ,PAR_DIR ,PARACD    ,PARC
     &          ,PAREL    ,PARNP      ,PRGC    ,QXYZ      ,RESCMAX
     &          ,RESCMAXOLD           ,RESIDMXTR          ,SOLUTION
     &          ,SOURCE   ,TABSOLUT   ,THETAT  ,TICAL     ,TICALAN
     &          ,TIME     ,TINC       ,TINCINI ,TINCLAST  ,TINTERVOBS
     &          ,TOLD     ,VD         ,WATVOL  ,WORK
     &          ,WSPECHEAT,WTHERMCON  ,XNORVD  ,ZEROT,ACTH,NPBTP)
********************************************************************************
* EXTERNAL VARIABLES: ARRAYS
*
*  ISOLEQ       Array containing the type of head/concentration
*                 solution desired for the user at each obs. time
*                 0. Transient
*                 1. Steady state.
*                 2. Read initial conditions
*                 3. Null initial conditions
*                 4. Not to solve.
*
*
* EXTERNAL VARIABLES: SCALAR
*
*  IENTRY       Indicates if is the first time access to the
*                 subroutine
*  IORDCH       Option for chain reactions
*  IORTS        Transport regime.
*                 0. Steady state transport.
*                 1. Transient trasnport with prescribed initial conditions.
*                 2. Transient trasnport with steady-state initial conditions.
*  IOVRWC       Equal to IOPTS(31).
*                 1.WATVOL calculated elementwise
*                 2.WATVOL calculated nodewise.
*                 DENSITY is calculate in the same way that WATVOL.
*  ISOLFL       Status of flow problem in current time step.
*                 0. Not solved.
*                 1. Steady flow solved.
*                 2. Transient flow solved.
*  ISOLTR       Status of transport problem in current time step.
*                 0. Not solved.
*                 1. Steady transport solved.
*                 2. Transient transport solved.
*  IREDTIMC     Indicator variable (=1 --> reduction time must be carried out
*
*
* INTERNAL VARIABLES: SCALARS.
*
*  IOCALMAT     If 1, ATRA and DTRA matrices must be computed.
*
********************************************************************************


      IMPLICIT NONE

      INTEGER*4::I,IOCALMAT,IOCONSRC,IOCONVT,INCLK,INTI,IONEWT,IDIMDENS,
     &           ITERTR, ITERTOTTR,IOTRLI,IODENS,IDIMFNT,INDFLTR,IOINV,
     &           IOPINITC,NCONVI,NINT,NPARALG,NPARNP,NTYPAR,NUMITER,
     &           NUMNP,NZPAR,ITERCHNGFL,ITERCHNGTR,ITERCHNGGL,
     &           ITERCONVFL,ITERCONVTR,ITERFL,ITERGL,INDSSTR,
     &           IOVRWC,LMXNDL,MAINF,IOFMLT,NFLAGS,
     &           NFNL,NPAREL,NPPEL,NUMEL,NPPNP,INCON,IDIMBB,
     &           NZPRG,IDIMQ,IDIMDQ,IDIMWORK,IODIM,IDIMDTRA,
     &           IDIMATRA,IORTS,MAXNB,MAXNBF,
     &           IATRADSC_ROWS,IATRADSC_COLS,IORDCH,IPROB,NTDMT,
     &           ITPTVAR,INARR,IDIMDFLU,IDIMCFLU,
     &           NBAND1,IODIRECT,IDESCTRA,NBAND,ITERM,
     &           ISPARSE,IRESCMAX,IDELCMAX,INDENDDT,IOWNR,
     &           NCONVITR,IREDTIMC,IUCAL,IOWAR,INPWR,IERROR,
     &           NROW,ISOLTR,ITYPATRA ,ITYPTRADSC ,ITYPDTRA,ITYPCFLU,
     &           IOITERTREND,IOTRS,IDELCGL,NUMDIVC,NPBTP

      INTEGER*4::LINMET(3,2),ISOLEQ(NINT,4),
     &           IBTCO(NUMNP),IXPARNP(NUMNP,NPARNP),NFTPAR(NZPAR),
     &           LNNDEL(NUMEL),IFLAGS(NFLAGS),IPAR_DIR(NPARALG),
     &           KXX(LMXNDL,NUMEL),LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR),
     &           NFNLPRG(8,NFNL),NFNLTIP(NFNL),NZONE_PAR(NTYPAR),
     &           LDIM(NUMEL),LTYPE(NUMEL),IADN_S(NUMNP),
     &           IAD_S(MAXNB, NUMNP),IADD_S(NUMNP),INORPAR(NTYPAR),
     &           IAFD_S(MAXNBF,NUMNP),IAFDD_S(NUMNP),IAFDN_S(NUMNP)

      REAL*8::THETAT,TINC,TINCINI,TINCLAST,TINTERVOBS,TICAL,TICALAN,
     &        EPSFLU,EPSTRA,DTIMET,TOLD,BETAC,
     &        CREF,DENSREF,RESCMAX,RESCMAXOLD,
     &        DELCMAX,DELCMAXOLD,DRELCMX,DCITMX,ZEROT,DABSMTR,
     &        DRELMXTR,RESIDMXTR,TABSOLUT,DELTACGL,WSPECHEAT,WTHERMCON

      REAL*8::BTRA(NUMNP),TIME(NINT),FNT(IDIMFNT,NINT),
     &        CFPARNP(NUMNP,NPARNP),CAUX1(NUMNP),CAUX2(NUMNP),
     &        CCALIT(NUMNP),CCALAN(NUMNP),CPREV1(NUMNP),CPREV2(NUMNP),
     &        PARC(NZPAR),AREA(NUMEL),
     &        CFPAREL(NUMEL,NPAREL),PARACD(3,NFNL),
     &        PAREL(NUMEL,NPPEL),PRGC(NZPRG),
     &        WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3),
     &        PARNP(NUMNP,NPPNP),VD(IODIM,NUMEL),
     &        ATRA(NUMEL, IDIMATRA),BIBI(IDIMBB,NUMEL),
     &        DAT_VD(IODIM,IDIMDQ,NUMEL),DVDH(LMXNDL,IODIM,NUMEL),
     &        DTRADFLU(NUMEL,LMXNDL*LMXNDL),
     &        DTRADTRA(NUMEL,LMXNDL*LMXNDL),
     &        DVDC(LMXNDL,IODIM,NUMEL),GRDFF(IODIM,LMXNDL,NUMEL),
     &        DTRA(NUMEL,IDIMDTRA),XNORVD(NUMEL),QXYZ(IDIMQ,NUMEL),
     &        CAUDAL(NUMNP),SOURCE(NUMNP),
     &        ATRADSC(IATRADSC_ROWS,IATRADSC_COLS),
     &        DQDFLU(NUMEL,LMXNDL*LMXNDL),
     &        ATRADSCF(MAXNBF,NUMNP),
     &        DQDTRA(NUMEL,LMXNDL*LMXNDL),
     &        CFLU(NUMEL,IDIMDFLU),SOLUTION(NUMNP),
     &        PAR_DIR(NPARALG),DPARELDC(NPPEL,NUMEL),
     &        DPARELDH(NPPEL,NUMEL),
     &        DENSITY(IDIMDENS),WORK(IDIMWORK),DELTAITER(NUMNP),
     &        DWDH(MAX(1,(IOVRWC-1)*2*LMXNDL),NUMEL),ACTH(NUMEL)

      CHARACTER*20 FILENAME (17)

C------------------------- Internal

      INTEGER*4::IOREGIMEN

      REAL*8::FACTOR

      REAL*8::XPARAM(8),DUMMY(1)

      REAL,ALLOCATABLE::DERB(:)


      IF (IFLAGS(3).NE.0) CALL IO_SUB('TRANSPORT',0)

C------------------------- Initilizes convergence and timestep reduction
C------------------------- indicators.
      IOCONVT = 0
      IREDTIMC= 0
      IOITERTREND = 0
      DELCMAX = HUGE(1D0)
      RESCMAX = HUGE(1D0)

      IF (IODENS.EQ.1) THEN
         DELTAITER = 0D0
      ENDIF

C------------------------- States the type of solution

      IF (INTI.EQ.0) THEN

          IF (IORTS.EQ.0 .OR. IORTS.EQ.2) THEN

              IOREGIMEN = 1 !Steady state or initial steady state
              INDSSTR = 0

          ELSE IF (IORTS.EQ.1) THEN

              IOREGIMEN = 4 !Initial conditions, but they have already
                            !been read during input data reading, so
                            !nothing is done.

          END IF !IORTS.EQ.0 ...

      ELSE

          IOREGIMEN = ISOLEQ(INTI,2)

          IF (IOREGIMEN.EQ.0) THEN
              INDSSTR = 1
          ELSE
              INDSSTR = 0
          END IF !IOREGIMEN.EQ.0

      END IF !INTI.EQ.0

C------------------------- Selects the type of solution for this timestep.

      SELECT CASE(IOREGIMEN)

C------------------------- Steady or transient solution
      CASE(0,1)

C------------------------- Current time-step iteration counter.
          ITERTR = 0

          IF (IOREGIMEN.EQ.1) THETAT = 1D0 ! To solve SS apropriately

C------------------------- If solving a non-linear problem,
C------------------------- state variable is initialized.
C------------------------- If solving a density-dependent problem,
C------------------------- inital value of has already been computed
C------------------------- while solving flow equation.
          IF (IOTRLI.NE.0 .AND. IODENS.EQ.0) THEN
              CALL STATE_VARIABLE_INIT
     &        (IDIMFNT  ,INDFLTR  ,INTI     ,IOINV    ,IOPINITC ,NCONVI
     &        ,NINT     ,NPARALG  ,NPARNP   ,NTYPAR   ,NUMITER  ,NUMNP
     &        ,NZPAR    ,THETAT   ,TICAL    ,TICALAN  ,TIME     ,TINC
     &        ,TINCINI  ,TINCLAST ,TINTERVOBS         ,BTRA     ,CFPARNP
     &        ,FNT      ,IBTCO    ,INORPAR  ,IXPARNP  ,NFTPAR   ,PARC
     &        ,CAUX1    ,CCALIT   ,CCALAN   ,CCALIT ,CPREV1   ,CPREV2)

              CALL UPDATE_VECTORS_AUX
     ;(NUMNP     ,TINC    ,THETAT    ,CCALAN   ,CCALIT   ,CAUX1  ,CAUX2)

          END IF

C------------------------- Iterates through transport equation solution while
C------------------------- convergence is not achieved or time step must not
C------------------------- be reduced.

          DO WHILE (IOCONVT.NE.1 .AND. IREDTIMC.NE.1
     &              .AND. IOITERTREND.NE.1)

C------------------------- Increments iteration counters.
              CALL INCREMENT_ITER
     &            (IOITERTREND ,ITERTOTTR ,ITERTR ,ITERCONVTR)


C------------------------- Choose linearization method.
              CALL CHOOSE_LINEARIZT_METHOD
     &            (IODENS    ,ITERCHNGFL ,ITERCHNGGL ,ITERCHNGTR
     &            ,ITERCONVFL,ITERCONVTR ,ITERFL     ,ITERGL
     &            ,ITERTR    ,LINMET)

C------------------------- Computes IONEWT.
              IF (LINMET(2,2).EQ.2 .OR. LINMET(3,2).EQ.2) THEN

                  IONEWT = 1
              ELSE

                  IONEWT = 0

              END IF

C------------------------- Computes element values of transport parameters.

              CALL COMP_PARAM_TRA
     &            (CCALAN   ,CCALIT   ,CFPAREL  ,CFPARNP  ,DPARELDC
     &            ,DPARELDH ,DTIMET   ,EPSTRA   ,FNT      ,IBTCO
     &            ,IDIMFNT  ,IFLAGS   ,INCLK    ,INCON    ,INDSSTR
     &            ,INORPAR  ,INTI     ,IODENS   ,IOFMLT   ,IOTRLI
     &            ,IPAR_DIR ,IXPARNP  ,KXX      ,LMXNDL   ,LNNDEL
     &            ,LXPAREL  ,MAINF    ,NFLAGS   ,NFNL     ,NFNLPAR
     &            ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT     ,NPARALG
     &            ,NPAREL   ,NPARNP   ,NPPEL    ,NPPNP    ,NTYPAR
     &            ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR    ,NZPRG
     &            ,PARACD   ,PARC     ,PAREL    ,PARNP    ,PRGC
     &            ,XPARAM   ,THETAT   ,TINC     ,TINTERVOBS)

C------------------------- Computes ATRA matrix if needed.

               IF (IOCALMAT.NE.0) THEN

C------------------------- Computes ATRA matrix.

                  CALL COMP_ATRA
     &      (AREA     ,ATRA     ,BETAC    ,BIBI     ,CAUX1    ,DAT_VD
     &      ,DENSITY  ,DENSREF  ,DTRADFLU ,DTRADTRA ,DVDH     ,DVDC
     &      ,DWDH     ,EPSFLU   ,EPSTRA   ,GRDFF    ,IDIMBB   ,IDIMDENS
     &      ,IDIMDQ   ,IDIMQ    ,IODENS   ,IODIM    ,IOINV    ,IOVRWC
     &      ,ITPTVAR  ,LINMET   ,LMXNDL   ,NPPEL    ,NTYPAR   ,NUMEL
     &      ,NUMNP    ,KXX      ,LDIM     ,LNNDEL   ,LTYPE    ,NZONE_PAR
     &      ,PAREL    ,QXYZ     ,THETAT   ,VD       ,WATVOL   ,WSPECHEAT
     &      ,WTHERMCON,XNORVD)

C------------------------- Adds areal recharge and first order decay terms
C------------------------- to ATRA matrix

                  CALL COMP_ATRA_REC_FOD
     &                (ACTH     ,AREA     ,ATRA     ,BETAC    ,CAUX1
     &                ,CREF     ,DENSREF  ,DTRADFLU ,DTRADTRA ,DWDH
     &                ,IODENS   ,IONEWT   ,IOVRWC   ,ITPTVAR  ,LINMET
     &                ,LMXNDL   ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP
     &                ,KXX      ,LNNDEL   ,NZONE_PAR,PAREL    ,THETAT
     &                ,WATVOL)

C------------------------- Initializes RHS (BTRA) to zero.

                  BTRA = 0D0

C------------------------- Initializes DTRA to zero if needed.
                  IF (NTDMT.NE.0 .OR. IOREGIMEN.EQ.0) THEN

                      DTRA = 0D0

                  END IF
C------------------------- Computes matrix diffusion contributions to BTRA
C------------------------- and DTRA, if any.

                  IF (NTDMT.NE.0)  THEN
                      CALL COMP_CNTR_DMT (BTRA, DTRA, NUMNP,  NUMEL,
     ;                     IDIMDTRA, KXX, LNNDEL, LMXNDL, IXPARNP(1,9))

C------------------------- Correction when solving heat transport

                      IF (ITPTVAR .EQ. 1) THEN

                          FACTOR = WSPECHEAT

                          IF (IODENS.EQ.0) THEN
C------------------------- If density is constant, water density appears
C------------------------- in this term (see TRANSDENS Guia rapida, 3.7).
                              FACTOR = DENSREF*FACTOR
                          END IF !IODENS.EQ.0

                          DTRA = DTRA/FACTOR
                          BTRA = BTRA/FACTOR

                      END IF !ITPTVAR.EQ.1


                  END IF !NTDMT.NE.0
C------------------------- Computes DTRA matrix if needed

                 IF (IOREGIMEN.EQ.0) THEN
                 
                   CALL COMP_DTRA
     &                (ACTH     ,AREA     ,BETAC    ,CAUX1    ,CAUX2
     &                ,CCALAN   ,CCALIT   ,CREF     ,DENSITY  ,DENSREF
     &                ,DTRA     ,DTRADFLU ,DTRADTRA ,DWDH     ,IDIMDTRA
     &                ,IODENS   ,IOINV    ,IOVRWC   ,ITPTVAR  ,KXX
     &                ,LINMET   ,LMXNDL   ,LNNDEL   ,NUMEL    ,NUMNP
     &                ,PAREL(1,14)       ,THETAT   ,WATVOL   ,WSPECHEAT)

                END IF !IOREGIMEN.EQ.0

              END IF  ! IOCALMAT.NE.0

C------------------------- If time increment has changed or
C------------------------- if solving transient transport or
C------------------------- matrices have changed, assembles transport
C------------------------- left hand side.

              IF (LINMET(2,2).NE.0) THEN
                  IF (DABS(TINC-TOLD).GE.1.D-6*TINC .OR. IOTRS.NE.0
     &                     .OR. IOCALMAT.EQ.1) THEN


C-------------------------  DBTRADTRA <=> PARNP(1,6) (only conc. leakage
C-------------------------  has derivatives w. r. t. conc.)


                      IF (IONEWT.GT.0) THEN

                          ALLOCATE (DERB(NUMNP))

                          CALL COMP_DER_BTRA
     &                        (CAUX1    ,DENSREF  ,DERB     ,IBTCO
     &                        ,IOCONSRC ,IODENS   ,ITPTVAR  ,NPPNP
     &                        ,NUMNP    ,PARNP    ,THETAT   ,WSPECHEAT)

                    END IF !IONEWT.GT.0

                    CALL ASSEMB_LHS
     &                  (IDIMATRA   ,NUMEL      ,IDIMATRA   ,NUMEL
     &                  ,1          ,NUMNP      ,IATRADSC_COLS
     &                  ,IATRADSC_ROWS          ,IDIMDTRA  ,NUMEL
     &                  ,INDSSTR    ,IONEWT     ,ITYPATRA  ,ITYPTRADSC
     &                  ,ITYPDTRA   ,4          ,1         ,LMXNDL
     &                  ,MAXNB      ,NBAND      ,NUMEL     ,NUMNP
     &                  ,THETAT     ,TINC       ,ATRA      ,ATRADSC
     &                  ,DTRA       ,DTRADTRA   ,DERB      ,IAD_S
     &                  ,IADD_S     ,IADN_S     ,KXX       ,LNNDEL
     &                  ,1D0)

                          
C------------------------- Sets mass flow boundary conditions
C------------------------- in LHS

                    CALL COMAT_BC
     &                  (ATRADSC  ,BETAC    ,CAUDAL   ,CAUX1    ,CREF
     &                  ,DENSREF  ,DQDFLU   ,DQDTRA   ,EPSTRA
     &                  ,IAD_S    ,IADD_S   ,IADN_S   ,IBTCO
     &                  ,IATRADSC_COLS      ,IATRADSC_ROWS      ,0
     &                  ,1        ,0        ,IODENS   ,IONEWT
     &                  ,ITYPTRADSC         ,KXX      ,LMXNDL   ,LNNDEL
     &                  ,MAXNB    ,NUMNP    ,NBAND    ,NPPNP    ,NUMEL
     &                  ,NUMNP    ,PARNP    ,THETAT   ,AREA
     &                  ,NZONE_PAR,PAREL    ,NPPEL    ,NTYPAR)

                    IDESCTRA = 2

                    IF (ALLOCATED(DERB)) THEN

                          DEALLOCATE(DERB)

                    END IF !ALLOCATED(DERB)

                  ELSE

                      IDESCTRA = 1

                  END IF !DABS(TINC-TOLD).GE.1.

              END IF !LINMET(2,2).NE.2


              IF (LINMET(2,2).NE.0) THEN !(if tpt is to be solved)

C------------------------- Sets conc. leakage boundary condition
C------------------------- in LHS

                  IF (NZONE_PAR(18).GT.0) THEN !NZCLK.GT.0

                      CALL PRESC_LEAK_BC
     &                    (ATRADSC    ,DUMMY      ,BETAC      ,CAUX1
     &                    ,CREF       ,DENSREF    ,IAD_S      ,IADD_S
     &                    ,IADN_S     ,IBTCO      ,IATRADSC_COLS
     &                    ,IATRADSC_ROWS          ,6          ,IOCALMAT
     &                    ,1          ,4          ,IOCONSRC   ,IODENS
     &                    ,IONEWT     ,ITPTVAR    ,ITYPTRADSC ,KXX
     &                    ,LMXNDL     ,LNNDEL     ,MAXNB      ,NBAND
     &                    ,NPPNP      ,NUMEL      ,NUMNP      ,PARNP
     &                    ,THETAT     ,CAUX1      ,WSPECHEAT)

                 END IF !NZONE_PAR(18).GT.0

C------------------------- Input mass

                  IF (IOCONSRC.EQ.1 .AND. ANY(IBTCO(:).EQ.4)) THEN

                      CALL INPUT_MASS
     &                    (ATRADSC  ,IAD_S    ,IADD_S   ,IADN_S
     &                    ,IATRADSC_COLS      ,IATRADSC_ROWS
     &                    ,IBTCO    ,ITYPTRADSC         ,KXX
     &                    ,LMXNDL   ,LNNDEL   ,MAXNB    ,NBAND
     &                    ,NPPNP    ,NUMEL    ,NUMNP    ,PARNP
     &                    ,THETAT)


                 END IF !IOCONCSRC.EQ.1 .AND. ANY(IBTCO(:).EQ.4)

             END IF !LINMET(2,2).NE.0


C------------------------- If there are radioactive chains and
C------------------------- current problem is not the first one,
C------------------------- contribution of its elder in radioactive
C------------------------- chain is added to the left hand side.

              IF (IORDCH.EQ.1 .AND. IPROB.NE.1) THEN

                  DO I=1,NUMNP
                      BTRA(I)=BTRA(I)+SOURCE(I)
                  END DO

              END IF

C------------------------- Assembles transport Right Hand Side
C------------------------- (prescribed concentration boundary
C------------------------- conditions NOT included).


              CALL ASSEMB_RHS
     &            (ATRA     ,BTRA     ,CAUX2    ,CFLU    ,DTRA
     &            ,iad_s    ,iadn_s   ,IDIMATRA
     &            ,IDIMCFLU ,IDIMDTRA
     &            ,1        ,INDSSTR  ,IODENS   ,IONEWT
     &            ,ITYPATRA ,ITYPDTRA ,ITYPCFLU
     &            ,KXX      ,LNNDEL   ,LMXNDL   ,NUMEL    ,NUMNP
     &            ,THETAT   ,TINC     ,CAUX1    ,CAUX2    ,CCALAN)

C------------------------- Adds areal recharge and boundary
C------------------------- conditions to RHS (prescribed
C------------------------- concentration included).

              CALL COMP_BTRA
     &            (AREA     ,BETAC    ,BTRA     ,CAUDAL   ,CAUX1
     &            ,CCALAN   ,CREF     ,DENSREF  ,IBTCO    ,INARR
     &            ,IOCONSRC ,IODENS   ,IONEWT   ,ITPTVAR  ,KXX
     &            ,LMXNDL   ,LNNDEL   ,LXPAREL  ,NPAREL   ,NPPEL
     &            ,NPPNP    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &            ,PAREL    ,PARNP    ,THETAT   ,WSPECHEAT)

C------------------------- If not solving flow and transport toghether...

              IF(LINMET(2,2).NE.0) THEN

C------------------------- ...Prescribed concentration boundary conditions
C------------------------- are  set to the system in LHS...

                  IF (ANY(IBTCO.EQ.1)) THEN
                      CALL PRESC_BC
     &            (ATRA          ,ATRADSC       ,BTRA          ,IDIMATRA
     &            ,IATRADSC_COLS ,IATRADSC_ROWS ,IADD_S        ,IADN_S
     &            ,IBTCO         ,1             ,1             ,INDSSTR
     &            ,IONEWT        ,ISPARSE       ,ITYPTRADSC    ,KXX
     &            ,LMXNDL        ,LNNDEL        ,NBAND1        ,NPPNP
     &            ,NUMEL         ,NUMNP         ,PARNP         ,THETAT
     &            ,CCALIT)
                  END IF

C------------------------- ...Computes maximum residue...


                      CALL COMPUTE_RESID_MAX
     &                    (IAD_S    ,IADN_S   ,1        ,0
     &                    ,IRESCMAX ,NUMNP    ,NUMNP
     &                    ,RESCMAX  ,RESCMAXOLD         ,IBTCO    ,BTRA
     &                    ,ATRADSC  ,KXX      ,LNNDEL   ,NUMEL
     &                    ,IATRADSC_COLS      ,IATRADSC_ROWS
     &                    ,ITYPTRADSC         ,CCALIT   ,SOLUTION
     &                    ,LMXNDL   ,IONEWT)


C------------------------- ...system is solved...

                  IF (IODIRECT.EQ.0) SOLUTION = CCALIT

                  IF (IFLAGS(26).GT.0) THEN

                      WRITE(622,*) 'dTdT en INTI = ', INTI
                      CALL WRITE_SQR_MATRIX
     &                    (ATRADSC  ,IATRADSC_COLS      ,IATRADSC_ROWS
     &                    ,IAD_S    ,IADN_S   ,ISPARSE  ,0
     &                    ,MAXNB    ,NUMNP    ,NBAND    ,622)

                      WRITE(622,*) ''

                      WRITE(332,*) 'BTRA en INTI = ', INTI
                      DO I=1,NUMNP

                          WRITE(332,*) BTRA(I)

                      END DO !I=1,NUMNP

                      WRITE(332,*)

                  END IF !IFLAGS(26).GT.0

                  CALL SOLVE
     &                (IATRADSC_ROWS      ,IATRADSC_COLS      ,NUMNP
     &                ,IDESCTRA ,IPAR_DIR(21)       ,IDIMWORK ,IODIRECT
     &                ,2        ,INTI     ,IPAR_DIR(18)
     &                ,IPAR_DIR(22)       ,0        ,ITERM
     &                ,IPAR_DIR(15),MAINF    ,NBAND1   ,MAXNBF
     &                ,IPAR_DIR(20)       ,IPAR_DIR(17)
     &                ,IPAR_DIR(16)       ,PAR_DIR(36)
     &                ,PAR_DIR(37)        ,PAR_DIR(38)        ,ATRADSC
     &                ,ATRADSCF ,BTRA     ,IAD_S    ,IADD_S   ,IADN_S
     &                ,IAFD_S   ,IAFDD_S  ,IAFDN_S  ,WORK     ,SOLUTION)

C------------------------- ...iteration updated and convergence
C------------------------- indicators computed...

                  CALL UPDATE_STATE_VARIABLE
     &           (DELCMAX  ,DELTACGL  ,DELCMAXOLD,DRELCMX  ,DCITMX
     &           ,1        ,IDELCGL   ,IDELCMAX  ,IOTRLI   ,0
     &           ,IODENS   ,IONEWT    ,NUMNP     ,ZEROT    ,DELTAITER
     &           ,SOLUTION ,CCALAN    ,CCALIT)

C------------------------- ...convergence and divergence checked...

                  CALL CHECK_CONVERGENCE
     &                (DABSMTR  ,DELCMAX ,DRELMXTR ,DRELCMX  ,IOCONVT
     &                ,INDENDDT ,IOWNR   ,NCONVITR ,RESIDMXTR,RESCMAX
     &                ,IOTRLI   ,IONEWT)




                  CALL CHECK_DIVERGENCE
     &        (DELCMAX  ,DELCMAXOLD,IONEWT,IREDTIMC,IPAR_DIR(27),NUMDIVC
     &        ,RESCMAX ,RESCMAXOLD)

C------------------------- If using Newton's method, write information.

                  IF (IONEWT.EQ.1) THEN
                      CALL WRITE_NONLIN_INFO
     &                    (DELCMAX  ,'TRANSPORT' ,IDELCMAX ,INDENDDT
     &                    ,IREDTIMC,IOWNR    ,IRESCMAX ,ITERTR  ,MAINF
     &                    ,NCONVI   ,RESCMAX  ,TABSOLUT ,TINC)
                  END IF

C------------------------- auxiliar vectors updated.

                  CALL UPDATE_VECTORS_AUX
     &                (NUMNP    ,TINC    ,THETAT    ,CCALAN   ,CCALIT
     &                ,CAUX1    ,CAUX2)

C------------------------- Updates Io and In functions as a consequence
C------------------------- of a concentration variation when
C------------------------- matrix diffusion is being considered. However
C------------------------- if transport inverse problem is being solved
C------------------------- this update MUST be done after the updating
C------------------------- of the derivatives of Io and In in JAC_C

       IF (NTDMT.NE.0 .AND. IOINV.LE.1)
     ;            CALL UPD_REC_DMT (CCALIT ,CCALAN ,NUMNP,IXPARNP(1,9))

              ELSE

                  IOCONVT = 1

              END IF !LINMET(2,2).NE.0


          END DO !WHILE (IOCONVT.NE.1 .AND. IREDTIMC.NE.1)

C------------------------- IOCALMAT must be set to zero for the next time step,
C------------------------- if density is constant and if we are in the last problem
C------------------------- (because in simultaneous problems in the same time step
C------------------------- matrices have to be computed for each problem)
c-parche-probs
c-no entiendo el comentario de arriba. Esto hace que no funcionen
c los problemas sucesivos. Lo comento.
c            IF (IODENS.NE.1 .AND. (IORDCH.EQ.0 .OR. IPROB.EQ.NPBTP) )
c     ;                                      IOCALMAT=0
c-fin-parche probs
C------------------------- The value of solved problem
C------------------------- indicator variable is set.

          IF (IOCONVT.EQ.1) THEN

               ISOLTR = INDSSTR + 1

C------------------------- In reactive chains, sources for the next
C------------------------- problem are computed

              IF (IORDCH.EQ.1 .AND. IPROB.NE.NPBTP)
     ;            CALL SOURCE_RDCH
     ; (LMXNDL   ,NPPEL    ,NUMEL    ,NUMNP    ,AREA     ,CAUX1
     ; ,KXX      ,LNNDEL   ,PAREL    ,SOURCE   ,WATVOL   ,IOVRWC)

            ELSE

                ISOLTR = 0

            END IF !IOCONVT.EQ.1

C------------------------- Read initial conditons.

      CASE(2)
C------------------------- If solving a density-dependent problem,
C------------------------- initial conditons have already been read
C------------------------- while solving flow equation.

          IF (IODENS.EQ.0) THEN
              CALL ENDATINICOND_AUX
     &            (NUMNP    ,2        ,IUCAL    ,MAINF    ,IOWAR
     &            ,INPWR    ,IERROR   ,NROW     ,FILENAME ,CCALIT)
              IOCONVT=1
          END IF

          ISOLTR = 0
C-------------------------Set initial conditions to zero.

      CASE(3)

          DO I=1,NUMNP
                 CCALIT(I)=0.D0
          ENDDO

          IOCONVT = 1
          ISOLTR = 0

C------------------------- Not solve.

      CASE(4)

          IOCONVT = 1
          ISOLTR = 0

      END SELECT !ISOLEQ(INTI,2)

      IF (IFLAGS(3).NE.0) CALL IO_SUB('TRANSPORT',1)

      END SUBROUTINE TRANSPORT
