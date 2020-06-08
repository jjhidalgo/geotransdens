      SUBROUTINE FLOW_EQN
     &(IAFLUDSC_COLS   ,IAFLUDSC_ROWS   ,IDIMAFLU
     &,IDIMBB          ,IDIMCFLU        ,IDIMDENS    ,IDIMDFLU
     &,IDIMFNT         ,IDIMGRAVEL  ,IDIMWORK
     &,IOCONVF         ,IODIRECT
     &,IENTRY          ,IERROR          ,INCLK       ,INCON    ,INDENDDT
     &,INTI            ,IOCONSRC    ,IODENS
     &,IODIM           ,IOFLLI          ,IOFLSAT     ,IOFMLF   ,IOFMLT
     &                 ,IOINV           ,IOITERFLEND ,ITPTVAR,IOTRLI     
     &,IOTRS
     &,IOVRWC          ,IOWAR           
     &,IREDTIMH        ,ISOLEQ          ,ISOLFL       ,ISOT         
     &,ISPARSE         ,ITERCHNGFL      ,ITERCHNGGL   ,ITERCHNGTR
     &,ITERCONVFL      ,ITERCONVTR      ,ITERFL       ,ITERGL
     &,ITERTOTFL       ,ITERTR          
     &,ITYPAFLU        ,ITYPAFLUDSC     ,ITYPCFLU     ,ITYPDFLU 
     &,LMXNDL          ,MAINF
     &,MAXNB           ,MAXNBF          ,NBAND        
     &,NBAND1                           ,NCONVI       ,NCONVIFL
     &                 ,NFLAGS          ,NFNL         ,NINT               
     &,NPAREL          ,NPARALG         ,NPARNP       ,NPBTP      
     &,NPPEL           ,NPPNP           ,NTRNRF       ,NTYPAR
     &,NUMEL           ,NUMITER         ,NUMNP        ,NWRITE
     &,NZTRA           ,NZPAR           ,NZPRG         
     &!REAL SCALARS
     &,BETAC          ,CREF             ,DENSREF      ,DHITMX
     &,DTIMEF         ,DTIMET        
     &,TABSOLUT       ,THETAF           ,THETAT
     &,TICAL          ,TICALAN
     &,TINC           ,TINCINI          ,TINCLAST     ,TINTERVOBS
     &,TOLD           ,VISCREF          ,VAR_REF      ,ZEROF
     &,FILENAME
     &!INTEGER ARRAYS
     &,IBCOD                                          ,IBTCO    ,IFLAGS
     &,INORPAR        ,IPAR_DIR         ,ISOZ         ,IXPARNP 
     &,KXX            ,LDIM             ,LINMET       ,LNNDEL
     &,LTYPE          ,LXPAREL          ,NFNLPAR      ,NFNLTIP
     &,NFNLPRG        ,NFTPAR           ,NROW         ,NZONE_PAR
     &,IAD_S          ,IADD_S           ,IADN_S       ,IOWRITE
     &,IAFD_S         ,IAFDD_S          ,IAFDN_S
     &!REAL ARRAYS
     &,AFLU             ,AFLUDSC      ,AFLUDSCF
     &,ALFA
     &,AREA           ,BFLU             ,BIBI
     &,CAUX1          ,CAUX2            ,COORD
     &,CCALAN         ,CCALIT           ,CFLU         ,CFPAREL
     &,CFPARNP        ,CPREV1           ,CPREV2       ,DENSITY
     &,DER_VISC       ,DNODALRH         ,DFLU
     &,DFLUDFLU       ,DFLUDTRA         ,DBFLUDFLU    ,DBFLUDTRA
     &,DPARELDH       ,DPARELDC
     &,FNT            ,GRAVEL           ,HAUX1
     &,HAUX2          ,HBASE            ,HCALAN       ,HCALIT
     &,HPREV1           ,HPREV2
     &,PARC
     &,PARACD         ,PAR_DIR          ,PAREL        ,PARNP
     &,PRGC
     &,SOLUTION       ,TIME             ,VISCOSITY    ,WATVOL       
     &!NUEVOS
     &,DELTAITER      ,DELTAHGL         ,IDELHGL      ,RESHMAX
     &,RESHMAXOLD     ,WORK
     &,BUOYANCY       ,DBUOYANCY        ,IPARTNER
     &,GRADLOC,GP_COORD,POINTWEIGHT,MAXPG,NTYPEL,grdff,DWDH,CONCFLOW
     &,ATRA,DTRA,IDIMATRA,IDIMDTRA,ITYPATRA,ITYPDTRA,CAUDAL,IORTS)


c     ******************************************************************
c     THE CODE PRINTED BELOW CORRESPONDS TO THE SUBROUTINE FLOW_EQUATION
c
c     NUEVAS VARIABLES
c     ----------------
c
c     IOTRS: indicates wether flow is stationary or not in the 
c                     current time step
C
C     QUEDA A HACER: ITERATION COUNTERS
C                    UPDATE_FILE_SCRATCH
C

c     ******************************************************************

C-------------------- step 0: initializing
C--------------------  the commands in the following block only need to be carried out 
C--------------------  once during direct problem solution
      IMPLICIT NONE

C-------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4 IAFLUDSC_COLS,IAFLUDSC_ROWS,IDIMDENS,IDIMAFLU
     &,IDIMBB,IDIMCFLU,IDIMDFLU,IDIMGRAVEL
     &,IDIMFNT,IDIMWORK,IDIMATRA,IDIMDTRA,IODIRECT,INCLK
     &,IENTRY,IERROR,INCON,INDENDDT,INTI,IOCONSRC
     &,IODENS,IODIM,IOFLLI,IOITERFLEND
     &,IOFLSAT,IOFMLF,IOFMLt,IOINV,IORTS,ITPTVAR,IOTRLI,IOTRS
     &,IOVRWC,IOWAR,IREDTIMH,ISOLFL
     &,ISOT,ISPARSE,ITERCHNGFL,ITERCHNGGL,ITERCHNGTR
     &,ITERCONVFL,ITERCONVTR,ITERFL,ITERGL,ITERTOTFL,ITERTR
     &,ITYPAFLU,ITYPAFLUDSC,ITYPDFLU,ITYPCFLU
     &,ITYPATRA,ITYPDTRA
     &,LMXNDL
     &,MAINF,MAXNB,NBAND1,NCONVI,NCONVIFL
     &,NFLAGS,NFNL,NINT 
     &,NPAREL,NPARALG,NPARNP,NPBTP,NPPEL,NPPNP,NTRNRF,NTYPAR,NUMEL
     &,NUMITER,NUMNP,NWRITE,NZTRA,NZPAR,NZPRG,MAXNBF,NBAND
     &,IDELHGL,MAXPG,NTYPEL



      REAL*8 BETAC,CREF,DENSREF,DHITMX,DTIMEF,DTIMET
     &,TABSOLUT,THETAF,THETAT,TICAL,TICALAN,TINC,TINCINI,TINCLAST
     &,TINTERVOBS,TOLD,VISCREF,DELHMAXOLD,VAR_REF,ZEROF,DELTAHGL
        
      CHARACTER FILENAME(20)*20

C-------------------- EXTERNAL VARIABLES: ARRAYS
      
      INTEGER*4 IBCOD(NUMNP)
     &,IFLAGS(NFLAGS),INORPAR(NTYPAR),IPAR_DIR(NPARALG) ,ISOZ(NZTRA)     
     &,IXPARNP(NUMNP,NPARNP) ,KXX(LMXNDL,NUMEL),LDIM(NUMEL)
     &,LINMET(3,2),LNNDEL(NUMEL),LTYPE(NUMEL)
     &,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR),NFNLTIP(NFNL)
     &,NFNLPRG(8,NFNL) 
     &,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)
     &,IAD_S(MAXNB*ISPARSE,NUMNP*ISPARSE),IADD_S(NUMNP*ISPARSE)
     &,IADN_S(NUMNP*ISPARSE)
     &,IAFD_S(MAXNBF,NUMNP),IAFDD_S(NUMNP)
     &,IAFDN_S(NUMNP),IOWRITE(NWRITE),ISOLEQ(NINT,4),IBTCO(NUMNP)
     &,IPARTNER(6,3,6)



      REAL*8 AFLU(NUMEL,IDIMAFLU)
     &,AFLUDSC(IAFLUDSC_ROWS,IAFLUDSC_COLS)
     &,AFLUDSCF(MAXNBF,NUMNP)
     &,ALFA(NUMNP),AREA(NUMEL)
     &,ATRA(NUMEL,IDIMATRA)
     &,BFLU(NUMNP),BIBI(IDIMBB,NUMEL)
     &,CAUX1(NUMNP),CAUX2(NUMNP),CCALAN(NUMNP),CCALIT(NUMNP)
     &,CFLU(NUMEL, IDIMCFLU),CONCFLOW(NUMNP),CPREV1(NUMNP),CPREV2(NUMNP)
     &,CFPAREL(NUMEL,NPAREL),CFPARNP(NUMNP,NPARNP)
     &,DENSITY(IDIMDENS),DER_VISC(NUMEL)
     &,DNODALRH(NUMNP,4)
     &,DFLU(NUMEL,IDIMDFLU),DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &,DFLUDTRA(NUMEL,LMXNDL*LMXNDL)  
     &,DBFLUDFLU(NUMNP),DBFLUDTRA(NUMNP)    
     &,DPARELDH(NPPEL,NUMEL)
     &,DPARELDC(NPPEL,NUMEL),DWDH(MAX(1,(IOVRWC-1)*2*LMXNDL),NUMEL)
     &,DTRA(NUMEL,IDIMDTRA)
     &,FNT(IDIMFNT,NINT)
     &,GRAVEL(IDIMGRAVEL,MIN(3*IDIMGRAVEL,3))
     &,HAUX1(NUMNP),HAUX2(NUMNP),HBASE(NUMEL)
     &,HCALAN(NUMNP),HCALIT(NUMNP),HPREV1(NUMNP),HPREV2(NUMNP)
     &,PARC(NZPAR),PARACD(3,NFNL) ,PAR_DIR(NPARALG),PAREL(NUMEL,NPPEL)
     &,PARNP(NUMNP,NPPNP),PRGC(NZPRG),TIME(NINT)
     &,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3,NPBTP),SOLUTION(NUMNP)
     &,VISCOSITY(NUMEL), WORK(IDIMWORK),DELTAITER(NUMNP)
     &,COORD(NUMNP,3)
     &,BUOYANCY(IODIM,LMXNDL,NUMEL),DBUOYANCY(IODIM,LMXNDL*LMXNDL,NUMEL)
     &,GRADLOC(IODIM,LMXNDL,MAXPG),GP_COORD(6,8,IODIM)
     &,POINTWEIGHT(NTYPEL),grdff(iodim,lmxndl,numel),CAUDAL(NUMNP)


C-------------------- INTERNAL VARIABLES

      REAL*8::XPARAM(8)

      INTEGER*4::I,IDESC, IOCALCDEV, IOCALCDEVF, IOCALCDEVT
     &,IDELHMAX,INDSSTR,IONEWT, IRESHMAX,ISYMETRIC,ITRAP,IORECATRA
     &,IOCONVF,INDCHANGES,INDSSTR_T,ITERM,N,NROW,NUMDIVH,L,IOREGIMEN
     &,IREGTRA

      REAL*8 DELHMAX,RESHMAX,DRELHMX,RESHMAXOLD


      IF (IFLAGS(3).NE.0) CALL IO_SUB('FLOW',0)

C------------------------- Initilizes  some internal variables

      IDESC = 0
      IOCALCDEV = 0
      IOCALCDEVF = 0
      IOCALCDEVT = 0
      IDELHMAX = 0
      INDSSTR = 0
      IONEWT = 0
      IRESHMAX = 0
      ISYMETRIC = 0
      ITRAP = 0
      IORECATRA = 0
      INDCHANGES = 0
      INDSSTR_T = 0
      ITERM = 0
      NROW = 0
      IREGTRA = 0
C------------------------- Initilizes convergence and timestep reduction
C------------------------- indicators.

      IOCONVF = 0
      IREDTIMH= 0
      IOITERFLEND = 0
      DELHMAX = HUGE(1D0) !Valores altos para que no se crea que
      RESHMAX = HUGE(1D0) ! el residuo crece en la primera iteración.
      NTRNRF = 0 

      IF (IODENS.EQ.1) THEN
         DELTAITER = 0D0
      END IF !IODENS.EQ.1

C------------------------- States the type of solution

      IF (INTI.EQ.0) THEN

          IF (IOTRS.EQ.0 .OR. IOTRS.EQ.2) THEN

              IOREGIMEN = 1 !Steady state or initial steady state

          ELSE IF (IOTRS.EQ.1) THEN

              IOREGIMEN = 4 

          END IF !IOTRS.EQ.0 ...

      ELSE

          IOREGIMEN = ISOLEQ(INTI,1)

      END IF !INTI.EQ.0

      INDSSTR = MIN(1,INTI)

C------------------------- Selects the type of solution for this timestep.

      SELECT CASE (IOREGIMEN)

C------------------------- Steady or transient solution

          CASE(0,1)

          IF (IOREGIMEN.EQ.1) THETAF = 1D0 ! To solve SS apropriately

C-------------------- if we are solving a nonlinear flow problem
C-------------------- or coupled flow and transport with newtons
C-------------------- method, we must make an initial estimate
C-------------------- for the heads in time k+1

c-parche-Entrar en esta subrutina da problemas.
          IF (IOFLLI.EQ.33333) THEN

              CALL STATE_VARIABLE_INIT
     &    (IDIMFNT  ,0        ,INTI        ,IOINV     ,IPAR_DIR(4)
     &    ,NCONVI   ,NINT     ,NPARALG     ,NPARNP    ,NTYPAR   ,NUMITER
     &    ,NUMNP    ,NZPAR    ,PAR_DIR(29) ,TICAL     ,TICALAN  ,TIME
     &    ,TINC     ,TINCINI  ,TINCLAST    ,TINTERVOBS,BFLU     ,CFPARNP
     &    ,FNT      ,IBCOD    ,INORPAR     ,IXPARNP   ,NFTPAR   ,PARC
     &    ,HAUX1    ,HCALIT   ,HCALAN      ,HCALIT  ,HPREV1   ,HPREV2)

              CALL UPDATE_VECTORS_AUX
     &            (NUMNP     ,TINC    ,THETAF    ,HCALAN   ,HCALIT
     &            ,HAUX1  ,HAUX2)

          END IF !IOFLLI.NE.0

C-------------------- Initializes transport variable if density is not
C-------------------- constant. Only in first global iteration
C-------------------- (in following ones transport has already been solved).

          IF (IODENS.EQ.1 .AND. ITERGL.EQ.0 .AND. INTI.GT.0) THEN
          
              IREGTRA = ISOLEQ(INTI,2)
c-parche- ULTRA PARCHE PORQUE SI NO, PONE CCALIT = INFINITO          
              IF (IREGTRA.EQ.58) THEN

                  CALL STATE_VARIABLE_INIT
     &(IDIMFNT  ,0        ,INTI        ,IOINV     ,IPAR_DIR(4)
     &,NCONVI   ,NINT     ,NPARALG     ,NPARNP    ,NTYPAR   ,NUMITER
     &,NUMNP    ,NZPAR    ,PAR_DIR(29) ,TICAL     ,TICALAN  ,TIME
     &,TINC     ,TINCINI  ,TINCLAST    ,TINTERVOBS,BFLU     ,CFPARNP
     &,FNT      ,IBCOD    ,INORPAR     ,IXPARNP   ,NFTPAR   ,PARC
     &,CAUX1    ,CCALIT   ,CCALAN      ,CCALIT  ,CPREV1   ,CPREV2)


C-------------------- If we are to read the concentration field from a file

              ELSE IF (IREGTRA.EQ.2) THEN

                  CALL ENDATINICOND_AUX
     &(NUMNP      ,1          ,15         ,MAINF      ,IOWAR
     &,0          ,IERROR     ,NROW       ,FILENAME  ,CCALIT)   

C-------------------- if all concentrations are to be put to zero

              ELSE IF (IREGTRA.EQ.3) THEN

                  CCALIT = 0D0

              END IF !IREGTRA.EQ.0,2,3

C-------------------- Allow to perturbate the concentration in a node
C-------------------- for validation of the derivatives

              IF (IFLAGS(27).GT.0) THEN

                  PRINT*,""
                  PRINT*, "NUDO A PERTURBAR (CONC): "
                  READ*,N
                  IF (N.NE.0) THEN
                      CCALIT(ABS(N)) = CCALIT(ABS(N)) + (N/ABS(N))*1D-5
                  END IF

              END IF !IFLAGS(27).GT.0


C-------------------- it might be necesary to calculate the auxiliary
C-------------------- vectors of transport

             IF (IREGTRA.EQ.0 .OR. IREGTRA.EQ.2 .OR. IREGTRA.EQ.3) THEN
                 CALL UPDATE_VECTORS_AUX
     &(NUMNP   ,TINC    ,PAR_DIR(30)  ,CCALAN   ,CCALIT  ,CAUX1  ,CAUX2)

             END IF !IREGTRA.EQ.0,2,3

          END IF !ISOLEQ(MAX(INTI,1),2).EQ.0 .AND. INTI.GT.0

C-------------------- Allow to perturbate the head in a node for
C-------------------- validation of the derivatives

          IF (IFLAGS(27).GT.0) THEN

              PRINT*,""
              PRINT*, "NUDO A PERTURBAR (NIVELES): "
              READ*,N
              IF (N.NE.0) THEN
                  HCALIT(ABS(N)) = HCALIT(ABS(N)) + (N/ABS(N))*1D-5
                  CALL UPDATE_VECTORS_AUX
     &                (NUMNP     ,TINC    ,THETAF    ,HCALAN   ,HCALIT
     &                ,HAUX1  ,HAUX2)
              END IF !N.NE.0

          END IF !IFLAGS(27).GT.0
             
C-------------------- initializing iteration counters

          ITERFL = 0
          NUMDIVH = 0


C-------------------- If there is variable density, transport parameters
C-------------------- have to be computed now.

          IF (IODENS.EQ.1) THEN

C-------------------- It is necesary to know transport regimen

              IF (INTI.EQ.0) THEN

                  IF (IORTS.EQ.0 .OR. IORTS.EQ.2) THEN

                      INDSSTR_T = 0

                  END IF !IORTS.EQ.0 ...

              ELSE

                  IF (ISOLEQ(INTI,2).EQ.0) THEN

                      INDSSTR_T = 1

                  ELSE

                      INDSSTR_T = 0

              END IF !ISOLEQ(INTI,2).EQ.0

              END IF !INTI.EQ.0

              CALL COMP_PARAM_TRA
     &            (CCALAN   ,CCALIT   ,CFPAREL  ,CFPARNP  ,DPARELDC
     &            ,DPARELDH ,DTIMET   ,PAR_DIR(32)        ,FNT
     &            ,IBTCO    ,IDIMFNT  ,IFLAGS   ,INCLK    ,INCON
     &            ,INDSSTR_T,INORPAR  ,INTI     ,IODENS   ,IOFMLT
     &            ,IOTRLI   ,IPAR_DIR ,IXPARNP  ,KXX      ,LMXNDL
     &            ,LNNDEL   ,LXPAREL  ,MAINF    ,NFLAGS   ,NFNL
     &            ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT
     &            ,NPARALG  ,NPAREL   ,NPARNP   ,NPPEL    ,NPPNP
     &            ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &            ,NZPRG    ,PARACD   ,PARC     ,PAREL    ,PARNP
     &            ,PRGC     ,XPARAM   ,THETAT   ,TINC     ,TINTERVOBS)

C------------------------- If solving a density dependent problem
C------------------------- fluid density and viscosity are computed.

              CALL CALC_DENS
     &            (BETAC    ,CAUX1    ,CREF    ,DENSREF
     &            ,DENSITY  ,IDIMDENS ,1       ,KXX
     &            ,LMXNDL   ,LNNDEL   ,NUMEL    ,NUMNP)

              CALL CALC_VISC
     &            (CAUX1    ,DER_VISC ,ITPTVAR  ,KXX     ,LMXNDL
     &            ,LNNDEL   ,NUMEL    ,NUMNP    ,VAR_REF
     &            ,VISCOSITY)


C------------------------- Checks if derivatives of buoyancy have to
C------------------------- be computed.

              IF (LINMET(2,2).EQ.2 .OR. LINMET(3,2).GT.1) THEN
                  IOCALCDEVT = 1
              END IF

              DO L=1,NUMEL

                  CALL COMP_BUOYANCY_CF
     &              (BETAC      ,BUOYANCY   ,CREF       ,CAUX1
     &              ,COORD      ,DBUOYANCY  ,LTYPE(L)   ,GRADLOC,GRAVEL
     &              ,LDIM(L)    ,IPARTNER   ,IOCALCDEVT ,KXX
     &              ,IODIM      ,L          ,LMXNDL     ,MAXPG
     &              ,LNNDEL(L)  ,NUMEL      ,NUMNP)

              END DO !L=1,NUMEL

          ELSE

              DENSITY = 1D0
              INDSSTR_T = 0

          END IF !IODENS.EQ.1

C-------------------- Loop until convergence,divergence or maximumiteration nr

          DO WHILE (IOCONVF.EQ.0 .AND. IREDTIMH.EQ.0 
     &                 .AND. ITERFL.LT.ITERCONVFL)

              CALL INCREMENT_ITER
     &            (IOITERFLEND,ITERTOTFL,ITERFL,ITERCONVFL) 

  
              CALL CHOOSE_LINEARIZT_METHOD
     &            (IODENS    ,ITERCHNGFL ,ITERCHNGGL ,ITERCHNGTR
     &            ,ITERCONVFL,ITERCONVTR ,ITERFL     ,ITERGL     ,ITERTR
     &            ,LINMET)  
         
C-------------------- to prevent getting stuck in Flow if linmet(1,2)=0

              IF (LINMET(1,2).EQ. 0) IOCONVF = 1

C-------------------- to define the kind of LHS desired

              IONEWT = 0

              IF (LINMET(1,2).EQ.0 .OR. LINMET(1,2).EQ.2
     &            .OR. LINMET(3,2).EQ.2) IONEWT = 1

              IOCALCDEVF = 0
              IOCALCDEVT = 0
      
C-------------------- if we are using newtons method, calculate DFLUDFLU

              IF (LINMET(1,2).GT.1 .OR. LINMET(3,2).EQ.2) IOCALCDEVF = 1
         
C-------------------- if we are using newton's method for transport and we have
C-------------------- variable density, then we need dfludtra to calculate dqdtra

              IF (IODENS.EQ.1 .AND. LINMET(2,2).EQ.2) IOCALCDEVT = 1
    

C-------------------- if we are using newtons method to solve flow & tpt
C-------------------- together, then calculate DFLUDTRA     

              IF (LINMET(3,2).GT.1) IOCALCDEVT = 1
      
      
C-------------------- if derivatives of some sort need to be calculated

              IOCALCDEV=MAX(IOCALCDEVF, IOCALCDEVT)

         
C-------------------- Step 5: Calculate flow parameters

              CALL  COMP_PARAM_FLOW 
     &            (CAUX1    ,COORD     ,DTIMEF    ,IDIMFNT
     &            ,INDSSTR  ,INTI      ,ITPTVAR   ,IODENS    ,IODIM
     &            ,IOFLLI   ,IOFLSAT   ,IOFMLF    ,ISOT      ,LMXNDL
     &            ,MAINF                          ,NFLAGS    ,NFNL
     &            ,NINT     ,NPARALG   ,NPAREL    ,NPARNP    ,NPPEL
     &            ,NPPNP    ,NTYPAR    ,NUMEL     ,NUMNP     ,NZPAR
     &            ,NZTRA    ,VISCREF              ,CFPAREL   ,CFPARNP
     &            ,DNODALRH ,DPARELDC  ,DPARELDH  ,FNT       ,GRAVEL
     &            ,HBASE    ,HCALAN    ,HCALIT    ,IBCOD
     &                      ,IFLAGS    ,INORPAR   ,IPAR_DIR  ,ISOZ
     &            ,IXPARNP  ,KXX       ,LNNDEL    ,LXPAREL   ,NFNLPAR
     &            ,NFNLTIP  ,NFNLPRG   ,NFTPAR    ,NZONE_PAR ,PARACD
     &            ,PARC     ,PAR_DIR   ,PAREL     ,PARNP     ,TINC
     &            ,TINTERVOBS,VAR_REF)

C------------------------- Concentration sources in flow equation.
C------------------------- (only at first iteration).

              IF (IOCONSRC.EQ.1 .AND. IODENS.EQ.1 .AND. ITERFL.EQ.1)THEN

C------------------------ States if ATRA contains recharge
C------------------------ (only when Picard method is used for transport)
C------------------------ (computed again, it may have changed inside TRANSPORT)

                  IF (LINMET(2,2).EQ.1 .OR. LINMET(3,2).EQ.1 )THEN
                      IORECATRA = 1
                  ELSE
                      IORECATRA = 0
                  END IF !LINMET(2,2).EQ.1 ...

                  CALL COMFLOW_CONC
     &                (AFLU     ,AREA     ,PAREL(1,8)         ,ATRA
     &                 ,BETAC    ,BUOYANCY ,CAUDAL   ,CAUX1    ,CAUX2
     &                ,CFLU     ,PAREL(1,15)        ,CONCFLOW ,COORD
     &                ,CREF     ,DENSITY  ,DENSREF  ,DFLU     ,DTRA
     &                ,GP_COORD ,GRADLOC  ,HAUX1    ,HAUX2    ,IAD_S
     &                ,IADN_S   ,IBCOD    ,IBTCO    ,IDIMAFLU ,IDIMATRA
     &                ,IDIMCFLU ,IDIMDFLU ,IDIMDTRA ,IODENS   ,IODIM
     &                ,IORECATRA,IOREGIMEN,IREGTRA  ,ISOZ     ,ITYPAFLU
     &                ,ITYPATRA ,ITYPCFLU ,ITYPDFLU ,ITYPDTRA ,KXX
     &                ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE ,LXPAREL(1,4)
     &                ,LXPAREL  ,MAXNB    ,MAXPG    ,NPAREL   ,NPPEL
     &                ,NPPNP    ,NTYPEL   ,NUMEL    ,NUMNP ,NZONE_PAR(3)
     &                ,NZTRA    ,PAREL    ,PARNP    ,POINTWEIGHT
     &                ,THETAT)

              END IF !IOCONSRC.EQ.1

C-------------------- If we are dealing with a linear flow problem
C-------------------- then the matrixes A and D should be calculated
C-------------------- only once during direct problem calculation.

              IF  (IOFLLI.NE.0 .OR. IODENS.EQ.1 .OR. IOREGIMEN.EQ.1
     &            .OR. (INTI.LE.1 .AND. IENTRY.LE.1))  THEN   

                  CALL COMP_AFLU
     &       (IDIMAFLU   ,IDIMBB     ,IOCALCDEVF ,IOCALCDEVT ,LMXNDL
     &       ,NPAREL     ,NPPEL      ,NUMEL      ,NUMNP      ,NZTRA
     &       ,AFLU       ,BETAC      ,BIBI       ,DENSITY    ,DFLUDTRA
     &       ,DFLUDFLU   ,DPARELDH   ,DPARELDC  ,PAR_DIR(31),PAR_DIR(32)
     &       ,HAUX1      ,ISOZ       ,KXX        ,LDIM       ,LNNDEL
     &       ,LXPAREL    ,PAREL      ,PAR_DIR(30))     

C-------------------- Storage matrixes C & D should only be
C-------------------- calculated for transient flow problems.

                  IF (IOREGIMEN.EQ.0) THEN

                      CALL COMP_DFLU
     &               (AREA     ,BETAC    ,CAUX1    ,CREF     ,DENSITY
     &               ,DENSREF  ,DFLU     ,DFLUDFLU ,DFLUDTRA ,DPARELDC
     &               ,DPARELDH ,PAR_DIR(31)        ,PAR_DIR(32)
     &               ,HAUX2    ,IDIMDFLU ,IFLAGS   ,IOCALCDEVF
     &               ,IOCALCDEVT         ,IODENS   ,IOVRWC   ,KXX
     &               ,LMXNDL   ,LNNDEL   ,LTYPE    ,MAINF    ,NFLAGS
     &               ,NPPEL    ,NUMEL    ,NUMNP    ,PAREL,PAR_DIR(30))

C-------------------- Storage matrix C should only be calculated
C-------------------- for density dependent flow.

                      IF (IODENS.EQ.1) THEN

                          CALL COMP_CFLU
     &                  (AREA     ,BETAC    ,CAUX1    ,CAUX2    ,CFLU
     &                  ,CREF     ,DENSITY  ,DENSREF  ,DFLUDFLU
     &                  ,DFLUDTRA ,DWDH     ,IDIMCFLU ,IOCALCDEV,IOVRWC
     &                  ,KXX      ,LMXNDL   ,LNNDEL   ,NUMEL    ,NUMNP
     &                  ,PAR_DIR(30)        ,WATVOL(1,1,1,1))

                      END IF !IODENS.EQ.1

                  END IF !IOTRS.GE.1

              END IF   !(IOFLLI.NE.0 .OR. (INTI.LE.1 .AND. IENTRY.LE.1))

              CALL COMP_BFLU
     &            (AREA     ,BETAC    ,BFLU     ,BUOYANCY ,CAUX1
     &            ,CONCFLOW ,COORD    ,CREF     ,DBFLUDFLU,DBFLUDTRA
     &            ,DBUOYANCY,DENSITY  ,DENSREF  ,DFLUDFLU ,DFLUDTRA
     &            ,DNODALRH ,DPARELDH ,PAR_DIR(31)        ,GP_COORD
     &            ,GRADLOC  ,HAUX1    ,HCALAN   ,IBCOD    ,IBTCO
     &            ,INDSSTR  ,IOCALCDEVF         ,IOCALCDEVT
     &            ,IOCONSRC ,IODENS   ,IODIM    ,IOFLLI   ,IONEWT
     &            ,ISOZ     ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL
     &            ,LTYPE    ,LXPAREL  ,MAXPG    ,NPAREL   ,NPPEL
     &            ,NPPNP    ,NTYPAR   ,NTYPEL   ,NUMEL    ,NUMNP
     &            ,NZONE_PAR,NZTRA    ,PAREL    ,PARNP    ,POINTWEIGHT
     &            ,THETAF   ,THETAT   ,gravel   ,grdff)

C-------------------- Assembling the flow right hand side


               IF (IOREGIMEN.EQ.0) THEN

                  CALL ASSEMB_RHS
     &                (AFLU     ,BFLU     ,CAUX2    ,CFLU     ,DFLU
     &                ,IAD_S    ,IADN_S   ,IDIMAFLU ,IDIMCFLU
     &                ,IDIMDFLU ,0        ,INDSSTR  ,IODENS  ,IONEWT
     &                ,ITYPAFLU ,ITYPDFLU ,ITYPCFLU ,KXX      ,LNNDEL
     &                ,LMXNDL   ,NUMEL    ,NUMNP    ,THETAF    ,TINC
     &                ,HAUX1    ,HAUX2    ,HCALAN)

              END IF !IOREGIMEN.EQ.0


              IF (LINMET(1,2).NE.0) THEN

C-------------------- To check out whether it is necessary 
C-------------------- to assemble the LHS term or whether 
C-------------------- the old will do.

                  CALL COMP_IND_CHANGE
     &(INORPAR(11) ,IENTRY      ,INDCHANGES  ,INTI       ,NFTPAR
     &,NZONE_PAR(6),NZPAR       ,TINC        ,TOLD       ,IOREGIMEN
     &,IODENS      ,IOFLLI)


C-------------------- Assembling the left hand side.

                  IF (INDCHANGES .NE. 0) THEN

                      CALL ASSEMB_LHS
     &(IDIMAFLU   ,NUMEL        ,LMXNDL*LMXNDL,NUMEL       ,1
     &,NUMNP      ,IAFLUDSC_COLS,IAFLUDSC_ROWS,IDIMDFLU  ,NUMEL
     &,1-IOREGIMEN,IONEWT
     &,ITYPAFLU   ,ITYPAFLUDSC  ,ITYPDFLU   ,4           ,1
     &,LMXNDL     ,MAXNB        ,NBAND      ,NUMEL
     &,NUMNP      ,THETAF       ,TINC       ,AFLU        ,AFLUDSC
     &,DFLU       ,DFLUDFLU     ,DBFLUDFLU  ,IAD_S      
     &,IADD_S     ,IADN_S       ,KXX        ,LNNDEL      ,1D0)

                  END IF


C-------------------- Establishing symmetry

                  IF (ITYPAFLUDSC.EQ.7) THEN
                      ISYMETRIC = 1
                  ELSE
                      ISYMETRIC = 0
                  END IF !ITYPAFLUDSC.EQ.7

C-------------------- Establishing wether the flow matrix must be decomposed 


C-------------------- the flow matrix must be decomposed if the
C-------------------- lhs has been calculated again.
C-------------------- the flow matrix must be decomposed if we
C-------------------- solve a lineal steady state pt

                  IDESC = 1

                  IF (INDCHANGES.NE.0 .OR. IOREGIMEN.EQ.1) THEN

                      IDESC = 2

                  END IF 

                  IF (LINMET(1,2).NE.0) THEN

C-------------------- Imposing Dirichlet boundary conditions (prescribed head).

                      IF (NZONE_PAR(4).NE.0) THEN

                          CALL PRESC_BC
     &            (AFLU          ,AFLUDSC       ,BFLU          ,IDIMAFLU
     &            ,IAFLUDSC_COLS ,IAFLUDSC_ROWS ,IADD_S        ,IADN_S
     &            ,IBCOD         ,INDCHANGES    ,0             ,INDSSTR
     &            ,IONEWT        ,ISPARSE       ,ITYPAFLUDSC   ,KXX
     &            ,LMXNDL        ,LNNDEL        ,NBAND1        ,NPPNP
     &            ,NUMEL         ,NUMNP         ,PARNP         ,THETAF
     &            ,HCALIT)

                      END IF !NZONE_PAR(4).NE.0 .AND. INDCHANGES.NE.0

C-------------------- Imposing Cauchy & Neumann boundary conditions.

                      IF (NZONE_PAR(6).NE.0 .AND. INDCHANGES.NE.0) THEN

                          CALL PRESC_LEAK_BC
     &                   (AFLUDSC    ,ALFA       ,BETAC      ,CAUX1
     &                   ,CREF       ,DENSREF    ,IAD_S      ,IADD_S
     &                   ,IADN_S     ,IBCOD      ,IAFLUDSC_COLS
     &                   ,IAFLUDSC_ROWS          ,3          ,INDCHANGES
     &                   ,0          ,1          ,IOCONSRC   ,IODENS
     &                   ,IONEWT     ,0          ,ITYPAFLUDSC,KXX
     &                   ,LMXNDL     ,LNNDEL     ,MAXNB      ,NBAND
     &                   ,NPPNP      ,NUMEL      ,NUMNP      ,PARNP
     &                   ,THETAF     ,HAUX1      ,1D0)

                      END IF  !if nzone_par(6).ne.0

C-------------------- Solving the system      

                      IF (IODIRECT.EQ.0) SOLUTION=HCALIT

C-------------------- Update the maximum residuo value.


                      CALL COMPUTE_RESID_MAX
     &            (iad_s    ,iadn_s     ,0        ,0        ,IRESHMAX 
     &            ,NUMNP    ,NUMNP
     &            ,RESHMAX  ,RESHMAXOLD         ,IBCOD    ,BFLU
     &            ,AFLUDSC  ,KXX      ,LNNDEL   ,NUMEL    ,IAFLUDSC_COLS
     &            ,IAFLUDSC_ROWS ,ITYPAFLUDSC   ,HCALIT   ,SOLUTION  
     &            ,LMXNDL  ,IONEWT)


C-------------------- Writes matrix for debugging purposes.

                      IF (IFLAGS(26).GT.0) THEN

                          WRITE(611,*) 'dFdF en INTI = ', INTI
                          CALL WRITE_SQR_MATRIX
     &                        (AFLUDSC  ,IAFLUDSC_COLS
     &                        ,IAFLUDSC_ROWS      ,IAD_S
     &                        ,IADN_S   ,ISPARSE  ,ISYMETRIC,MAXNB
     &                        ,NUMNP    ,NBAND    ,611)

                          WRITE(611,*) ''

                          WRITE(331,*) 'BFLU en INTI = ', INTI
                          DO I=1,NUMNP

                              WRITE(331,*) BFLU(I)

                          END DO !I=1,NUMNP
                          WRITE(331,*) ''

                      END IF !IFLAGS(26).GT.0


                      CALL SOLVE   
     &(IAFLUDSC_ROWS      ,IAFLUDSC_COLS      ,NUMNP    ,IDESC
     &,IPAR_DIR(21)       ,IDIMWORK ,IODIRECT ,1        ,INTI
     &,IPAR_DIR(18)       ,IPAR_DIR(22)       ,ISYMETRIC         
     &,ITERM ,IPAR_DIR(15),MAINF    ,NBAND1   ,MAXNBF
     &,IPAR_DIR(20)
     &,IPAR_DIR(17)       ,IPAR_DIR(16)       ,PAR_DIR(36)   
     &,PAR_DIR(37)        ,PAR_DIR(38)        ,AFLUDSC  ,AFLUDSCF
     &,BFLU     ,IAD_S    ,IADD_S   ,IADN_S   ,IAFD_S   ,IAFDD_S
     &,IAFDN_S  ,WORK     ,SOLUTION)


C-------------------- Updating state variables. 
C-------------------- in the lineal case hcal contains hcalit
C-------------------- in the newton case the h-increment.

                      CALL UPDATE_STATE_VARIABLE
     &(DELHMAX  ,DELTAHGL  ,DELHMAXOLD,DRELHMX  ,DHITMX   
     &,0        ,IDELHGL   ,IDELHMAX  ,IOFLLI    ,0
     &,IODENS   ,IONEWT    ,NUMNP     ,ZEROF    ,DELTAITER
     &,SOLUTION ,HCALAN    ,HCALIT)                          


C------------------------- auxiliar vectors updated.

                      CALL  UPDATE_VECTORS_AUX
     &(NUMNP     ,TINC    ,THETAF    ,HCALAN   ,HCALIT   ,HAUX1  ,HAUX2)

                      IF (IOFLLI.NE.0) THEN

                          CALL WRITE_NONLIN_INFO
     &(DELHMAX  ,'   FLOW   ',IDELHMAX ,INDENDDT ,IREDTIMH ,IOWRITE(13)
     &,IRESHMAX ,ITERFL     ,MAINF    ,NCONVI   ,RESHMAX  ,TABSOLUT   
     &,TINC)

                      ENDIF  ! IOFLLI.NE.0

                      IF(IODENS.EQ.1 .OR. IOVRWC.GT.0) THEN

                          CALL CALC_WATVOL
     &                  (HCALAN   ,HCALIT   ,INTI     ,IOVRWC   ,KXX
     &                  ,LMXNDL   ,LNNDEL   ,NPBTP    ,NPPEL    ,NUMEL
     &                  ,NUMNP    ,PAREL    ,THETAT   ,WATVOL)

                          IF (IOCALCDEVF.EQ.1) THEN

                              CALL DER_WATVOL_H
     &                         (DPARELDH  ,DWDH      ,PAR_DIR(32),HCALAN
     &                         ,HCALIT    ,IOVRWC    ,KXX        ,LMXNDL
     &                         ,LNNDEL    ,NPPEL     ,NUMEL      ,NUMNP
     &                         ,PAREL)

                          END IF !IOCALCDEV.EQ.1

                      END IF !IODENS.EQ.1 .OR. IOVRWC.GT.0

C-------------------- Writes information about the process
C-------------------- in the main output file

                      CALL CHECK_CONVERGENCE
     &                   (PAR_DIR(5) ,DELHMAX    ,PAR_DIR(4) ,DRELHMX
     &                   ,IOCONVF    ,INDENDDT   ,IOWRITE(13),NCONVI
     &                   ,PAR_DIR(6) ,RESHMAX    ,IOFLLI     ,IONEWT)

                      CALL CHECK_DIVERGENCE
     &(DELHMAX,DELHMAXOLD,IONEWT,IREDTIMH,IPAR_DIR(26),NUMDIVH,RESHMAX
     &,RESHMAXOLD)
       
C-------------------- Checks if timestep reduction is needed.
c-parche
      nconvifl=nconvi
      if(ioflli.ne.0 .and. iterfl.ge.iterconvfl) then
          iredtimh=1
      end if
c-fin-parche

                  END IF !LINMET(1,2).NE.0
              END IF     !LINMET(1,2).NE.0

          END DO !Do While IOCONVF.EQ.0 .AND. IREDTIMH.EQ.0 .AND. ITERFL.LT.ITERCONVFL

C------------------------- The value of solved problem
C------------------------- indicator variable is set.
 
          IF (IOCONVF.EQ.1) THEN

              ISOLFL = INDSSTR + 1

          ELSE

              ISOLFL = 0

          END IF !IOCONVT.EQ.1

C------------------------- Read initial conditons.

          CASE(2)


              CALL ENDATINICOND_AUX
     &          (NUMNP      ,1          ,15         ,MAINF      ,IOWAR
     &          ,0          ,IERROR     ,NROW       ,FILENAME   ,HCALIT)
        
              IOCONVF = 1
              ISOLFL = 0

C-------------------------Set initial conditions to zero.

          CASE(3)

              HCALIT = 0D0
              IOCONVF = 1
              ISOLFL = 0

C------------------------- Not solve.

          CASE(4)

              IOCONVF = 1
              ISOLFL = 0
      
      END SELECT !IOREGIMEN

      IF (INTI.EQ.0) INDENDDT=1

      IF (IFLAGS(3).NE.0) CALL IO_SUB('FLOW',1)

      END SUBROUTINE FLOW_EQN
