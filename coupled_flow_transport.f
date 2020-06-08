      SUBROUTINE COUPLED_FLOW_TRANSPORT
     &          (A_COUPL_DSC   ,A_COUPL_DSCF  ,AFLU          ,AREA
     &          ,ATRA          ,BCOUPLED      ,BETAC         ,BFLU
     &          ,BTRA          ,CAUDAL        ,CAUX1         ,CAUX2
     &          ,CCALAN        ,CCALIT        ,CFLU          ,CREF
     &          ,DBFLUDFLU     ,DBFLUDTRA     ,DELTAC_SJ     ,DELTAC_TR
     &          ,DELTACOLD     ,DELTAH_FL     ,DELTAH_SJ     ,DELTAHOLD
     &          ,DELTAITER     ,DENSREF       ,DFLU          ,DFLUDFLU
     &          ,DFLUDTRA      ,DQDFLU        ,DQDTRA        ,DRELCMX
     &          ,DRELHMX       ,DTRA          ,DTRADFLU      ,DTRADTRA
     &          ,HAUX1         ,HAUX2         ,HCALAN        ,HCALIT
     &          ,IA_COUPLED_DSC_COLS          ,IA_COUPLED_DSC_ROWS
     &          ,IAD_D         ,IAD_S         ,IADD_D        ,IADN_D
     &          ,IADN_S        ,IAFD_D        ,IAFDD_D       ,IAFDN_D
     &          ,IBCOD         ,IBTCO         ,IDELCGL       ,IDELHGL
     &          ,IDIMAFLU      ,IDIMATRA      ,IDIMCFLU      ,IDIMDFLU
     &          ,IDIMDTRA      ,IDIMWORK      ,IFLAGS        ,INDENDDT
     &          ,INDSSTR       ,INTI          ,IOCONSRC      ,IOCONVF
     &          ,IOCONVGL      ,IOCONVT       ,IODENS        ,IODIRECT
     &          ,IOFLLI        ,IOITERGLEND   ,IOPTS         ,IORTS
     &          ,IOTRLI        ,IOTRS         ,IOWRITE       ,IPAR_DIR
     &          ,IREDTIMC      ,IREDTIMGL     ,IREDTIMH      ,ISOLFL
     &          ,ISOLTR        ,ISPARSE       ,ITERGL        ,ITERGLMX
     &          ,ITERM         ,ITERTOTGL     ,ITPTVAR
     &          ,ITYPACOUPLDSC ,ITYPCFLU      ,ITYPDFLU      ,KXX
     &          ,LINMET        ,LMXNDL        ,LNNDEL        ,MAINF
     &          ,MAXNB         ,MAXNBF        ,MAXNN         ,NBAND
     &          ,NBAND1        ,NCONVIFL      ,NCONVITR      ,NFLAGS
     &          ,NOPTS         ,NPARALG       ,NPBFL         ,NPBMX
     &          ,NPBTP         ,NPPEL         ,NPPNP         ,NTYPAR
     &          ,NUMDIVCGL     ,NUMDIVHGL     ,NUMEL         ,NUMNP
     &          ,NWRITE        ,NZONE_PAR     ,PAR_DIR       ,PAREL
     &          ,PARNP         ,RESCMAX       ,RESCMAXGL
     &          ,RESCMAXGLOLD  ,RESCMAXOLD    ,RESHMAX       ,RESHMAXGL
     &          ,RESHMAXGLOLD  ,RESHMAXOLD    ,SOLUTION      ,SOURCE
     &          ,TINC          ,WATVOL        ,WORK          ,WSPECHEAT)


      IMPLICIT NONE

      INTEGER*4::IA_COUPLED_DSC_COLS,IA_COUPLED_DSC_ROWS
     &          ,IDELCGL  ,IDELHGL ,IDIMAFLU
     &          ,IDIMATRA ,IDIMCFLU ,IDIMDFLU ,IDIMDTRA
     &          ,IDIMWORK ,INDENDDT ,INDSSTR  ,INTI    ,IOCONSRC
     &          ,IOCONVF  ,IOCONVGL ,IOCONVT  ,IODENS   ,IODIRECT
     &          ,IOFLLI   ,IOITERGLEND        ,IORTS    ,IOTRLI
     &          ,IOTRS    ,IREDTIMC ,IREDTIMH ,ISOLFL   ,ISOLTR
     &          ,ISPARSE  ,ITERGL   ,ITERGLMX ,ITERM    ,ITERTOTGL
     &          ,ITPTVAR  ,ITYPACOUPLDSC      ,ITYPCFLU ,ITYPDFLU
     &          ,LMXNDL   ,MAINF    ,MAXNB    ,MAXNBF   ,MAXNN
     &          ,NBAND    ,NBAND1   ,NCONVIFL ,NCONVITR ,NFLAGS
     &          ,NOPTS    ,NPARALG  ,NPBFL    ,NPBMX    ,NPBTP
     &          ,NPPEL    ,NPPNP    ,NTYPAR   ,NUMEL    ,NUMNP
     &          ,NWRITE

      INTEGER*4::IAD_D(2*MAXNB,2*NUMNP)       ,IAD_S(MAXNB,NUMNP)
     &          ,IADD_D(2*NUMNP)              ,IADN_D(2*NUMNP)
     &          ,IADN_S(NUMNP)                ,IAFD_D(2*MAXNBF,NUMNP*2)
     &          ,IAFDD_D(NUMNP*2)             ,IAFDN_D(NUMNP*2)
     &          ,IBCOD(NUMNP,NPBFL)           ,IBTCO(NUMNP,NPBTP)
     &          ,IFLAGS(NFLAGS)               ,IOPTS(NOPTS)
     &          ,IOWRITE(NWRITE)              ,IPAR_DIR(NPARALG)
     &          ,KXX(LMXNDL,NUMEL)            ,LINMET(3,2)
     &          ,LNNDEL(NUMEL)                ,NZONE_PAR(NTYPAR)

      REAL*8::BETAC       ,CREF        ,DELTAC_SJ   ,DELTAH_SJ
     &       ,DELTACOLD   ,DELTAHOLD   ,DELTAC_TR   ,DELTAH_FL
     &       ,DENSREF     ,DRELCMX     ,DRELHMX     ,RESCMAX
     &       ,RESCMAXOLD  ,RESHMAX     ,RESHMAXOLD  ,RESHMAXGL
     &       ,RESCMAXGL   ,RESHMAXGLOLD,RESCMAXGLOLD,TINC
     &       ,WSPECHEAT
     
     
      REAL*8::A_COUPL_DSC(IA_COUPLED_DSC_ROWS,IA_COUPLED_DSC_COLS)
     &       ,A_COUPL_DSCF(MAXNB,2*NUMNP),AFLU(NUMEL,IDIMAFLU,NPBFL)
     &       ,AREA(NUMEL)                ,ATRA(NUMEL,IDIMATRA,NPBTP)
     &       ,BCOUPLED(2*NUMNP)          ,BFLU(NUMNP,NPBFL)
     &       ,BTRA(NUMNP,NPBTP)          ,CAUDAL(NUMNP,NPBFL)
     &       ,CAUX1(NUMNP,NPBTP)         ,CAUX2(NUMNP,NPBTP)
     &       ,CCALAN(NUMNP,NPBTP)        ,CCALIT(NUMNP,NPBTP)
     &       ,CFLU(NUMEL,IDIMCFLU)       ,DBFLUDFLU(NUMNP,NPBFL)
     &       ,DBFLUDTRA(NUMNP)           ,DELTAITER(NUMNP)
     &       ,DFLU(NUMEL,IDIMDFLU,NPBFL)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL,NPBFL)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL) 
     &       ,DQDFLU(NUMEL,LMXNDL*LMXNDL),DQDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DTRA(NUMEL,IDIMDTRA,NPBTP) ,DTRADFLU(NUMEL,LMXNDL*LMXNDL)    
     &       ,DTRADTRA(NUMEL,LMXNDL*LMXNDL,NPBTP)
     &       ,HAUX1(NUMNP,NPBFL)         ,HAUX2(NUMNP,NPBFL)
     &       ,HCALAN(NUMNP,NPBFL)        ,HCALIT(NUMNP,NPBFL)
     &       ,PAR_DIR(NPARALG)           ,PAREL(NUMEL,NPPEL,NPBMX)
     &       ,PARNP(NUMNP,NPPNP,NPBMX)   ,SOLUTION(2*NUMNP)
     &       ,SOURCE(NUMNP)
     &       ,WATVOL(MAX(1,(IOPTS(31)-1)*LMXNDL),NUMEL,3,NPBTP)
     &       ,WORK(IDIMWORK)
 
      
      INTEGER*4 IREDTIMGL,NUMDIVCGL,NUMDIVHGL


C------------------------- Internal

      INTEGER*4:: IDELCMAX,IDELHMAX, IONEWT, IRESCMAXGL, IRESHMAXGL
      REAL*8:: DUMMY, DELTAC, DELTAH, THETAF, THETAT

      LOGICAL::IOINITSS

      REAl*8,ALLOCATABLE::DBTRADTRA(:)

      integer*4::i

      IF (IFLAGS(3).NE.0) CALL IO_SUB('COUPLED_FLOW_TRANSPORT',0)

C------------------------- States if the steady state is being solved at
C------------------------- INTI = 0.

      !IOINITSS = INTI.EQ.0 .AND. IOTRS.EQ.0 .AND. IORTS.EQ.0
      IOINITSS = INTI.EQ.0 .AND. IOTRS.NE.1 .AND. IORTS.NE.1

C------------------------- If solving the initial steady state or the
C------------------------- 'normal' time-line is running...

      IF (IOINITSS .OR.  INTI.GT.0) THEN

C------------------------- ..."global" iteration counters are incremented,

          CALL INCREMENT_ITER(IOITERGLEND ,ITERTOTGL,ITERGL ,ITERGLMX)
          PRINT*,"Iteracion GLOBAL: ", ITERGL

	    THETAF = PAR_DIR(29)
	    THETAT = PAR_DIR(30)


C------------------------- ...if we use Newton's method, one that
C------------------------- solves flow and transport toghether...
      
          IF (LINMET(3,2).EQ.2) THEN

              ALLOCATE (DBTRADTRA(NUMNP))
              DBTRADTRA = 0D0

              CALL COMP_DER_BTRA
     &            (CAUX1    ,DENSREF  ,DBTRADTRA,IBTCO(1,1) ,IOCONSRC
     &            ,IODENS   ,ITPTVAR  ,NPPNP    ,NUMNP    ,PARNP
     &            ,THETAT   ,WSPECHEAT)

C------------------------- the system is assembled and solved,

              CALL ASSEMBLE_LHS_COUPLED
     &          (A_COUPL_DSC    ,AFLU     ,ATRA     ,CFLU     ,dbfludflu
     &          ,DBFLUDTRA,DBTRADTRA
     &          ,DFLU     ,DFLUDFLU ,DFLUDTRA
     &          ,DTRA    ,DTRADFLU 
     &          ,DTRADTRA       ,ITYPACOUPLDSC
     &          ,ITYPCFLU       ,ITYPDFLU ,IA_COUPLED_DSC_COLS 
     &          ,IA_COUPLED_DSC_ROWS
     &          ,IAD_D          ,IADD_D   ,IADN_D   ,IDIMAFLU ,IDIMATRA
     &          ,IDIMCFLU       ,IDIMDFLU ,IDIMDTRA ,INDSSTR
     &          ,ISPARSE        ,KXX      ,LMXNDL   ,LNNDEL   ,2*MAXNB    
     &          ,2*NUMNP        ,NBAND    ,NUMEL    ,NUMNP    ,TINC     
     &          ,THETAF         ,THETAT   ,1D0      ,1D0)


              IF (ALLOCATED(DBTRADTRA)) THEN

                  DEALLOCATE(DBTRADTRA)

              END IF !ALLOCATED(DBTRADTRA)

              CALL ASSEMBLE_RHS_COUPLED(BCOUPLED,BFLU,BTRA,NUMNP)

C------------------------- Flow boundary conditions and contribution of
C------------------------- nodal fluxes to transport.

              CALL COMAT_BC
     &       (A_COUPL_DSC   ,BETAC         ,CAUDAL        ,CAUX1
     &       ,CREF          ,DENSREF       ,DQDFLU        ,DQDTRA
     &       ,PAR_DIR(32)   ,IAD_D         ,IADD_D        ,IADN_D
     &       ,IBTCO         ,IA_COUPLED_DSC_COLS
     &       ,IA_COUPLED_DSC_ROWS          ,1             ,1
     &       ,0             ,IODENS        ,1             ,ITYPACOUPLDSC
     &       ,KXX           ,LMXNDL        ,LNNDEL        ,MAXNB
     &       ,MAXNN         ,NBAND         ,NPPNP         ,NUMEL
     &       ,NUMNP         ,PARNP         ,THETAT
     &       ,AREA,NZONE_PAR,PAREL,NPPEL,NTYPAR)


C------------------------- Prescribed head and concentration boundary conditions.

              IF (ANY(IBCOD.EQ.1 .OR. IBTCO.EQ.1)) THEN

                  CALL PRESC_BC_COUPLED
     &             (IA_COUPLED_DSC_ROWS  ,IA_COUPLED_DSC_COLS  ,ISPARSE   
     &             ,MAXNN      ,NBAND      ,NPPNP      ,NUMNP      
     &             ,A_COUPL_DSC,BCOUPLED             ,CCALIT     ,HCALIT
     &             ,IADD_D     ,IADN_D               ,IBCOD      ,IBTCO
     &             ,PARNP(1,1,1))
               
              END IF !ANY(IBCOD.EQ.1 .OR. IBTCO.EQ.1)


              IF (IFLAGS(26).GT.0) THEN

                  CALL UNSHUFFLE_LHS (INTI,A_COUPL_DSC,NUMNP,NBAND)

                  WRITE(337,*) INTI
                  WRITE(331,*) 'BFLU en INTI = ', INTI
                  WRITE(332,*) 'BTRA en INTI = ', INTI

                 DO I=1,2*NUMNP

                      WRITE(337,*) BCOUPLED(I)

                      IF (MOD(I,2).NE.0) THEN

                          WRITE(331,*) BCOUPLED(I)

	                ELSE

                          WRITE(332,*) BCOUPLED(I)

                      END IF !MOD(I,2).NE.0

                 END DO !1,2*NUMNP


              END IF !IFLAGS(26).GT.0


              CALL SOLVE
     &            (IA_COUPLED_DSC_ROWS,IA_COUPLED_DSC_COLS,2*NUMNP
     &            ,2             ,IPAR_DIR(21)  ,IDIMWORK
     &            ,IODIRECT      ,3             ,INTI
     &            ,IPAR_DIR(18)  ,IPAR_DIR(22)  ,0
     &            ,ITERM         ,IPAR_DIR(15)  ,MAINF
     &            ,2*NBAND1      ,2*IPAR_DIR(24),IPAR_DIR(20)
     &            ,IPAR_DIR(17)  ,IPAR_DIR(16)  ,PAR_DIR(36)
     &            ,PAR_DIR(37)   ,PAR_DIR(38)   ,A_COUPL_DSC
     &            ,A_COUPL_DSCF  ,BCOUPLED      ,IAD_D
     &            ,IADD_D        ,IADN_D        ,IAFD_D
     &            ,IAFDD_D       ,IAFDN_D       ,WORK
     &            ,SOLUTION)


C------------------------- convergence and divergence for flow
C------------------------- and transport is checked,

C------------------------- Updates flow variable and computes increments.

              CALL UPDATE_STATE_VARIABLE
     &        (DELTAH_SJ  ,DUMMY     ,DELTAHOLD ,DRELHMX    ,PAR_DIR(11) 
     &        ,0          ,IDELHGL   ,IDELHMAX  ,IOFLLI     ,1
     &        ,IODENS     ,1         ,NUMNP     ,PAR_DIR(8) ,DELTAITER
     &        ,SOLUTION   ,HCALAN    ,HCALIT)

C------------------------- Updates transport variable and computes
C------------------------- increments.

              CALL UPDATE_STATE_VARIABLE
     &        (DELTAC_SJ  ,DUMMY     ,DELTACOLD ,DRELCMX    ,PAR_DIR(12)  
     &        ,1          ,IDELCGL   ,IDELCMAX  ,IOTRLI     ,1
     &        ,IODENS     ,1         ,NUMNP     ,PAR_DIR(9) ,DELTAITER
     &        ,SOLUTION   ,CCALAN    ,CCALIT)    

              CALL  UPDATE_VECTORS_AUX
     &(NUMNP     ,TINC    ,THETAF    ,HCALAN   ,HCALIT   ,HAUX1  ,HAUX2)

              CALL  UPDATE_VECTORS_AUX
     &(NUMNP     ,TINC    ,THETAT    ,CCALAN   ,CCALIT   ,CAUX1  ,CAUX2)

C------------------------- Computes maximum resid. (flow variable)

              CALL COMPUTE_RESID_MAX
     &            (IAD_S    ,IADN_S   ,0        ,1        ,IRESHMAXGL
     &            ,NUMNP    ,2*NUMNP  ,RESHMAXGL,RESHMAXGLOLD
     &            ,IBCOD    ,BCOUPLED ,A_COUPL_DSC        ,KXX
     &            ,LNNDEL   ,NUMEL    ,IA_COUPLED_DSC_COLS
     &            ,IA_COUPLED_DSC_ROWS,ITYPACOUPLDSC      ,HCALIT
     &            ,SOLUTION ,LMXNDL   ,1)
 

C------------------------- Computes maximum resid. (tpt. variable)

              CALL COMPUTE_RESID_MAX
     &            (IAD_S    ,IADN_S   ,1        ,1        ,IRESCMAXGL
     &            ,NUMNP    ,2*NUMNP  ,RESCMAXGL,RESCMAXGLOLD
     &            ,IBTCO    ,BCOUPLED ,A_COUPL_DSC        ,KXX
     &            ,LNNDEL   ,NUMEL    ,IA_COUPLED_DSC_COLS
     &            ,IA_COUPLED_DSC_ROWS,ITYPACOUPLDSC      ,HCALIT
     &            ,SOLUTION ,LMXNDL   ,1) 

          END IF ! (LINMET(3,2).EQ.2)


C------------------------- If we use picard, 
C------------------------- set the global residuals equal to the local ones

          IF (LINMET(3,2).EQ.0) THEN

              RESHMAXGL = RESHMAX
              RESCMAXGL = RESCMAX
              RESHMAXGLOLD = RESHMAXOLD
              RESCMAXGLOLD = RESCMAXOLD
              DELTAH = DELTAH_FL
              DELTAC = DELTAC_TR

          ELSE

              DELTAH = DELTAH_SJ
              DELTAC = DELTAC_SJ 

          END IF !LINMET(3,2).EQ.0


          IONEWT = MAX(LINMET(3,2)-1,0)

          CALL CHECK_CONVERGENCE
     &        (PAR_DIR(5) ,DELTAH     ,PAR_DIR(4) ,DRELHMX
     &        ,IOCONVF    ,INDENDDT   ,IOWRITE(13),NCONVIFL
     &        ,PAR_DIR(6) ,RESHMAXGL  ,1          ,IONEWT)

          CALL CHECK_CONVERGENCE
     &        (PAR_DIR(35),DELTAC     ,PAR_DIR(34),DRELCMX
     &        ,IOCONVT    ,INDENDDT   ,IOWRITE(17),NCONVITR
     &        ,PAR_DIR(7) ,RESCMAXGL  ,1          ,IONEWT)

          CALL CHECK_DIVERGENCE
     &        (DELTAH     ,DELTAHOLD  ,IONEWT     ,IREDTIMH
     &        ,IPAR_DIR(26)           ,NUMDIVHGL  ,RESHMAXGL
     &        ,RESHMAXGLOLD)

          CALL CHECK_DIVERGENCE
     &        (DELTAC     ,DELTACOLD  ,IONEWT     ,IREDTIMC
     &        ,IPAR_DIR(27)           ,NUMDIVCGL  ,RESCMAXGL
     &        ,RESCMAXGLOLD)

          DELTAHOLD = DELTAH
          DELTACOLD = DELTAC

C------------------------- Writes non-lineal info (flow)
c-parche  comentado para no liar el res.out
c         CALL WRITE_NONLIN_INFO
c     &       (DELTAH  ,'COUPLED - FLOW'     ,IDELHMAX
c     &       ,INDENDDT ,IREDTIMH ,IOWRITE(13),IRESHMAX
c     &       ,ITERGL   ,MAINF    ,NCONVI     ,RESHMAX
c     &       ,TABSOLUT ,TINC)
c
C------------------------- Writes non-lineal info (tpt.)
c
c         CALL WRITE_NONLIN_INFO
c     &       (DELTAC  ,'COUPLED - TRANSPORT',IDELCMAX
c     &       ,INDENDDT ,IREDTIMC ,IOWRITE(17),IRESCMAX
c     &       ,ITERGL   ,MAINF    ,NCONVI     ,RESCMAX
c     &       ,TABSOLUT ,TINC)
c- fin parche

          WRITE(*,10) RESHMAXGL,PAR_DIR(6)
          WRITE(*,20) RESCMAXGL,PAR_DIR(7)

          WRITE(*,30) DELTAH,PAR_DIR(5)
          WRITE(*,40) DELTAC,PAR_DIR(35)

   10 FORMAT(1X,'Res. max. Flujo:',1X,E10.3E3,4X,'Criterio:',1X,E10.3E3)
   20 FORMAT(1X,'Res. max. Trans:',1X,E10.3E3,4X,'Criterio:',1X,E10.3E3)
   30 FORMAT(1X,'Delt.max. Flujo:',1X,E10.3E3,4X,'Criterio:',1X,E10.3E3)
   40 FORMAT(1X,'Delt.max. Trans:',1X,E10.3E3,4X,'Criterio:',1X,E10.3E3)

      END IF !IOINITSS .OR.  INTI.GT.0

C------------------------- Global convergence.

      IOCONVGL = IOCONVF*IOCONVT

C------------------------- Tiem-step reduction.

      IREDTIMGL = MAX(IREDTIMH,IREDTIMC)


C------------------------- sources for simultaneous problems computed
C------------------------- if transport and flow have been solved succesfully.


      IF (IOCONVGL.EQ.1) THEN

C------------------------------ update indicators of having solved flow
C------------------------------ and transport for the inital steady
C------------------------------ state or the 'normal', time line.
C------------------------------ Otherwiese keeps the values from
C------------------------------ FLOW_EQN and TRANSPORT subroutines.

          IF (IOINITSS .OR.  INTI.GT.0) THEN

C-parche- Esto deería estar contralado por ISOLEQ y no por IOTRS y IORTS
              ISOLFL = 1
              ISOLTR = 1

              IF (INTI .GT. 0 .AND. IOTRS .GT. 0) ISOLFL = 2

              IF (INTI .GT. 0 .AND. IORTS .GT. 0) ISOLTR = 2

C-fin-parche

              CALL SOURCE_RDCH
     &            (LMXNDL ,NPPEL  ,NUMEL  ,NUMNP
     &            ,AREA   ,CAUX1  ,KXX    ,LNNDEL
     &            ,PAREL  ,SOURCE ,WATVOL,IOPTS(31))

              PRINT*,"He convergido!"
              PRINT*,""
          END IF !IOINITSS .OR.  INTI.GT.0


      ELSE

          PRINT*,"No he convergido"
          PRINT*,""
          IOCONVF = 0
          IOCONVT = 0
          ISOLFL = 0
          ISOLTR = 0

      END IF !IOCONVGL.EQ.1

      IF (IFLAGS(3).NE.0) CALL IO_SUB('COUPLED_FLOW_TRANSPORT',0)

      END SUBROUTINE COUPLED_FLOW_TRANSPORT
