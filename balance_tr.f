      SUBROUTINE BALANCE_TR
     &          (ACTH     ,AFLU     ,AREA     ,ARRC     ,ATRA
     &          ,BETAC    ,BM_ND_TT ,BM_ZN_TT ,BUOYANCY ,CAUDAL
     &          ,CAUX1    ,CAUX2    ,CCALAN   ,CCALIT   ,CFLU
     &          ,CLKCF    ,COECEL   ,COECNP   ,COORD    ,CRDC
     &          ,CREF     ,DELTAT   ,DENSITY  ,DENSREF  ,DFLU
     &          ,DIVQ     ,DTRA     ,FODC     ,GP_COORD ,GRADLOC
     &          ,HAUX1    ,HAUX2    ,IAD_S    ,IADN_S   ,IBCOD
     &          ,IBTCO    ,IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU
     &          ,IDIMDTRA ,IFLAGS   ,INDINIT  ,INTI     ,IOBALDC
     &          ,IOBALGC  ,IODENS   ,IODIM    ,IORECATRA,IOVRWC
     &          ,ISOLFL   ,ISOLTR   ,ISOZ     ,ITPTVAR  ,ITYPAFLU
     &          ,ITYPATRA ,ITYPCFLU ,ITYPDFLU ,ITYPDTRA ,IXCLK
     &          ,IXCON    ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL
     &          ,LTYPE    ,LXARR    ,LXCOE    ,LXFOD    ,LXPAREL
     &          ,LXPOR    ,LXZOR    ,MAXNB    ,MAXPG    ,NFLAGS
     &          ,NMAXT    ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL
     &          ,NUMEL    ,NUMNP    ,NZARR    ,NZCLK    ,NZCOE
     &          ,NZFOD    ,NZPOR    ,NZTRA    ,NZZOR    ,PAREL
     &          ,PARNP    ,POINTWEIGHT        ,THETAT   ,WATVOL
     &          ,WSPECHEAT,ZORC)

********************************************************************************
*
* PURPOSE Computes transport mass balance at k-th time step
*
* DESCRIPTION This subroutine computes the temporal transport mass balance 
*             related to four types of processes: 1) Storage capacity
*                                                 2) First and zero order decay
*                                                    reactions
*                                                 3) Inflow/outflow masses
*                                                 4) Lateral mass fluxes
*             The arrays calculated at this subr. are organized accordingly:
*             1) BM_ZN_TT (NMAXT,2), where the first index indicates the 
*                summation over types of processes of the number of zones 
*                corresponding to each process type. That is:
*                NMAXT=NZPOR+NZFOD+NZZOR+NZCOE*6+NZCLK
*                The second index takes two values: =1: temporal t.m.b. at k-th 
*                                                       time step
*                                                   =2: contribution to global
*                                                        t.m.b. (not used here)
*             2) BM_ND_TT (NUMNP,4,2), where the first index relates mesh nodes.
*                The second one is organized per type of process (same 
*                numeration as above) and the third indicates temporal or global
*                t.m.b. (also as above). 
*           
*             Each process is calculated over elements (storage capacity, decay
*             reactions and recharge with concentration. The inflow/outflow mass 
*             is calculated over mesh nodes.
*             Finally, contribution to global mass balance of the k-th time step
*             is computed.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  ARRC                   Computed areal recharge zonal parameter.              
*  BM_ND_FL               Array containing nodal flow mass balance information
*  BM_ND_TT               Array containing nodal tpt. mass balance information  
*  BM_ZN_TT               Array containing zonal tpt. mass balance information  
*  CAUDAL                 Input/output flow at every node.                      
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  CAUX2                  Array containing diference of concentrations in two   
*                         consecutives times, related to time increment         
*  CCALIT                 Computed concentration at every node                  
*  CCALAN                 Computed concentrations in the previous time step.    
*  COECEL                 External concentration defined over elements  
*  COECNP                 External concentration defined over nodal points
*  CRDC                   Computed retardation coefficient zonal parameter.     
*  FODC                   Computed first order decay zonal values               
*  IBTCO                  Transport boundary condition index                    
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  IXCON                  Ext. concent. (steady) zone number at a given node    
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LXARR                  Areal recharge (steady) zone number at a given element
*  LXCOE                  External concentration zone number at a given element 
*  LXFOD                  First order decay zone number at a given element      
*  LXPOR                  Porosity zone number at a given element               
*  LXZOR                  Zero order reaction number at a given element      
*  WATVOL                 Array containing water volume at mesh elements 
*  ZORC                   Zero order reaction parameter at  mesh elements
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  DELTAT                 Current time increment
*  IOBALDC                Option of computation of nodal flow and/or transport 
*                         mass balance
*  INDBALC                Transport mass balance option:
*                           - 0 = Nothing is done
*                           - 1 = Computation only
*                           - 2 = Computation and writing
*  IOBALGC                Integrated mass balance option.
*                           - 0 = Does not compute gobal m.b.(gmb) and nodal and
*                                 zonal m.b. are instantaneous (i.e., not
*                                 integrated).
*                           - 1 = Comp. gmb, instantaneous zonal and nodal m.b.
*                           - 2 = Does not comp. gmb, integrated zonal and nodal.
*                           - 3 = Comp. gmb, integ.zonal and nodal   
*  ISOLTR                  Transport regime                                      
*  ISOLFL                  Flow regime                                           
*  LMXNDL                 Maximum number of nodes per element                   
*  NFLAGS                 Maximum number of allowed flags                       
*  NMAXT                  Maximum number of processes.
*                         As maximum: NZPOR+NZFOD+NZZOR+NZCOE
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZARR                  Number of areal recharge zones                        
*  NZCOE                  Number of external concentration zones                
*  NZFOD                  Number of zones of first order decay                  
*  NZPOR                  Number of porosity zones                              
*  NZZOR                  Number of zero order reaction zone
*  THETAF                 Time ponderation factor for transport equation
*
* INTERNAL VARIABLES: SCALARS
*
*  ALMN                   Summation variable related to nodal storage capacity
*  ARRN                   Summation variable related to nodal recharge*conc.
*  DIVQ                   Divergence of Darcy's velocities (nodal point)    
*  FODN                   Summation variable related to nodal F.O.D.        
*  INOD                   Dummy counter: loop over nodes of a given element
*  IPOINT                 Denotes the pointer to the beginning of next process
*                         type mass balance information in array BM_ZN_TT
*  IZONE                  Dummy counter: loop over transport parameters'zones
*  L                      Dummy counter: loop over elements
*  NNOD                   Node number: global connectivity matrix
*  NNUD                   Number of nodes of the current element                
*  NPRO                   Dummy counter: loop over mass balance processes
*  NZONEARR               Number of areal recharge zone of a given element
*  NZONECOE               Number of external conc. zone of a given element      
*  NZONEFOD               Number of F.O.D. reaction zone of a given element
*  NZONEPOR               Number of porosity zone of a given element
*  NZONEZOR               Number of zero order reaction zone of a given element
*  SUMALM                 Summation dummy variable related to storage capacity
*                         (element)
*  SUMFOD                 Summation dummy variable related to F.O.D. (element) 
*  WTV_K1                 WATVOL at k+1 associated to a node
*  WTV_K                  WATVOL at k associated to a node
*  WTV_KT                 WATVOL at k+theta associated to a node
*  ZOREL                  Dummy variable related to Z.O.R. 
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  IO_SUB     Checks if subroutine has been entered and passed through
*             successfully
*
* HISTORY: JLF: First coding (Nov-2000)
*          AAR: First coding (Nov-2000)
*
********************************************************************************

      IMPLICIT NONE

C-------------------- External

      INTEGER*4::IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU ,IDIMDTRA
     &          ,INDINIT  ,INTI     ,IOBALDC  ,IOBALGC  ,IODENS
     &          ,IODIM    ,IORECATRA,IOVRWC   ,ISOLFL   ,ISOLTR
     &          ,ITPTVAR  ,ITYPAFLU ,ITYPATRA ,ITYPCFLU ,ITYPDFLU
     &          ,ITYPDTRA ,LMXNDL   ,MAXNB    ,MAXPG    ,NFLAGS
     &          ,NMAXT    ,NPAREL   ,NPPEL    ,NPPNP    ,NUMEL
     &          ,NTYPEL   ,NUMNP    ,NZARR    ,NZCLK    ,NZCOE
     &          ,NZFOD    ,NZTRA    ,NZPOR    ,NZZOR


      REAL*8::BETAC    ,CREF     ,DELTAT   ,DENSREF  ,THETAT
      REAL*8::DENS

      INTEGER*4::IAD_S(MAXNB,NUMNP)  ,IADN_S(NUMNP)     ,IBCOD(NUMNP)
     &          ,IBTCO(NUMNP)        ,IFLAGS(NFLAGS)    ,ISOZ(NZTRA)
     &          ,IXCLK(NUMNP)        ,IXCON(NUMNP)    ,KXX(LMXNDL,NUMEL)
     &          ,LDIM(NUMEL)         ,LNNDEL(NUMEL)     ,LTYPE(NUMEL)
     &          ,LXARR(NUMEL)        ,LXCOE(NUMEL)      ,LXFOD(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL)                  ,LXPOR(NUMEL)
     &          ,LXZOR(NUMEL)       


      REAL*8::ACTH(NUMEL)                  ,AFLU(NUMEL,IDIMAFLU)
     &       ,AREA(NUMEL)                  ,ARRC(NUMEL)
     &       ,ATRA(NUMEL,IDIMATRA)         ,BM_ND_TT(NUMNP,12,2)
     &       ,BM_ZN_TT(NMAXT,2)            ,BUOYANCY(IODIM,LMXNDL,NUMEL)
     &       ,CAUDAL(NUMNP)                ,CAUX1(NUMNP)
     &       ,CAUX2(NUMNP)                 ,CCALAN(NUMNP)
     &       ,CCALIT(NUMNP)                ,CFLU(NUMEL,IDIMCFLU)
     &       ,CLKCF(NUMNP)                 ,COECEL(NUMEL)
     &       ,COECNP(NUMNP)                ,COORD(NUMNP,3)
     &       ,CRDC(NUMEL)                  ,DENSITY(NUMEL)
     &       ,DIVQ(NUMNP)                  ,DTRA(NUMEL,IDIMDTRA)
     &       ,FODC(NUMEL)                  ,GP_COORD(6,8,IODIM)
     &       ,GRADLOC(IODIM,LMXNDL,MAXPG)  ,HAUX1(NUMNP)
     &       ,PAREL(NUMEL,NPPEL)           ,PARNP(NUMNP,NPPNP)
     &       ,POINTWEIGHT(MAXPG,NTYPEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
     &       ,ZORC(NUMEL),dflu(numel,idimdflu),haux2(numnp)

C-------------------- Internal

      INTEGER*4::I        ,IBF      ,IBT      ,INODE    ,IPCLK
     &          ,IPBF     ,IPFLX    ,IPFOD
     &          ,IPOS     ,IPZOD    ,L        ,NNUD
     &          ,NPRO     ,NZONEARR ,NZONECLK ,NZONECOE ,NZONEFOD
     &          ,NZONEPOR ,NZONEZOR


      REAL*8::ALMN     ,AREAL    ,AREALN   ,ARREL    ,ARRELN   ,CAUD
     &       ,COE      ,CONC     ,DENSK    ,DENSK1   ,DENSKTH  ,DENSREC
     &       ,DTMNT    ,FLOW     ,FLOWLEAK ,FLUXDENS ,FODN     ,MASSFLUX
     &       ,RECHRG   ,RETARD   ,RHOAREAL ,RHOAREALN,SUMALM   ,SUMFOD
     &       ,THT1     ,THTINV   ,WTV_KT_AVG         ,WSPECHEAT,ZOREL
     &       ,ZORI

      REAL*8,ALLOCATABLE::PRESCFLUX(:),WTV_K1(:),WTV_K(:),WTV_KT(:)

C-------------------- Indexes for zonal mass balance.

      INTEGER*4,PARAMETER::
     &      IHEAD = 1,
     &      IQ = 2,
     &      ILEAK = 3,
     &      IREC = 4,
     &      ICONC = 5,
     &      IMASS = 6

C-------------------- Indexes for nodal mass balance.

      INTEGER*4,PARAMETER::
     &       INDPOR  = 1      ,INDMASS  =  7
     &      ,INDH    = 2      ,INDCLK   =  8
     &      ,INDQ    = 3      ,INDFOD   =  9
     &      ,INDLK   = 4      ,INDZOD   = 10 
     &      ,INDREC  = 5      ,INDMTDIF = 11
     &      ,INDCON  = 6      ,INDLTFL  = 12


C-------------------- First executable statement.

C-------------------- Checks the succesfull entry to current subr.

      IF(IFLAGS(3).EQ.1) CALL IO_SUB('BALANCE_TR',0)

C-------------------- Initialises the temporal mass balance arrays 
C-------------------- (zones).If so desired, initialises the temporal
C-------------------- mass balance arrays (nodal points)

      IF (INDINIT.NE.0) THEN

          BM_ZN_TT(:,1) = 0D0
          IF (IOBALDC.NE.0) BM_ND_TT(:,:,1) = 0D0

      END IF !INDINIT.NE.0

C-------------------- Define a false DT (DTMNT) which equals 1.0 for
C-------------------- instantaneous mass balance, and the true DT for
C-------------------- integrated mass balance computations

      IF (IOBALGC.LE.1) THEN

          DTMNT = 1D0

      ELSE

          DTMNT = DELTAT

      END IF !IOBALGC.LE.1


********************************************************************************
*        Storage capacity, F.O.D. reactions. and Z.O.D                         *
********************************************************************************

      IPFLX = NZPOR
      IPCLK = IPFLX + NZCOE*6
      IPFOD = IPCLK + NZCLK
      IPZOD = IPFOD + NZFOD
*     IPMATDIF = IPZOD + NZZOR Not used yet      


C-------------------- Notice that computations related to storage
C-------------------- capacity are only done while dealing transient
C-------------------- transport cases.

      THTINV = 1D0/THETAT
      THT1 = THETAT - 1D0

      ALLOCATE (WTV_K1(LMXNDL),WTV_K(LMXNDL),WTV_KT(LMXNDL))

      WTV_K1 = 0D0
      WTV_K = 0D0
      WTV_KT = 0D0

      DO L=1,NUMEL

C-------------------- Initialize element variables

          NNUD = LNNDEL(L)
          AREAL = AREA(L)
          AREALN = AREAL/NNUD
          RHOAREAL = DENSITY(L)*AREA(L)
          RHOAREALN = RHOAREAL/NNUD

          NZONEPOR = LXPOR(L)           ! Porosity zone of current elem.
          NZONEFOD = LXFOD(L)           ! F.O.D. zone of current elem.


          RETARD = CRDC(L)*ACTH(L)

          SUMFOD = 0D0
          SUMALM = 0D0

C-------------------- Computes the water content at current time(k+1),
C-------------------- on basis of water content at time k+THETAF 
C-------------------- (WATVOL (L,2)) and at time k (WATVOL(L,1)). Only 
C-------------------- done when one is solving transient flow, as in 
C-------------------- this case the water content can vary, and when it 
C-------------------- is going to be used (transient tpt.)
C-------------------- Notice that water content is computed at each
C-------------------- node. Then it is not necesary to divide by the
C-------------------- number of nodes when computing mass balance (the
C-------------------- usual AREAL/NNUD factor must be only AREAL).

          IF (IOVRWC.LT.2) THEN !Elementwise


              IF (ISOLTR.EQ.2.AND.ISOLFL.EQ.2.AND.THETAT.GT.0.D0) THEN


                  WTV_K(1:NNUD) = WATVOL(1,L,1)

                  WTV_KT(1:NNUD) = WATVOL(1,L,2)

                  WTV_K1(1:NNUD) = 
     &                      (WTV_KT(1:NNUD)+(THT1*WTV_K(1:NNUD)))*THTINV


              ELSE

                  WTV_K(1:NNUD) = WATVOL(1,L,2)
                  WTV_K1(1:NNUD) = WTV_K(1:NNUD)
                  WTV_KT(1:NNUD) = WTV_K(1:NNUD)

              END IF !ISOLTR.EQ.2.AND.ISOLFL.EQ.2.AND.THETAF.GT.0.D0

                  WTV_KT_AVG = WATVOL(1,L,2)

          ELSE !Nodewise

              IF (ISOLTR.EQ.2.AND.ISOLFL.EQ.2.AND.THETAT.GT.0.D0) THEN

                  WTV_K(1:NNUD) = WATVOL(1:NNUD,L,1)

                  WTV_KT(1:NNUD) = WATVOL(1:NNUD,L,2)

                  WTV_K1(1:NNUD) =
     &                      (WTV_KT(1:NNUD)+(THT1*WTV_K(1:NNUD)))*THTINV


              ELSE

                  WTV_K(1:NNUD) =  WATVOL(1:NNUD,L,2)
                  WTV_K1(1:NNUD) = WTV_K(1:NNUD)
                  WTV_KT(1:NNUD) = WTV_K(1:NNUD)

              END IF !ISOLTR.EQ.2.AND.ISOLFL.EQ.2.AND.THETAF.GT.0.D0

              WTV_KT_AVG = 0D0

              DO I=1,NNUD

                  WTV_KT_AVG = WTV_KT_AVG + WTV_KT(I)

              END DO !I=1,NNUD

                  WTV_KT_AVG = WTV_KT_AVG/NNUD

          END IF !IOVRWC.LT.2


          DO I=1,NNUD

              INODE = KXX(I,L)

C-------------------- Computation of nodal mass balance related to 
C-------------------- F.O.D. reactions

              IF (NZONEFOD.NE.0) THEN

                  IF (IODENS.EQ.1) THEN

                      DENSKTH = DENS(DENSREF,BETAC,CAUX1(INODE),CREF)

                  ELSE

                      DENSKTH = 1d0

                  END IF !IODENS.EQ.1

                  FODN = FODC(L)*(WTV_KT(I)+RETARD)*DENSKTH*AREALN
     &                    *CAUX1(INODE)

                  IF (IOBALDC.NE.0) THEN

                      BM_ND_TT(INODE,INDFOD,1)=BM_ND_TT(INODE,INDFOD,1)
     &                                        - FODN*DTMNT

                      BM_ND_TT(INODE,INDLTFL,1) = 
     &                            BM_ND_TT(INODE,INDLTFL,1) + FODN*DTMNT

                  END IF !IOBALDC.NE.0

C--------------------  Updates dummy variables

                  SUMFOD = SUMFOD - FODN

              END IF !NZONEFOD.NE.0

C--------------------  Computation of nodal mass balance related to 
C--------------------  storage capacity

              IF (INTI.NE.0 .AND. ISOLTR.GT.1) THEN !Transient transport

                  IF (IODENS.EQ.1) THEN

                      DENSK = DENS(DENSREF,BETAC,CCALAN(INODE),CREF)
                      DENSK1 = DENS(DENSREF,BETAC,CCALIT(INODE),CREF)

                      IF (IOVRWC.LT.2) THEN
                          DENSKTH = DENSITY(L)
                      ELSE
                         DENSKTH = DENS(DENSREF,BETAC,CAUX1(INODE),CREF)
                      END IF !IOVRWC.LT.2

                  ELSE

                      DENSK = 1D0
                      DENSK1 = 1D0
                      DENSKTH = 1D0

                  END IF !IODENS.EQ.1

                  IF (IFLAGS(28).LT.1) THEN

                      IF (ITPTVAR.EQ.0) THEN !Solute transport

                          ALMN = (
     &                        (WTV_K1(I) + RETARD)*CCALIT(INODE)*DENSK1
     &                       -(WTV_K(I) + RETARD)*CCALAN(INODE)*DENSK
     &                         )*AREALN/DELTAT

                      ELSE IF (ITPTVAR.EQ.1) THEN !heat transport

                          ALMN = (
     &               (DENSK1*WTV_K1(I)*WSPECHEAT + RETARD)*CCALIT(INODE)
     &              -(DENSK*WTV_K(I)*WSPECHEAT   + RETARD)*CCALAN(INODE)
     &                            )*AREALN/DELTAT

                      END IF !ITPTVAR.EQ.0,1

                  ELSE


                      ALMN = DENSKTH*(WTV_KT(I) + RETARD)*CAUX2(INODE)
     &                       *AREALN

                  END IF !IFLAGS(28).LT.1

C--------------------  nodal mass balance?

                  IF (IOBALDC.NE.0) THEN

                      BM_ND_TT(INODE,INDPOR,1)= BM_ND_TT(INODE,INDPOR,1)
     &                                         + ALMN*DTMNT

                  END IF !IOBALDC.NE.0

C--------------------  Updates dummy variables

                  SUMALM = SUMALM + ALMN

              END IF  ! INTI.NE.0 .AND. ISOLTR.GT.1

          END DO !I=1,NNUD

C--------------------  Computes (only sense in transient transport) the
C--------------------  zonal mass balance due to storage capacity

          IF (ISOLTR.EQ.2) THEN     ! Transient transport
              
              BM_ZN_TT(NZONEPOR,1) = BM_ZN_TT(NZONEPOR,1) + SUMALM*DTMNT

          END IF !ISOLTR.EQ.2

C-------------------- Finally, stores the information of F.O.D. m.b

          IF (NZONEFOD.NE.0)  THEN  ! If that information exists...

              IPOS = IPFOD+NZONEFOD
              BM_ZN_TT(IPOS,1) = BM_ZN_TT(IPOS,1) + SUMFOD*DTMNT

          END IF !NZONEFOD.NE.0


C--------------------  Z.O.D.

          IF (NZZOR.NE.0) THEN

              NZONEZOR = LXZOR(L)

              IF (NZONEZOR.NE.0) THEN ! A zero O.D.R. is ocurring

                  ZOREL = ZORC(L)*(WTV_KT_AVG+RETARD)*RHOAREAL

                  IPOS = IPZOD + NZONEZOR
                  BM_ZN_TT(IPOS,1)= BM_ZN_TT(IPOS,1) + ZOREL*DTMNT

                  IF (IOBALDC.NE.0) THEN  ! Nodal mass balance

                      NNUD = LNNDEL(L)

                        DO I=1,NNUD

                          INODE = KXX(I,L)

                          ZORI = ZORC(L)*(WTV_KT(I)+RETARD)*RHOAREAL

                          BM_ND_TT(INODE,INDZOD,1) = 
     &                             BM_ND_TT(INODE,INDZOD,1) + ZORI*DTMNT

                      END DO !I=1,NNUD

                  END IF !IOBALDC.NE.0

              END IF !NZONEZOR.NE.0
          END IF !NZZOR.NE.0)

*******************************************************************************     
*                      Mass fluxes from/to outside
*******************************************************************************

C--------------------  1st. subprocess: recharge*ext. concentration

          IF (NZARR.NE.0.AND.NZCOE.NE.0) THEN

              NZONEARR = LXARR(L)
              NZONECOE = LXCOE(L)
              RECHRG = ARRC(L)

              IF (NZONEARR.NE.0.AND.NZONECOE.NE.0.AND.RECHRG.GT.0) THEN

                  CONC = COECEL(L)


                  IF (IODENS.EQ.1) THEN

                      DENSREC = DENS(DENSREF,BETAC,CONC,CREF)

                  ELSE

                      DENSREC = 1D0

                  END IF !IODENS.EQ.0


                  ARREL = RECHRG*CONC*DENSREC*AREAL*DTMNT

                  IF (ITPTVAR.EQ.1) THEN

                      ARREL = ARREL*WSPECHEAT

                  END IF !ITPTVAR.EQ.1

                  IPOS = IPFLX + (NZONECOE-1)*6 + IREC
                  BM_ZN_TT(IPOS,1) = BM_ZN_TT(IPOS,1)+ARREL

                  IF (IOBALDC.NE.0) THEN

                      ARRELN = ARREL/NNUD

                      DO I=1,NNUD

                          INODE = KXX(I,L)

                          BM_ND_TT(INODE,INDREC,1) = 
     &                                 BM_ND_TT(INODE,INDREC,1) + ARRELN

                          IF (IORECATRA.EQ.1) THEN

                              BM_ND_TT(INODE,INDLTFL,1) =
     &                            BM_ND_TT(INODE,INDLTFL,1) + ARRELN

                          END IF !IONEWT.EQ.0

                      END DO !I=1,NNUD

                  END IF !IOBALDC.NE.0

              END IF !NZONEARR.NE.0.AND.NZONECOE.NE.0.AND.RECHRG.GT.0

          END IF !NZARR.NE.0.AND.NZCOE.NE.0

      END DO !L=1,NUMEL

      DEALLOCATE(WTV_K,WTV_K1,WTV_KT)

C--------------------  2nd. subprocess: inflow/outflow masses

      ALLOCATE(PRESCFLUX(NUMNP))
      PRESCFLUX(:)= 0D0

      CALL COMFLOW_PRESC_CONC
     &    (AFLU     ,AREA     ,ARRC     ,ATRA     ,BETAC
     &    ,BUOYANCY ,CAUDAL   ,CAUX1    ,CAUX2    ,CFLU
     &    ,COECEL   ,COORD    ,CREF     ,DFLU     ,DENSITY
     &    ,DENSREF  ,DTRA     ,PRESCFLUX,GP_COORD ,GRADLOC
     &    ,HAUX1    ,HAUX2    ,IAD_S    ,IADN_S   ,IBCOD
     &    ,IBTCO    ,IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU
     &    ,IDIMDTRA ,IODENS   ,IODIM    ,IORECATRA,ISOLFL
     &    ,ISOLTR   ,ISOZ     ,ITYPAFLU ,ITYPATRA ,ITYPCFLU
     &    ,ITYPDFLU ,ITYPDTRA ,KXX      ,LDIM     ,LMXNDL
     &    ,LNNDEL   ,LTYPE    ,LXARR    ,LXPAREL  ,MAXNB
     &    ,MAXPG    ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL
     &    ,NUMEL    ,NUMNP    ,NZARR    ,NZTRA    ,PAREL
     &    ,PARNP    ,POINTWEIGHT        ,THETAT)
     
      DO I=1,NUMNP

          IBT = IBTCO(I)
          IBF = IBCOD(I) 

          NZONECOE = IXCON(I) ! External concentration zone
          COE = COECNP(I)

          SELECT CASE(IBT)

C-------------------- Prescribed Concentration

          CASE(1)

              MASSFLUX = PRESCFLUX(I)*DTMNT

              IF (ITPTVAR.EQ.1) THEN !Heat transport

                  MASSFLUX = MASSFLUX*WSPECHEAT

              END IF !ITPTVAR.EQ.1

              IPOS = IPFLX + (NZONECOE-1)*6 + ICONC

              BM_ZN_TT(IPOS,1) = BM_ZN_TT(IPOS,1) + MASSFLUX

              IF (IOBALDC.NE.0) THEN

                  BM_ND_TT(I,INDCON,1) = BM_ND_TT(I,INDCON,1) + MASSFLUX

              END IF !IOBALDC.NE.0

C-------------------- Mass flow

          CASE(2,3)

              SELECT CASE(IBF)

              CASE(1,2,3,4) !Head, flow, leakage, mixed

                  IPBF = IBF
                  IF (IBF.EQ.4) IPBF = 2 !mixed written in flow pos.

                  CAUD = CAUDAL(I)

                  IPOS = IPFLX + (NZONECOE-1)*6 + IPBF

                  IF (CAUD.GE.0D0) THEN  ! Input mass flux

                      CONC = COE

                  ELSE ! Output mass flux

                      CONC = CAUX1(I)

                  END IF !CAUDAL(I).GE.0D0

                  IF (IODENS.EQ.1) THEN

                      FLUXDENS = DENS(DENSREF,BETAC,CONC,CREF)

                  ELSE

                      FLUXDENS = 1D0

                  END IF !IODENS.EQ.1

                  MASSFLUX = FLUXDENS*CONC*CAUD*DTMNT

                  IF (ITPTVAR.EQ.1) THEN !Heat transport

                      MASSFLUX = MASSFLUX*WSPECHEAT

                  END IF !ITPTVAR.EQ.1

                  BM_ZN_TT(IPOS,1) = BM_ZN_TT(IPOS,1) + MASSFLUX

                  IF (IOBALDC.NE.0) THEN

                      IPOS = INDPOR + IBF

                      BM_ND_TT(I,IPOS,1) = BM_ND_TT(I,IPOS,1) + MASSFLUX

                  END IF !IOBALDC.NE.0

              END SELECT !IBF

C-------------------- Prescribed mass

          CASE(4)

              IPOS = IPFLX + (NZONECOE-1)*6 + IMASS

              BM_ZN_TT(IPOS,1) = BM_ZN_TT(IPOS,1) + COE*DTMNT

              IF (IOBALDC.NE.0) THEN

                  BM_ND_TT(I,INDMASS,1)=BM_ND_TT(I,INDMASS,1)+COE*DTMNT

              END IF !IOBALDC.NE.0

C-------------------- Concentration leakage

          CASE(5)

              NZONECLK = IXCLK(I)

              IPOS = IPCLK + NZONECLK

              MASSFLUX = CLKCF(I)*(COE - CAUX1(I))*DTMNT

              BM_ZN_TT(IPOS,1) = BM_ZN_TT(IPOS,1) + MASSFLUX

              IF (IOBALDC.NE.0) THEN

                  BM_ND_TT(I,INDCLK,1) = BM_ND_TT(I,INDCLK,1) + MASSFLUX

              END IF !IOBALDC.NE.0

          END SELECT !IBT
 
      END DO !I=1,NUMNP

      DEALLOCATE(PRESCFLUX)

********************************************************************************
*        Lateral fluxes.
********************************************************************************

C--------------------  Lateral fluxes are equal to
C--------------------  -ATRA*C -c*DIVq + FOD*C - Rec*C
C--------------------  F.O.D. and recharge contribution have been already
C--------------------  computed (recharge only when Picard method is used).

      IF (IOBALDC.NE.0) THEN  ! Nodal mass balance

C--------------------  First, computes  -(ATRA*ck+theta)*DeltaT

          CALL PROD_MAT_VEC
     &        (-DTMNT     ,IAD_S    ,IADN_S   ,IDIMATRA ,NUMEL    ,NUMNP
     &        ,1        ,ITYPATRA ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &        ,KXX      ,LNNDEL   ,ATRA     ,BM_ND_TT(1,INDLTFL,1)
     &        ,CAUX1)     

C-------------------- Substracts c*div(q)

          BM_ND_TT(:,INDLTFL,1) = BM_ND_TT(:,INDLTFL,1)
     &                          - DIVQ(:)*CAUX1(:)*DTMNT


      END IF !IOBALDC.NE.0

********************************************************************************
*            Global mass balance
********************************************************************************

      IF (ISOLTR.EQ.2.AND.(IOBALGC.EQ.1 .OR. IOBALGC.EQ.3)) THEN

C--------------------  Integration factor assigned, if needed.
C--------------------  If instantaneous m. b. has been already computed
C--------------------  DTMNT = 1, in order not to integrate again.
C--------------------  Otherwise, DTMNT = DELTAT, to integrate the m. b.

          IF (IOBALGC.EQ.1) THEN  

              DTMNT = DELTAT

          ELSE IF (IOBALGC.EQ.3) THEN

              DTMNT=1D0

          END IF !IOBALGC.EQ.1,3


          BM_ZN_TT(1:NMAXT,2) = BM_ZN_TT(1:NMAXT,2)
     &                         + BM_ZN_TT(1:NMAXT,1)*DTMNT

          IF (IOBALDC.NE.0) THEN
                                 
              DO NPRO=1,12 ! Starts loop over process

                  BM_ND_TT(1:NUMNP,NPRO,2) = BM_ND_TT(1:NUMNP,NPRO,2)
     &                                  + BM_ND_TT(1:NUMNP,NPRO,1)*DTMNT

              END DO !NPRO=1,12

          END IF !IOBALDC.NE.0

      END IF !ISOLTR.EQ.2.AND.(IOBALGC.EQ.1 .OR. IOBALGC.EQ.3

C--------------------  Checks the succesfull entry to current subr.

      IF(IFLAGS(3).EQ.1) CALL IO_SUB('BALANCE_TR',1)


      END SUBROUTINE BALANCE_TR
