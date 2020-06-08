      SUBROUTINE BALANCE_FL
     &          (AFLU     ,ALFC     ,AREA     ,ARRC     ,ATRA
     &          ,BETAC    ,BM_ND_FL ,BM_ZN_FL ,BUOYANCY ,CAUDAL
     &          ,CAUX1    ,CAUX2    ,CFLU     ,CHPC     ,CLK
     &          ,COECEL   ,COECNP   ,COORD    ,CREF     ,DBUOYANCY
     &          ,DELTAT   ,DENSITY  ,DENSREF  ,DFLU     ,DFLUDFLU
     &          ,DFLUDTRA ,DPARELDH ,DTRA     ,GP_COORD ,GRADLOC
     &          ,GRAVEL   ,GRDFF    ,HCALAN   ,HCALIT   ,HAUX1
     &          ,HAUX2    ,IAD_S    ,IADN_S   ,IBCOD    ,IBTCO
     &          ,IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU ,IDIMDTRA
     &          ,IFLAGS   ,INTI     ,IOBALDH  ,IOBALGH  ,IOCONSRC
     &          ,IODENS   ,IODIM    ,IOEQT    ,IOFLLI   ,IORECATRA
     &          ,IOVRWC   ,ISOLFL   ,ISOLTR   ,ISOZ     ,ITYPAFLU
     &          ,ITYPATRA ,ITYPCFLU ,ITYPDFLU ,ITYPDTRA ,IXALF
     &          ,IXCHP    ,IXCONC   ,IXQQP    ,KXX      ,LDIM
     &          ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXARR    ,LXPAREL
     &          ,LXSTG    ,MAXNB    ,MAXPG    ,NFLAGS   ,NMAXF
     &          ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL   ,NUMEL
     &          ,NUMNP    ,NZALF    ,NZARR    ,NZCHP    ,NZQQP
     &          ,NZSTG    ,NZTRA    ,PAREL    ,PARNP    ,POINTWEIGHT
     &          ,QQPC     ,STGC     ,THETAT   ,WATVOL)

********************************************************************************
*     
*     PURPOSE Computes flow mass balances: time to time f.m.b. and its contribution
*     to the global one (if so desired). Also computes de divergence of 
*     Darcy's fluxes, to be used for tpt. m.b. purposes.
*     
*     
*     DESCRIPTION Loop over types of flow parameters
*     
*     EXTERNAL VARIABLES: ARRAYS
*     
*     AFLU                   Matrix of finite elements equations for flow problem  
*     No boundary conditions are included on it.            
*     ALFC                   Computed leakage zonal parameter                      
*     AREA                   Element size (length for 1-D elem, area for 2-D,      
*                            volume for 3-D)
*     ARRC                   Computed areal recharge zonal parameter.              
*     BM_ND_FL               Nodal flow mass balance_array
*     BM_ZN_FL               Zonal flow mass balance_array                         
*     CHPC                   Computed prescribed head zonal parameter              
*     DFLU                   Matrix of finite elements equations for flow          
*                            problem related to storage term.
*     HAUX1                  Array containing HEADS, ponderated by THETAF          
*     time factor                                           
*     HAUX2                  Array containing diference of HEADS in two            
*                            consecutives times.
*     IBCOD                  Flow boundary condition index                         
*     IFLAGS                 Array with different writing options. Used mainly for 
*                            debugging.
*     IXALF                  Leakage (steady) zone number at a given node          
*     IXCHP                  Presc. head (steady) zone number at a given node      
*     IXQQP                  Presc. flow (steady) zone number at a given node      
*     KXX                    Node numbers of every element (counterclockwise       
*                            order).
*     LNNDEL                 Number of nodes at every element                      
*     LXARR                  Areal recharge (steady) zone number at a given element
*     LXSTG                  Storage coefficient zone number at a given element    
*     QQPC                   Computed presc. flow zonal values                     
*     STGC                   Computed storage zonal values                         
*     
*     INTERNAL VARIABLES: ARRAYS
*     
*     
*     EXTERNAL VARIABLES: SCALARS
*     
*     DELTAT                 Current time increment
*     IOBALDH                Detailed flow mass balance option (nodal f.m.b.)
*     IOBALGH                Global flow mass balance (integrated in time)         
*     IOCNSF                 Scheme for storage term in flow problem               
*     ISOLFL                 Indicates if flow eq. has been solved in the current time
*                            step. (0 not solved, 1 steady solved, 2 transient solved)
*     ISOLTR                 Indicates if tpt. eq. has been solved in the current time
*                            step. (0 not solved, 1 steady solved, 2 transient solved)
*     LMXNDL                 Maximum number of nodes per element                   
*     NBAND                  Half Bandwith (maximum difference between the         
*                            numbers of two nodes belonging to the same element)
*     NBAND1                 Used to dimension. It is equal to NBAND+1             
*     NFLAGS                 Maximum number of allowed flags                       
*     NMAXF                  Maximum number of flow parameter zones
*     NUMEL                  Number of elements                                    
*     NUMNP                  Number of nodes                                       
*     NZARR                  Number of areal recharge zones                        
*     NZCHP                  Number of prescribed head zones                       
*     NZQQP                  Number of prescribed flow zones                       
*     NZSTG                  Number of storage Coefficient zones                   
*     
*     INTERNAL VARIABLES: SCALARS
*     
*     CAP                    Nodal capacity
*     CONTRIB                Contribution of recharge to f.m.b (zonal or nodal)
*     FLOW                   Lateral flows    
*     IGROUP                 Grouping index
*     INOD                   Current node                                     
*     IPOINT                 Pointer to last position of previous flow event
*     IZONE                  Number of event zone       
*     KNOD                   Auxiliar variable. Similar to INODE
*     L                      Number of current element                   
*     NNOD                   Node number (Glbal connectivity)
*     NNUD                   Number of nodes of the current element                
*     NZONEALF               Current leakage zone                      
*     NZONEARR               Current recharge zone
*     NZONECHP               Current presc. head level zone
*     NZONEQQP               Current presc. flow zone                        
*     NZONESTO               Current storavity zone
*     SUM                    Auxiliar variable (storage)
*     
*     FUNCTIONS AND SUBROUTINES REFERENCED
*     
*     IO_SUB                                                                       
*     MUL_MATBAN_VEC                                                               
*     ZERO_ARRAY                                                                   
*     
*     HISTORY: AAR: First coding (Nov-2000)
********************************************************************************

      IMPLICIT NONE

C-------------------- External integer variables: scalars 

      INTEGER*4::IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU ,IDIMDTRA
     &          ,INTI     ,IOBALDH  ,IOBALGH  ,IOCONSRC ,IODENS
     &          ,IODIM    ,IOEQT    ,IOFLLI   ,IORECATRA,IOVRWC
     &          ,ISOLFL   ,ISOLTR   ,ITYPAFLU ,ITYPATRA ,ITYPCFLU
     &          ,ITYPDFLU ,ITYPDTRA ,LMXNDL   ,MAXNB    ,MAXPG
     &          ,NFLAGS   ,NMAXF    ,NPAREL   ,NPPEL    ,NPPNP
     &          ,NTYPEL   ,NUMEL    ,NUMNP    ,NZALF    ,NZARR
     &          ,NZCHP    ,NZQQP    ,NZSTG    ,NZTRA

      REAL*8::BETAC    ,CREF     ,DELTAT   ,DENSREF  ,THETAT
      REAL*8::DENS

C-------------------- External integer variables: arrays

      INTEGER*4::IAD_S(MAXNB,NUMNP),IADN_S(NUMNP)        ,IBCOD(NUMNP)
     &          ,IBTCO(NUMNP)
     &          ,IFLAGS(NFLAGS)    ,ISOZ(NZTRA)          ,IXALF(NUMNP)
     &          ,IXCHP(NUMNP)      ,IXCONC(NUMNP)        ,IXQQP(NUMNP)
     &          ,KXX(LMXNDL,NUMEL) ,LDIM(NUMEL)          ,LNNDEL(NUMEL)
     &          ,LTYPE(NUMEL)      ,LXARR(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL)                   ,LXSTG(NUMEL)


C-------------------- External real variables

      REAL*8::AFLU(NUMEL,IDIMAFLU) ,ALFC(NUMNP)          ,AREA(NUMEL)
     &       ,ARRC(NUMEL)          ,ATRA(NUMEL,IDIMATRA)
     &       ,BM_ND_FL(NUMNP,8,2)  ,BM_ZN_FL(NMAXF,2)    ,BUOYANCY(*)
     &       ,CAUDAL(NUMNP)        ,CAUX1(NUMNP)         ,CAUX2(NUMNP)
     &       ,CFLU(NUMEL,IDIMCFLU) ,CHPC(NUMNP)          ,CLK(NUMNP)
     &       ,COECEL(NUMEL)        ,COECNP(NUMNP)        ,COORD(NUMNP,3)
     &       ,DBUOYANCY(*)         ,DENSITY(NUMEL)  
     &       ,DFLU(NUMNP,IDIMDFLU) ,DFLUDFLU(*)          ,DFLUDTRA(*)
     &       ,DPARELDH(NPPEL,NUMEL),DTRA(NUMEL,IDIMDTRA)
     &       ,GP_COORD(6,8,IODIM)  ,GRADLOC(*)
     &       ,GRAVEL(NUMEL,3)      ,GRDFF(IODIM,LMXNDL,NUMEL)
     &       ,HCALAN(NUMNP)        ,HCALIT(NUMNP)
     &       ,HAUX1(NUMNP)         ,HAUX2(NUMNP)
     &       ,PAREL(NUMEL,NPPEL)   ,PARNP(NUMNP,NPPNP)   ,POINTWEIGHT(*)
     &       ,QQPC(NUMNP)          ,STGC(NUMEL)   
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)


C-------------------- Internal variables: integers

      INTEGER*4::I        ,IB       ,IBT      ,IGROUP   ,INODE
     &          ,IPOINT   ,IPOS     ,IPOSC    ,IPOSH    ,IPOSL
     &          ,IPOSNOD  ,IPOSQ    ,L        ,NNUD
     &          ,NZONEALF ,NZONEARR ,NZONECEXT,NZONECHP,NZONECONC
     &          ,NZONEQQP ,NZONESTO ,NZSTGW

C-------------------- Internal variables: real*8

      REAL*8::AREALN,CAP      ,CEXT     ,CONTRIB  ,DENSREC  ,FLOW
     &       ,MASS            ,SUMW     ,SUMWI    ,THETAT1  ,THETATINV
     &       ,DENSNODE        ,SUMHI    ,SUMH

      REAL*8,ALLOCATABLE::FLOWAUX(:),WTVK(:),WTVK1(:),WTVKTH(:)



C-------------------- Checks the succesfull entry to current subr.

      IF(IFLAGS(3).EQ.1) CALL IO_SUB('BALANCE_FL',0)

C-------------------- Initializes temporal components of f.m.b. arrays

      CALL ZERO_ARRAY(BM_ZN_FL,NMAXF) ! Zonal mass b.

      IF (IOBALDH.NE.0) THEN    ! Nodal mass.b.

          DO I=1,8
              CALL ZERO_ARRAY (BM_ND_FL(1,I,1),NUMNP)
          ENDDO

      END IF !IOBALDH.NE.0

********************************************************************************
***                      STORAGE CAPACITY                                    ***
********************************************************************************

      IF (IODENS.EQ.1) THEN

          ALLOCATE (WTVK(LMXNDL),WTVK1(LMXNDL),WTVKTH(LMXNDL))
          WTVK = 0D0
          WTVK1 = 0D0
          WTVKTH = 0D0

      END IF !IODENSEQ.1

C-------------------- Notice that computations related to storage capacity 
C-------------------- are only done while dealing transient cases.

      IF (INTI.NE.0 .AND. ISOLFL.GT.1) THEN !Only sense in transient state

          THETATINV = 1D0/THETAT
          THETAT1 = THETAT - 1D0

          DO L=1,NUMEL          ! Starts loop over elements

              SUMH = 0D0         ! Dummy variable

              NZONESTO = LXSTG(L) ! Storage zone

              NNUD = LNNDEL(L)    ! Number of nodes of current element

              AREALN = AREA(L)/NNUD

              CAP = STGC(L)*AREALN


              IF(IODENS.EQ.1) THEN

                  SUMW = 0D0

                  SELECT CASE(IOVRWC)

                      CASE(0) !constant

                      WTVK(:) = WATVOL(1,L,1)
                      WTVK1(:) = WTVK(:)
                      WTVKTH(:) = WTVK(:)

                      CASE(1) !variable elementwise

C-------------------- W(k+1) = 1/theta * (W(k+theta) +(theta-1)*W(k))

                        WTVK(:) = WATVOL(1,L,1)
                        WTVKTH(:) = WATVOL(1,L,2)
                        WTVK1(:) = THETATINV*(WATVOL(1,L,2)
     &                                          + THETAT1*WATVOL(1,L,1))

                      CASE(2) !variable nodewise

                          DO I=1,NNUD

                            WTVK(I) = WATVOL(I,L,1)
                            WTVKTH(I) = WATVOL(I,L,2)
                            WTVK1(I) = THETATINV*(WATVOL(I,L,2)
     &                                          + THETAT1*WATVOL(I,L,1))
                          END DO !I=1,NNUD

                  END SELECT !IOVRWC

              END IF !IODENS.EQ:1

              DO I=1,NNUD

                  INODE = KXX(I,L) 

                  IF (IODENS.EQ.1) THEN

                      IF (IOVRWC.LT.2) THEN

                          DENSNODE = DENSITY(L)

                      ELSE

                          DENSNODE=DENS(DENSREF,BETAC,CAUX1(INODE),CREF)

                      END IF !IOVRWC.LT.2
c-opción-1
c                     SUMWI = DENSITY(L)*AREALN*BETAC*
c     &                                 (WTVK1(I)*CCALIT(INODE)
c     &                                   - WTVK(I)*CCALAN(INODE))/DELTAT

c-opción-2

                      SUMWI = DENSNODE*AREALN*BETAC
     &                            *WTVKTH(I)*CAUX2(INODE)

c-opción-3
c                   DENSK = DENS(DENSREF,BETAC,CCALAN(INODE),CREF)
c                   DENSK1 = DENS(DENSREF,BETAC,CCALIT(INODE),CREF)
c                   SUMWI = WTVKTH(I)*(DENSK1-DENSK)*AREALN/DELTAT

                      SUMW = SUMW + SUMWI

                  ELSE

                      DENSNODE = 1D0

                  END IF ! IODENS.EQ.1

                      SUMHI = DENSNODE*(HCALIT(INODE)
     &                                 - HCALAN(INODE))/DELTAT

                      SUMH = SUMH + SUMHI

C-------------------- Nodal flow mass balance. Storage capacity contr.

                  IF (IOBALDH.NE.0) THEN

                      BM_ND_FL(INODE,1,1) = BM_ND_FL(INODE,1,1)
     &                                      + CAP*SUMHI

                      IF (IODENS.EQ.1) THEN

                          BM_ND_FL(INODE,2,1) = BM_ND_FL(INODE,2,1)
     &                                          + SUMWI

                      END IF !IODENS.EQ.0

                  END IF !IOBALDH.NE.0

              END DO !I=1,NNUD

C-------------------- Zonal flow mass balance. Storage capacity contr.

              BM_ZN_FL(NZONESTO,1) = BM_ZN_FL(NZONESTO,1) + CAP*SUMH

              IF (IODENS.EQ.1) THEN

                  NZSTGW = NZONESTO+NZSTG
                  BM_ZN_FL(NZSTGW,1) = BM_ZN_FL(NZSTGW,1) + SUMW

              END IF !IODENS.EQ.0

          END DO !L=1,NUMEL

      END IF !INTI.NE.0 .AND. ISOLFL.GT.1


      IF (IODENS.EQ.1) THEN

          IF (ALLOCATED(WTVK1)) THEN

              DEALLOCATE (WTVK1,WTVK,WTVKTH)

          END IF  !ALLOCATED(WTVK1)

      END IF !IODENSEQ.1

********************************************************************************
***                            AREAL RECHARGE                                ***
********************************************************************************

C-------------------- Initializes pointer IPOINT to last position of
C-------------------- storativity terms (twice NZSTG if density dependent flow
C-------------------- is being solved).

      IPOINT = NZSTG*(IODENS + 1)

      IF (NZARR.NE.0) THEN

          DO L=1,NUMEL

              NZONEARR = LXARR(L)

              IF (NZONEARR.NE.0) THEN

                  NNUD = LNNDEL(L)

C-------------------- Zonal flow mass balance. Areal recharge contr.

                  NZONECEXT = LXPAREL(L,10)

                  IF (NZONECEXT.GT.0 .AND. ARRC(L).GT.0) THEN

                      CEXT = PAREL(L,15)

                  ELSE

                      CEXT = 0D0

                  END IF !NZONECEXT.GT.0 .AND. RECHARGE.GT.0

                  IF (IODENS.EQ.1) THEN

                      DENSREC =  DENS(DENSREF,BETAC,CEXT,CREF)

                  ELSE

                      DENSREC = 1.D0

                  END IF !IODENS.EQ.1

                  CONTRIB = DENSREC*ARRC(L)*AREA(L) ! Recharge*size(l)

                  BM_ZN_FL(NZONEARR+IPOINT,1) =
     &                 BM_ZN_FL(NZONEARR+IPOINT,1) + CONTRIB

                  CONTRIB = CONTRIB/NNUD ! Nodal contrib.

                  DO I=1,NNUD ! Loop over element nodes

                      INODE=KXX(I,L) ! Node number

C-------------------- Nodal flow mass balance. Areal recharge contr.

                      IF (IOBALDH.NE.0) THEN

                          BM_ND_FL(INODE,3,1) = BM_ND_FL(INODE,3,1)
     &                                          + CONTRIB

                      END IF !IOBALDH

                  END DO ! I=1,NNUD

              END IF ! NZONEARR.NE.0

          END DO ! L=1,NUMEL


      END IF !NZARR.NE.0

********************************************************************************
***                            BOUNDARY CONDITIONS                           ***
********************************************************************************

C-------------------- Initializes pointers to mass balance vector.

      IPOSH = IPOINT + NZARR
      IPOSQ = IPOSH + NZCHP
      IPOSL = IPOSQ + NZQQP
      IPOSC = IPOSL + NZALF

C-------------------- Fluxes at prescribed head nodes.

      IF (NZCHP.NE.0 .AND. IOEQT.EQ.1) THEN

C-------------------- If not already computed for transport equation nodal
C-------------------- fluxes in nodes with prescribed head are calculated.

              CAUDAL = 0D0

              CALL COMFLOW_PRESC_HEAD
     &          (AFLU     ,AREA     ,BETAC    ,BUOYANCY ,CAUDAL
     &          ,CAUX2    ,CFLU     ,COORD    ,CREF     ,DBUOYANCY
     &          ,DENSITY  ,DENSREF  ,DFLU     ,DFLUDFLU ,DFLUDTRA
     &          ,DPARELDH ,GP_COORD ,GRADLOC  ,GRAVEL   ,GRDFF
     &          ,HAUX1    ,HAUX2    ,IBCOD    ,IDIMAFLU ,IDIMDFLU
     &          ,IODENS   ,IODIM    ,IOFLLI   ,ISOLFL   ,ISOZ
     &          ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &          ,LXARR    ,LXPAREL  ,MAXPG    ,NPAREL   ,NPPEL
     &          ,NTYPEL   ,NUMEL    ,NUMNP    ,NZTRA    ,PAREL
     &          ,POINTWEIGHT        ,THETAT)

      END IF !NZCHP.NE.0

C-------------------- Concentration sources at prescribed concentration nodes.

      IF (IODENS.EQ.1 .AND. IOCONSRC.EQ.1) THEN

          ALLOCATE(FLOWAUX(NUMNP))

          FLOWAUX = 0D0

          CALL COMFLOW_PRESC_CONC
     &        (AFLU     ,AREA     ,ARRC     ,ATRA     ,BETAC
     &        ,BUOYANCY ,CAUDAL   ,CAUX1    ,CAUX2    ,CFLU
     &        ,COECEL   ,COORD    ,CREF     ,DFLU     ,DENSITY
     &        ,DENSREF  ,DTRA     ,FLOWAUX  ,GP_COORD ,GRADLOC
     &        ,HAUX1    ,HAUX2    ,IAD_S    ,IADN_S   ,IBCOD
     &        ,IBTCO    ,IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU
     &        ,IDIMDTRA ,IODENS   ,IODIM    ,IORECATRA,ISOLFL
     &        ,ISOLTR   ,ISOZ     ,ITYPAFLU ,ITYPATRA ,ITYPCFLU
     &        ,ITYPDFLU ,ITYPDTRA ,KXX      ,LDIM     ,LMXNDL
     &        ,LNNDEL   ,LTYPE    ,LXARR    ,LXPAREL  ,MAXNB
     &        ,MAXPG    ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL
     &        ,NUMEL    ,NUMNP    ,NZARR    ,NZTRA    ,PAREL
     &        ,PARNP    ,POINTWEIGHT        ,THETAT)

      END IF !IODENS.EQ.1 .AND. IOCONSRC.EQ.1



      DO I=1,NUMNP

          IB = IBCOD(I)
          IBT = IBTCO(I)

C-------------------- Prescribed head

          IF (IB.EQ.1) THEN

              NZONECHP = IXCHP(I)
              IPOS = IPOSH + NZONECHP
              IPOSNOD = 4
              FLOW = CAUDAL(I)

              IF (IOEQT.NE.1) THEN

                  CALL CALC_FLOWDENS
     &                (BETAC    ,CAUX1    ,CREF     ,DENSREF  ,FLOW
     &                ,I        ,IODENS   ,NUMNP    ,NPPNP    ,PARNP)

              END IF !IOEQT.NE.1

              BM_ZN_FL(IPOS,1) = BM_ZN_FL(IPOS,1) + FLOW

              IF (IOBALDH.NE.0) THEN

                  BM_ND_FL(I,IPOSNOD,1) = BM_ND_FL(I,IPOSNOD,1) + FLOW

              END IF !IOBALDH.NE.0

          END IF !IB.EQ.1

C-------------------- Prescribed flow and mixed


          IF(IB.EQ.2 .OR. IB.EQ.4) THEN

              NZONEQQP = IXQQP(I)
              IPOS = IPOSQ+NZONEQQP
              IPOSNOD = 5
              FLOW = QQPC(I)

              IF (IB.EQ.4) THEN

                  FLOW = FLOW + ALFC(I)*(CHPC(I)-HAUX1(I))

              END IF
              CALL CALC_FLOWDENS
     &            (BETAC    ,CAUX1    ,CREF     ,DENSREF  ,FLOW
     &            ,I        ,IODENS   ,NUMNP    ,NPPNP    ,PARNP)

              BM_ZN_FL(IPOS,1) = BM_ZN_FL(IPOS,1) + FLOW

              IF (IOBALDH.NE.0) THEN

                  BM_ND_FL(I,IPOSNOD,1) = BM_ND_FL(I,IPOSNOD,1) + FLOW

              END IF !IOBALDH.NE.0

          END IF !IB.EQ.2 .OR. IB.EQ.4

C-------------------- Leakage

          IF (IB.EQ.3) THEN

              NZONEALF = IXALF(I)
              IPOS = IPOSL+NZONEALF
              IPOSNOD = 6

              FLOW = ALFC(I)*(CHPC(I)-HAUX1(I))

              CALL CALC_FLOWDENS
     &            (BETAC    ,CAUX1    ,CREF     ,DENSREF  ,FLOW
     &            ,I        ,IODENS   ,NUMNP    ,NPPNP    ,PARNP)



              BM_ZN_FL(IPOS,1) = BM_ZN_FL(IPOS,1) + FLOW

              IF (IOBALDH.NE.0) THEN

                  BM_ND_FL(I,IPOSNOD,1) = BM_ND_FL(I,IPOSNOD,1) + FLOW

              END IF !IOBALDH.NE.0

          END IF !IB.GE.3

C-------------------- Mass sources (only density dependent flow).

          IF (IODENS.EQ.1 .AND. IOCONSRC.EQ.1) THEN

              NZONECONC = IXCONC(I)
              IPOS = IPOSC + NZONECONC
              IPOSNOD = 8
              MASS = 0D0

              IF (IBT.EQ.1) THEN

                  MASS = FLOWAUX(I)

              ELSE IF (IBT.EQ.4) THEN

                  MASS = COECNP(I)

              ELSE IF(IBT.EQ.5) THEN

                  MASS = CLK(I)*(COECNP(I) - CAUX1(I))

              END IF !IBT.EQ.1,4,5

              BM_ZN_FL(IPOS,1) = BM_ZN_FL(IPOS,1) + MASS

              IF (IOBALDH.NE.0) THEN

                  BM_ND_FL(I,IPOSNOD,1) = BM_ND_FL(I,IPOSNOD,1) + MASS

              END IF !IOBALDH.NE.0

          END IF !IODENS.EQ.1 .AND. IOCONSRC.EQ.1


      END DO !I=1,NUMNP

********************************************************************************
***                            DIVERGENCE OF DARCY'S FLOWS                   ***
********************************************************************************

C-------------------- Divergence of Darcy's velocities are computed
C-------------------- when the nodal mass balance is required

      IF (IOBALDH.NE.0) THEN

          IF (.NOT. ALLOCATED(FLOWAUX)) THEN

              ALLOCATE(FLOWAUX(NUMNP))

          END IF !.NOT. ALLOCATED(FLOWAUX)

          FLOWAUX = 0D0

          CALL COMP_DIV_Q
     &        (AFLU     ,AREA     ,BUOYANCY ,COORD    ,DENSITY
     &        ,FLOWAUX  ,GP_COORD ,GRADLOC  ,HAUX1    ,IAD_S
     &        ,IADN_S   ,IDIMAFLU ,IODENS   ,IODIM    ,ISOZ
     &        ,ITYPAFLU ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL
     &        ,LTYPE    ,LXPAREL  ,MAXNB    ,MAXPG    ,NPAREL
     &        ,NPPEL    ,NTYPEL   ,NUMEL    ,NUMNP    ,NZTRA
     &        ,PAREL    ,POINTWEIGHT        ,THETAT)


          BM_ND_FL(:,7,1) = BM_ND_FL(:,7,1) + FLOWAUX(:)


      END IF ! IOBALDH.NE.0

********************************************************************************
***                            GLOBAL FLOW MASS BALANCE                      ***
********************************************************************************

C-------------------- Finally, computes and stores the contribution
C-------------------- of temporal flow mass balance to global f.m.b.

      IF (ISOLFL.EQ.2) THEN     ! Only sense when transient flow has been solved

          IF (IOBALGH.NE.0) THEN ! Global mass balance needed

C-------------------- Zonal

              BM_ZN_FL(:,2) = BM_ZN_FL(:,2) + BM_ZN_FL(:,1)*DELTAT

C-------------------- Nodal

              IF (IOBALDH.NE.0) THEN ! Nodal

                  DO IGROUP=1,7 ! Loop over flow parameter types

                      BM_ND_FL(:,IGROUP,2) = BM_ND_FL(:,IGROUP,2)
     &                                  + BM_ND_FL(:,IGROUP,1)*DELTAT 

                  END DO !IGROUP=1,7 ! Next type of flow parameter

              END IF ! IOBALDH.NE.0

          END IF  !IOBALGH.NE.0

      END IF ! SS/Transient

C--------------------Checks the pass through of current subr.

      IF(IFLAGS(3).EQ.1) CALL IO_SUB('BALANCE_FL',1)

      IF (ALLOCATED(FLOWAUX)) THEN

          DEALLOCATE(FLOWAUX)

      END IF !ALLOCATED(FLOWAUX)



      END SUBROUTINE BALANCE_FL
