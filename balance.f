      SUBROUTINE BALANCE
     &          (ACTH     ,AFLU     ,AREA     ,ATRA     ,BETAC
     &          ,BM_ND_FL ,BM_ND_TT ,BM_ZN_FL ,BM_ZN_TT ,BUOYANCY
     &          ,CAUDAL   ,CAUX1    ,CAUX2    ,CCALAN   ,CCALIT
     &          ,CFLU     ,COORD    ,CREF     ,DBUOYANCY,DELTAT
     &          ,DENSITY  ,DENSREF  ,DFLU     ,DFLUDFLU ,DFLUDTRA
     &          ,DPARELDH ,DTRA     ,GP_COORD ,GRADLOC  ,GRAVEL
     &          ,GRDFF    ,HCALAN   ,HCALIT   ,HAUX1    ,HAUX2
     &          ,I_REC    ,IAD_S    ,IADN_S   ,IBCOD    ,IBTCO
     &          ,IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDTRA ,IDIMDFLU
     &          ,IENTRY   ,IFLAGS   ,INDENDDT ,INTI     ,IOBALC
     &          ,IOBALDC  ,IOBALDH  ,IOBALGC  ,IOBALGH  ,IOBALH
     &          ,IOCONSRC ,IODENS   ,IODIM    ,IOEQT    ,IOFLLI
     &          ,IOSMFL   ,IOSMTP   ,IORECATRA,IOVRWC   ,ISOLEQ
     &          ,ISOLFL   ,ISOLTR   ,ISOZ     ,ITPTVAR  ,ITYPAFLU
     &          ,ITYPATRA ,ITYPCFLU ,ITYPDFLU ,ITYPDTRA ,IXPARNP
     &          ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &          ,LXPAREL  ,MAXNB    ,MAXPG    ,NFLAGS   ,NINT
     &          ,NMAXF    ,NMAXT    ,NPARALG  ,NPAREL   ,NPARNP
     &          ,NPBFL    ,NPBMX    ,NPBTP    ,NPPEL    ,NPPNP
     &          ,NTYPAR   ,NTYPEL   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &          ,NZTRA    ,PAR_DIR  ,PAREL    ,PARNP    ,POINTWEIGHT
     &          ,TABSOLUT ,TIME     ,WSPECHEAT,WATVOL)

********************************************************************************
*
* PURPOSE Computes flow and/or transport mass balance 
*
* DESCRIPTION This subroutine consists of 3 basic steps
*             1) Decision of whether or not flow and/or transport mass balance
*                (temporal and/or global) must be computed and/or written. 
*                Actually, transport mass balance can be point in time defined
*                or integrated over a time interval.
*             2) Computes flow and/or tpt. mass bal.(global mass balances are 
*                computed at last time step at particular subroutines 
*                BALANCE_FL and BALANCE_TR)
*             3) Writing of mass balance results 
*
* EXTERNAL VARIABLES: ARRAYS
*
*  ATRA                   Matrix of finite elements equations for transport     
*                         problem. Only mass flow boundary conditions included. 
*  BM_ND_FL                                                                     
*  BM_ND_TT                                                                     
*  BM_ZN_FL                                                                     
*  BM_ZN_TT                                                                     
*  CFPARNP                Array containing node coefinient of node j            
*                         corresponding to INpar index parameter zone.          
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  IXPARNP                Array containing zone number of each node j,          
*                         corresponding to INpar index parameter zone.          
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*  PARNP                  Parameter values at every node and current time for   
*                         all nodal parameters (each value is computed as the   
*                         product of up to four terms:                          
*                           nodal coeff*zonal value*time funct.*nonl. funct. )  
*  PAR_DIR                Array containing all real direct problem              
*                         parameters                                            
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  AFLU                   Matrix of finite elements equations for flow problem  
*                         No boundary conditions are included on it.            
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CAUDAL                 Input/output flow at every node.                      
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  CCAL                   Computed concentration at every node                  
*  CCALAN                 Computed concentrations in the previous time step.    
*  DELTAT                                                                       
*  DFLU                   Matrix of finite elements equations for flow          
*                         problem related to storage term.                      
*  HAUX1                  Array containing HEADS, ponderated by THETAF          
*                         time factor                                           
*  HAUX2                  Array containing diference of HEADS in two            
*                         consecutives times.                                   
*  IBCOD                  Flow boundary condition index                         
*  IBTCO                  Transport boundary condition index                    
*  IENTRY                                                                       
*  INDENDDT               Controls the coincidence of the next computing        
*                         time with the end of the current time observation     
*                         interval                                              
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOBALC                 Zonal transport mass balance at observation times     
*  IOBALDC                                                                      
*  IOBALDH                                                                      
*  IOBALGC                Global transport  mass balance (integrated in time)   
*  IOBALGH                Global flow mass balance (integrated in time)         
*  IOBALH                 Zonal flow mass balance at observation times          
*  IOCNSF                 Scheme for storage term in flow problem               
*  IORTS                  Transport regime                                      
*  IOTRS                  Flow regime
*  ISOLFL                 If 0, no flow has been solved at current time.
*                         If 1, steady flow has been solved at current time.
*                         If 2, transient flow has been solved at current time
*  ISOLTR                 If 0, no transport has been solved at current time.
*                         If 1, steady transport has been solved at current time
*                         If 2, trans. transport has been solved at current time
*  I_REC                                                                        
*  KINT                   Number of solution time increments, between           
*                         successive observation times.                         
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LMXNDL                 Maximum number of nodes per element                   
*  LNNDEL                 Number of nodes at every element                      
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NBAND                  Half Bandwith (maximum difference between the         
*                         numbers of two nodes belonging to the same element)   
*  NBAND1                 Used to dimension. It is equal to NBAND+1             
*  NBAND2                 Used to dimension. It is equal to 2*NBAND+1           
*  NFLAGS                 Maximum number of allowed flags                       
*  NINT                   Number of observation times                           
*  NMAXF                                                                        
*  NMAXT                                                                        
*  NPARALG                Maximum number of algorithm parameters, including     
*                         minimization and simulation (used for dimensioning)   
*  NPAREL                 Number of element parameters in current problem       
*  NPARNP                 Number of nodal parameters in current problem         
*  NPBTP                  Number of simultaneous transport problems             
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NPPNP                  Total number of parameters by nodes (not confuse      
*                         with NPARNP, because in this casethere is no          
*                         difference between a given parameter in steady or tr.)
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  TABSOLUT               Current absolut computation time                      
*  TIME                   Observation times.                                    
*  TIMEINT                                                                      
*  WATVOL                                                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Dummy counter variable
*  INDBALC                Transport mass balance option:
*                           - 0= Nothing is done
*                           - 1= Computation only
*                           - 2= Computation and writing
*  INDBALH                Flow mass balance option (same as INDBALC)
*  INDINIT                If non zero, zonal and nodal t.m.b. arrays have to be
*                         initialised. Conditions: first computation, after 
*                         writing tmb output or if instantaneous t.m.b. is 
*                         requested
*  J                      Dummy counter variable
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  BALANCE_FL     Computes and writes flow mass balance
*  BALANCE_TR     Computes and writes transport mass balance
*  BALANCE_WRITE  Writes the computation results
*  IO_SUB         Checks if subroutine has been entered and passed through 
*                 successfully 
*  
* HISTORY: JLF: First coding (Oct-2000)
*          AAR: First coding (Oct-2000)
*          JLF,JCR,AAR: Revision (June-2001)
*          AMS: Modification to allow succesive/simultaneous problems (7-2003)
*
********************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::I_REC    ,IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU
     &          ,IDIMDTRA ,IENTRY   ,INDENDDT ,INTI     ,IOBALC
     &          ,IOBALDC  ,IOBALDH  ,IOBALGC  ,IOBALGH  ,IOBALH
     &          ,IOCONSRC ,IODENS   ,IODIM    ,IOEQT    ,IOFLLI
     &          ,IOSMFL   ,IOSMTP   ,IORECATRA,IOVRWC   ,ISOLFL
     &          ,ISOLTR   ,ITPTVAR  ,ITYPAFLU ,ITYPATRA ,ITYPCFLU
     &          ,ITYPDFLU ,ITYPDTRA ,LMXNDL   ,MAXNB    ,MAXPG
     &          ,NFLAGS   ,NINT     ,NMAXF    ,NMAXT    ,NPARALG
     &          ,NPAREL   ,NPARNP   ,NPBFL    ,NPBMX    ,NPBTP
     &          ,NPPEL    ,NPPNP    ,NTYPAR   ,NTYPEL   ,NUMEL
     &          ,NUMNP    ,NZTRA
     

      REAL*8::BETAC    ,CREF     ,DELTAT   ,DENSREF  ,TABSOLUT
     &       ,WSPECHEAT


      INTEGER*4::IAD_S(MAXNB,NUMNP)         ,IADN_S(NUMNP)
     &          ,IBCOD(NUMNP,NPBFL)         ,IBTCO(NUMNP,NPBTP)
     &          ,IFLAGS(NFLAGS)             ,ISOLEQ(NINT,4)
     &          ,ISOZ(NZTRA)                ,IXPARNP(NUMNP,NPARNP,NPBMX)
     &          ,KXX(LMXNDL,NUMEL)          ,LDIM(NUMEL)
     &          ,LNNDEL(NUMEL)              ,LTYPE(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL,NPBMX),NZONE_PAR(NTYPAR)



      REAL*8::ACTH(NUMEL)                  ,AFLU(NUMEL,IDIMAFLU,NPBFL)
     &       ,AREA(NUMEL)                  ,ATRA(NUMEL,IDIMATRA,NPBTP)
     &       ,BM_ND_FL(NUMNP,8,2)          ,BM_ND_TT(NUMNP,12,2)
     &       ,BM_ZN_FL(NMAXF,2,NPBFL)      ,BM_ZN_TT(NMAXT,2,NPBTP)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL) ,CAUDAL(NUMNP)
     &       ,CAUX2(NUMNP,NPBTP)           ,CAUX1(NUMNP,NPBTP)
     &       ,CCALAN(NUMNP,NPBTP)          ,CCALIT(NUMNP,NPBTP)
     &       ,CFLU(NUMNP,IDIMDFLU,NPBTP)   ,COORD(NUMNP,3)
     &       ,DBUOYANCY(IODIM,LMXNDL*LMXNDL,NUMEL)
     &       ,DENSITY(NUMEL)               ,DFLU(NUMNP,IDIMDFLU,NPBFL)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL),DPARELDH(NPPEL,NUMEL)
     &       ,DTRA(NUMEL,IDIMDTRA,NPBTP)   ,GP_COORD(6,8,IODIM)
     &       ,GRADLOC(IODIM,LMXNDL,MAXPG)  ,GRAVEL(NUMEL,3)
     &       ,GRDFF(IODIM,LMXNDL,NUMEL)    ,HCALAN(NUMNP,NPBFL)
     &       ,HCALIT(NUMNP,NPBFL)          ,HAUX1(NUMNP,NPBFL)
     &       ,HAUX2(NUMNP,NPBFL)           ,PAR_DIR(NPARALG)
     &       ,PAREL(NUMEL,NPPEL,NPBMX)     ,PARNP(NUMNP,NPPNP,NPBMX)
     &       ,POINTWEIGHT(MAXPG,NTYPEL)    ,TIME(NINT)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)

C------------------------- Internal     

      INTEGER*4::IPROB,INDBALC,INDBALH,INDINIT,IPBFL,IPBTP,IUBALC,IUBALH

C______________________________ Checks the succesfull entry to current subr.

      IF(IFLAGS(3).EQ.1) CALL IO_SUB('BALANCE',0)

C______________________________ Step 1. Control of flow mass balance options

      INDBALH=0

      IF (IOBALGH.NE.0) THEN              ! Global flow mass balance is required
 
        INDBALH=1                                            ! Computation only.
   
      END IF

      IF (IOBALH.LT.0) THEN                       ! F. m. b. at simulation times

        IF (MOD(IENTRY,IOBALH).EQ.0) INDBALH=2       ! Computation and writting 
                                                     ! PROVISIONAL (IENTRY)

      ELSE IF (IOBALH.GT.0.AND.INDENDDT.EQ.1) THEN

                                                  ! At first or last obs. time,
                                                  ! or if corresponds, f.m.b. is
                                                  ! computed and written

        IF (MOD(INTI+1,IOBALH).EQ.0.OR.INTI.EQ.0.OR.INTI.EQ.NINT)
     ;    INDBALH=2           

      END IF

C______________________________Step 1. Control of transport mass balance options


      INDBALC=0
                                     ! Global transport mass balance is required
      IF (IOBALGC.NE.0) THEN 

        INDBALC=1                                             ! Only computation

      END IF !IOBALGC.NE.0

      IF (IOBALC.LT.0) THEN                                   ! Simulation times

        IF (MOD(IENTRY,IOBALC).EQ.0) INDBALC=2     
 
      ELSE IF (IOBALC.GT.0.AND.INDENDDT.EQ.1) THEN                  ! Obs. times

                                                  ! At first or last obs. time,
                                                  ! or if corresponds, f.m.b. is
                                                  ! computed and written

        IF (MOD(INTI+1,IOBALC).EQ.0.OR.INTI.EQ.0.OR.INTI.EQ.NINT)
     ;    INDBALC=2           

      END IF

C______________________________ Conditions to compute  f.m.b. at current, as 
C______________________________ detailed nodal tpt.  m. b. is requested

      IF (INDBALH.EQ.0.AND.IOBALDC.NE.0.AND.INDBALC.NE.0) THEN

        IF (ISOLFL.NE.0) INDBALH=1

      END IF

C______________________________ Definition of INDINIT (controls initializing of
C______________________________ tmb arrays)

      IF (IOBALGC.LE.1.OR.INTI.EQ.0) INDINIT=1

C______________________________ Step 2. Flow mass balance computations


      IF (ISOLFL.EQ.0) INDBALH=0        ! Flow equation is not computed

      IF (INDBALH.NE.0) THEN            ! Temporal flow mass balance is required

         DO IPROB=1,MAX(1,IOSMFL*NPBFL)

C------------------------- IPBFL is the number of the current flow problem

            IF (IOSMFL.EQ.1) THEN
               IPBFL=IPROB
            ELSE
               IPBFL=MAX(1,ISOLEQ(MAX(1,INTI),3))
            ENDIF

C-------------------- NZONE_PAR(1)               !NZTRA
C-------------------- NZONE_PAR(6)               !NZALF        
C-------------------- NZONE_PAR(3)               !NZARR
C-------------------- NZONE_PAR(4)               !NZCHP
C-------------------- NZONE_PAR(5)               !NZQQP
C-------------------- NZONE_PAR(2)               !NZSTG
C-------------------- PAR_DIR(30)                !THETAT
C-------------------- PARNP(1,3,IPROB)           !ALFC
C-------------------- PAREL(1,8,IPROB)           !ARRC
C-------------------- PARNP(1,1,IPROB)           !CHPC      
C-------------------- IXPARNP(1,4+ISOLFL,IPBFL)  !IXALF
C-------------------- IXPARNP(1,ISOLFL,IPBFL)    !IXCHP
C-------------------- IXPARNP(1,8,IPROB)         !IXCONC
C-------------------- IXPARNP(1,2+ISOLFL,IPBFL)  !IXQQP
C-------------------- LXPAREL(1,2+ISOLFL,IPBFL)  !LXARR
C-------------------- LXPAREL(1,2,IPBFL)         !LXSTG
C-------------------- PARNP(1,2,IPROB)           !QQPC
C-------------------- PAREL(1,7,IPROB)           !STGC
C-------------------- PARNP(1,6,IPROB)           !CLKCF
C-------------------- PAREL(1,15,IPROB)          !COECEL
C-------------------- PARNP(1,4,IPROB)           !COECNP

            CALL BALANCE_FL
     &          (AFLU(1,1,IPROB)    ,PARNP(1,3,IPROB)   ,AREA
     &          ,PAREL(1,8,IPROB)   ,ATRA     ,BETAC    ,BM_ND_FL
     &          ,BM_ZN_FL(1,1,IPROB),BUOYANCY ,CAUDAL   ,CAUX1  ,CAUX2
     &          ,CFLU     ,PARNP(1,1,IPROB)
     &          ,PARNP(1,6,IPROB)   ,PAREL(1,15,IPROB)
     &          ,PARNP(1,4,IPROB)   ,COORD    ,CREF     ,DBUOYANCY
     &          ,DELTAT   ,DENSITY  ,DENSREF  ,DFLU     ,DFLUDFLU
     &          ,DFLUDTRA ,DPARELDH ,DTRA     ,GP_COORD ,GRADLOC
     &          ,GRAVEL   ,GRDFF    ,HCALAN(1,IPROB)    ,HCALIT(1,IPROB)
     &          ,HAUX1(1,IPROB)     ,HAUX2(1,IPROB)    ,IAD_S    ,IADN_S
     &          ,IBCOD(1,IPBFL),IBTCO    ,IDIMAFLU ,IDIMATRA ,IDIMCFLU
     &          ,IDIMDFLU ,IDIMDTRA ,IFLAGS   ,INTI     ,IOBALDH
     &          ,IOBALGH  ,IOCONSRC ,IODENS   ,IODIM    ,IOEQT
     &          ,IOFLLI   ,IORECATRA,IOVRWC   ,ISOLFL   ,ISOLTR ,ISOZ
     &          ,ITYPAFLU ,ITYPATRA ,ITYPCFLU ,ITYPDFLU ,ITYPDTRA
     &          ,IXPARNP(1,4+ISOLFL,IPBFL)    ,IXPARNP(1,ISOLFL,IPBFL)
     &          ,IXPARNP(1,8,IPBFL) ,IXPARNP(1,2+ISOLFL,IPBFL)
     &          ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &          ,LXPAREL(1,2+ISOLFL,IPBFL)    ,LXPAREL
     &          ,LXPAREL(1,2,IPBFL) ,MAXNB    ,MAXPG    ,NFLAGS
     &          ,NMAXF    ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL
     &          ,NUMEL    ,NUMNP    ,NZONE_PAR(6)       ,NZONE_PAR(3)
     &          ,NZONE_PAR(4)       ,NZONE_PAR(5)       ,NZONE_PAR(2)
     &          ,NZONE_PAR(1)   ,PAREL(1,1,IPROB)   ,PARNP(1,1,IPROB)
     &          ,POINTWEIGHT        ,PARNP(1,2,IPROB)
     &          ,PAREL(1,7,IPROB)   ,PAR_DIR(30)     ,WATVOL)

 
         END DO !IPROB=1,MAX(1,IOSMFL*NPBFL)
      END IF !INDBALH.NE.0

C______________________________ Step 2. Transport mass balance computations

      IF (ISOLTR.EQ.0) INDBALC=0          ! Transport equation is not computed

      IF (INDBALC.NE.0) THEN              !  Transport mass balance is requested

         DO IPROB=1,MAX(1,IOSMTP*NPBTP)

C------------------------- IPBTP is the number of the current transport problem

            IF (IOSMTP.EQ.1) THEN
               IPBTP=IPROB
               IPBFL=1
            ELSE
               IPBTP=MAX(1,ISOLEQ(MAX(1,INTI),4))
               IPBFL=MAX(1,ISOLEQ(MAX(1,INTI),3))
            ENDIF

C-------------------- NZONE_PAR(3)              --> NZARR
C-------------------- NZONE_PAR(13)             --> NZCOE
C-------------------- NZONE_PAR(11)             --> NZFOD
C-------------------- NZONE_PAR(10)             --> NZPOR
C-------------------- NZONE_PAR(15)             --> NZZOR
C-------------------- NZONE_PAR(18)             --> NZCLK
C-------------------- PAR_DIR(30)               --> THETAT
C-------------------- PAREL(1,8,IPROB)          --> ARRC
C-------------------- PARNP(1,6,IPROB)          --> CLKCF
C-------------------- PAREL(1,15,IPROB)         --> COECEL
C-------------------- PARNP(1,4,IPROB)          --> COECNP
C-------------------- PAREL(1,14,IPROB)         --> CRDC
C-------------------- PAREL(1,13,IPROB)         --> FODC
C-------------------- IXPARNP(1,6+ISOLTR,IPBTP) --> IXCON
C-------------------- LXPAREL(1,2+ISOLTR,IPBTP) --> LXARR
C-------------------- LXPAREL(1,10,IPBTP)       --> LXCOE
C-------------------- LXPAREL(1,8,IPBTP)        --> LXFOD
C-------------------- LXPAREL(1,7,IPBTP)        --> LXPOR
C-------------------- LXPAREL(1,11,IPBTP)       --> LXZOR
C-------------------- PAREL(1,16,IPROB)         --> ZORC   
C-------------------- BM_ND_FL(1,6,1)           --> DIVQ

            CALL BALANCE_TR
     &          (ACTH     ,AFLU     ,AREA     ,PAREL(1,8,IPROB)
     &          ,ATRA(1,1,IPROB)    ,BETAC    ,BM_ND_TT 
     ;          ,BM_ZN_TT(1,1,IPROB)
     &          ,BUOYANCY ,CAUDAL   ,CAUX1(1,IPROB)     ,CAUX2(1,IPROB)
     &          ,CCALAN(1,IPROB)    ,CCALIT(1,IPROB )   ,CFLU
     &          ,PARNP(1,6,IPROB)   ,PAREL(1,15,IPROB) ,PARNP(1,4,IPROB)
     &          ,COORD    ,PAREL(1,14,IPROB)  ,CREF     ,DELTAT
     &          ,DENSITY  ,DENSREF  ,DFLU     ,BM_ND_FL(1,7,1)
     &          ,DTRA     ,PAREL(1,13,IPROB)  ,GP_COORD ,GRADLOC
     &          ,HAUX1    ,HAUX2    ,IAD_S    ,IADN_S   ,IBCOD(1,IPBFL)
     &          ,IBTCO(1,IPBTP) ,IDIMAFLU ,IDIMATRA ,IDIMDFLU ,IDIMDFLU
     &          ,IDIMDTRA ,IFLAGS   ,INDINIT  ,INTI     ,IOBALDC
     &          ,IOBALGC  ,IODENS   ,IODIM    ,IORECATRA,IOVRWC
     &          ,ISOLFL   ,ISOLTR   ,ISOZ     ,ITPTVAR  ,ITYPAFLU
     &          ,ITYPATRA ,ITYPCFLU ,ITYPDFLU ,ITYPDTRA
     &          ,IXPARNP(1,10,IPBTP),IXPARNP(1,6+ISOLTR,IPBTP)
     &          ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &          ,LXPAREL(1,2+ISOLTR,IPBTP)    ,LXPAREL(1,10,IPBTP)
     &          ,LXPAREL(1,8,IPBTP) ,LXPAREL  ,LXPAREL(1,7,IPBTP)
     &          ,LXPAREL(1,11,IPBTP),MAXNB    ,MAXPG    ,NFLAGS
     &          ,NMAXT    ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL
     &          ,NUMEL    ,NUMNP    ,NZONE_PAR(3)       ,NZONE_PAR(18)
     &          ,NZONE_PAR(13)      ,NZONE_PAR(11)      ,NZONE_PAR(10)
     &          ,NZTRA    ,NZONE_PAR(15)      ,PAREL(1,1,IPROB)
     ;          ,PARNP(1,1,IPROB)
     &          ,POINTWEIGHT        ,PAR_DIR(30)        ,WATVOL
     &          ,WSPECHEAT,PAREL(1,16,IPROB))

         END DO !IPROB=1,MAX(1,IOSMTP*NPBTP)
      END IF !INDBALC.NE.0

C------------------------- Writes flow mass balance

      IF (INDBALH.GE.1) THEN

	    IUBALH = 35

         DO IPROB=1,MAX(1,IOSMFL*NPBFL)

C------------------------- IPBFL is the number of the current flow problem

            IF (IOSMFL.EQ.1) THEN
               IPBFL=IPROB
               IPBTP=1        ! only 1 tpt problem for simultaneous flow probs
            ELSE
               IPBFL=MAX(1,ISOLEQ(MAX(1,INTI),3))
               IPBTP=MAX(1,ISOLEQ(MAX(1,INTI),4))
            ENDIF

C------------------------- Temporal mass balance

            IF (INDBALH.EQ.2)
     &         CALL BALANCE_WRITE_FL
     &             (BM_ND_FL(1,1,1)    ,BM_ZN_FL(1,1,IPROB),IPBFL
     &             ,I_REC    ,ISOLFL-1 ,INTI     ,IOBALH   ,IOBALDH
     &             ,IOCONSRC ,IODENS   ,IUBALH   ,NINT     ,8
     &             ,NMAXF    ,NUMNP    ,NZONE_PAR(6)       ,NZONE_PAR(3)
     &             ,NZONE_PAR(4)       ,NZONE_PAR(13)      ,NZONE_PAR(5)
     &             ,NZONE_PAR(2)       ,TABSOLUT ,TIME) 

C------------------------- If current time step is the last simulation time
C------------------------- of the current flow problem, global mass balance
C------------------------- is written (only in transient flow ISOLFL=2)

         IF (IOBALGH.NE.0 .AND. INDENDDT.EQ.1 .AND. ISOLFL.EQ.2 .AND.
     ;        ISOLEQ(INTI,3).NE.ISOLEQ(INTI+1,3) ) THEN

            CALL BALANCE_WRITE_FL
     &          (BM_ND_FL(1,1,2)    ,BM_ZN_FL(1,2,IPROB),IPBFL
     &          ,I_REC    ,2        ,INTI     ,IOBALH   ,IOBALDH
     &          ,IOCONSRC ,IODENS   ,IUBALH   ,NINT     ,8
     &          ,NMAXF    ,NUMNP    ,NZONE_PAR(6)       ,NZONE_PAR(3)
     &          ,NZONE_PAR(4)       ,NZONE_PAR(13)      ,NZONE_PAR(5)
     &          ,NZONE_PAR(2)       ,TABSOLUT ,TIME) 

C------------------------- Global balance arrays must be initialized to zero for
C------------------------- the next problem

              IF (IOBALH.NE.0 .OR. IOBALGH.NE.0)
     ;           CALL ZERO_ARRAY (BM_ZN_FL(1,2,IPROB),NMAXF)
              IF (IOBALDH.NE.0) CALL ZERO_ARRAY(BM_ND_FL(1,1,2),6*NUMNP)

           ENDIF
         ENDDO   ! IPROB=1,MAX(1,IOSMFL*NPBFL)
      ENDIF      ! INDBALH.GE.1

C------------------------- Writes transport mass balance

      IF (INDBALC.GE.1) THEN

	    IUBALC = 36

         DO IPROB=1,MAX(1,IOSMTP*NPBTP)

C------------------------- IPBTP is the number of the current transport problem

            IF (IOSMTP.EQ.1) THEN
               IPBTP=IPROB
               IPBFL=1
            ELSE
               IPBTP=MAX(1,ISOLEQ(MAX(1,INTI),4))
               IPBFL=MAX(1,ISOLEQ(MAX(1,INTI),3))
            ENDIF

            IF (INDBALC.EQ.2) THEN

C------------------------- Temporal mass balance

               CALL BALANCE_WRITE_TR
     &             (BM_ND_TT(1,1,1)    ,BM_ZN_TT(1,1,IPROB),IPBTP
     &             ,I_REC    ,IBTCO(1,IPBTP)     ,ISOLTR-1 ,INTI
     &             ,IOBALC   ,IOBALDC  ,IUBALC   ,NINT     ,NMAXT
     &             ,NUMNP    ,NZONE_PAR(18)      ,NZONE_PAR(13)
     &             ,NZONE_PAR(11)      ,NZONE_PAR(16)
     &             ,NZONE_PAR(10)      ,NZONE_PAR(17)      ,TABSOLUT
     &             ,TIME)       
        
               INDINIT=1 
 
            ELSE

               INDINIT=0

            ENDIF

C------------------------- If current time step is the last simulation time
C------------------------- of the current transp. problem, global mass balance
C------------------------- is written (only in transient transp. ISOLTR=2)

            IF (IOBALGC.EQ.1.OR.IOBALGC.EQ.3) THEN
               IF (INDENDDT.EQ.1.AND .ISOLTR.EQ.2 .AND.
     ;                        ISOLEQ(INTI,4).NE.ISOLEQ(INTI+1,4) ) THEN


	         CALL BALANCE_WRITE_TR
     &             (BM_ND_TT(1,1,2)    ,BM_ZN_TT(1,2,IPROB),IPBTP
     &             ,I_REC    ,IBTCO(1,IPBTP)     ,2        ,INTI
     &             ,IOBALC   ,IOBALDC  ,IUBALC   ,NINT     ,NMAXT
     &             ,NUMNP    ,NZONE_PAR(18)      ,NZONE_PAR(13)
     &             ,NZONE_PAR(11)      ,NZONE_PAR(16)
     &             ,NZONE_PAR(10)      ,NZONE_PAR(17)      ,TABSOLUT
     &             ,TIME)       

C------------------------- Global balance arrays must be initialized to zero for
C------------------------- the next problem

                  IF (IOBALC.NE.0 .OR. IOBALGC.NE.0)
     ;               CALL ZERO_ARRAY (BM_ZN_TT(1,2,IPROB),NMAXT)
                  IF (IOBALDC.NE.0) 
     ;               CALL ZERO_ARRAY(BM_ND_TT(1,1,2),4*NUMNP)

               ENDIF
            ENDIF


         ENDDO   ! IPROB=1,MAX(1,IOSMTP*NPBTP)
      ENDIF      ! INDBALC.GE.1
      IF(IFLAGS(3).EQ.1) CALL IO_SUB('BALANCE',1)
      RETURN
      END
