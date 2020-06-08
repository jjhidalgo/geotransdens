      SUBROUTINE COMP_PARAM_FLOW
     &          (CAUX1    ,COORD    ,DTIMEF   ,IDIMFNT
     &          ,INDSSTR  ,INTI     ,ITPTVAR  ,IODENS   ,IODIM
     &          ,IOFLLI   ,IOFLSAT  ,IOFMLF   ,ISOT     ,LMXNDL
     &          ,MAINF                        ,NFLAGS   ,NFNL
     &          ,NINT     ,NPARALG  ,NPAREL   ,NPARNP   ,NPPEL
     &          ,NPPNP    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZPAR
     &          ,NZTRA    ,VISCREF            ,CFPAREL  ,CFPARNP
     &          ,DNODALRH ,DPARELDC ,DPARELDH ,FNT      ,GRAVEL
     &          ,HBASE    ,HCALAN   ,HCALIT   ,IBCOD
     &                    ,IFLAGS   ,INORPAR  ,IPAR_DIR ,ISOZ
     &          ,IXPARNP  ,KXX      ,LNNDEL   ,LXPAREL  ,NFNLPAR
     &          ,NFNLTIP  ,NFNLPRG  ,NFTPAR   ,NZONE_PAR,PARACD
     &          ,PARC     ,PAR_DIR  ,PAREL    ,PARNP    ,TINC
     &          ,TINTERVOBS,VAR_REF)
   

********************************************************************************
*
* PURPOSE Computes the values of flow parameters 
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  CFPARNP                Array containing node coefinient of node j            
*                         corresponding to INpar index parameter zone.          
*  COORD                  Nodal coordinates.
*  DERSTGH                Array contanining the derivative of storativity with 
*                         respect to head. Used to assemble matrix DERADFLU
*  DNODALRH               1st column: Derivative of presc. head w.r.t. head 
*  (NUMNP,4)              2nd column: Derivative of Qb w.r.t. head 
*                         3rd column: Derivative of leakage coeff. w.r.t. head 
*                         4th column: Value of Qb. This array is used to 
*                         assemble matrix DERBFLU and AFLU and array BFLU
*  DTRH                   Array contanining the derivative of transmissivity 
*                         wrt. head level. Will be used to assemble matrix
*                         DERADFLU                          
*  FNT                    Array containing time functions values                
*  HBASE                  Aquifer bottom. Onlu used for non linear problems
*  HCALAN                 Head level at previous time                           
*  HCALIT                 Computed heads in last iteration                      
*  IBCOD                  Flow boundary condition index                         
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IPAR_DIR               Array containing all integer direct problem           
*                         parameters                                            
*  ISOZ                   Anisotropy of every transmissivity zone               
*  IXPARNP                Array containing zone number of each node j,          
*                         corresponding to INpar index parameter zone.          
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  NFNLPAR                Vector containing non-linear function order           
*                         afecting every parameter at each zone.                
*  NFNLPRG                Generic parameter zone number for every nonlinear     
*                         function                                              
*  NFNLTIP                Type of non-linear function                           
*  NFTPAR                 Vector containing time function number at every       
*                         parameter zone                                        
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PARACD                 Agreement parameters                                  
*  PARC                   Vector containing calculated values for all           
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
*  DERVEC                 Array contaning derivatives wrt. head level 
*  XPARAM                 Auxiliar array for non linear calculations          
*
* EXTERNAL VARIABLES: SCALARS
*
*  DTIMEF                 Current time for the computation of flow time         
*                         functions (counted since the beginning of the         
*                         current observation interval, not since the beginning 
*                         of the problem) divided by the length of the          
*                         observation interval (interval between two            
*                         consecutive observation times). Used only to make a   
*                         linear interpolation of time function values.         
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  INDSSTR                Current problem state. If 1, transient state,         
*                         if 0, steady-state                                    
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOFLSAT                Indicates the possibility that one part of the domain 
*                         reaches unsaturated state.                            
*  IOFMLF                 Flow Formulation number                               
*  ISOT                   Maximum hydraulic conductivity tensor anisotropy      
*                         degree in the problem.                                
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFLAGS                 Maximum number of allowed flags                       
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NINT                   Number of observation times                           
*  NPARALG                Maximum number of algorithm parameters, including     
*                         minimization and simulation (used for dimensioning)   
*  NPAREL                 Number of element parameters in current problem       
*  NPARNP                 Number of nodal parameters in current problem         
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
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  NZTRA                  Number of transmissivity zones                        
*
* INTERNAL VARIABLES: SCALARS
*
*  ALFC                   Computed leakage zonal parameter                      
*  ARRC                   Computed areal recharge zonal parameter.              
*  CFALF                  Leakage nodal coefficient.                            
*  CFARR                  Steady-state recharge element coefficient.            
*  CFCHP                  Steady state prescribed head nodal coefficient.       
*  CFPOR                  Porosity element coefficient.                         
*  CFQQP                  Steady state prescribed flow nodal coefficient.       
*  CFSTG                  Storage element coefficient.                          
*  CFTRA                  Transmissivity element coefficient.                   
*  CHPC                   Computed prescribed head zonal parameter              
*  DERSC                  Derivative of a parameter wrt. head level
*  DYDXOLD                Dummy ???             
*  EPSFLU                 Time weighting parameter for nonlinear flow problems  
*  FT1                    Time function evaluated at time INTI
*  FT2                    Time function evaluated at time INTI+1
*  IB                     Boundary condition of a given node                
*  IBLOCK                 Block at which a given element belongs
*  INDEX                  Pointer to PARC, NFTPAR,NFNLPAR                
*  INPORC                 Location of first porosity zone in array              
*                         variables PARC, PARM, STPAR ... minus 1               
*  INPRGC                 Location of first generic parameter zone in array     
*                         variables PARC, PARM, STPAR ... minus 1               
*  IPIPO                  A specific pilot point parameterising a given block
*  IPRESC                 Prescribed head zone
*  IS                     Dummy counter of degrees of anisotropy
*  IST                    Anisotropy degree of a given transmissivity zone
*  J                      Dummy counter      
*  L                      Dummy counter                                         
*  N                      Dummy counter                                         
*  NFNLALF                Leakage nonlinear function number                     
*  NFNLARR                Recharge nonlinear function number                    
*  NFNLCHP                Presc. head nonlinear function number                 
*  NFNLQQP                Presc. flow nonlinear function number                 
*  NFNLSTG                Storage nonlinear function number                     
*  NFNLTRA                Transmissivity nonlinear function number              
*  NFTALF                 Leakage time function number                          
*  NFTARR                 Recharge time function number                         
*  NFTCHP                 Presc. head time function number                      
*  NFTIP                  Non linear function type
*  NFTPOR                 Porosity time function number                         
*  NFTQQP                 Presc. flow time function number                      
*  NFTSTG                 Storage time function number                          
*  NFTTRA                 Transmissivity time function number                   
*  NNUD                   Number of nodes of the current element                
*  NUMPIPO                Number of pilot points parameterising a given block
*  NZONEALF               Leakage zone to which a particular node belongs to   
*  NZONEARR               Recharge zone to which a particular element belongs to
*  NZONEPOR               Porosity zone to which a particular element belongs to
*  NZONESTO               Storavity zone to which a part. element belongs to
*  NZONTRA                Transmissivity zone to which a part. element belongs to
*  NZPRG                  Total number of generic parameter zones               
*  PIPO                   Number of zones of generic parameters
*  PORC                   Computed porosity zonal values                        
*  QQPC                   Computed presc. flow zonal values                     
*  STGC                   Computed storage zonal values                         
*  STOR                   Zonal parameter                 
*  SUM                    Dummy counter                                       
*  TEMPCOEF               Contribution of time function
*  TRAC                   Computed Transmissivity zonal values                  
*  YOLD                   Dummy variable                         
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  PARAM_VALUE            Computes the total value of a given parameter, either
*                         nodal or elemental
*
* HISTORY: AAR    First coding   Feb-2000
*
********************************************************************************

C-------------------------  Step 0: Declaration of variables

      IMPLICIT NONE

C------------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4::IDIMFNT  ,INDSSTR  ,INTI     ,IODENS   ,IODIM
     &          ,IOFLLI   ,IOFLSAT  ,IOFMLF   ,ISOT
     &          ,ITPTVAR  ,LMXNDL   ,MAINF
     &          ,NFLAGS   ,NFNL     ,NINT     ,NPARALG  ,NPAREL
     &          ,NPARNP   ,NPPEL    ,NPPNP    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZPAR    ,NZTRA

      REAL*8::DERVISCDTRA,DTIMEF,TINC,TINTERVOBS,VAR_REF,VISC,VISCREF
      

C------------------------- EXTERNAL VARIABLES: ARRAYS

      INTEGER*4  IBCOD(NUMNP)           ,IFLAGS(NFLAGS)
     &          ,INORPAR(NTYPAR)        ,IPAR_DIR(NPARALG)
     &          ,ISOZ(NZTRA)            ,IXPARNP(NUMNP,NPARNP)
     &          ,KXX(LMXNDL,NUMEL)      ,LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL)  ,NFNLPAR(NZPAR)
     &          ,NFNLPRG(8,NFNL)        ,NFNLTIP(NFNL)
     &          ,NFTPAR(NZPAR)          ,NZONE_PAR(NTYPAR)


      REAL*8  DNODALRH(NUMNP,4)
     &       ,CAUX1(NUMNP)            ,CFPAREL(NUMEL,NPAREL)
     &       ,CFPARNP(NUMNP,NPARNP)   ,COORD(NUMNP,3)
     &       ,DPARELDH(NPPEL,NUMEL)   ,DPARELDC(NPPEL,NUMEL)
     &       ,FNT(IDIMFNT,NINT)       ,GRAVEL(NUMEL)
     &       ,HBASE(NUMEL)            ,HCALAN(NUMNP)
     &       ,HCALIT(NUMNP)           ,PAR_DIR(NPARALG)
     &       ,PARACD(3,NFNL)          ,PARC(NZPAR)
     &       ,PAREL(NUMEL,NPPEL)      ,PARNP(NUMNP,NPPNP)

C------------------------- INTERNAL VARIABLES: SCALARS

      INTEGER*4::I        ,IB       ,INDEX    ,INODE    ,INPORC
     &          ,INPRGC   ,IPRESC   ,IS       ,IST      ,J
     &          ,JNODE    ,L        ,NFNLALF  ,NFNLARR  ,NFNLCHP
     &          ,NFNLQQP  ,NFNLSTG  ,NFNLTRA  ,NFTALF   ,NFTARR
     &          ,NFTCHP   ,NFTIP    ,NFTPOR   ,NFTQQP   ,NFTSTG
     &          ,NFTTRA   ,NNUD     ,NZONEALF,NZONEARR  ,NZONEPOR
     &          ,NZONESTO ,NZONTRA  ,NZPRG
      
      REAL*8::ALFC     ,ARRC     ,CFALF    ,CFARR    ,CFCHP    ,CFPOR
     &       ,CFQQP    ,CFSTG    ,CFTRA    ,CHPC     ,CMED     ,DERSCC
     &       ,DERSCH   ,DERVIS   ,DTIMAUX  ,EPSFLU   ,FT1      ,FT2
     &       ,GRAVMOD  ,PORC     ,QQPC     ,STGC     ,STOR     ,TEMPCOEF
     &      ,TERM     ,TRAC     ,VISCO    ,ZAVG

C-------------------------  INTERNAL VARIABLES: ARRAYS
      
      REAL*8 XPARAM(8)
   


C------------------------- Step 1: Identifies some useful variables

      NZPRG = NZONE_PAR(14)
      EPSFLU = PAR_DIR(31)
      INPRGC = INORPAR(19) + 1

C-------------------------  Step 2: Loop over mesh elements

      DO L=1,NUMEL

          NNUD = LNNDEL(L)

C-------------------------  Step 2.1: Calculates element transmissivity

              NZONTRA = LXPAREL(L,1)
              CFTRA = CFPAREL(L,1)  
              IST = ISOZ(NZONTRA)

C-------------------------  If transport state variable is temperature
C-------------------------  and density and viscosity are variable the
C-------------------------  mean temperature is need to compute
C-------------------------  viscosity.
              IF (IODENS.EQ.1 .AND. ITPTVAR.EQ.1) THEN

                  CMED = 0D0

                  DO I=1,NNUD

                      INODE = KXX(I,L)
                      CMED = CMED + CAUX1(INODE)

                  END DO !I=1,NNUD

                  CMED = CMED/NNUD

              END IF !IODENS.EQ.1 .AND. ITPTVAR.EQ.1


              DO IS=1,IST

                  INDEX = INORPAR(IS) + NZONTRA
                  TRAC = PARC(INDEX)
                  NFTTRA = NFTPAR(INDEX)

                  IF (IOFLLI.NE.0.AND.INDSSTR.NE.0) THEN

                      NFNLTRA = NFNLPAR(INDEX)
                      NFTIP = 0
                      IF (NFNLTRA.NE.0) NFTIP = NFNLTIP(NFNLTRA)
                      IF (NFTIP.EQ.1.OR.NFTIP.EQ.2) XPARAM(1) = HBASE(L)

C------------------------- Elevation average for Van Genuchten retention.
C------------------------- Only in IS=1 no to recompute it IST-times.
	                IF (NFTIP.EQ.4 .AND. IS.EQ.1) THEN

                          ZAVG = 0
	                    GRAVMOD = 0D0
	                    DO I=1,3

                              GRAVMOD = GRAVMOD + GRAVEL(I)*GRAVEL(I)

                              DO J=1,NNUD
	                          JNODE = KXX(J,L)
	                          ZAVG = ZAVG + GRAVEL(I)*COORD(JNODE,I)
	                      END DO


	                    END DO !I=1,3

	                    ZAVG = ZAVG/(NNUD*DSQRT(GRAVMOD))

	                    XPARAM(1) = ZAVG

	                END IF !NFNLTIP(NFNLTRA).EQ.4

                  ELSE

                      NFNLTRA = 0
                      NFTIP = 0

                  END IF !IOFLLI.NE.0.AND.INDSSTR.NE.0

                  CALL PARAM_VALUE
     &                (CFTRA    ,DERSCC   ,DERSCH  ,DTIMEF   ,EPSFLU
     &                ,IDIMFNT  ,INDSSTR  ,INTI    ,IOFMLF   ,L
     &                ,LMXNDL   ,NFLAGS   ,NFNL    ,NFNLTRA  ,NFTIP
     &                ,NFTTRA   ,NINT     ,NNUD    ,NPARALG  ,NUMEL
     &                ,NUMNP    ,NZPRG    ,TRAC    ,PAREL(L,IS)
     &                ,PARC(INPRGC)       ,FNT     ,HCALAN   ,HCALIT
     &                ,IFLAGS   ,IPAR_DIR ,KXX     ,NFNLPRG  ,PARACD
     &                ,XPARAM)  


C------------------------- Viscosity variation effect on transmissivity.

                  IF (IODENS.EQ.1 .AND. ITPTVAR.EQ.1) THEN

                      VISCO = VISC(ITPTVAR,CMED,VAR_REF)
                      TERM = VISCREF/VISCO
                      PAREL(L,IS) = PAREL(L,IS)*TERM

                  END IF !IODENS.EQ.1 .AND. ITPTVAR.EQ.1

C------------------------- Assigns derivatives

                  IF (IOFLLI.NE.0) THEN

                      DPARELDH(IS,L) = DERSCH

                  END IF !IOFLLI.NE.0

                  IF (IODENS.EQ.1) THEN

                      DPARELDC(IS,L) = DERSCC

                      IF (ITPTVAR.EQ.1) THEN

                          DPARELDH(IS,L)=DERSCH*TERM

                          DERVIS = DERVISCDTRA(ITPTVAR,CMED,VAR_REF)

                          DPARELDC(IS,L) = (DPARELDC(IS,L)
     &                                       - DERVIS/VISCO)*TERM
                      
                      END IF !ITPTVAR.EQ.1

                  END IF !IODENS.EQ.1

              END DO !IS=1,IST

C------------------------- Fill the repeated components

              IF (IST.LT.IODIM) THEN

                  DO I=IST+1,IODIM

                      PAREL(L,I) = PAREL(L,1)

                  END DO !I=IST+1,IODIM

              END IF !IST.LT.IODIM
                
C------------------------- Step 2.2: Calculates storavity

          IF (INDSSTR.NE.0.AND.NZONE_PAR(2).NE.0) THEN

              NZONESTO = LXPAREL(L,2)
              INDEX = INORPAR(7) + NZONESTO
              STGC = PARC(INDEX)
              CFSTG = CFPAREL(L,2)
              NFTSTG = NFTPAR(INDEX)
        
              IF (IOFLLI.NE.0) THEN

                  NFNLSTG = NFNLPAR(INDEX)
                  NFTIP = 0
                  IF (NFNLSTG.NE.0) NFTIP = NFNLTIP(NFNLSTG)

              ELSE

                  NFNLSTG = 0
                  NFTIP = 0

              END IF !IOFLLI.NE.0   

              IF (IOFLSAT.NE.0) THEN

                  STOR = CFSTG*STGC

                  IF (NFTSTG.NE.0) THEN

                      FT1 = FNT(NFTSTG,INTI)
                      FT2 = FNT(NFTSTG,INTI+1)
                      TEMPCOEF = FT2*DTIMEF+(1D0-DTIMEF)*FT1
                      STOR = STOR*TEMPCOEF

                  END IF !IOFLSAT.NE.0

                  INPORC = INORPAR(15) + NZONEPOR  
                  NZONEPOR = LXPAREL(L,7)
                  PORC = PARC(INPORC)
                  CFPOR = CFPAREL(L,7)
                  NFTPOR = NFTPAR(INPORC)
                  PORC = PORC*CFPOR

                  IF (NFTPOR.NE.0) THEN

                      FT1 = FNT(NFTPOR,INTI)
                      FT2 = FNT(NFTPOR,INTI+1)
                      TEMPCOEF = FT2*DTIMEF+(1D0-DTIMEF)*FT1
                      PORC = PORC*TEMPCOEF

                  END IF !NFTPOR.NE.0

                  XPARAM(1) = PORC/STOR
                  XPARAM(2) = STOR

              ELSE

                  XPARAM(1) = 0D0
                  XPARAM(2) = 0D0

              END IF !IOFLSAT.NE.0

              CALL PARAM_VALUE
     &            (CFSTG    ,DERSCC   ,DERSCH   ,DTIMEF   ,EPSFLU
     &            ,IDIMFNT  ,INDSSTR  ,INTI     ,IOFMLF   ,L
     &            ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLSTG  ,NFTIP
     &            ,NFTSTG   ,NINT     ,NNUD     ,NPARALG  ,NUMEL
     &            ,NUMNP    ,NZPRG    ,STGC     ,PAREL(L,7)
     &            ,PARC(INPRGC)       ,FNT      ,HCALAN   ,HCALIT
     &            ,IFLAGS             ,IPAR_DIR ,KXX      ,NFNLPRG     
     &            ,PARACD    ,XPARAM)

              IF (IOFLLI.NE.0) THEN

                  DPARELDH(7,L) = DERSCH
                  IF (IODENS.EQ.1) THEN
                  
                      DPARELDC(7,L) = DERSCC

                  END IF !IODENS.EQ.1

              END IF !IOFLLI.NE.0

          END IF !INDSSTR.NE.0.AND.NZONE_PAR(2).NE.0
 
C------------------------- Step 2.3: Areal recharge

          IF (NZONE_PAR(3).NE.0) THEN

              NZONEARR=LXPAREL(L,3+INDSSTR)

              IF (NZONEARR.NE.0) THEN

                  INDEX = INORPAR(8) + NZONEARR
                  CFARR = CFPAREL(L,3+INDSSTR)
                  ARRC = PARC(INDEX)
                  NFTARR = NFTPAR(INDEX)
            
                  IF (IOFLLI.NE.0) THEN

                      NFNLARR = NFNLPAR(INDEX)
                      NFTIP = 0

                      IF (NFNLARR.NE.0) NFTIP = NFNLTIP(NFNLARR)

                  ELSE

                      NFNLARR = 0
                      NFTIP = 0

                  END IF !IOFLLI.NE.0

                  CALL PARAM_VALUE
     &                (CFARR    ,DERSCC   ,DERSCH   ,DTIMEF   ,EPSFLU
     &                ,IDIMFNT  ,INDSSTR  ,INTI     ,IOFMLF   ,L
     &                ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLARR  ,NFTIP
     &                ,NFTARR   ,NINT     ,NNUD     ,NPARALG  ,NUMEL
     &                ,NUMNP    ,NZPRG    ,ARRC     ,PAREL(L,8)
     &                ,PARC(INPRGC)       ,FNT      ,HCALAN   ,HCALIT
     &                ,IFLAGS   ,IPAR_DIR ,KXX      ,NFNLPRG  ,PARACD
     &                ,XPARAM)

              END IF !NZONEARR.NE.0

          END IF !NZONE_PAR(3).NE.0

      END DO !L=1,NUMEL

C------------------------- Step 3: Loop over mesh nodal points

      DO I=1,NUMNP

          IB = IBCOD(I)

          IF (IB.NE.0) THEN                          ! Non-zero boundary condition

C------------------------- Step 3.1: Prescribed head

              IF (IB.EQ.1 .OR. IB.GE.3) THEN

                  IPRESC = IXPARNP(I,1+INDSSTR)
                  INDEX = INORPAR(9) + IPRESC
                  CFCHP = CFPARNP(I,1+INDSSTR)
                  CHPC = PARC(INDEX)
                  NFTCHP = NFTPAR(INDEX)
  
                  IF (IOFLLI.NE.0) THEN

                      NFNLCHP = NFNLPAR(INDEX)
                      NFTIP = 0

                      IF (NFNLCHP.NE.0) NFTIP = NFNLTIP(NFNLCHP)

                  ELSE

                      NFNLCHP = 0
                      NFTIP = 0

                  END IF !IOFLLI.NE.0

C------------------------ Prescribed head boundary condition must be evaluated
C------------------------ at k+1 (contrarily to k+theta used for the rest of BC)

              DTIMAUX = DTIMEF
              IF (IB.EQ.1 .AND. INDSSTR.NE.0) THEN
              
                  DTIMAUX = DTIMAUX + (1D0-PAR_DIR(29))*TINC/TINTERVOBS

	        END IF !IB.EQ.1 .AND. INDSSTR.NE.0

                  CALL PARAM_VALUE
     &                (CFCHP    ,DERSCC   ,DERSCH   ,DTIMAUX  ,EPSFLU
     &                ,IDIMFNT  ,INDSSTR  ,INTI     ,1        ,I
     &                ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLCHP  ,NFTIP
     &                ,NFTCHP   ,NINT     ,1        ,NPARALG  ,NUMEL
     &                ,NUMNP    ,NZPRG    ,CHPC     ,PARNP(I,1)
     &                ,PARC(INPRGC)       ,FNT      ,HCALAN   ,HCALIT
     &                ,IFLAGS   ,IPAR_DIR ,KXX      ,NFNLPRG  ,PARACD
     &                ,XPARAM)


                  IF (IOFLLI.NE.0) DNODALRH(I,1) = DERSCH  

              END IF ! IBCOD(I).EQ.1 .OR. IBCOD(I).GE.3

C------------------------- Step 3.2: Prescribed flow

              IF (IB.EQ.2 .OR. IB.EQ.4) THEN

                  INDEX = INORPAR(10) + IXPARNP(I,3+INDSSTR)
                  CFQQP = CFPARNP(I,3+INDSSTR)
                  QQPC = PARC(INDEX)
                  NFTQQP = NFTPAR(INDEX)

                  IF (IOFLLI.NE.0) THEN

                    NFNLQQP = NFNLPAR(INDEX)
                    NFTIP = 0

                    IF (NFNLQQP.NE.0) NFTIP = NFNLTIP(NFNLQQP)

                  ELSE

                      NFNLQQP = 0
                      NFTIP = 0

                  END IF !IOFLLI.NE.0

                  CALL PARAM_VALUE
     &                (CFQQP    ,DERSCH   ,DERSCC   ,DTIMEF   ,EPSFLU
     &                ,IDIMFNT  ,INDSSTR  ,INTI     ,1        ,I
     &                ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLQQP  ,NFTIP
     &                ,NFTQQP   ,NINT     ,1        ,NPARALG  ,NUMEL
     &                ,NUMNP    ,NZPRG    ,QQPC     ,PARNP(I,2)
     &                ,PARC(INPRGC)       ,FNT      ,HCALAN   ,HCALIT
     &                ,IFLAGS   ,IPAR_DIR ,KXX      ,NFNLPRG  ,PARACD
     &                ,XPARAM)
            
              END IF ! IBCOD(I).EQ.2 .OR. IBCOD(I).EQ.4

C------------------------- Step 3.3: Leakage coefficient

              IF (IB.GE.3) THEN
            
                  NZONEALF = IXPARNP(I,5+INDSSTR)
                  INDEX = INORPAR(11) + NZONEALF
                  CFALF = CFPARNP(I,5+INDSSTR)
                  ALFC = PARC(INDEX)
                  NFTALF = NFTPAR(INDEX)

                  IF (IOFLLI.NE.0.AND.INDSSTR.NE.0) THEN

                      NFNLALF = NFNLPAR(INDEX)

                      IF (NFNLALF.NE.0) THEN

                          NFTIP = NFNLTIP(NFNLALF)
                          XPARAM(1) = PARNP(I,1)
                          XPARAM(2) = DNODALRH(I,1)
                          XPARAM(3) = CFALF*ALFC

                          IF (NFTALF.NE.0) THEN

                              FT1 = FNT(NFTALF,INTI)
                              FT2 = FNT(NFTALF,INTI+1)
                              TEMPCOEF = DTIMEF*FT2 + (1D0-DTIMEF)*FT1
                              XPARAM(3) = XPARAM(3)*TEMPCOEF

                          END IF !NFTALF.NE.0

                      END IF !NFNLALF.NE.0

                  ELSE

                      NFNLALF = 0
                      NFTIP = 0
                      XPARAM(1) = 0D0
                      XPARAM(2) = 0D0
                      XPARAM(3) = 0D0

                  END IF !IOFLLI.NE.0.AND.INDSSTR.NE.0

                  CALL PARAM_VALUE
     &                (CFALF    ,DERSCC   ,DERSCH   ,DTIMEF   ,EPSFLU
     &                ,IDIMFNT  ,INDSSTR  ,INTI     ,1        ,I
     &                ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLALF  ,NFTIP
     &                ,NFTALF   ,NINT     ,1        ,NPARALG  ,NUMEL
     &                ,NUMNP    ,NZPRG    ,ALFC     ,PARNP(I,3)
     &                ,PARC(INPRGC)       ,FNT      ,HCALAN   ,HCALIT
     &                ,IFLAGS   ,IPAR_DIR ,KXX      ,NFNLPRG  ,PARACD
     &                ,XPARAM)


                  IF (IOFLLI.NE.0.AND.INDSSTR.NE.0) THEN

                      DNODALRH(I,2) = XPARAM(5) ! Der. of Qb wrt. head level
                      DNODALRH(I,3) = DERSCH    ! Der. of ALFA wrt. head level
                      DNODALRH(I,4) = XPARAM(4) ! Value of Qb

                  END IF !IOFLLI.NE.0.AND.INDSSTR.NE.0

              END IF ! IBCOD(I).GE.3

          END IF ! IBCOD(I).NE.0

      END DO !I=1,NUMNP

C------------------------- Writes PAREL,PARNP in MAINF

      IF (IFLAGS(30).NE.0) THEN

          WRITE(MAINF,1000)
 1000     FORMAT(//,10X,'PAREL ARRAY',/,10X,'===== =====',/,
     &            5X,'TRA1',9X,'TRA2',9X,'TRA3',9X,'TRA4',9X,
     &               'TRA5',9X,'TRA6',9X,'STG',10X,'ARR')

          DO L=1,NUMEL

              WRITE(MAINF,1100) L,(PAREL(L,J),J=1,MAX(ISOT,IODIM)),
     &                          (PAREL(L,J),J=MAX(ISOT,IODIM)+1,8)
          END DO !L=1,NUMEL

          WRITE(MAINF,1200)
 1200     FORMAT(//,10X,'PARNP ARRAY',/,10X,'===== =====',/,
     &            5X,'CHP',10X,'QQP',10X,'ALF')

          DO I=1,NUMNP

              WRITE(MAINF,1100) I,(PARNP(I,J),J=1,3)

          END DO !I=1,NUMNP

 1100     FORMAT(I5,8G10.3)

      END IF !IFLAGS(30).NE.0

      END SUBROUTINE COMP_PARAM_FLOW
