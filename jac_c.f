      SUBROUTINE JAC_C
     &          (ACTH     ,AFLU     ,AREA     ,ATRA     ,ATRADSC
     &          ,ATRADSCF ,BETAC    ,BIBI     ,BUOYANCY ,CAUDAL
     &          ,CAUX1    ,CAUX2    ,CCALAN   ,CCALIT   ,CFPAREL
     &          ,CFPARNP  ,COORD    ,CREF     ,DAT_VD   ,DENSITY
     &          ,DENSREF  ,DERC     ,DERCS    ,DERH     ,DFLU
     &          ,DTIM     ,DTRA     ,DTRADTRA ,DVDP     ,EPSFLU
     &          ,EPSTRA   ,FNT      ,GP_COORD ,GRADLOC  ,GRDFF
     &          ,HAUX1    ,HAUX2    ,HBASE    ,HCALAN   ,HCALIT
     &          ,HINI     ,IAD_S    ,IADD_S   ,IADN_S   ,IADSC_COLS
     &          ,IADSC_ROWS         ,IAFD_S   ,IAFDD_S  ,IAFDN_S
     &          ,IATRADSC_COLS      ,IATRADSC_ROWS      ,IBCOD
     &          ,IBTCO    ,IDESC    ,IDIMAFLU ,IDIMATRA ,IDIMBB
     &          ,IDIMDERC ,IDIMDERH ,IDIMDFLU ,IDIMDQ   ,IDIMDTRA
     &          ,IDIMFNT  ,IDIMQ    ,IDIMWORK ,IDIRECT  ,IFLAGS
     &          ,INEW     ,INEWT    ,INORPAR  ,INTI     ,IOCAP
     &          ,IOCTRA   ,IODENS   ,IODIM    ,IOFLLI   ,IOFMLF
     &          ,IOFMLT   ,IOINV    ,IOLD     ,IOLDT    ,IOLG_PAR
     &          ,IORDCH   ,IOTRLI   ,IOVRWC   ,IPAR_DIR ,IPROB
     &          ,ISOLFL   ,ISOLTR   ,ISOZ     ,ITERM    ,ITPTVAR
     &          ,ITYPAFLU ,ITYPATRA ,ITYPDFLU ,ITYPDTRA ,IVPAR
     &          ,IXPARNP  ,IXPARFL  ,KXX      ,LDIM     ,LMXNDL,LNNDEL
     &          ,LTYPE    ,LXPAREL  ,LXPARFL  ,MAINF    ,MAXNB
     &          ,MAXNBF   ,MAXPG    ,NBAND1   ,NFLAGS   ,NFNL
     &          ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT
     &          ,NPAR     ,NPARALG  ,NPAREL   ,NPARF    ,NPARNP
     &          ,NPBTP    ,NPPEL    ,NPPNP    ,NTYPAR   ,NTYPEL
     &          ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR    ,NZTRA
     &          ,PARACD   ,PARC     ,PAREL    ,PARNP    ,POINTWEIGHT
     &          ,QXYZ     ,SOLUTION ,THETAF   ,THETAT   ,TINC
     &          ,TINTERVOBS         ,VD       ,WATVOL   ,WORK
     &          ,WSPECHEAT,WTHERMCON,DERH1,DERH2
     ;          ,IDIMWGT   ,WGT_PAR  ,IPNT_PAR ,IPOS    ,DERIV
     &          ,IOCONSRC  ,DTRADFLU ,CFLU     ,DFLUDFLU,DFLUDTRA
     &          ,DNODALRH  ,DQDFLU   ,DQDTRA   ,LINMET  ,DPARELDH)

*******************************************************************************
*
* PURPOSE
*
*     Computes derivatives of concentration w.r.t. estimated parameters
*
* DESCRIPTION
*
*     Computes derivatives of concentration w.r.t. estimated parameters. 
*     First computes the RHS of the derivatives and then solves it.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AFLU                   Matrix of finite elements equations for flow problem  
*                         No boundary conditions are included on it.            
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  ATRA                   Matrix of finite elements equations for transport     
*                         problem. Only mass flow boundary conditions included. 
*  BIBI                   Array containing the product of interpolation         
*                         functions gradient, for a given element               
*  CAUDAL                 Input/output flow at every node.                      
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  CAUX2                  Array containing diference of concentrations in two   
*                         consecutives times, related to time increment         
*  CCALIT                 Computed concentration at every node                  
*  CCALAN                 Computed concentrations in the previous time step.    
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  CNST                   CNST(i,j,k) is the integral of the product of         
*                         interpolation functions i and j in an element of      
*                         type k divided by the AREA of the element. It is used 
*                         only in consistent scheme.                            
*  COORD                  Nodal coordinates                                     
*  DAT_VD                 Derivatives of ATRA w.r.t. the components of Darcy's  
*                         velocity (qx, qy, qz).                                
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
*  DERH                   Nodal head derivatives with respect to estimated      
*                         flow parameters.                                      
*  DVDP                   Derivatives of Darcy's velocity with respect to       
*                         flow parameters. In some cases it is dimensioned      
*                         as (IODIM,NUMEL) to reduce storage.                   
*  DFLU                   Matrix of finite elements equations for flow          
*                         problem related to storage term.                      
*  DTRA                   Matrix of finite elements equations for transport     
*                         problem related to storage.                           
*  FNT                    Array containing time functions values                
*  GRAVEL                 Projection of gravity at each element                 
*  GRDFF                  Array containing the product between interpolation    
*                         functions integrals and interp. functions gradient    
*  HAUX1                  Array containing heads, ponderated by THETAF time     
*                         weight. Is equal to THETAF*HCAL+(1-THETAF)*HCALAN     
*  HAUX2                  Array containing difference of heads in two           
*                         consecutives times. Is equal to HCAL-HCALAN/TIME STEP 
*  HCAL                   Computed heads at every node                          
*  HCALAN                 Head level at previous time                           
*  HEAD                   Auxiliar array used to store nodal heads              
*  IBCOD                  Flow boundary condition index                         
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  ISOZ                   Anisotropy of every transmissivity zone               
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  IXPARNP                Array containing zone number of each node j,          
*                         corresponding to INpar index parameter zone.          
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  LXPAREL                Array containing zone number of each element j,       
*                         corresponding to INpar index parameter zone.          
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
*
* INTERNAL VARIABLES: ARRAYS
*
*  IAUX                   If zero, there is not estimated any parameter of
*                         type J. It is used locally to make some additional
*                         computations needed due to the explicit dependence
*                         of some expresions w.r.t some parameters 
*
* EXTERNAL VARIABLES: SCALARS
*
*  DERH1                                                                        
*  DERH2                                                                        
*  DTIM                                                                         
*  EPSTRA                 Time weighting parameter for nonlinear transport      
*                         problems                                              
*  IDIMAFLU               Used to dimension array AFLU                          
*  IDIMBB                 Used to dimension array BIBI. Is equal to IDIMQ times 
*                         the maximum possible anisotropy of the problem        
*  IDIMDQ                 Used to dimension array DAT_VD (second dimension). It 
*                         is the number of different terms (ATRA)ij varying     
*                         i and  j in the local numeration of one element.      
*                         Is equal to LMXNDL*(LMXNDL-1)/2                       
*  IDIMFNT                First dimension of array FNT, it coincides with       
*  IDIMQ                  Used to dimension array QXYZ. It is equal to          
*                         IODIM*(IODIM+1)/2                                     
*  INEW                   Third dimension of array DERH used to store           
*                         derivatives of H w.r.t. param. at current time step.
*                         However, at this point it stores the derivatives at 
*                         previous time step
*  INEWT                  Third dimension of array DERC used to store           
*                         derivatives of C w.r.t. parameters at next time step  
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOCNSF                 Scheme for storage term in flow problem               
*  IOCNST                 Scheme for mass storage term in transport problem     
*  IODIM                  Maximum dimension of the problem                      
*  IOFMLT                 Transport formulation number                          
*  IOLD                   Third dimension of array DERH used to store           
*                         derivatives of H w.r.t. parameters at previous time   
*                         step. However, at this point it stores the 
*                         derivatives at current time step                                    
*  IOLDT                  Third dimension of array DERC used to store           
*                         derivatives of C w.r.t. parameters at previous time   
*                         step                                                  
*  IOLG_PAR               Array containing all logarithmic options of           
*                         estimated parameters                                  
*  IOPRHED                Indicates whether the flow state variable state is    
*                         preasure (set to 1) or head (set to 0)                
*  IOTRLI                 Idem to IOFLLI, in the case of transport.             
*  ISOLFL                 If 0, no flow has been solved at current time.
*                         If 1, steady flow has been solved at current time.
*                         If 2, transient flow has been solved at current time
*  ISOLTR                 If 0, no transport has been solved at current time.
*                         If 1, steady transport has been solved at current time
*                         If 2, trans. transport has been solved at current time
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NBAND                  Half Bandwith (maximum difference between the         
*                         numbers of two nodes belonging to the same element)   
*  NBAND1                 Used to dimension. It is equal to NBAND+1             
*  NBAND2                 Used to dimension. It is equal to 2*NBAND+1           
*  NFLAGS                 Maximum number of allowed flags                       
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NINT                   Number of observation times                           
*  NPAR                   Total number of parameters to be estimated            
*  NPAREL                 Number of element parameters in current problem       
*  NPARF                  Number of transient parameters to be estimated        
*  NPARNP                 Number of nodal parameters in current problem         
*  NPBFL                  Number of simultaneous flow problems                  
*  NPBMX                  Max(NPBFL,NPBTP). Used to dimension some arrays       
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
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  PARNP                  Parameter values at every node and current time for   
*                         all nodal parameters (each value is computed as the   
*                         product of up to four terms:                          
*                           nodal coeff*zonal value*time funct.*nonl. funct. )  
*  QXYZ                   Products between the different components of          
*                         Darcy's velocity divided by its norm                  
*  THETAF                 Time weighting parameter for flow problems            
*  THETAT                 Time weighting parameter for transport problems       
*  TINC                   Current time increment                                
*  WATVOL                 Array containing the water content of every element   
*  WORK                   Workspace array used in some subroutines.             
*  XNORVD                 Euclidean norm of Darcy's velocity                    
*
* INTERNAL VARIABLES: SCALARS
*
*  ATRADSC                Coefficient matrix of transport system. Most often in 
*                         its decomposed form.                                  
*  IDIMDFLU               Used to dimension array DFLU                          
*  INARR                  Index for areal recharge                              
*                         in array variables (LXPAREL and CFPAREL)              
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOFMLF                 Flow Formulation number                               
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  DERCOE                                                                       
*  DERCRD                                                                       
*  DERDFM                                                                       
*  DERDSL                                                                       
*  DERDST                                                                       
*  DERFOD                                                                       
*  DERIVATIVE_VD                                                                
*  DERPOR                                                                       
*  DERQ_ARR                                                                     
*  DERQ_GEN                                                                     
*  DERQ_NUD                                                                     
*  DERQ_STG                                                                     
*  DER_ATRA_FP_VD                                                               
*  DER_ATRA_VD                                                                  
*  DER_VD_TRA                                                                   
*  ENS_IND_T3                                                                   
*  ENS_IND_T4                                                                   
*  LEQT1B                                                                       
*
* HISTORY
*
*     AMS      3-2002     First coding
*
*******************************************************************************

* OJO CON CADENAS DERIVADAS DE CRD, FOD Y  WATVOL 

      IMPLICIT NONE


C------------------------- External

      INTEGER*4::IADSC_COLS         ,IADSC_ROWS         ,IATRADSC_COLS
     &          ,IATRADSC_ROWS      ,IDESC    ,IDIMAFLU ,IDIMATRA
     &          ,IDIMBB   ,IDIMDERC ,IDIMDERH ,IDIMDFLU ,IDIMDQ,IDIMDTRA
     &          ,IDIMFNT  ,IDIMQ    ,IDIMWORK ,IDIRECT  ,INEW
     &          ,INEWT    ,INTI     ,IOCAP    ,IOCTRA   ,IODENS
     &          ,IODIM    ,IOFLLI   ,IOFMLF   ,IOFMLT   ,IOINV
     &          ,IOLD     ,IOLDT    ,IORDCH   ,IOTRLI   ,IOVRWC
     &          ,IPROB    ,ISOLFL   ,ISOLTR   ,ITERM    ,ITPTVAR
     &          ,ITYPAFLU ,ITYPATRA ,ITYPDFLU ,ITYPDTRA ,LMXNDL
     &          ,MAINF    ,MAXNB    ,MAXNBF   ,MAXPG    ,NBAND1
     &          ,NFLAGS   ,NFNL     ,NINT     ,NPAR     ,NPARALG
     &          ,NPAREL   ,NPARF    ,NPARNP   ,NPBTP    ,NPPEL
     &          ,NPPNP    ,NTYPAR   ,NTYPEL   ,NUMEL    ,NUMNP
     &          ,NZPAR    ,NZTRA
     ;          ,IDIMWGT  ,IPNT_PAR ,IPOS(NPAR)
     &          ,IOCONSRC


      REAL*8::BETAC    ,CREF     ,DENSREF  ,DTIM     ,EPSFLU   ,EPSTRA
     &       ,THETAF   ,THETAT   ,TINC     ,TINTERVOBS
     &       ,WSPECHEAT,WTHERMCON,DERIV(NPAR)    ,WGT_PAR


      INTEGER*4::IAD_S(MAXNB,NUMNP),IADD_S(NUMNP),IADN_S(NUMNP)
     &          ,IAFD_S(MAXNBF,NUMNP),IAFDD_S(NUMNP),IAFDN_S(NUMNP)
     &          ,IBCOD(NUMNP),IBTCO(NUMNP),IFLAGS(NFLAGS)
     &          ,INORPAR(NTYPAR),IOLG_PAR(NTYPAR,2),IPAR_DIR(NPARALG)
     &          ,ISOZ(NZTRA),IVPAR(NZPAR),IXPARNP(NUMNP,NPARNP)
     &          ,KXX(LMXNDL,NUMEL),LDIM(NUMEL),IXPARFL(NUMNP,NPARNP)
     &          ,LNNDEL(NUMEL),LTYPE(NUMEL),LXPAREL(NUMEL,NPAREL)
     &          ,LXPARFL(NUMEL,NPAREL),NFNLPAR(NZPAR),NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL),NFTPAR(NZPAR),NZONE_PAR(NTYPAR)
     &          ,LINMET(*)

      REAL*8::ACTH(NUMEL),AFLU(NUMEL,IDIMAFLU),AREA(NUMEL)
     &       ,ATRA(NUMEL,IDIMATRA),ATRADSC(IATRADSC_ROWS,IATRADSC_COLS)
     &       ,ATRADSCF(MAXNBF,NUMNP),BIBI(IDIMBB,NUMEL)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL),CAUDAL(NUMNP),CAUX1(NUMNP)
     &       ,CAUX2(NUMNP),CCALAN(NUMNP),CCALIT(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)
     &       ,CFPARNP(NUMNP,NPARNP),COORD(NUMNP,3)
     &       ,DAT_VD(IODIM,IDIMDQ,NUMEL),DENSITY(NUMEL)
     &       ,DERC(NUMNP,NPAR,IDIMDERC),DERCS(NUMNP,NPAR,2)
     &       ,DERH(NUMNP,NPARF,IDIMDERH)
     &       ,DFLU(NUMNP),DTRA(NUMNP)
     &       ,DTRADTRA(NUMEL,LMXNDL*LMXNDL),DVDP(IODIM,NPARF,NUMEL)
     &       ,FNT(IDIMFNT,NINT),GP_COORD(6,8,IODIM)
     &       ,GRADLOC(IODIM,LMXNDL,MAXPG),GRDFF(IODIM,LMXNDL,NUMEL)
     &       ,HAUX1(NUMNP),HAUX2(NUMNP),HBASE(NUMEL),HCALAN(NUMNP)
     &       ,HCALIT(NUMNP),HINI(NUMNP),PARACD(3,NFNL)
     &       ,PARC(NZPAR),PAREL(NUMEL,NPPEL),PARNP(NUMNP,NPPNP)
     &       ,POINTWEIGHT(6,8),QXYZ(IDIMQ,NUMEL),SOLUTION(NUMNP)
     &       ,VD(IODIM,NUMEL),WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
     &       ,WORK(IDIMWORK),DERH1(NUMNP),DERH2(NUMNP)
     &       ,DTRADFLU(NUMEL,LMXNDL*LMXNDL),CFLU(NUMEL,IDIMDFLU)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DNODALRH(NUMNP,4)
     &       ,DQDFLU(NUMEL,LMXNDL*LMXNDL),DQDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DPARELDH(NPPEL,NUMEL)

C------------------------- Internal

      INTEGER*4::I        ,ICMPDERS ,IFP      ,IP       ,ISSFL    ,J
     &          ,NP

      REAL*8::EPST1,FACTOR,THETAF1,THETAT1,TINCINV

      CHARACTER::strFmt1*20

      INTEGER*4::IAUX(20)

      REAL*8::DUMMY(1)

      REAL*8,ALLOCATABLE::DERB(:)
	
C------------------------- First Executable Statment

      IF(IFLAGS(3).EQ.1) CALL IO_SUB('JAC_C',0)
	  
      strFmt1 = ''
      WRITE (strFmt1,*) NPAR
      strFmt1 = '(I5,'//Trim(AdjustL(strFmt1))//'F15.10)'

C------------------------- Initializes to zero the RHS of the derivatives
C------------------------- (father (the first time) and the "son" if necessary)
C------------------------- Only if density is constant, otherwise, it is already
C------------------------- initialized in jac_h.

      IF (IPROB.EQ.1 .OR. IORDCH.EQ.0 .AND. IODENS.EQ.0) THEN
          DERC(:,:,INEWT) = 0D0
      END IF !IPROB.EQ.1 .OR. IORDCH.EQ.0


      IF (IPROB.NE.NPBTP .AND. IORDCH.NE.0) THEN
          DERCS(:,:,INEWT) = 0D0
      END IF !IPROB.NE.NPBTP .AND. IORDCH.NE.0

C------------------------- Sets variable ICMPDERS

      IF (IPROB.NE.NPBTP .AND. IORDCH.NE.0) THEN
          ICMPDERS = 1
      ELSE
          ICMPDERS = 0
      END IF !IPROB.NE.NPBTP .AND. IORDCH.NE.0

C------------------------- Adds the nonlinear part to inverse problem
C------------------------- RHS, related to the previous time step.

      IF (IOTRLI.NE.0 .OR. IODENS.EQ.1) THEN

          FACTOR = (EPSTRA - 1D0)/EPSTRA

          ALLOCATE (DERB(NUMNP))
          DERB = 0D0

          CALL COMP_DER_BTRA
     &        (CAUX1    ,DENSREF  ,DERB     ,IBTCO    ,IOCONSRC
     &        ,IODENS   ,ITPTVAR  ,NPPNP    ,NUMNP
     &        ,PARNP    ,THETAT   ,WSPECHEAT)

          CALL RHS_NOLI_IN 
     &        (DTRADTRA ,DERB     ,FACTOR   ,-1        ,KXX
     &        ,LMXNDL   ,LNNDEL   ,NPARF    ,NUMEL    ,NUMNP
     &        ,DERC(1,1,INEWT)    ,DERC(1,1,IOLDT))

C------------------------- Corrects contribution of DTRA due to the 
C------------------------- backward integration in time scheme.

          EPST1 = 1D0 - EPSTRA

          IF ( (IODENS.GT.0 .AND. EPST1.NE.EPSTRA) .OR.
     &         (EPST1.NE.EPSTRA .OR. IOFLLI.GT.0)) THEN

              CALL  CORRECT_DDTRA_RHS
     &            (AREA     ,BETAC    ,CAUX2    ,CCALAN    ,CCALIT
     &            ,CREF     ,DENSREF  ,DERC     ,DERH      ,DPARELDH
     &            ,DTRA     ,EPSTRA   ,HCALAN   ,HINI      ,IDIMDTRA
     &            ,INEWT    ,IODENS   ,IOFLLI   ,IOLD      ,IOLDT
     &            ,IOVRWC   ,ITPTVAR  ,KXX      ,LMXNDL    ,LNNDEL
     &            ,NPAR     ,NPPEL    ,NUMEL    ,NUMNP     ,PAREL
     &            ,THETAT   ,WATVOL)

          END IF !(IODENS.GT.0 .AND. EPST1.NE.EPSTRA) ...

          IF (IODENS.EQ.1) THEN

              CALL RHS_NOLI_IN
     &          (DTRADFLU ,DUMMY    ,FACTOR   ,0        ,KXX
     &          ,LMXNDL   ,LNNDEL   ,NPARF    ,NUMEL    ,NUMNP
     &          ,DERC(1,1,INEWT)    ,DERH(1,1,IOLD))

C------------------------- For nodal fluxes when Newton method.

              CALL COMP_DER_FLOW
     &            (AFLU     ,BETAC    ,CAUDAL   ,CFLU     ,CREF
     &            ,DENSREF  ,DFLU     ,DFLUDFLU ,DFLUDTRA ,DNODALRH
     &            ,DQDFLU   ,DQDTRA   ,EPSFLU   ,EPSTRA   ,HAUX1
     &            ,IBCOD    ,IDIMAFLU ,IDIMDFLU ,1        ,IODENS
     &            ,ISOLFL   ,KXX      ,LINMET   ,LMXNDL   ,LNNDEL
     &            ,NPPNP    ,NUMEL    ,NUMNP    ,PARNP    ,TINC)

              CALL COMAT_BC
     &            (DUMMY    ,BETAC    ,CAUDAL   ,CAUX1    ,CREF
     &            ,DENSREF  ,DQDFLU   ,DQDTRA   ,EPSTRA
     &            ,DUMMY    ,DUMMY    ,DUMMY    ,IBTCO    ,1
     &            ,1        ,1        ,0        ,1        ,IODENS
     &            ,1        ,1        ,KXX      ,LMXNDL   ,LNNDEL
     &            ,1        ,1        ,1        ,NPPNP    ,NUMEL
     &            ,NUMNP    ,PARNP    ,THETAT
     &            ,AREA,NZONE_PAR,PAREL,NPPEL,NTYPAR)

             CALL RHS_NOLI_IN
     &           (DQDFLU   ,DUMMY    ,1D0     ,0        ,KXX
     &           ,LMXNDL   ,LNNDEL   ,NPAR     ,NUMEL    ,NUMNP
     &           ,DERC(1,1,INEWT)    ,DERH(1,1,IOLD))

             CALL RHS_NOLI_IN
     &           (DQDTRA   ,DUMMY    ,1D0     ,0        ,KXX
     &           ,LMXNDL   ,LNNDEL   ,NPAR     ,NUMEL    ,NUMNP
     &           ,DERC(1,1,INEWT)    ,DERC(1,1,IOLDT))

          END IF !IODENS.EQ.1

      END IF !IOFLLI.NE.0 .OR. IODENS.EQ.1

      IF (ISOLTR.EQ.2) THEN

          THETAT1 = THETAT - 1D0
          TINCINV = 1D0/TINC

          DO IP=1,NPAR

C------------------------- (THETAT-1)*ATRA*DERC(OLD)

              IF (THETAT.LT.1D0) THEN

                  CALL PROD_MAT_VEC
     &                (THETAT1  ,IAD_S   ,IADN_S
     &                ,IDIMATRA ,NUMEL    ,NUMNP    ,1
     &                ,ITYPATRA ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &                ,KXX      ,LNNDEL   ,ATRA     ,DERC(1,IP,INEWT)
     &                ,DERC(1,IP,IOLDT))

              END IF !THETAT.LT.1D0

C------------------------- 1/dt *DTRA*DERC(OLD)

              CALL PROD_MAT_VEC
     &            (TINCINV  ,IAD_S    ,IADN_S     
     &            ,IDIMDTRA ,NUMEL    ,NUMNP        ,1
     &            ,ITYPDTRA ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &            ,KXX      ,LNNDEL   ,DTRA     ,DERC(1,IP,INEWT)
     &            ,DERC(1,IP,IOLDT))

C------------------------- Nodal flow in nodes with mass flow bnd. cond.
C------------------------- (THETAT - 1)*rho*Q*DERC(OLD)
C------------------------- Only when Picard iterations (i.e. constant density).

              IF (THETAT.LT.1D0 .AND. IODENS.EQ.0) THEN

                  CALL RHS_NODE_FLUX_INV
     &                (BETAC    ,CAUDAL   ,CAUX1    ,CREF     ,DENSREF
     &                ,DERC     ,IBTCO    ,IDIMDERC ,INEWT    ,IODENS
     &                ,IOLDT    ,NUMNP    ,NPAR     ,THETAT)

              END IF !THETAT.LT.1D0 .AND. IODENS.EQ.0

C------------------------- Concentration leakage

              IF (NZONE_PAR(18).NE.0) THEN

C------------------------- This contribution is equal to (TH-1)*LEAK*DERC(K)

c                      DERC(1:NUMNP,IP,INEWT) = DERC(1:NUMNP,IP,INEWT)
c     &                 + THETAT1*PARNP(1:NUMNP,6)*DERC(1:NUMNP,IP,IOLDT)

              END IF !(NZONE_PAR(18).NE.0)

          END DO !IP=1,NPAR

      END IF !ISOLTR.EQ.2



      IF (IOINV.EQ.3) THEN !If any flow parameter is being estimate.

C------------------------- Assigns auxiliar variables. IAUX is used to 
C------------------------- make a unique CALL to subroutines 
C-------------------------     DER_VD_TRA, DERQ_TRA, DERQ_STG,
C-------------------------     DERQ_ARR,   DERQ_NUD, DER_WTV_STG

          DO I=1,6 !IOLG(:,2)=1 => se estima alguna zona. 0 => no.

              IF (IOLG_PAR(I,2).NE.0) THEN
                  IAUX(I)=1
              ELSE
                  IAUX(I)=0
              END IF !IOLG_PAR(I,2).NE.0

          END DO !I=1,6

C------------------------- Assigns ISSFL (0 for steady flow and 1 for 
C------------------------- transient) variable when no flow has been solved at
C------------------------- the current time step

          IF (ISOLFL.EQ.0 .OR. ISOLFL.EQ.1) THEN
             ISSFL=0          ! Steady flow
          ELSE
             ISSFL=1          ! Transient flow
          ENDIF !ISOLFL.EQ.0 .OR. ISOLFL.EQ.1


C------------------------- Loop over flow parameters

          DO IFP=1,NPARF

C------------------------- Computes auxiliar arrays DERH1 and DERH2. At this 
C------------------------- point, IOLD location contains current time head 
C------------------------- derivative and location INEW contains previous time
C------------------------- head derivative
C------------------------- Only for constant density.

              IF (IODENS.EQ.0) THEN

                  IF (ISOLFL.EQ.2) THEN  ! Transient flow

                      THETAF1 = 1D0 - THETAF
                      THETAT1 = 1D0 - THETAT

                      DERH1(:) = THETAF*DERH(:,IFP,IOLD) +
     &                               THETAF1*DERH(:,IFP,INEW)

                      DERH2(:) = (DERH(:,IFP,IOLD) - DERH(:,IFP,INEW))
     &                           /TINC

                  END IF !ISOLFL.EQ.2

C------------------------- Computes derivatives of VD w.r.t. flow parameter 
C------------------------- IFP (only common part, i.e. derivatives w.r.t
C------------------------- the explicit dependence with transmissivity not 
C------------------------- included).
C------------------------  Only for constant density.

                  CALL DER_VD
     &   (IODIM    ,LMXNDL   ,MAINF    ,NFLAGS   ,NPAREL   ,NPPEL
     &   ,NUMEL    ,NUMNP    ,NZONE_PAR(1)       ,DERH1,DVDP
     &   ,GRDFF    ,IFLAGS   ,ISOZ     ,KXX      ,LDIM     ,LNNDEL
     &   ,LXPARFL  ,PAREL    ,VD)

C------------------------- Assembles derivative of ATRA w.r.t. VD and 
C------------------------- derivatives of VD w.r.t flow parameter IFP 
C------------------------- Only for constant density.

                  CALL DER_ATRA_FP_VD
     &           (IDIMDQ   ,IODIM    ,LMXNDL   ,NUMEL    ,NUMNP    ,AREA
     &           ,CAUX1    ,DAT_VD   ,DVDP     ,GRDFF    ,KXX      ,LDIM
     &           ,LNNDEL   ,DERC(1,IFP,INEWT))

              END IF !IODENS.EQ.0

C------------------------- Computes derivatives of VD w.r.t. transmissivity 
C------------------------- (only explicit dependence). As it is computed 
C------------------------- simultaneously for all T's, the call to DER_VD_TRA
C------------------------- has to be done only once


              IF (IAUX(1).NE.0) THEN

                  CALL DER_VD_TRA
     &                (AREA     ,BUOYANCY ,CAUX1    ,CFPAREL  ,COORD
     &                ,DAT_VD   ,DENSITY  ,DVDP     ,DTIM     ,EPSFLU
     &                ,FNT      ,GP_COORD ,GRADLOC  ,GRDFF    ,HAUX1
     &                ,HBASE    ,HCALAN   ,HCALIT   ,IDIMDQ   ,IDIMFNT
     &                ,IFLAGS   ,ISSFL    ,INORPAR  ,INTI     ,IOCAP
     &                ,0        ,IODENS   ,IODIM    ,IOFLLI   ,IOFMLF
     &                ,ISOZ     ,IVPAR    ,KXX      ,LDIM     ,LMXNDL
     &                ,LNNDEL   ,LTYPE    ,LXPARFL  ,MAINF    ,MAXPG
     &                ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &                ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NTYPAR
     &                ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR    ,NZTRA
     &                ,PARACD   ,PARC     ,POINTWEIGHT
     &                ,DERC(1,IFP,INEWT)
     ;                ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)

              END IF !IAUX(1).NE.0

C------------------------- Computes derivatives of ATRA w.r.t flow parameters
C------------------------- through the dependence with boundary flow (only 
C------------------------- indirect dependence, i.e., through heads)
C------------------------- Only for constant density.

              IF (IODENS.EQ.0) THEN

                  CALL DERQ_GEN
     & (IAD_S    ,IADN_S   ,IDIMDFLU ,ISSFL    ,IFP      ,NFLAGS
     & ,NPAR     ,NPPNP    ,NUMNP    ,AFLU     ,CAUX1   ,DERC(1,1,INEWT)
     & ,DERH(1,1,IOLD)     ,DERH1    ,DERH2    ,DFLU     ,IBCOD   ,IBTCO
     & ,IFLAGS   ,PARNP    ,CAUDAL   ,IDIMAFLU ,ITYPAFLU ,ITYPDFLU,NUMEL   
     & ,LMXNDL   ,KXX   ,LNNDEL)

              END IF !IODENS.EQ.0
      
C------------------------- Computes derivatives of ATRA w.r.t flow parameters
C------------------------- through the dependence with boundary flow (direct
C------------------------- dependence)

              IF (IAUX(1).NE.0) THEN

                  CALL DERQ_TRA
     &          (AREA     ,BIBI     ,BUOYANCY ,CAUDAL   ,CAUX1
     &          ,CFPAREL  ,COORD    ,DENSITY  ,DERC(1,1,INEWT)
     &          ,DTIM     ,EPSFLU   ,FNT      ,GP_COORD ,GRADLOC
     &          ,HAUX1    ,HBASE    ,HCALAN   ,HCALIT   ,IBCOD
     &          ,IDIMBB   ,IDIMFNT  ,IFLAGS   ,ISOLTR-1 ,INORPAR
     &          ,INTI     ,0        ,IOCTRA   ,IODENS   ,IODIM
     &          ,IOFLLI   ,IOFMLF   ,ISOZ     ,IVPAR    ,KXX
     &          ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPARFL
     &          ,MAXPG    ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &          ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR     ,NPAREL
     &          ,NPPNP    ,NTYPAR   ,NTYPEL   ,NUMEL    ,NUMNP
     &          ,NZONE_PAR,NZPAR    ,NZONE_PAR(1)       ,PARACD
     &          ,PARC     ,PARNP    ,POINTWEIGHT        ,IDIMWGT
     &          ,WGT_PAR  ,IPNT_PAR ,IPOS               ,DERIV)
    

                IAUX(1) = 0

              END IF !IAUX(1).NE.0

              IF (IAUX(2).NE.0) THEN

                  CALL DERQ_STG
     &          (AREA     ,BETAC    ,CAUDAL   ,CAUX1    ,CAUX2
     &          ,CFPAREL  ,CREF     ,DENSITY  ,DENSREF  ,DERC
     &          ,DTIM     ,EPSFLU   ,FNT      ,HAUX1    ,HAUX2
     &          ,HCALAN   ,HCALIT   ,HINI     ,IBCOD    ,IBTCO
     &          ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,ISOLTR-1 ,INEWT
     &          ,INORPAR  ,INTI     ,0        ,IPAR_DIR(9)
     &          ,IODENS   ,IOFLLI   ,IOFMLF   ,IOVRWC   ,IVPAR
     &          ,KXX      ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPARFL(1,2)
     &          ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPPNP
     &          ,1        ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &          ,NZPAR    ,PARACD   ,PARC     ,PARNP    ,IDIMWGT
     &          ,WGT_PAR  ,IPNT_PAR ,IPOS    ,DERIV)


                  IF (IOVRWC.EQ.0) THEN
                      IAUX(2) = 0
                  END IF !IOVRWC.EQ.0
                  
              END IF !IAUX(2).NE.0

              IF (IAUX(3).NE.0) THEN


                  CALL DERQ_ARR
     &                (AREA     ,BETAC    ,CAUDAL   ,CAUX1    ,CFPAREL
     &                ,CREF     ,DENSREF  ,DERC     ,DTIM     ,EPSFLU
     &                ,FNT      ,HCALAN   ,HCALIT   ,IBCOD    ,IBTCO
     &                ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,ISOLTR-1 ,INEWT
     &                ,INORPAR  ,INTI     ,IOCAP    ,IODENS   ,IOFLLI
     &                ,IOFMLF   ,IVPAR    ,KXX      ,LMXNDL   ,LNNDEL
     &                ,LXPARFL(1,3+ISSFL) ,NFLAGS   ,NFNL     ,NFNLPAR
     &                ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR
     &                ,NPAREL   ,NPPEL    ,NPPNP    ,1        ,NTYPAR
     &                ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD
     &                ,PARC     ,PAREL    ,PARNP    ,IDIMWGT  ,WGT_PAR
     &                ,IPNT_PAR ,IPOS    ,DERIV)

                  IAUX(3) = 0

              END IF !IAUX(3).NE.0

              IF (IAUX(4).NE.0 .OR. IAUX(5).NE.0 .OR. IAUX(6).NE.0) THEN

                  CALL DERQ_NUD
     &               (BETAC    ,CAUDAL   ,CAUX1    ,CCALAN   ,CCALIT
     &               ,CFPARNP  ,CREF     ,DENSREF  ,DERC     ,DTIM
     &               ,EPSFLU   ,FNT      ,HAUX1    ,IBCOD    ,IBTCO
     &               ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,ISOLTR-1 ,INEWT
     &               ,INORPAR  ,INTI     ,0        ,IODENS   ,IOFLLI
     &               ,IOFMLF   ,IVPAR    ,IXPARFL  ,KXX      ,LMXNDL
     &               ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &               ,NFTPAR   ,NINT     ,NPAR     ,NPARNP   ,NPPNP
     7               ,1        ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &               ,NZPAR    ,PARACD   ,PARC     ,PARNP    ,IDIMWGT
     &               ,WGT_PAR  ,IPNT_PAR ,IPOS    ,DERIV)

                  IAUX(4) = 0
                  IAUX(5) = 0
                  IAUX(6) = 0

              END IF !IAUX(4).NE.0 .OR. IAUX(5).NE.0 .OR. IAUX(6).NE.0

              IF (IOVRWC.GE.1) THEN

                  IF (IODENS.EQ.0) THEN

                      CALL DER_WATVOL_GEN
     &                    (IDIMBB   ,IDIMDERC ,IDIMDERH ,INEWT
     &                    ,IOVRWC   ,IFP      ,LMXNDL   ,NPAR
     &                    ,NPARF    ,NPPEL    ,NUMEL    ,NUMNP
     &                    ,THETAT   ,AREA     ,BIBI     ,CAUX1
     &                    ,CAUX2    ,DERC     ,DERH     ,KXX
     &                    ,LDIM     ,LNNDEL   ,PAREL    ,INEW
     &                    ,IOLD)

                  END IF !IODENS.EQ.0

                  IF (IAUX(2).NE.0) THEN
                      
                      CALL DER_WTV_STG
     &                (AREA     ,BETAC    ,BIBI     ,CAUX1    ,CAUX2
     &                ,CFPAREL  ,CREF     ,DENSITY  ,DENSREF  ,DERC
     &                ,DTIM     ,EPSFLU   ,HCALAN   ,HCALIT   ,HINI
     &                ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,ISOLTR-1 ,INEWT
     &                ,INORPAR  ,INTI     ,0        ,IODENS   ,IOFLLI
     &                ,IOFMLF   ,IOVRWC   ,IVPAR    ,KXX      ,LDIM
     &                ,LMXNDL   ,LNNDEL   ,LXPAREL  ,NFLAGS   ,NFNL
     &                ,NINT     ,NPAR     ,NPAREL   ,NPPEL    ,1
     &                ,NTYPAR   ,NUMEL    ,NUMNP    ,NZPAR    ,PARC
     &                ,PAREL    ,THETAT   ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR)

 
                  IAUX(2)=0

                  END IF !IAUX(2).NE.0

              END IF !IOVRWC.GE.1

          END DO ! IFP = 1, NPARF

      END IF     ! IOINV.EQ.3

C------------------------- Derivatives w.r.t areal recharge, explicit dependence
C------------------------- in transport equation

      IF (IOLG_PAR(3,2).NE.0) THEN

          CALL DERARR_TPT         !(inarr)
     &        (DTIM     ,EPSFLU   ,IDIMDERC ,IDIMFNT  ,ISSFL+3  ,ISSFL
     &        ,INEWT    ,INTI     ,IOCAP    ,IOFLLI   ,IOFMLF   ,LMXNDL
     &        ,NFLAGS   ,NFNL     ,NINT     ,NPAR     ,NPAREL   ,NPPEL
     &        ,NTYPAR   ,NUMEL    ,NUMNP    ,NZPAR    ,AREA     ,BETAC
     &        ,CAUX1    ,CFPAREL  ,DERC     ,FNT      ,HCALIT   ,HCALAN
     &        ,IFLAGS   ,INORPAR  ,IVPAR    ,KXX      ,LDIM     ,LNNDEL
     &        ,LXPARFL  ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR
     &        ,NZONE_PAR,PARACD   ,PARC     ,PAREL    ,IODENS   ,DENSREF
     &        ,CREF     ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)

      END IF ! IOLG_PAR(3,2).NE.0 (Derivatives w.r.t recharge)


C------------------------- Derivatives w.r.t transport parameters

C------------------------- Derivatives with respect to dispersivity

      IF (IOLG_PAR(7,2).NE.0) THEN

C------------------------- Derivatives with respect to longitudinal disp.

          CALL DERDSL
     &        (DTIM     ,EPSTRA   ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IDIMQ
     &        ,ISOLTR-1 ,INEWT    ,INTI     ,IOCAP    ,IODENS   ,IOFMLT
     &        ,IOTRLI   ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT     ,NPAR
     &        ,NPAREL   ,NTYPAR   ,NUMEL    ,NUMNP    ,NZPAR    ,BIBI
     &        ,CAUX1    ,CCALIT   ,CCALAN   ,CFPAREL  ,DENSITY  ,DERC
     &        ,FNT      ,IFLAGS   ,INORPAR  ,IVPAR    ,KXX      ,LDIM
     &        ,LNNDEL   ,LXPAREL  ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR
     &        ,NZONE_PAR,PARACD   ,PARC     ,QXYZ     ,IDIMWGT  ,WGT_PAR
     &        ,IPNT_PAR ,IPOS     ,DERIV)

     
      END IF !IOLG_PAR(7,2).NE.0

C------------------------- Derivatives with respect to transversal disp.

      IF (IOLG_PAR(8,2).NE.0) THEN

          CALL DERDST
     &       (DTIM     ,EPSTRA   ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IDIMQ
     &       ,ISOLTR-1 ,INEWT    ,INTI     ,IOCAP    ,IODENS   ,IOFMLT
     &       ,IOTRLI   ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT     ,NPAR
     &       ,NPAREL   ,NTYPAR   ,NUMEL    ,NUMNP    ,NZPAR    ,BIBI
     &       ,CAUX1    ,CCALIT   ,CCALAN   ,CFPAREL  ,DENSITY  ,DERC
     &       ,FNT      ,IFLAGS   ,INORPAR  ,IVPAR    ,KXX      ,LDIM
     &       ,LNNDEL   ,LXPAREL  ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR
     &       ,NZONE_PAR,PARACD   ,PARC     ,QXYZ     ,IDIMWGT  ,WGT_PAR
     &       ,IPNT_PAR ,IPOS     ,DERIV)

      END IF !IOLG_PAR(8,2).NE.0

C------------------------- Derivatives with respect to molecular diffusion

      IF (IOLG_PAR(9,2).NE.0) THEN

          CALL DERDFM
     &        (BIBI     ,CAUX1    ,CCALAN   ,CCALIT   ,CFPAREL
     &        ,DENSITY  ,DERC     ,DTIM     ,EPSTRA   ,FNT
     &        ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,ISOLTR-1
     &        ,INEWT    ,INORPAR  ,INTI     ,0        ,IODENS
     &        ,IOFMLT   ,IOTRLI   ,IOVRWC   ,ITPTVAR  ,IVPAR
     &        ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LXPAREL
     &        ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &        ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPPEL
     &        ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &        ,PARACD   ,PARC     ,PAREL    ,WATVOL   ,WSPECHEAT
     ;        ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR)


      END IF !IOLG_PAR(9,2).NE.0

C------------------------- Derivatives with respect to porosity

      IF (IOLG_PAR(10,2).NE.0) THEN

          CALL DERPOR
     &        (ACTH     ,AREA     ,BETAC    ,BIBI     ,CAUX1
     &        ,CAUX2    ,CCALAN   ,CCALIT   ,CFPAREL  ,CREF
     &        ,DENSITY  ,DENSREF  ,DERC     ,DERCS    ,DTIM
     &        ,EPSTRA   ,FNT      ,ICMPDERS ,IDIMBB   ,IDIMDERC
     &        ,IDIMFNT  ,IFLAGS   ,ISOLTR-1 ,INEWT    ,INORPAR
     &        ,INTI     ,IOCAP    ,IODENS   ,IOFMLT   ,IOTRLI
     &        ,IOVRWC   ,ITPTVAR  ,IVPAR    ,KXX      ,LDIM
     &        ,LMXNDL   ,LNNDEL   ,LXPAREL  ,NFLAGS   ,NFNL
     &        ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT
     &        ,NPAR     ,NPAREL   ,NPPEL    ,NTYPAR   ,NUMEL
     7        ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &        ,PAREL    ,THETAT   ,WSPECHEAT,WTHERMCON,IDIMWGT
     &        ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)

	    IF (IODENS.EQ.1) THEN

              CALL DERQ_POR
     &            (ACTH     ,AREA     ,BETAC    ,CAUDAL   ,CAUX1
     &            ,CAUX2    ,CCALAN   ,CCALIT   ,CFPAREL  ,CREF
     &            ,DENSITY  ,DENSREF  ,DERC     ,DTIM     ,EPSTRA
     &            ,FNT      ,IBCOD    ,IBTCO    ,IDIMDERC ,IDIMFNT
     &            ,IFLAGS   ,ISOLTR-1 ,INEWT    ,INORPAR  ,INTI
     &            ,IOCAP    ,IOFMLT   ,IOTRLI   ,IOVRWC   ,IVPAR
     &            ,KXX      ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPAREL(1,7)
     &            ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &            ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPPNP
     &            ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &            ,PARACD   ,PARC     ,PARNP    ,IDIMWGT  ,WGT_PAR
     &            ,IPNT_PAR ,IPOS     ,DERIV)

          END IF !IODENS.EQ.1

      END IF !IOLG_PAR(10,2).NE.0

C------------------------- Derivatives with respect to first order decay 

      IF (IOLG_PAR(11,2).NE.0) THEN

          CALL DERFOD
     &        (ACTH     ,AREA     ,BETAC    ,CAUX1    ,CCALAN
     &        ,CCALIT   ,CFPAREL  ,CREF     ,DENSREF  ,DERC
     &        ,DERCS    ,DTIM     ,EPSTRA   ,FNT      ,ICMPDERS
     &        ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,ISOLTR-1 ,INEWT
     &        ,INORPAR  ,INTI     ,0        ,IODENS   ,IOFMLT
     &        ,IOTRLI   ,IOVRWC   ,IVPAR    ,KXX      ,LDIM
     &        ,LMXNDL   ,LNNDEL   ,LXPAREL  ,NFLAGS   ,NFNL
     &        ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT
     &        ,NPAR     ,NPAREL   ,NPPEL    ,NTYPAR   ,NUMEL
     &        ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &        ,PAREL    ,WATVOL   ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR
     &        ,IPOS     ,DERIV)

      END IF !IOLG_PAR(11,2).NE.0

C------------------------- Derivatives with respect to retardation

      IF (IOLG_PAR(12,2).NE.0) THEN

          CALL DERCRD
     &          (ACTH     ,AREA     ,BETAC    ,CAUX1    ,CAUX2
     &          ,CCALAN   ,CCALIT   ,CFPAREL  ,CREF     ,DENSREF
     &          ,DERC     ,DERCS    ,DTIM     ,EPSTRA   ,FNT
     &          ,ICMPDERS ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,ISOLTR-1
     &          ,INEWT    ,INORPAR  ,INTI     ,IOCAP    ,IODENS
     &          ,IOFMLT   ,IOTRLI   ,IOVRWC   ,ITPTVAR  ,IVPAR
     &          ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LXPAREL
     &          ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPPEL
     &          ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &          ,PARACD   ,PARC     ,PAREL    ,THETAT   ,WSPECHEAT
     &          ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)


      END IF !IOLG_PAR(12,2).NE.0

C------------------------- Derivatives with respect to concentration leakage

      IF (IOLG_PAR(18,2).NE.0) THEN

          CALL DERCLK
     &        (CAUX1    ,CCALAN   ,CCALIT   ,CFPARNP  ,DERC
     &        ,DTIM     ,EPSTRA   ,FNT      ,IDIMDERC ,IDIMFNT
     &        ,IFLAGS   ,ISOLTR-1 ,INEWT    ,INORPAR  ,INTI
     &        ,IOCAP    ,IOFMLT   ,IOTRLI   ,IXPARNP  ,IVPAR
     &        ,KXX      ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLPAR
     &        ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR
     &        ,NPARNP   ,NPPNP    ,NTYPAR   ,NUMEL    ,NUMNP
     &        ,NZONE_PAR,NZPAR    ,PARACD   ,PARC     ,PARNP
     ;        ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV
     &        ,ITPTVAR  ,WSPECHEAT)

      END IF !IOLG_PAR(18,2).NE.0

C------------------------- Derivatives with respect to external concentration

C------------------------- Derivatives in nodes with prescribed concentration
C------------------------- must be set to zero.

      DO I=1,NUMNP

          IF (IBTCO(I).EQ.1) THEN

              DERC(I,1:NPAR,INEWT) = 0D0

          END IF !IBTCOD(I).EQ.1

      END DO !I=1,NUMNP

      IF (IOLG_PAR(13,2).NE.0) THEN

C------------------------- INCON = 6+ISOLTR
C------------------------- INARR = ISSFL+3 

          CALL DERCOE
     &        (AREA     ,BETAC    ,CAUDAL   ,CAUX1    ,CCALAN
     &        ,CCALIT   ,CFPAREL  ,CFPARNP  ,CREF     ,DENSREF
     &        ,DERC     ,DTIM     ,EPSTRA   ,FNT      ,IBCOD
     &        ,IBTCO    ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,ISSFL+3
     &        ,6+ISOLTR ,ISOLTR-1 ,INEWT    ,INORPAR  ,INTI
     &        ,IOCAP    ,IODENS   ,IOFMLT   ,IOTRLI   ,IVPAR
     &        ,IXPARNP  ,KXX      ,LMXNDL   ,LNNDEL   ,LXPAREL
     &        ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &        ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPARNP
     &        ,NPPEL    ,NPPNP    ,1        ,NTYPAR   ,NUMEL
     &        ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &        ,PAREL    ,PARNP    ,THETAT   ,TINC     ,TINTERVOBS
     ;        ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV
     &        ,ITPTVAR  ,WSPECHEAT)


       END IF !IOLG_PAR(13,2).NE.0

C------------------------- Writes the RHS of the inverse problem (sensitivity
C------------------------- equations)

       IF (IFLAGS(25).GT.0 .AND. IODENS.NE.1) THEN

          WRITE(746,*) 'DERC BEFORE INTI = ',INTI
          DO I=1,NUMNP
              WRITE(746,strFmt1) I,(DERC(I,J,INEWT),J=1,NPAR)
          END DO !I=1,NUMNP

       END IF !IFLAGS(25).GT.0

      IF (IFLAGS(10).NE.0) THEN

          DO NP=1,NPAR

              WRITE(MAINF,*)'RHS TPT PARAM. ESTIM. NUM.',NP,'INTI=',INTI

                  DO I=1,NUMNP
                      WRITE(MAINF,'(2I5,G20.13)') NP,I,DERC(I,NP,INEWT)
                  END DO !I=1,NUMNP

          END DO !NP=1,NPAR

      END IF !IFLAGS(10).NE.0

C------------------------- Computed concentration derivatives at current time

      IF (IODENS.EQ.0) THEN

          DO NP=1,NPAR

              CALL SOLVE
     &            (IADSC_ROWS    ,IADSC_COLS    ,NUMNP         ,IDESC         
     &            ,IPAR_DIR(21)  ,IDIMWORK      ,IDIRECT       ,1        
     &            ,INTI          ,IPAR_DIR(18)  ,IPAR_DIR(22)
     &            ,0             ,ITERM         ,IPAR_DIR(15)  ,MAINF
     &            ,NBAND1        ,IPAR_DIR(24)
     &            ,IPAR_DIR(20)  ,IPAR_DIR(17)  ,IPAR_DIR(16)
     &            ,IPAR_DIR(36)  ,IPAR_DIR(37)  ,IPAR_DIR(38)  ,ATRADSC
     &            ,ATRADSCF      ,DERC(1,NP,INEWT)             ,IAD_S
     &            ,IADD_S        ,IADN_S        ,IAFD_S        ,IAFDD_S
     &            ,IAFDN_S       ,WORK          ,SOLUTION)

              DERC(1:NUMNP,NP,INEWT) = SOLUTION(1:NUMNP)

          END DO !NP=1,NPAR

          IF (IFLAGS(25).GT.0) THEN

              WRITE(746,*) 'DERC AFTER INTI = ',INTI
              DO I=1,NUMNP
                  WRITE(746,strFmt1) I,(DERC(I,J,INEWT),J=1,NPAR)
              END DO !I=1,NUMNP

c    1             FORMAT(I5,<NPAR>F15.10)

          END IF !IFLAGS(25).GT.0

          IF (IFLAGS(4).GT.0) THEN

              WRITE(750,*) 'DERC AFTER INTI = ',INTI

              DO I=1,NUMNP

                  WRITE(750,strFmt1) I,(DERC(I,J,INEWT),J=1,NPAR)

              END DO !I=1,NUMNP

          END IF !IFLAGS(4).GT.0

C------------------------- Writes the derivative of nodal conc. with respect
C------------------------- to each parameter

          IF (IFLAGS(13).NE.0) THEN

              DO NP=1,NPAR

                  WRITE(MAINF,10) NP,INTI
   10             FORMAT(' DERIV. PARAM. ESTIM. NUM.',I5,'INTI=',I5)

                  DO I=1,NUMNP

                      WRITE(MAINF,20) NP,I,DERC(I,NP,INEWT)
   20                 FORMAT(2I5,G20.13)

                  END DO !I=1,NUMNP

              END DO !NP=1,NPAR

          END IF !IFLAGS(13).NE.0

C------------------------- The contribution of derivatives of the RHS of the 
C------------------------- "son" w.r.t the parameters of the "father" must be
C------------------------- updated

          IF (IPROB.NE.NPBTP .AND. IORDCH.NE.0)  THEN

              CALL MODIF_JAC_C_RDCH
     &            (IDIMDERC ,INEWT    ,IOLDT    ,LMXNDL   ,NPAR
     &            ,NPPEL    ,NUMEL    ,NUMNP    ,THETAT   ,AREA
     &            ,DERC     ,DERCS    ,KXX      ,LNNDEL   ,PAREL
     &            ,WATVOL)

          END IF !IPROB.NE.NPBTP .AND. IORDCH.NE.0

      END IF !IODENS.EQ,0

	  IF(IFLAGS(3).EQ.1) CALL IO_SUB('JAC_C',1)
      END SUBROUTINE JAC_C
