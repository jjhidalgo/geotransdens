      SUBROUTINE JAC_H
     &          (ACTH    ,AFLU     ,AFLUDSC  ,AFLUDSCF ,ALFA
     &         ,AREA     ,BETAC    ,BIBI     ,BUOYANCY
     &         ,CAUX1    ,CAUX2    ,CCALAN   ,CCALIT   ,CFLU
     &         ,CFPAREL  ,CFPARNP  ,COORD    ,CREF     ,DBFLUDFLU
     &         ,DBFLUDTRA,DENSITY  ,DENSREF  ,DERC     ,DERH
     &         ,DFLU     ,DFLUDFLU ,DFLUDTRA ,DTIM     ,EPSFLU
     &         ,EPSTRA   ,FNT      ,GP_COORD ,GRADLOC  ,HAUX1
     &         ,HAUX2    ,HBASE    ,HCALAN   ,HCALIT   ,IAD_S
     &         ,IADD_S   ,IADN_S   ,IAFLUDSC_COLS      ,IAFLUDSC_ROWS
     &         ,IAFD_S   ,IAFDD_S  ,IAFDN_S  ,IBCOD
     &         ,IDESC    ,IDIMAFLU ,IDIMBB   ,IDIMCFLU
     &         ,IDIMDENS ,IDIMDERC ,IDIMDERH ,IDIMDFLU ,IDIMFNT
     &         ,IDIMWORK ,IDIRECT  ,IFLAGS   ,INALF    ,INARR
     &         ,INCHP    ,INDSSTR  ,INEW     ,INEWT    ,INORPAR
     &         ,INQQP    ,INTI     ,IOCAP    ,IODENS   ,IODIM
     &         ,IOFLLI   ,IOFMLF   ,IOFMLT   ,IOINV
     &         ,IOLD     ,IOLDT    ,IOTRLI   ,IOVRWC   ,IPAR_DIR
     &         ,ISOZ     ,ISPARSE  ,ISYMETRIC,ITYPAFLU ,ITYPCFLU
     &         ,ITYPDFLU ,IVPAR    ,IXPARNP  ,KXX      ,LDIM
     &         ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPAREL  ,MAINF
     &         ,MAXNB    ,MAXNBF   ,MAXPG    ,NBAND1
     &         ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &         ,NFNLTIP  ,NFTPAR   ,NINT     ,NPARALG  ,NPAREL
     &         ,NPARF    ,NPARNP   ,NPART    ,NPPEL    ,NPPNP
     &         ,NTYPAR   ,NTYPEL   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &         ,NZPAR    ,NZTRA    ,PARACD   ,PARC     ,PAREL
     &         ,PARNP    ,POINTWEIGHT        ,SOLUTION ,THETAF
     &         ,THETAT   ,TINC     ,TINTERVOBS         ,WORK
     ;         ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV   ,NPAR)

********************************************************************************
*
* PURPOSE
*
*     Computes the RHS of the flow inverse problem
*
* DESCRIPTION
*
*     Computes the RHS of the flow inverse problem, ie., the RHS of heads
*     (pressures) with respect to all flow estimated parameters
*     The RHS is stored in DERH(.,.,INEW) and at the end it is also stored 
*     the actual derivatives at time k+1. DERH(.,.,IOLD) stores the previous 
*     computed derivatives at time k.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AFLU                   Matrix of finite elements equations for flow problem  
*                         No boundary conditions are included on it.            
*  AFLUDSC                Coefficient matrix of flow system (2.15), most        
*                         often in its decomposed form.                         
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  CFPARNP                Array containing node coefinient of node j            
*                         corresponding to INpar index parameter zone.          
*  DERH                   Nodal head derivatives with respect to estimated      
*                         flow parameters.                                      
*  HAUX1                  Array containing HEADS, ponderated by THETAF          
*                         time factor                                           
*  HAUX2                  Array containing diference of HEADS in two            
*                         consecutives times.                                   
*  HBASE                  Bottom level of the aquifer
*  HCAL                   Computed heads at every node                          
*  HCALAN                 Head level at previous time                           
*  IBCOD                  Flow boundary condition index                         
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  IXPARNP                Array containing zone number of each node j,          
*                         corresponding to INpar index parameter zone.          
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
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
*
* INTERNAL VARIABLES: ARRAYS
*
*  DERIV                  Derivatives of the current parameter with respect
*                         to its zonal value plus its generic parameters
*  IPOS                   Location order for estimated parameters
*  XPARAM                 Array used to store auxiliar values for the
*                         computation of the non-linear function
*
* EXTERNAL VARIABLES: SCALARS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  BIBI                   Array containing the product of interpolation         
*                         functions gradient, for a given element               
*  CNST                   Interpolation functions gradient for a given element  
*                         nodes                                                 
*  DERADFLU               Contains the derivatives of AFLU and DFLU matrices    
*                         with respect to head. This term appears on the left   
*                         hand side of both flow direct and inverse problems,   
*                         as well as on the right hand side of flow inverse     
*                         problem.                                              
*  DFLU                   Matrix of finite elements equations for flow          
*                         problem related to storage term.                      
*  DTIM                   Value to compute the time function value at the
*                         current time (k+theta)
*  EPSFLU                 Time weighting parameter for nonlinear flow problems  
*  FNT                    Array containing time functions values                
*  GRAV                   Gravity array direction                               
*  GRAVEL                 Gravity array direction at every element
*  GRDFF                  Array containing the product between interpolation    
*                         functions integrals and interp. functions gradient    
*  IDIMAFLU               Used to dimension array AFLU                          
*  IDIMBB                 Used to dimension array BIBI                          
*  IDIMDADFLU             Second dimension of array DERADFLU
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INALF                  Index for leakage                                     
*                         in array variables (IXPARNP and CFPARNP)              
*  INARR                  Index for areal recharge                              
*                         in array variables (LXPAREL and CFPAREL)              
*  INCHP                  Index for prescribed head                             
*                         in array variables (IXPARNP and CFPARNP)              
*  INDSSTR                Current problem state. If 1, transient state,         
*                         if 0, steady-state                                    
*  INEW                   Index to locate the RHS of the derivatives of head
*                         at current time with respect to parameters inside
*                         array DERH
*  INQQP                  Index for prescribed flow                             
*                         in array variables (IXPARNP and CFPARNP)              
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOCNSF                 Scheme for storage term in flow problem               
*  IODIM                  Maximum dimension of any element included             
*                         in the problem                                        
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOFMLF                 Flow Formulation number                               
*  IOLD                   Index to locate the derivatives of head with 
*                         respect to parameters at the previous time step
*  IOPRHED                Indicates whether the flow state variable state is    
*                         preasure (set to 1) or head (set to 0)                
*  IOWRITE                Array containing all output options                   
*  IPAR_DIR               Array containing all integer direct problem           
*                         parameters                                            
*  ISOZ                   Anisotropy of every transmissivity zone               
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NBAND                  Half Bandwith (maximum difference between the         
*                         numbers of two nodes belonging to the same element)   
*  NBAND1                 Used to dimension. It is equal to NBAND+1             
*  NFLAGS                 Maximum number of allowed flags                       
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NINT                   Number of observation times                           
*  NPARALG                Maximum number of algorithm parameters, including     
*                         minimization and simulation (used for dimensioning)   
*  NPAREL                 Number of element parameters in current problem       
*  NPARF                  Number of transient parameters to be estimated        
*  NPARNP                 Number of nodal parameters in current problem         
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NWRITE                 Number of output options (used for dimensioning)      
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  THETAF                 Time weighting parameter for flow problems            
*  TINC                   Current time increment                                
*  WORK                   Workspace array used in some subroutines.             
*
* INTERNAL VARIABLES: SCALARS
*
*  DTH                    1-THETAF
*  FACTOR                 1-EPSFLU
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  DERARR                                                                       
*  DERSTG                                                                       
*  DERTRA                                                                       
*  DER_PARAM                                                                    
*  ENS_IND_T10                                                                  
*  ENS_IND_T9                                                                   
*  FLOW_GRAV_IN                                                                 
*  LEQT1B                                                                       
*  LUELPB                                                                       
*  RHS_IN_CHP                                                                   
*  RHS_IN_CHP_TRCNS                                                             
*  RHS_IN_CHP_TRLUM                                                             
*  RHS_NOLI_IN                                                                  
*  RT_ERROR                                                                     
*
* HISTORY
*
*     SCR      6-1997     First coding
*     AMS     11-1998     Full modification
*     AMS      2-1999     Revision of derivatives in the nonlinear case.
*                         Inclusion of derivatives with respect to nl leakage
*
********************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IAFLUDSC_COLS      ,IAFLUDSC_ROWS      ,IDESC
     &          ,IDIMAFLU ,IDIMBB   ,IDIMCFLU ,IDIMDENS ,IDIMDERC
     &          ,IDIMDERH ,IDIMDFLU ,IDIMFNT  ,IDIMWORK ,IDIRECT
     &          ,INALF    ,INARR    ,INCHP    ,INDSSTR  ,INEW
     &          ,INEWT    ,INQQP    ,INTI     ,IOCAP    ,IODENS
     &          ,IODIM    ,IOFLLI   ,IOFMLF   ,IOFMLT
     &          ,IOINV    ,IOLD     ,IOLDT    ,IOTRLI   ,IOVRWC
     &          ,ISPARSE  ,ISYMETRIC,ITYPAFLU ,ITYPCFLU ,ITYPDFLU
     &          ,LMXNDL   ,MAINF    ,MAXNB    ,MAXNBF   ,MAXPG
     &          ,NBAND1   ,NFLAGS   ,NFNL     ,NPAR
     &          ,NINT     ,NPARALG  ,NPAREL   ,NPARF    ,NPARNP
     &          ,NPART    ,NPPEL    ,NPPNP    ,NTYPAR   ,NTYPEL
     &          ,NUMEL    ,NUMNP    ,NZPAR    ,NZTRA
     ;          ,IDIMWGT  ,IPNT_PAR ,IPOS(NPAR)
       

      REAL*8::BETAC    ,CREF     ,DENSREF  ,DTIM     ,DTIMAUX
     &       ,EPSFLU   ,EPSTRA   ,THETAF   ,THETAT   ,TINC
     &       ,TINTERVOBS         ,DERIV(NPAR)    ,WGT_PAR



      INTEGER*4::IAD_S(MAXNB,NUMNP),IADD_S(NUMNP),IADN_S(NUMNP)
     &          ,IAFD_S(MAXNBF,NUMNP),IAFDD_S(NUMNP),IAFDN_S(NUMNP)
     &          ,IBCOD(NUMNP)
     &          ,IFLAGS(NFLAGS),INORPAR(NTYPAR),IPAR_DIR(NPARALG)
     &          ,ISOZ(NZTRA),IVPAR(NZPAR),IXPARNP(NUMNP,NPARNP)
     &          ,KXX(LMXNDL,NUMEL)
     &          ,LDIM(NUMEL),LNNDEL(NUMEL),LTYPE(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR),NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL),NFTPAR(NZPAR),NZONE_PAR(NTYPAR)

      REAL*8::ACTH(NUMEL),AFLU(NUMEL,IDIMAFLU)
     &       ,AFLUDSC(IAFLUDSC_ROWS,IAFLUDSC_COLS)
     &       ,AFLUDSCF(MAXNBF,NUMNP),ALFA(NUMNP),AREA(NUMEL)
     &       ,BIBI(IDIMBB,NUMEL)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL),CAUX1(NUMNP),CAUX2(NUMNP)
     &       ,CCALAN(NUMNP),CCALIT(NUMNP),CFLU(NUMEL,IDIMCFLU)
     &       ,CFPAREL(NUMEL,NPAREL),CFPARNP(NUMNP,NPARNP),COORD(NUMNP,3)
     &       ,DBFLUDFLU(NUMNP),DBFLUDTRA(NUMNP),DENSITY(IDIMDENS)
     &       ,DERC(NUMNP,NPART,IDIMDERC),DERH(NUMNP,NPARF,IDIMDERH)
     &       ,DFLU(NUMEL,IDIMDFLU),DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL),FNT(IDIMFNT,NINT)
     &       ,GP_COORD(6,8,IODIM),GRADLOC(IODIM,LMXNDL,MAXPG)
     &       ,HAUX1(NUMNP),HAUX2(NUMNP),HBASE(NUMEL)
     &       ,HCALAN(NUMNP),HCALIT(NUMNP)
     &       ,PARACD(3,NFNL),PARC(NZPAR),PAREL(NUMEL,NPPEL)
     &       ,PARNP(NUMNP,NPPNP),POINTWEIGHT(MAXPG,NTYPEL)
     &       ,SOLUTION(NUMNP),WORK(IDIMWORK)

C------------------------- Internal

      INTEGER*4::I      ,IB     ,IP     ,ITERM
     &          ,NP     ,j

      REAL*8::FACTOR

      CHARACTER::strFmt1*20

C------------------------- First executable statement

      IF(IFLAGS(3).EQ.1) CALL IO_SUB('JAC_H',0)
	  
      strFmt1 = ''
      WRITE (strFmt1,*) NPAR
      strFmt1 = '(I5,'//Trim(AdjustL(strFmt1))//'F15.10)'


      DERH(:,:,INEW) = 0D0

      IF (IODENS.EQ.1) THEN
          DERC(:,:,INEWT) = 0D0
      END IF !IODENS.EQ.1


C------------------------- Adds the nonlinear part to inverse problem
C------------------------- RHS, related to the previous time step.


      IF (IOFLLI.NE.0 .OR. IODENS.EQ.1) THEN

          FACTOR = (EPSFLU - 1D0)/EPSFLU

          CALL RHS_NOLI_IN 
     &        (DFLUDFLU ,DBFLUDFLU,FACTOR   ,1        ,KXX
     &        ,LMXNDL   ,LNNDEL   ,NPARF    ,NUMEL    ,NUMNP
     &        ,DERH(1,1,INEW)     ,DERH(1,1,IOLD))


          IF (IODENS.EQ.1) THEN

              CALL RHS_NOLI_IN
     &          (DFLUDTRA ,DBFLUDTRA,FACTOR   ,1        ,KXX
     &          ,LMXNDL   ,LNNDEL   ,NPARF    ,NUMEL    ,NUMNP
     &          ,DERH(1,1,INEW)     ,DERC(1,1,IOLDT))



          END IF !IODENS.EQ.1

      END IF !IOFLLI.NE.0 .OR. IODENS.EQ.1

C------------------------- Adds the contribution of previous time step to 
C------------------------- inverse problem  RHS (DERH(.,.,INEW)


      IF (INDSSTR.NE.0) THEN

          DO IP=1,NPARF

C------------------------- (THETAF-1)*AFLU*DERH(OLD)

              CALL PROD_MAT_VEC
     &            (THETAF-1D0,IAD_S   ,IADN_S   
     &            ,IDIMAFLU ,NUMEL    ,NUMNP        ,1
     &            ,ITYPAFLU ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &            ,KXX      ,LNNDEL   ,AFLU     ,DERH(1,IP,INEW)
     &            ,DERH(1,IP,IOLD))

C------------------------- 1/dt *DFLU*DERH(OLD) (only transient flow)

              IF (INDSSTR.EQ.1) THEN

                  CALL PROD_MAT_VEC
     &                (1D0/TINC  ,IAD_S   ,IADN_S
     &                ,IDIMDFLU ,NUMEL    ,NUMNP        ,1
     &                ,ITYPDFLU ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &                ,KXX      ,LNNDEL   ,DFLU     ,DERH(1,IP,INEW)
     &                ,DERH(1,IP,IOLD))

              END IF !INDSSTR.EQ.1

C------------------------- Leakage

              IF (NZONE_PAR(6).NE.0 .AND. IODENS.EQ.0) THEN

C------------------------- Contributrion of leakage is added to RHS.
C------------------------- This contribution is equal to (TH-1)*RHO*LEAK*DERH(K)
C------------------------- When received ALFA storages ALFA = TH*RHO*LEAK
C------------------------- So, it must be multiplied by FACTOR=(TH-1)/TH
C------------------------- to obtain ALFA = (TH-1)*RHO*LEAK_CF
C------------------------- This contribution is already taken into account by
C------------------------- DBFLUDLU when solving by Newton method.

                  FACTOR = (THETAF-1D0)/THETAF
                  ALFA(1:NUMNP) = FACTOR*ALFA(1:NUMNP)
                  

                  DERH(1:NUMNP,IP,INEW) = DERH(1:NUMNP,IP,INEW)
     &                             + ALFA(1:NUMNP)*DERH(1:NUMNP,IP,IOLD)

              END IF !(NZONE_PAR(6).NE.0)


              IF (IODENS.EQ.1) THEN

C------------------------- 1/dt *CFLU*DERC(OLD)

                  CALL PROD_MAT_VEC
     &                (1D0/TINC ,IAD_S   ,IADN_S   
     &                ,IDIMCFLU ,NUMEL    ,NUMNP    ,1
     &                ,ITYPCFLU ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &                ,KXX      ,LNNDEL   ,CFLU     ,DERH(1,IP,INEW)
     &                ,DERC(1,IP,IOLDT))


              END IF !IODENS.EQ.1

          END DO !IP=1,NPARF

      END IF !(INDSSTR.NE.0)

C------------------------- Derivatives with respect to transmissivity

          CALL DERTRA
     &          (AREA     ,BIBI     ,BUOYANCY ,CFPAREL  ,COORD
     &          ,DENSITY  ,DERH     ,DTIM     ,EPSFLU   ,FNT
     &          ,GP_COORD ,GRADLOC  ,HAUX1    ,HBASE    ,HCALAN
     &          ,HCALIT   ,IDIMBB   ,IDIMDENS ,IDIMDERH ,IDIMFNT
     &          ,IFLAGS   ,INDSSTR  ,INEW     ,INORPAR  ,INTI
     &          ,IOCAP    ,2        ,IODENS   ,IODIM    ,IOFLLI
     &          ,IOFMLF   ,ISOZ     ,IVPAR    ,KXX      ,LDIM
     &          ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPAREL  ,MAXPG
     &          ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NINT     ,NPAREL   ,NPARF    ,NTYPAR
     &          ,NTYPEL   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &          ,NZTRA    ,PARACD   ,PARC     ,POINTWEIGHT
     &          ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,NPAR     ,IPOS
     &          ,DERIV)

C------------------------- Derivatives with respect to areal recharge 

      IF (NZONE_PAR(3).NE.0) THEN

          CALL DERARR
     &        (AREA     ,BETAC    ,CFPAREL  ,CREF     ,DENSREF
     &        ,DERH     ,DTIM     ,EPSFLU   ,FNT      ,HCALIT
     &        ,HCALAN   ,IDIMDERH ,IDIMFNT  ,IFLAGS   ,INARR
     &        ,INDSSTR  ,INEW     ,INORPAR  ,INTI     ,IODENS
     &        ,IOFLLI   ,IOFMLF   ,IVPAR    ,KXX      ,LMXNDL
     &        ,LNNDEL   ,LXPAREL  ,NFLAGS   ,NFNL     ,NFNLPAR
     &        ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAREL
     &        ,NPARF    ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP
     &        ,NZONE_PAR,NZPAR    ,PARACD   ,PARC     ,PAREL
     &        ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,NPAR     ,IPOS
     &        ,DERIV)

C------------------------- Derivatives with respect external concentration of
C------------------------- areal recharge

          IF (IODENS.EQ.1) THEN

	        CALL DERCOE_ARR_FLU
     &            (AREA     ,BETAC    ,CCALAN   ,CCALIT   ,CFPAREL
     &            ,CREF     ,DENSREF  ,DERH     ,DTIM     ,EPSTRA
     &            ,FNT      ,IBCOD    ,IDIMDERH ,IDIMFNT  ,IFLAGS
     &            ,INARR    ,INDSSTR  ,INEW     ,INORPAR  ,INTI
     &            ,IOCAP    ,IOFMLT   ,IOTRLI   ,IVPAR    ,KXX
     &            ,LMXNDL   ,LNNDEL   ,LXPAREL  ,NFLAGS   ,NFNL
     &            ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT
     &            ,NPARF    ,NPAREL   ,NPPEL    ,1        ,NTYPAR
     &            ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD
     &            ,PARC     ,PAREL    ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR
     &            ,IPOS     ,DERIV)

	    END IF !IODENS.EQ.1
     
      END IF !NZONE_PAR(3).NE.0

C------------------------- Derivatives with respect to storage
       

      IF (INDSSTR.NE.0) THEN

          CALL DERSTG
     &        (AREA     ,BETAC    ,CAUX1    ,CAUX2    ,CFPAREL
     &        ,CREF     ,DENSITY  ,DENSREF  ,DERH     ,DTIM
     &        ,EPSFLU   ,FNT      ,HAUX2    ,HBASE    ,HCALAN
     &        ,HCALIT   ,IDIMDERH ,IDIMFNT  ,IFLAGS   ,INDSSTR
     &        ,INEW     ,INORPAR  ,INTI     ,IPAR_DIR(9)
     &        ,IODENS   ,IOFLLI   ,IOFMLF   ,IOVRWC   ,IVPAR
     &        ,KXX      ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPAREL
     &        ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &        ,NFTPAR   ,NINT     ,NPAREL   ,NPARF    ,NTYPAR
     &        ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD
     &        ,PARC     ,THETAT   ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR
     &        ,NPAR     ,IPOS     ,DERIV)

      END IF !INDSSTR.NE.0

C------------------------- Derivatives with respect to transport parameters,
C------------------------- namely porosity, when solving variable density problems

      IF (IOINV.EQ.3 .AND. IODENS.EQ.1) THEN

	    CALL DERPOR_FLU
     &        (ACTH     ,AREA     ,BETAC    ,CAUX1    ,CAUX2
     &        ,CCALAN   ,CCALIT   ,CFPAREL  ,CREF     ,DENSITY
     &        ,DENSREF  ,DERH     ,DTIM     ,EPSTRA   ,FNT
     &        ,IDIMDERH ,IDIMFNT  ,IFLAGS   ,INDSSTR  ,INEW
     &        ,INORPAR  ,INTI     ,IOCAP    ,IOFMLT   ,IOTRLI
     &        ,IOVRWC   ,IVPAR    ,KXX      ,LMXNDL   ,LNNDEL
     &        ,LXPAREL  ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &        ,NFNLTIP  ,NFTPAR   ,NINT     ,NPARF    ,NPAREL
     &        ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &        ,PARACD   ,PARC     ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR
     &        ,IPOS     ,DERIV)

      END IF !IOINV.EQ.3 .AND. IODENS.EQ.1


C------------------------- Correction for the derivatives with respect 
C------------------------- to presc. head parameter (due to the  
C------------------------- symmetrization of flow system matrix in 
C------------------------- boundary cond. of type 1).
C------------------------- It has no sense in two cases:
C-------------------------    i) Nonlinear flow, because the flow system
C-------------------------        matrix is non-symmetric.
C-------------------------   ii) When system is solved with a sparse
C-------------------------       matrix (is non-symetric too).

      IF (IODENS.EQ.0 .AND. ISPARSE.EQ.0 .AND. IOFLLI.EQ.0
     &   .AND. ISYMETRIC.EQ.1 .AND. NZONE_PAR(4).NE.0) THEN

          
          IF (INDSSTR.EQ.0) THEN          ! Steady state flow

              CALL RHS_IN_CHP
     &            (AFLU     ,CFPARNP(1,INCHP)   ,DERH     ,IBCOD
     &            ,IDIMAFLU ,IDIMDERH ,INEW     ,IVPAR(INORPAR(9)+1)
     &            ,IXPARNP(1,INCHP)   ,KXX      ,LMXNDL   ,LNNDEL
     &            ,NPARF    ,NUMEL    ,NUMNP    ,NZONE_PAR(4)
     &            ,IDIMWGT  ,NZPAR    ,IPNT_PAR)


          ELSE                            ! Transient, lumped

* DTIMAUX, para que se calcule en k+1

              DTIMAUX = DTIM + (1D0-THETAF)*TINC/TINTERVOBS      

              CALL RHS_IN_CHP_TRLUM
     &            (AFLU     ,CFPARNP(1,INCHP)   ,DERH     ,DFLU
     &            ,DTIMAUX  ,FNT      ,IBCOD    ,IDIMAFLU ,IDIMDERH
     &            ,IDIMDFLU,IDIMFNT  ,INEW     ,INTI,IVPAR(INORPAR(9)+1)
     &            ,IXPARNP(1,INCHP)   ,KXX      ,LMXNDL   ,LNNDEL
     &            ,NFTPAR(INORPAR(9)+1)         ,NINT     ,NPARF
     &            ,NUMEL    ,NUMNP    ,NZONE_PAR(4)       ,THETAF
     &            ,TINC     ,IDIMWGT  ,NZPAR    ,IPNT_PAR)

          END IF !INDSSTR.EQ.0

      END IF !IODENS.EQ.0 .AND. ISPARSE.EQ.0 .AND...

C------------------------- Derivatives with respect to any parameter (except 
C------------------------- prescribed head that will be set later) have to be 
C------------------------- set equal to zero at nodes with presc. bound. cond.
       
      DO NP=1,NPARF

          DO I=1,NUMNP

              IF (IBCOD(I).EQ.1) THEN

                  DERH(I,NP,INEW) = 0D0

	        END IF !IBCOD(I).EQ.1

          END DO !I=1,NUMNP

      END DO !NP=1,NPARF


C------------------------- Derivatives with respect to nodal parameters

C------------------------- Cross all nodes

      DO I=1,NUMNP

          IB = IBCOD(I)

C------------------------- Derivative with respect to prescribed head

          IF (IB.EQ.1) THEN

              CALL DER_PRESC_HEAD
     &            (CFPARNP  ,DERH     ,DTIM     ,EPSFLU   ,FNT
     &            ,HCALAN   ,HCALIT   ,IDIMDERH ,IDIMFNT  ,IFLAGS
     &            ,INCHP    ,INDSSTR  ,INEW     ,INORPAR  ,INTI
     &            ,IOCAP    ,IOFLLI   ,IOFMLF   ,IVPAR    ,IXPARNP
     &            ,KXX      ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLPAR
     &            ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT     ,I
     &            ,NPARF    ,NPARNP   ,1        ,NTYPAR   ,NUMEL
     &            ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &            ,THETAF   ,TINC     ,TINTERVOBS         ,IDIMWGT
     &            ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)


              
          END IF !IB.EQ.1 <==> Prescribed head condition

C------------------------- Derivative with respect to boundary flow

          IF (IB.EQ.2 .OR. IB.EQ.4) THEN

	        CALL DER_PRESC_FLOW
     &          (BETAC    ,CAUX1    ,CFPARNP  ,CREF     ,DENSREF
     &          ,DERH     ,DTIM     ,EPSFLU   ,FNT      ,HCALAN
     &          ,HCALIT   ,IDIMDERH ,IDIMFNT  ,IFLAGS   ,INDSSTR
     &          ,INEW     ,INORPAR  ,INQQP    ,INTI     ,IOCAP
     &          ,IODENS   ,IOFLLI   ,IOFMLF   ,IVPAR    ,IXPARNP
     &          ,KXX      ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLPAR
     &          ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT     ,I
     &          ,NPARF    ,NPARNP   ,NPPNP    ,1        ,NTYPAR
     &          ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD
     &          ,PARC     ,PARNP    ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR
     &          ,IPOS     ,DERIV)


          END IF !IB.EQ.2 .OR. IB.EQ.4


C------------------------- Derivatives of leakage condition

          IF (IB.EQ.3 .OR. IB.EQ.4) THEN

              CALL DER_LEAK
     &            (BETAC    ,CAUX1    ,CFPARNP  ,CREF     ,DENSREF
     &            ,DERH     ,DTIM     ,EPSFLU   ,FNT      ,HAUX1
     &            ,HCALAN   ,HCALIT   ,IDIMDERH ,IDIMFNT  ,IFLAGS
     &            ,INALF    ,INCHP    ,INDSSTR  ,INEW     ,INORPAR
     &            ,INTI     ,IOCAP    ,IODENS   ,IOFLLI   ,IOFMLF
     &            ,IVPAR    ,IXPARNP  ,KXX      ,LMXNDL   ,NFLAGS
     &            ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR
     &            ,NINT     ,I        ,NPARF    ,NPARNP   ,NPPNP
     &            ,1        ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &            ,NZPAR    ,PARACD   ,PARC     ,PARNP    ,IDIMWGT
     &            ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)

          END IF !IB.EQ.3 .OR. IB.EQ.4

	    IF (IODENS.EQ.1) THEN

C------------------------- Derivatives of external concentration associated to
C------------------------- nodal fluxes.

	        CALL DERCOE_NOD_FLU
     &            (BETAC    ,CCALAN   ,CCALIT   ,CFPARNP  ,CREF
     &            ,DENSREF  ,DERH     ,DTIM     ,EPSTRA   ,FNT
     &            ,HAUX1    ,IBCOD    ,IDIMDERH ,IDIMFNT  ,IFLAGS
     &            ,INDSSTR  ,INEW     ,INORPAR  ,INTI     ,IOCAP
     &            ,IOFMLT   ,IOTRLI   ,IVPAR    ,IXPARNP  ,KXX
     &            ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &            ,NFNLTIP  ,NFTPAR   ,NINT     ,I        ,NPARF
     &            ,NPARNP   ,NPPNP    ,1        ,NTYPAR   ,NUMEL
     &            ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &            ,PARNP    ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS
     &            ,DERIV)


	    END IF !IODENS.EQ.1

       END DO        ! Nodes

C------------------------- Writes the RHS of the inverse problem (sensitivity 
C------------------------- equations) only if density is constant.

       IF (IFLAGS(10).NE.0 .AND. IODENS.EQ.0) THEN

          DO NP=1,NPARF

             WRITE(MAINF,*)'RHS PARAMETRO A ESTIMAR NUM., inti',NP,INTI

             DO I=1,NUMNP

                 WRITE(MAINF,10) NP,I,DERH(I,NP,INEW)
   10            FORMAT (2I5,G20.13)

             END DO !I=1,NUMNP

          END DO !NP=1,NPARF

       END IF !IFLAGS(10).NE.0


C------------------------- Solves flow inverse problem linear system 

      IF (IODENS.EQ.0) THEN

C-------------------- Writes DERH before solving

          IF (IFLAGS(25).GT.0) THEN

              IF(IODENS.EQ.0) THEN
                  WRITE(747,*) 'DERH BEFORE INTI = ',INTI
	            DO I=1,NUMNP
	                WRITE(747,strFmt1) I,(DERH(I,J,INEW),J=1,NPARF)
	            END DO !I=1,NUMNP
	        END IF !IODENS.EQ.0
c    1         FORMAT(I5,<NPARF>F15.10)

          END IF !IFLAGS(25).GT.0

C------------------------- RHS stored in DERH(.,NP,INEW) on the call to the 
C------------------------- linear solvers. On output the actual derivatives 
C------------------------- are stored on the same variable
          DO NP=1,NPARF

              CALL SOLVE
     &            (IAFLUDSC_ROWS ,IAFLUDSC_COLS ,NUMNP         ,IDESC         
     &            ,IPAR_DIR(21)  ,IDIMWORK      ,IDIRECT       ,1        
     &            ,INTI          ,IPAR_DIR(18)  ,IPAR_DIR(22)
     &            ,ISYMETRIC     ,ITERM         ,IPAR_DIR(15)  ,MAINF
     &            ,NBAND1        ,IPAR_DIR(24)
     &            ,IPAR_DIR(20)  ,IPAR_DIR(17)  ,IPAR_DIR(16)
     &            ,IPAR_DIR(36)  ,IPAR_DIR(37)  ,IPAR_DIR(38)  ,AFLUDSC
     &            ,AFLUDSCF      ,DERH(1,NP,INEW)              ,IAD_S
     &            ,IADD_S        ,IADN_S        ,IAFD_S        ,IAFDD_S
     &            ,IAFDN_S       ,WORK          ,SOLUTION)

 
              DERH(:,NP,INEW) = SOLUTION(:)

          END DO !NP=1,NPARF

	    IF (IFLAGS(25).GT.0) THEN

              WRITE(747,*) 'DERH AFTER INTI = ',INTI
              DO I=1,NUMNP
	            WRITE(747,strFmt1) I,(DERH(I,J,INEW),J=1,NPARF)
              END DO !I=1,NUMNP

          END IF !IFLAGS(25).GT.0

	    IF (IFLAGS(4).GT.0) THEN

              WRITE(749,*) 'DERH AFTER INTI = ',INTI

              DO I=1,NUMNP

                  WRITE(749,strFmt1) I,(DERH(I,J,INEW),J=1,NPAR)

              END DO !I=1,NUMNP

          END IF !IFLAGS(4).GT.0

      END IF !IODENS.EQ.0


C------------------------- Writes the derivative of nodal head with respect 
C------------------------- to every parameter

       IF (IFLAGS(13).NE.0 .AND. IODENS.EQ.0) THEN

          DO NP=1,NPARF

             WRITE(MAINF,*)'DERIVADA DEL PARAMETRO A ESTIMAR NUM.',NP

             DO I=1,NUMNP
                WRITE(MAINF,'(2I5,G20.13)') NP,I,DERH(I,NP,INEW)
             END DO !I=1,NUMNP

          END DO !NP=1,NPARF

       END IF !IFLAGS(13).NE.0 .AND. IODENS.EQ.0
	   
	  IF(IFLAGS(3).EQ.1) CALL IO_SUB('JAC_H',1)

      END SUBROUTINE JAC_H
