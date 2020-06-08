      SUBROUTINE COMP_BTRA
     &          (AREA     ,BETAC    ,BTRA     ,CAUDAL   ,CAUX1
     &          ,CCALAN   ,CREF     ,DENSREF  ,IBTCO    ,INARR
     &          ,IOCONSRC ,IODENS   ,IONEWT   ,ITPTVAR  ,KXX
     &          ,LMXNDL   ,LNNDEL   ,LXPAREL  ,NPAREL   ,NPPEL
     &          ,NPPNP    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &          ,PAREL    ,PARNP    ,THETAT   ,WSPECHEAT)

********************************************************************************
*
* PURPOSE
*
*  Manages the computation of BTRA vector with boundary conditions.
*  Notice that the derivatives of BTRA w. r. t. head and mass fraction are null.
*
* DESCRIPTION
*
*      Includes areal concentration contribution and boundary conditions
*      to transport RHS
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*
*  BTRA                   Right hand side of transport discretized equation     
*  CAUDAL                 Input/output flow at every node.                      
*  CCALAN                 Computed concentrations in the previous time step.    
*  CCALIT                 Computed concentration in last iteration              
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  DENSITY                Array containig density of each element.
*  IBTCO                  Transport boundary condition index                    
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR
*  IPAR_DIR               Array containing all integer direct problem           
*                         parameters                                            
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LINMET                   Linearization method
*                         (see CHOOSE_LINEARITZ_METHOD for details).
*  LNNDEL                 Number of nodes at every element                      
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  NFNLPAR                Vector containing non-linear function order           
*                         afecting every parameter at each zone.                
*  NFNLPRG                Generic parameter zone number for every nonlinear     
*                         function                                              
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
*  XPARAM                                                                       
*
* EXTERNAL VARIABLES: SCALARS
*
*  BETAC                  Derivative of the density to C divided by density
*  DENSREF                Density of reference (at C=CREF )
*  DTIMBTRA               Relative time at which transport right hand side is   
*                         evaluated. If 1, rhs is evaluated at the              
*                         next computing time, if 0, rhs is evaluated           
*                         at the last computing time.                           
*  EPSTRA                 Time weighting parameter for nonlinear transport      
*                         problems                                              
*  FNT                    Array containing time functions values                
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  INARR                  Index for areal recharge                              
*                         in array variables (LXPAREL and CFPAREL)              
*  INCOE                  Index for external concentration (elements)           
*                         in array variables (LXPAREL and CFPAREL)              
*  INDSSTR                Current problem state. If 1, transient state,         
*                         if 0, steady-state                                    
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOFMLT                 Transport formulation number
*  IOINV                  Inverse problem option
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFLAGS                 Maximum number of allowed flags                       
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NINT                   Number of observation times                           
*  NPARALG                Maximum number of algorithm parameters, including     
*                         minimization and simulation (used for dimensioning)   
*  NPAREL                 Number of element parameters in current problem       
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  NZPRG                  Total number of generic parameter zones               
*  PRGC                   Generic parameters
*  CREF                   Mass fraction o reference.
*
*  WSPECHEAT              Water specific heat.
*
* INTERNAL VARIABLES: SCALARS
*
*  NNUD                   Number of nodes of the current element                
*  NZA                    Areal recharge zone at current element
*  NZC                    External concentration zone at current element
*  DENRATIO               Ratio between external density and density of the
*                         element.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  PARAM_VALUE                                                                  
*  ZERO_ARRAY                                                                   
*  DENS
*
* HISTORY
*
*     AMS        1988     First coding
*     AMS      3-1999     Revision: Inclusion of nonlinear calls, addition 
*                         of comments and common elimination
*     JHG      5-2003     Revision: Inclusion of density dependent 
*                         flow an transport.
********************************************************************************

      IMPLICIT NONE

C------------------------- External
      
      INTEGER*4::INARR    ,IOCONSRC,IODENS   ,IONEWT   ,ITPTVAR
     &          ,LMXNDL   ,NPAREL   ,NPPEL   ,NPPNP    ,NTYPAR
     &          ,NUMEL    ,NUMNP

      INTEGER*4::IBTCO(NUMNP)          ,KXX(LMXNDL,NUMEL) ,LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL) ,NZONE_PAR(NTYPAR)

      REAL*8::BETAC    ,CREF     ,DENSREF  ,THETAT ,WSPECHEAT

      REAL*8::DENS

      REAL*8::AREA(NUMEL)         ,BTRA(NUMNP)    ,CAUDAL(NUMNP)
     &       ,CAUX1(NUMNP)        ,CCALAN(NUMNP)  ,PARNP(NUMNP,NPPNP)
     &       ,PAREL(NUMEL,NPPEL)

C------------------------- Internal

      INTEGER*4::I    ,IB   ,K    ,KNODE,L    ,NNUD ,NZA

      REAL*8::ALFAX    ,AREALN   ,CAUD     ,CAUDNOD
     &       ,CONC     ,CNODE    ,COE      ,DENSEXT  ,DENSRECH ,FACTOR
     &       ,RECN     ,RECHRG   ,THTT1

C------------------------- First executable statement


C------------------------- Performed only if areal recharge zones are present

      IF (NZONE_PAR(3).NE.0) THEN

          DO L=1,NUMEL

              NZA = LXPAREL(L,INARR)

              IF (NZA.NE.0) THEN

                  NNUD = LNNDEL(L)
                  AREALN = AREA(L)/NNUD

                  RECHRG = PAREL(L,8)
                  RECN = RECHRG*AREALN

                  IF (RECHRG.GT.0) THEN

                      COE = PAREL(L,15)

                  ELSE !Evaporation

                      COE = 0D0

                  END IF !RECHRG.GT.0

C------------------------- Density dependence effect on recharge.

                  IF(IODENS.EQ.1) THEN

                      DENSRECH = DENS(DENSREF,BETAC,COE,CREF)

                  ELSE

                      DENSRECH = 1D0

                  END IF !IODENS.EQ.1


                  DO K=1,NNUD

                      KNODE = KXX(K,L)

                      IF (IONEWT.EQ.0) THEN !Picard

                          BTRA(KNODE) = BTRA(KNODE) + RECN*COE*DENSRECH

                      ELSE !NEWTON

                          CNODE = CAUX1(KNODE)

                          BTRA(KNODE) = BTRA(KNODE)
     &                                 + RECN*DENSRECH*(COE - CNODE)

                      END IF !(IONEWT.EQ.0)
                             
                  END DO !K=1,NNUD

              END IF ! NZA .NE. 0

          END DO ! L=1,NUMEL

      END IF ! NZONE_PAR(3).NE.0

C------------------------- Boundary conditions

      THTT1 = THETAT - 1D0

      DO I=1,NUMNP

          IB = IBTCO(I)

          SELECT CASE (IB)

C------------------------- Mass flow.

          CASE(2,3)

              CAUDNOD = CAUDAL(I)

              IF (CAUDNOD.GT.0) THEN

                  COE = PARNP(I,4)

                  IF (IODENS.EQ.1) THEN

                      DENSEXT = DENS(DENSREF,BETAC,COE,CREF)

                  ELSE

                      DENSEXT = 1D0

                  END IF !IODENS.EQ.1


                  IF (IONEWT.EQ.1) THEN

                      CAUD = DENSEXT*CAUDNOD*(COE - CAUX1(I))

                  ELSE

                      CAUD = DENSEXT*CAUDNOD*(COE + THTT1*CCALAN(I))

                  END IF !IONEWT.EQ.1

                  BTRA(I) = BTRA(I) + CAUD

              END IF !CAUDNOD.GT.0

C------------------------- Prescribed concentration.

          CASE(1)

                  BTRA(I) = PARNP(I,4)

C------------------------- Input mass

          CASE(4)

              IF (IOCONSRC.EQ.1) THEN

                  FACTOR = (1D0 + THTT1*CCALAN(I))

              ELSE

                  FACTOR = 1D0

              END IF !IOCONSRC.EQ.1

C------------------------- Whe solving heat transport the heat flow
C------------------------- has to be dived by the water specific heat
C------------------------- to be consistent with the rest of the equation.

              IF (ITPTVAR.EQ.1) THEN

                  IF (IODENS.EQ.0) THEN

C------------------------- If density is constant, water density appears
C------------------------- in this term (see TRANSDENS Guia rapida, 3.7).
                      FACTOR = FACTOR/(DENSREF*WSPECHEAT)

                  ELSE

                      FACTOR = FACTOR/WSPECHEAT

                  END IF !IODENS.EQ.0
             END IF !ITPTVAR.EQ.1

              BTRA(I) = BTRA(I) + PARNP(I,4)*FACTOR

C--------------------------- Conc. leakage

          CASE(5)

              ALFAX = PARNP(I,6)
              CONC = PARNP(I,4)

C------------------------- Whe solving heat transport the heat flow
C------------------------- has to be dived by the water specific heat
C------------------------- to be consistent with the rest of the equation.

              IF (ITPTVAR.EQ.1) THEN

                  IF (IODENS.EQ.0) THEN

C------------------------- If density is constant, water density appears
C------------------------- in this term (see TRANDES Guia rapida, 3.7).
                      ALFAX = ALFAX/(DENSREF*WSPECHEAT)

                  ELSE

                      ALFAX = ALFAX/WSPECHEAT

                  END IF !IODENS.EQ.0

             END IF !ITPTVAR.EQ.1

C--------------------------- Constant density or heat transport or
C--------------------------- solute transport with variable density and
C--------------------------- no concentration sources [alf*(w*-w)].

              IF (IODENS.EQ.0 .OR. ITPTVAR.EQ.1 .OR.
     &           ((IODENS.EQ.1 .AND. ITPTVAR.EQ.0 .AND. IOCONSRC.EQ.0)
     &           )) THEN

                  IF (IONEWT.EQ.0) THEN

                      BTRA(I) = BTRA(I) + ALFAX*(CONC + THTT1*CCALAN(I))

                  ELSE

                      BTRA(I) = BTRA(I) + ALFAX*(CONC-CAUX1(I))

                  END IF !IONEWT.EQ.0

C--------------------------- Solute transport with variable density and
C--------------------------- concentration sources [alf*(w*-w)*(1-w)].

	        ELSE IF (IODENS.EQ.1 .AND. ITPTVAR.EQ.0
     &                .AND. IOCONSRC.EQ.1) THEN

                  IF (IONEWT.EQ.0) THEN

                      BTRA(I) = BTRA(I) 
     &                 + ALFAX*(1D0-CAUX1(I))*(CONC + THTT1*CCALAN(I))


                  ELSE

                      BTRA(I) = BTRA(I)
     &                         + ALFAX*(CONC-CAUX1(I))*(1D0-CAUX1(I))

                  END IF !IONEWT.EQ.0

              END IF !IODENS.EQ.0 ...

          END SELECT !IB

      END DO ! I=1,NUMNP

      END SUBROUTINE COMP_BTRA
