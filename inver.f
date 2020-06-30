       SUBROUTINE INVER
     &           (A_COUPL_DSC        ,A_COUPL_DSCF       ,ACTH
     &           ,AFLU     ,AFLUDSC  ,AFLUDSCF ,ALFA     ,AREA
     &           ,ATRA     ,ATRADSC  ,ATRADSCF ,BCOUPLED ,BETAC    ,BFLU
     &           ,BIBI     ,BM_ND_FL ,BM_ND_TT ,BM_ZN_FL
     &           ,BM_ZN_TT ,BTRA     ,BUOYANCY ,DBUOYANCY
     &           ,CAUDAL   ,CAUX1    ,CAUX2
     &           ,CCALAN   ,CCALIT   ,CFLU
     &           ,CFPAREL  ,CFPARNP  ,CONCFLOW ,COORD    ,COVINV
     &           ,COVPAR   ,CPREV1   ,CPREV2   ,CREF     
     &           ,DAT_VD   ,DBFLUDFLU,DBFLUDTRA
     &           ,DELTAITER,DENSITY  ,DENSREF
     &           ,DERVISC  
     &           ,DERC     ,DERH
     &           ,DFLU     ,DFLUDFLU ,DFLUDTRA ,DNODALRH ,DPARELDC
     &           ,DPARELDH ,DQDFLU   ,DQDTRA   
     &           ,DTMXDS   ,DTPREVINV,DTRA     ,DTRADFLU
     &           ,DTRADTRA ,DVDC     ,DVDH     ,DVDP     ,DVOBS  ,DWDH
     &           ,FILENAME ,FNT
     &           ,FOBJ_WGT ,GP_COORD ,GRAD     ,GRAVEL   
     &           ,GRDFF    ,HAUX1    ,HAUX2    ,HBASE    ,HCALAN
     &           ,HCALIT   ,HESS     ,HESSAUX  ,HPREV1
     &           ,HPREV2
     &           ,IA_COUPLED_DSC_COLS,IA_COUPLED_DSC_ROWS,IAD_D
     &           ,IAD_S    ,IADD_D   ,IADD_S   ,IADN_D
     &           ,IADN_S   ,IAFD_D   ,IAFD_S   ,IAFDD_D  ,IAFDD_S
     &           ,IAFDN_D  ,IAFDN_S  ,IAFLUDSC_COLS      ,IAFLUDSC_ROWS
     &           ,IATRADSC_COLS      ,IATRADSC_ROWS      ,IBCOD
     &           ,IBTCO    ,IDIMAFLU
     &           ,IDIMATRA ,IDIMBB   ,IDIMCFLU ,IDIMCOV  ,IDIMDENS
     &           ,IDIMDERH ,IDIMDERC ,IDIMDFLU ,IDIMDQ   ,IDIMDTRA
     &           ,IDIMFNT  ,IDIMGRAVEL,IDIMHESS,IDIMQ
     &           ,IDIMWORK
     &           ,IFLAGS   ,INDEXNOD ,INDPAR   ,INORPAR
     &           ,IOCONSRC ,IOCRITRAP,IODENS_INI
     &           ,IODEVICE ,IODIM    ,IODIRECT ,IOEQT    ,IOFLSAT
     &           ,IOFLLI   ,IOFMLF   ,IOFMLT   ,IOINV
     &           ,IOLG_PAR ,IOPINITC ,IOPINITH ,IOPINVDT,IOPTS
     &           ,IORTS    ,IOTRLI   ,IOTRS    ,IOWRITE
     &           ,IPAR_DIR ,IPAR_INV ,IPARTNER           ,ISOLEQ
     &           ,ISOT     ,ISOZ     ,ISPARSE  ,ITERGLMX
     &           ,ITPTVAR  ,ITYPACOUPLDSC      ,ITYPAFLU
     &           ,ITYPATRA ,ITYPCFLU
     &           ,ITYPDFLU ,ITYPDTRA
     &           ,ITYPFLUDSC,ITYPTRADSC        ,IVPAR    ,IXPARNP
     &           ,KINT     ,KXX      ,LDIM
     &           ,LINMET   ,LMXNDL   ,LNNDEL   ,GRADLOC  ,LTYPE
     &           ,LXPAREL  ,MAINF    ,MAXNB    ,MAXNBF
     &           ,MAXNN    ,MAXPG    ,MEASTYP  ,NBAND
     &           ,NBAND1   ,NBANDCOV ,NDEVS
     &           ,NFLAGS   ,NFNL
     &           ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR
     &           ,NINT     ,NMAXF    ,NMAXT    ,NOOBSIT  ,NOPTS
     &           ,NPAR     ,NPARALG
     &           ,NPAREL   ,NPARF    ,NPARNP
     &           ,NPBFL    ,NPBMX    ,NPBTP
     &           ,NPPEL    ,NPPNP    ,NROW     ,NSTAT    
     &           ,NTDMT    ,NTYPAR   ,NTYPEL   ,NUMEL    ,NUMNP
     &           ,NUMTIT   ,NUMTNOD  ,NUMTOBS  ,NWRITE
     &           ,NZONE_PAR,NZPAR    ,NZPRG    ,NZTRA    ,DLT_PAR
     &           ,PAR_DIR  ,PAR_INV  ,PAR_WGT  ,PARACD   ,PARAUX
     &           ,PARZ     ,PAREL    ,PARM     ,PARNP    ,POINTWEIGHT
     &           ,PRGC     ,QXYZ     ,SOLUTION
     &           ,SOURCE             ,TIME     ,TIT      ,TOBS
     &           ,VAR_REF  ,VD       ,VISCOSITY,VISCREF  ,VJAC
     &           ,VOBS     ,VOBSC    ,WATVOL   ,WORK
     &           ,WTOBSN   ,WTOBST   ,XNORVD   ,HINI     ,WSPECHEAT
     &           ,WTHERMCON
     ; ,IOINV_GS    ,IDIMDATASC_GS
     ; ,MXGRPZN,NGROUP_ZN,MXMEASPP_GS,IDIMVAR_GS,MXNVAR_GS
     ; ,IDIMCROSS_GS,MXKRIG_GS,IDIMIVARIO_GS,MXNPRIM_GS
     ; ,MXNZON_GS,MXCLOSE_GS,MXNPP_GS,IDIMZONPP_GS,MXSB_GS
     ; ,MXDISC_GS,MXROT_GS,MXSAM_GS,IOPT_GS,IO_KG_GS
     ; ,IZN_PP_GS,IZN_NPP_GS,IPOLDRIFT_GS,IVARIO_GS,NUMSB_GS,ISUPBL_GS
     ; ,ICROSSCOV_GS,POSMEAS_GS,VMEAS_GS,VSTATS_GS,SEARCH_GS
     ; ,TRIM_GS,SUPBL_GS,KRISYS_GS,CLOSESAM_GS,ZNWGT_GS,KRISOL_GS
     ; ,VARIO_GS,ROTMAT_GS,POSZN_GS,KRIGAUX_GS,CROSSCOV_GS
     ; ,ESTKRIG_GS,POSDISAUX_GS,POSDIS_GS,EXDRZN_GS,MXZONPP_GS
     ; ,DATASC_GS,IDIMWGT,IPNT_PAR,WGT_PAR,NPARDET,ICHECK_GS,LDIM_GS
     ; ,COORDGR_GS,PARC,WGT_UNK,IPOS,DERIV,PARGOOD,IACTSIMU
     ; ,PARC_GS,IFLAG_SIMUL,COVPAR_GR, PARNAME)

*****************************************************************************
*
* PURPOSE
*       Solves the inverse problem
*
* DESCRIPTION
*
*       Manages all processes related o inverse problem: Simulation and
*       jacobiens computation, Minimization process, Statistical results, etc.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CCAL                   Computed concentration at every node                  
*  CCALIT                 Computed concentration in last iteration              
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  CFPARNP                Array containing node coefinient of node j            
*                         corresponding to INpar index parameter zone.          
*  CJAC                   Jacobian matrix of concentrations with respect        
*                         to estimated parameters at the observation points.    
*  COBS                   Measured concentrations at the observation points.    
*  COBSC                  Computed concentrations at the observation points.    
*  COORD                  Nodal coordinates                                     
*  COORD_OBS              Observation points coordinates                        
*  DERH                   Nodal head derivatives with respect to estimated      
*                         flow parameters.                                      
*  FOBJ_WGT               Array containing all objective function weights for   
*                         state variables (heads, concentrations, etc)          
*  GRAD                   Vector containing objective function gradient         
*  GRDFF                  Array containing the product between interpolation    
*                         functions integrals and interp. functions gradient    
*  HCALIT                 Computed heads in last iteration                      
*  HEAD                                                                         
*  HESS                   Hessian matrix of objective function.                 
*  HJAC                   Jacobian matrix of heads with respect to estimated    
*                         parameters at observation points                      
*  HOBS                   Measured heads at observation points                  
*  HOBSC                  Computed heads at observation points                  
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IOLG_PAR               Array containing all logarithmic options of           
*                         estimated parameters                                  
*  IOPTS                                                                        
*  IOWRITE                Array containing all output options                   
*  IPAR_DIR               Array containing all integer direct problem           
*                         parameters                                            
*  IPAR_INV               Array containing all integer inverse problem          
*                         parameters                                            
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  IXPARNP                Array containing zone number of each node j,          
*                         corresponding to INpar index parameter zone.          
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  NFNLPAR                Vector containing non-linear function order           
*                         afecting every parameter at each zone.                
*  NFNLTIP                Type of non-linear function                           
*  NFTPAR                 Vector containing time function number at every       
*                         parameter zone                                        
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PAR                    Parameter's increment at every iteration.             
*  PARAUX                 Auxiliary variable to store the best computed         
*                         parameters                                            
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*  PARM                   Vector containing measured values for all             
*                         parameters                                            
*  PARNP                  Parameter values at every node and current time for   
*                         all nodal parameters (each value is computed as the   
*                         product of up to four terms:                          
*                           nodal coeff*zonal value*time funct.*nonl. funct. )  
*  PAR_DIR                Array containing all real direct problem              
*                         parameters                                            
*  PAR_INV                Array containing all real inverse problem             
*                         parameters                                            
*  PAR_WGT                Array containing objective function weights for       
*                         all estimated parameters                              
*  PNAME                  Name of observation points                            
*  QXYZ                   Products between the different components of          
*                         Darcy's velocity divided by its norm                  
*  STCOBS                 Standard deviation of concentration data              
*  STHOBS                 Standard deviation of head data                       
*  STPAR                  Vector containing standard deviation errors of        
*                         all parameters prioo information                      
*  TIME                   Observation times.                                    
*  WORK                   Workspace array used in some subroutines.             
*
* INTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*
* EXTERNAL VARIABLES: SCALARS
*
*  ACTH                   Aquifer thickness of every element. Cross sectional   
*                         area for 1-D elements, thickness for 2-D elements.    
*  AFLU                   Matrix of finite elements equations for flow problem  
*                         No boundary conditions are included on it.            
*  AFLUDSC                Coefficient matrix of flow system (2.15), most        
*                         often in its decomposed form.                         
*  ALFA                                                                         
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  ATRA                   Matrix of finite elements equations for transport     
*                         problem. Only mass flow boundary conditions included. 
*  ATRADSC                Coefficient matrix of transport system. Most often in 
*                         its decomposed form.                                  
*  BFLU                   Right hand side of flow discretized equation.         
*  BFOBS                  Basis functions values at observation points.         
*  BIBI                   Array containing the product of interpolation         
*                         functions gradient, for a given element               
*  BM_ND_FL                                                                     
*  BM_ND_TT                                                                     
*  BM_ZN_FL                                                                     
*  BM_ZN_TT                                                                     
*  BTRA                   Right hand side of transport discretized equation     
*  CAUDAL                 Input/output flow at every node.                      
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  CAUX2                  Array containing diference of concentrations in two   
*                         consecutives times, related to time increment         
*  CBAS                   Concentration base line                               
*  CBASE                  Base concentration                                    
*  CCALAN                 Computed concentrations in the previous time step.    
*  CNST                   CNST(i,j,k) is the integral of the product of         
*                         interpolation functions i and j in an element of      
*                         type k divided by the AREA of the element. It is used 
*                         only in consistent scheme.                            
*  COVTRA                 Covariance matrix for transmissivity                  
*  CRECOV                 Concentration computed at k-2 time step               
*  CREL                                                                         
*  DERADFLU               Contains the derivatives of AFLU and DFLU matrices    
*                       ) with respect to head. This term appears on the left   
*                         hand side of both flow direct and inverse problems,   
*                         as well as on the right hand side of flow inverse     
*                         problem.                                              
*  DERADTRA               Contains the derivatives of ATRA and DTRA matrices    
*                         with respect to concentration. This term appears on   
*                         the left hand side of both transp. direct and inverse 
*                         problems, as well as on the right hand side of        
*                         transport inverse problem.                            
*  DERBFLU                                                                      
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
*  DERCAN                                                                       
*  DERCDT                 Derivative of concentrations respect to time          
*  DFLU                   Matrix of finite elements equations for flow          
*                         problem related to storage term.                      
*  DTMXDS                 Maximum allowed time increment for nonlinear          
*                         problems                                              
*  DTPREVINV                                                                    
*  DTRA                   Matrix of finite elements equations for transport     
*                         problem related to storage.                           
*  DVDP                                                                         
*  DYDXOLD                                                                      
*  FLOWGRAV                                                                     
*  FNT                    Array containing time functions values                
*  GRAV                   Gravity array direction                               
*  GRAVEL                 Projection of gravity at each element                 
*  HAUX1                  Array containing HEADS, ponderated by THETAF          
*                         time factor                                           
*  HAUX2                  Array containing diference of HEADS in two            
*                         consecutives times.                                   
*  HBAS                   Head measurements treshold. Measurements under HBAS   
*                         are not taken into account.                           
*  HBASE                                                                        
*  HCALAN                 Head level at previous time                           
*  HESSAUX                                                                      
*  HPREV1                                                                       
*  HPREV2                                                                       
*  IBCOD                  Flow boundary condition index                         
*  IBTCO                  Transport boundary condition index                    
*  IDIMAFLU               Used to dimension array AFLU                          
*  IDIMATRA                                                                     
*  IDIMBB                 Used to dimension array BIBI                          
*  IDIMDADF               Second dimension of matrix DERADFLU.                  
*  IDIMDADFLU                                                                   
*  IDIMDADTRA                                                                   
*  IDIMDFLU               Used to dimension array DFLU                          
*  IDIMDQ                 Used to dimension array DAT_VD (second dimension)       
*  IDIMDTRA               Used to dimension array DTRA (second dimension)       
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  IDIMHESS               Used to dimension array HESS. It is equal to          
*                         NPAR*(NPAR+1)/2                                       
*  IDIMQ                  Used to dimension array QXYZ                          
*  INDPAR                 Array of 0's and 1's. 0 means that the parameter is   
*                         estimated aritmethically, otherwise logarithmically.  
*  INTAUX                                                                       
*  IOALFT                                                                       
*  IOCNSF                 Scheme for storage term in flow problem               
*  IOCNST                 Scheme for mass storage term in transport problem     
*  IOCRITRAP              Option of treatement on the direct problem            
*                         convergence criteria                                  
*  IODIM                  Maximum dimension of any element included             
*                         in the problem                                        
*  IOEQT                  Type of problem to be solved                          
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOFLSAT                Indicates the possibility that one part of the domain 
*                         reaches unsaturated state.                            
*  IOFMLF                 Flow Formulation number                               
*  IOFMLT                 Transport formulation number                          
*  IOFOBJ                 If 0, the code solves completely the direct problem   
*                         at each iteration of the inverse problem. If sets     
*                         1 it computes completely the simulation only for      
*                         succesfull iterations of the inverse problem.         
*  IOINV                  Inverse problem option                                
*  IOOBS                  If it equals 1, all nodal points are observation      
*                         points. Otherwise they are a set of specified points  
*                         (usually 0)                                           
*  IOPINITC               Option for the extrapolation of concentrations        
*                         at the next time step in the Newton process.          
*  IOPINITH               Option for the extrapolation of heads or pressures    
*                         at the next time step in the Newton process.          
*  IOPINVDT                                                                     
*  IOPRHED                Indicates whether the flow state variable state is    
*                         preasure (set to 1) or head (set to 0)                
*  IORTS                  Transport regime                                      
*  IOSTC                  If it equals 1, standard deviation of measured        
*                         concentration errors is taken constant (usually 0)    
*  IOSTH                  If it equals 1, standard deviation of measured head   
*                         errors is taken constant (usually 0)                  
*  IOSUCHUM               Indicates if measures are given in terms of pressure  
*                         or piezometric heads (set to 0), or in terms of       
*                         saturation degree                                     
*  IOTRLI                 Idem to IOFLLI, in the case of transport.             
*  IOTRS                  Flow regime                                           
*  IPARTRA                Used for dimensioning, it equals                      
*                         MAX(1,NPARTRA*(NPARTRA-1)/2)                          
*  ISOLEQ                                                                       
*  ISOT                   Maximum hydraulic conductivity tensor anisotropy      
*                         degree in the problem.                                
*  ISOZ                   Anisotropy of every transmissivity zone               
*  KINT                   Number of solution time increments, between           
*                         successive observation times.                         
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LMXNDL                 Maximum number of nodes per element                   
*  LNNDEL                 Number of nodes at every element                      
*  LOBS                   Element number to which observation points belong to  
*  LTYPE                  Vector containing the type of each element            
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MAXNEOP                Used to reserve some space. It is equal to            
*                         MAX (NUMEL,NUMNP,NUOBS,NPAR)                          
*  NBAND                  Half Bandwith (maximum difference between the         
*                         numbers of two nodes belonging to the same element)   
*  NBAND1                 Used to dimension. It is equal to NBAND+1             
*  NBAND2                 Used to dimension. It is equal to 2*NBAND+1           
*  NFLAGS                 Maximum number of allowed flags                       
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NFNLPRG                Generic parameter zone number for every nonlinear     
*                         function                                              
*  NFNT                   Number of time functions used for describing time     
*                         dependence of all transient parameters                
*  NINT                   Number of observation times                           
*  NOPTS                                                                        
*  NPAR                   Total number of parameters to be estimated            
*  NPARALG                Maximum number of algorithm parameters, including     
*                         minimization and simulation (used for dimensioning)   
*  NPAREL                 Number of element parameters in current problem       
*  NPARF                  Number of transient parameters to be estimated        
*  NPARFPRG               Number of uncertain generic parameter zones involved  
*                         in the non-linear flow inverse problem                
*  NPARNP                 Number of nodal parameters in current problem         
*  NPARPRG                Total number of uncertain generic parameter zones     
*                         involved in the inverse problem                       
*  NPARTRA                Number of transmissivity zones to be estimated        
*                         with non-diagonal prior information covariance matrix 
*  NPBFL                  Number of simultaneous flow problems                  
*  NPBTP                  Number of simultaneous transport problems             
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NPPNP                  Total number of parameters by nodes (not confuse      
*                         with NPARNP, because in this casethere is no          
*                         difference between a given parameter in steady or tr.)
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)               
*  NTDMT                  If diferent form zero, number of terms in matrix      
*                         diffusion(if zero, no diffusion)                      
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NUOBS                  Number of observation points                          
*  NWRITE                 Number of output options (used for dimensioning)      
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  PARACD                 Agreement parameters                                  
*  PARDIF                                                                       
*  PRGC                                                                         
*  SOURCE                                                                       
*  STATVAR                                                                      
*  VD                     Darcy's velocity                                      
*  VOLNOD                                                                       
*  WATVOL                                                                       
*  XNORVD                 Euclidean norm of Darcy's velocity                    
*  YOLD                                                                         
*
* INTERNAL VARIABLES: SCALARS
*
*  ALF                                                                          
*  DABSMX                 Absolute convergence criterion                        
*  DABSMXF                                                                      
*  DABSMXT                                                                      
*  DRELMX                 Relative convergence criterion                        
*  DRELMXF                                                                      
*  DRELMXT                                                                      
*  FNEW                   Objective function value in the current iteration     
*  FOLD                   Objective function computed value in last iteration   
*  GNORM                                                                        
*  GNORM1                                                                       
*  I                                                                            
*  IDIMDADT                                                                     
*  IDIMTRAC                                                                     
*  INDCOS                                                                       
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOMIN                  AlgoritHm minimization option                         
*  ISUMFO                                                                       
*  ITACUMF                                                                      
*  ITTOTV                                                                                                                                              
*  MIN_STOP               Control variable for stopping minimization process    
*  NDCON                                                                        
*  NDHED                                                                        
*  NFORM                                                                        
*  NITERF1                                                                      
*  NITERF2                                                                      
*  NLOCMAX                                                                      
*  NMAXF                                                                        
*  NMAXT                                                                        
*  NUMITER                Current iteration in inverse problem process          
*  OBJCON                 Concentration contribution to objective function      
*  OBJGEO                 Geophysics observations to objective function               
*  OBJHED                 Head contribution to objective function               
*  OBJPAR                 Parameters contribution to objective function         
*  XMAXIM                                                                       
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ASSIGN                 Assigns new computed parameters to PARC variable      
*  FUN_PAR                Computes the objective parameters function            
*  INIT_VAR               Inicialises state variables for Newton process        
*  IO_SUB                                                                       
*  MAX_PAR                                                                      
*  MINIMIZATION           Selects minimization algoristm                        
*  MODIF_MAXITER          Reads maximum iteration allwed for non-linear flow    
*  MOD_RAPSON_CRITERIA    Modifies the Newton-raphson's convergence criteria    
*                         as function of Mardquart convergence process          
*  RESID                                                                        
*  SIM_JAC                Solves numerical state equations (Direct problem)     
*  STORE_PARAM            Stores last good parameters in PARAUX variable        
*  UPD_PAR_DMT                                                                  
*  WRIPAR_IT                                                                    
*  ZERO_ARRAY                                                                   
*
*****************************************************************************
       USE forpy_mod      
       IMPLICIT REAL*8 (A-H,O-Z)

       CHARACTER FILENAME(20)*20

C     EXTERNAL VARIABLES: SCALARS


C     EXTERNAL VARIABLES: ARRAYS
       DIMENSION IAD_D(*)                   ,IAD_S(*)
     &          ,IADD_D(*)                  ,IADD_S(*)
     &          ,IADN_D(*)                  ,IADN_S(*)
     &          ,IAFD_D(*)                  ,IAFD_S(*)
     &          ,IAFDD_D(*)                 ,IAFDD_S(*)
     &          ,IAFDN_D(*)                 ,IAFDN_S(*)
     &          ,IBCOD(NUMNP)               ,IBTCO(NUMNP)
     &          ,IFLAGS(NFLAGS)                  ,INDEXNOD(NUMTNOD)
     &          ,INDPAR(NPAR)
     &          ,INORPAR(NTYPAR)
     &          ,IODEVICE(NDEVS+1,10)            ,IOLG_PAR(NTYPAR,2)
     &          ,IOPTS(NOPTS)                    ,IOWRITE(NWRITE)
     &          ,KINT(NINT)                      ,IPAR_DIR(NPARALG)
     &          ,IPAR_INV(NPARALG)               ,IPARTNER(6,3,6)
     &          ,ISOLEQ(NINT,4)                  ,ISOZ(NZTRA)
     &          ,IVPAR(NZPAR,4)                    
     &          ,IXPARNP(NUMNP,NPARNP,NPBMX)     ,KXX(LMXNDL,NUMEL)
     &          ,LDIM(NUMEL)
     &          ,LINMET(3,2)                     ,LNNDEL(NUMEL)
     &          ,LTYPE(NUMEL)                    
     &          ,LXPAREL(NUMEL,NPAREL)           ,MEASTYP(NUMTOBS)
     &          ,NFNLPAR(NZPAR)                  ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)                   ,NFTPAR(NZPAR)
     &          ,NOOBSIT(NUMTIT)                 ,NZONE_PAR(NTYPAR)
              
       DIMENSION   A_COUPL_DSC(IA_COUPLED_DSC_ROWS,IA_COUPLED_DSC_COLS)
     &          ,A_COUPL_DSCF(MAXNB,2*NUMNP)
     &          ,ACTH(NUMEL)                ,AFLU(NUMEL,IDIMAFLU)
     &          ,AFLUDSC(IAFLUDSC_ROWS,IAFLUDSC_COLS)
     &          ,AFLUDSCF(*)
     &          ,ALFA(NUMNP)                ,AREA(NUMEL)
     &          ,ATRA(NUMEL,IDIMATRA)       
     &          ,ATRADSC(IATRADSC_ROWS,IATRADSC_COLS)
     &          ,ATRADSCF(*),BFLU(NUMNP)    ,BCOUPLED(2*NUMNP)
     &          ,BIBI(IDIMBB,NUMEL)         
     &          ,BM_ND_FL(NUMNP,8,2)        ,BM_ND_TT(NUMNP,12,2)
     &          ,BM_ZN_FL(NMAXF,2,NPBFL),BM_ZN_TT(NMAXT,2,NPBTP)
     &          ,BTRA(NUMNP)            ,BUOYANCY(IODIM,LMXNDL,NUMEL)
     &          ,DBUOYANCY(IODIM,LMXNDL*LMXNDL,NUMEL)
     &          ,CAUDAL(NUMNP)              ,CAUX1(NUMNP)
     &          ,CAUX2(NUMNP)
     &          ,CCALAN(NUMNP)
     &          ,CCALIT(NUMNP)              ,CFLU(NUMEL,IDIMCFLU)
     &          ,CFPAREL(NUMEL,NPAREL)      ,CFPARNP(NUMNP,NPARNP)
     &          ,CONCFLOW(NUMNP)            ,COORD(NUMNP,3)
     &          ,COVINV(IDIMCOV)
     &          ,CPREV1(NUMNP)
     &          ,CPREV2(NUMNP),DAT_VD(IODIM,IDIMDQ,NUMEL)
     &          ,DBFLUDFLU(NUMNP,NPBFL)     ,DBFLUDTRA(NUMNP)
     &          ,DELTAITER(NUMNP)
     &          ,DENSITY(IDIMDENS)
     &          ,DERVISC(IDIMDENS)
     &          ,DERC(NUMNP,NPAR,IDIMDERC,NPBTP)
     &          ,DERH(NUMNP,NPARF,IDIMDERH,NPBFL)
     &          ,DFLU(NUMEL,IDIMDFLU)     ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &          ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL)     ,DNODALRH(NUMNP,4)
     &          ,DPARELDC(NPPEL,NUMEL)      ,DPARELDH(NPPEL,NUMEL)
     &          ,DQDFLU(NUMEL,LMXNDL*LMXNDL),DQDTRA(NUMEL,LMXNDL*LMXNDL)
     &          ,DTMXDS(NINT)
     &          ,DTRA(NUMEL,IDIMDTRA)     ,DTRADFLU(NUMEL,LMXNDL*LMXNDL)
     &          ,DTRADTRA(NUMEL,LMXNDL*LMXNDL,NPBTP),DVDC(NUMEL)
     &          ,DVDH(LMXNDL,IODIM,NUMEL)   ,DVDP(IODIM,NPARF,NUMEL)
     &          ,DVOBS(NPAR,3)
     &          ,DWDH(MAX(1,(IOPTS(31)-1)*2*LMXNDL),NUMEL,NPBTP)
     &          
     &          ,FNT(IDIMFNT,NINT)          ,FOBJ_WGT(NSTAT)
     &          ,GP_COORD (6,8,3)
     &          ,GRAD(NPAR)         
     &          ,GRAVEL(IDIMGRAVEL,MIN(3*IDIMGRAVEL,3))
     &          ,GRDFF(IODIM,LMXNDL,NUMEL)  ,HAUX1(NUMNP)
     &          ,HAUX2(NUMNP),HBASE(NUMEL)  
     &          ,HCALAN(NUMNP)              ,HCALIT(NUMNP)
     &          ,HESS(IDIMHESS)             ,HESSAUX(IDIMHESS)
     &          ,HINI(NUMNP)
     &          ,HPREV1(NUMNP)
     &          ,HPREV2(NUMNP)
     &          ,GRADLOC(IODIM,LMXNDL,MAXPG)
     &          ,DLT_PAR(NPAR)              ,PAR_DIR(NPARALG)
     &          ,PAR_INV(NPARALG)           ,PAR_WGT(NTYPAR)
     &          ,PARACD(3,NFNL)             ,PARAUX(NPAR)
     &          ,PARC(NPAR)                 ,PAREL(NUMEL,NPPEL)
     &          ,PARM(NZPAR),PARNP(NUMNP,NPPNP,NPBMX)
     &          ,POINTWEIGHT(MAXPG,NTYPEL)
     &          ,PRGC(NZPRG)                ,QXYZ(IDIMQ,NUMEL)
     &          ,SOLUTION(NUMNP)
     &          ,SOURCE(NUMNP)              ,PARZ(NZPAR)
     &          ,TIME(NINT)                 ,TIT(NUMTIT)
     &          ,TOBS(NUMTOBS,2)            ,VD(IODIM,NUMEL)
     &          ,VISCOSITY(IDIMDENS)        ,VJAC(NUMTOBS,NPAR)
     &          ,VOBS(NUMTOBS)              ,VOBSC(NUMTOBS+NDEVS)
     &          ,WATVOL(MAX(1,(IOPTS(31)-1)*LMXNDL),NUMEL,3,NPBTP)
     &          ,WORK(IDIMWORK)             ,WTOBSN(NUMTNOD)
     &          ,WTOBST(NUMTOBS),XNORVD(NUMEL)
     ;   ,COVPAR(NPAR*(NPAR+1)/2),IPNT_PAR(NZPAR*IDIMWGT)
     ;   ,WGT_PAR(NZPAR*IDIMWGT),IOPT_GS(MXGRPZN,20)
     ;   ,POSMEAS_GS(MXMEASPP_GS,3,2,NGROUP_ZN),COORDGR_GS(6,NGROUP_ZN)
     ;   ,IO_KG_GS(MXGRPZN,16),EXDRZN_GS(MXNZON_GS,4,NGROUP_ZN)
     ;   ,VMEAS_GS(MXMEASPP_GS,IDIMVAR_GS,2,NGROUP_ZN)
     ;   ,TRIM_GS(8,NGROUP_ZN),WGT_UNK(NPAR)
     ;   ,POSZN_GS(MXNZON_GS,3,NGROUP_ZN),PARGOOD(NZPAR)
     ;   ,PARC_GS(MXNPP_GS,NGROUP_ZN)

       CHARACTER::PARNAME(NPAR)*4 

       IF(IFLAGS(3).EQ.1) CALL IO_SUB('INVER',0)

C------------------------- Rewinds inverse problem related files

       IF (IOWRITE(15).NE.0) REWIND(70)
       IF (IOWRITE(14).NE.0) REWIND(71)

C------------------------- Sets up variable IOBMCMP to control the computation
C------------------------- of mass balance

       IF (IOINV.GT.0) THEN
          IOBMCMP=0
       ELSE
          IOBMCMP=1
       ENDIF


C__________OJO!!!!!!!!!! ESTO VA EN LA ENTRADA DE DATOS

       IOMIN=1        ! MARQUARDT UNICO ALGORITMO DE MOMENTO!!!!!!

C_______________________Inicialization for control of convergence

       NUMITER=0
       NITERF1=0
       NITERF2=0
       FNEW=1.D80
       ITTOTV=0        ! Cumulative N.R. iterations (non-linear problems)
       MIN_STOP=0
       GNORM=0.D0
       GNORM1=0.D0
       XMAXIM=0.D0
       PHI = 0.D0

c-parche
      IOCTRA = 0 !se pasa a der_param, pero alli sólo funciona la opción 0
c-fin-parche
C_______________________Starts inverse problem convergence loop

         DO WHILE (MIN_STOP.EQ.0)

            NUMITER=NUMITER+1

C___________________________ GEOSTATISTICS CALCULATIONS

          IF (IOINV_GS.GT.0) THEN
             CALL GEO_CALC
     ;(IDIMCROSS_GS   ,IDIMIVARIO_GS ,IDIMVAR_GS
     ;,IDIMWGT        ,IDIMWORK      ,IDIMZONPP_GS    ,IERROR          
     ;,IOINV          ,LMXNDL        ,MAINF
     ;,MXCLOSE_GS     ,MXDISC_GS     ,MXGRPZN         ,MXKRIG_GS       
     ;,MXMEASPP_GS    ,MXNPP_GS      ,MXNPRIM_GS      ,MXNVAR_GS       
     ;,MXNZON_GS      ,MXROT_GS      ,MXSAM_GS        ,MXSB_GS         
     ;,MXZONPP_GS     ,NFLAGS        ,NGROUP_ZN       ,NPAR            
     ;,NPARDET        ,NPAREL        ,NTYPAR          ,NUMEL           
     ;,NUMNP          ,NUMITER       ,NZPAR           ,NWRITE          
     ;,AREA           ,CLOSESAM_GS   ,COORD           ,COORDGR_GS      
     ;,COVPAR         ,CROSSCOV_GS   ,DATASC_GS       ,ESTKRIG_GS      
     ;,EXDRZN_GS      ,FILENAME      ,ICHECK_GS       ,ICROSSCOV_GS    
     ;,IFLAGS         ,INDPAR        ,INORPAR         ,IO_KG_GS        
     ;,IOLG_PAR       ,IOPT_GS       ,IOWRITE         ,IPNT_PAR        
     ;,IPOLDRIFT_GS   ,ISOZ          ,ISUPBL_GS       ,IVARIO_GS       
     ;,IVPAR          ,IZN_NPP_GS    ,IZN_PP_GS       ,KRIGAUX_GS      
     ;,KRISOL_GS      ,KRISYS_GS     ,KXX             ,LDIM_GS         
     ;,NUMSB_GS       ,NZONE_PAR     ,PAR_WGT         ,PARC            
     ;,PARM           ,PARZ          ,POSDIS_GS       ,POSDISAUX_GS    
     ;,POSMEAS_GS     ,POSZN_GS      ,ROTMAT_GS       ,SEARCH_GS       
     ;,SUPBL_GS       ,TRIM_GS       ,VARIO_GS        ,VMEAS_GS        
     ;,VSTATS_GS      ,WGT_PAR       ,WGT_UNK         ,WORK            
     ;,ZNWGT_GS       ,IACTSIMU      ,IDIMDATASC_GS   ,PARC_GS
     ;,IFLAG_SIMUL    ,COVPAR_GR)


          ENDIF ! IOINV_GS.NE.0

C------------------------- Assigns IOLG_PAR(I,2) to 1 if any zone of 
C------------------------- parameter type I is estimated, otherwise zero

          IF (IOINV.GT.0) THEN
             IPOSIC=1
             DO J=1,19
                IOLG_PAR(J,2)=0
                IF (J.EQ.1 .AND. NZONE_PAR(1).NE.0) THEN
                   DO I=1,MAX (ISOT,IODIM)*NZONE_PAR(1)
                      IF (IVPAR(I,1).NE.0) IOLG_PAR(J,2)=1
                   ENDDO
                   IPOSIC=IPOSIC+MAX (ISOT,IODIM)*NZONE_PAR(1)
                ELSE
                   IF (NZONE_PAR(J).NE.0) THEN
                      DO I=IPOSIC,IPOSIC+NZONE_PAR(J)-1
                         IF (IVPAR(I,1).NE.0) IOLG_PAR(J,2)=1
                      ENDDO
                      IPOSIC=IPOSIC+NZONE_PAR(J)
                   ENDIF
                ENDIF
             ENDDO
          END IF ! Inverse problem ?

C___________________________ Saves first set of parameters

          IF (NUMITER.EQ.1 .AND. IOINV.GT.0) 
     ;        CALL EQUAL_ARRAY(PARGOOD,PARZ,NZPAR)

ccc modificacion; escribimos aqui el PARZ y PARC

          IF (NUMITER == 1 .AND. IOINV.GT.0 .AND. IOINV_GS > 0) THEN
             WRITE(778,*) 'Prior information at Pilot point'
             WRITE(778,*) PARC
          END IF

ccc fin de la modificacion

C_________________________ Initializes array ALFA, containing leakage coeff. 

            CALL ZERO_ARRAY(ALFA,NUMNP)

C_________________________ Rewinds auxiliar file where state variable or 
C_________________________ any other dependent one at all obs. devices and 
C_________________________ at all simulation times are stored

          IF (NUMITER.GT.1.AND. (IOWRITE(5).EQ.2.OR.IOWRITE(6).EQ.2))
     ;        REWIND(81)

          IF (IOINV.GT.0) THEN
             CALL MODIF_MAXITER (IPAR_INV (4))
             IF (IOINV.NE.2) CALL ZERO_ARRAY
     ;                  (DERH,2*NUMNP*NPARF*MAX(1,IOPTS(28)*NPBFL))
             IF (IOINV.NE.1) CALL ZERO_ARRAY
     ;                  (DERC,2*NUMNP*NPAR*MAX(1,IOPTS(29)*NPBTP))
          ENDIF
               
C_______________________Builds some file-variables inherents to the inverse 
C_______________________problem for non-linear case. Inicialization of state
C_______________________variable for Newton-Raphson process

C__________________________________Flow


            CALL INIT_VAR
     ;(  IOINV     ,IOPINITH     ,IOFLLI     ,ISUMFO     ,IPAR_INV(4)
     ;  ,NUMITER   ,NUMNP        ,HCALIT )
C__________________________________Transport

            CALL INIT_VAR
     ;(  IOINV     ,IOPINITC     ,IOTRLI     ,ISUMFO     ,IPAR_INV(4)
     ;  ,NUMITER   ,NUMNP        ,CCALIT )

C_______________________Modifies the Newton-raphson's convergence criteria 
C_______________________as function of Mardquart's convergence process

C________________________________Flow criteria

            IF (ISUMFO.EQ.0. AND .NUMITER.GT.1. AND .IOCRITRAP.NE.0)THEN

              IF (IOFLLI.NE.0) THEN
                CALL MOD_RAPSON_CRITERIA 
     ;(  PAR_DIR(5)  ,PAR_DIR(15)   ,PAR_DIR(19)   ,PAR_DIR(4)
     ;  ,PAR_DIR(16) ,PAR_DIR(20)   ,OBJHED        ,PAR_DIR(13)
     ;  ,PAR_DIR(17) ,PAR_DIR(6)    ,PAR_DIR(14)   ,PAR_DIR(18) )

                DRELMX=PAR_DIR(4)
                DABSMX=PAR_DIR(5)
              ENDIF          
        
C________________________________Transport criteria


              IF (IOTRLI.NE.0) THEN
                CALL MOD_RAPSON_CRITERIA 
     ;(  PAR_DIR(35) ,PAR_DIR(23)   ,PAR_DIR(27)   ,PAR_DIR(34)
     ;  ,PAR_DIR(24) ,PAR_DIR(28)   ,OBJCON        ,PAR_DIR(21)
     ;  ,PAR_DIR(25) ,PAR_DIR(7)    ,PAR_DIR(22)   ,PAR_DIR(26))

C____________whether there is non-linear flow, state
C____________variable change criteria(drelmx dabsmx)
C____________will be assumed in the similar flow's segment

                IF (IOFLLI.EQ.0)THEN
                  DRELMX=PAR_DIR(34)
                  DABSMX=PAR_DIR(35)
                ENDIF
              ENDIF
            ENDIF

C_______________________Computes objective parameters function

          IF (IOINV.GT.0) CALL FUN_PAR
     ;(NPAR     ,OBJPAR   ,COVPAR   ,PARC     ,PARM
     ;,WORK     ,WGT_UNK)

C_________________Renovates objective function
      
           FOLD=FNEW
           FNEW=OBJPAR
           OBJHED=0D0
           OBJCON=0D0
           OBJGEO=0D0

C________________________Opens file:hcal and derh

           IF (IFLAGS(4).EQ.1 .AND. NUMITER.EQ.1) THEN
c-parche
             OPEN(749,FILE='DerH.out',STATUS='UNKNOWN')
	       OPEN(750,FILE='DerC.out',STATUS='UNKNOWN')
c             WRITE(6,*)'P (1) O DELTAP(2)'
c             READ(5,*) I
c             IF (I.EQ.1) THEN
c                OPEN(69,FILE='P.OUT',STATUS='NEW')
c             ELSE IF (I.EQ.2) THEN
c                OPEN(69,FILE='DELTAP.OUT',STATUS='NEW')
c             ENDIF
c-fin-parche
           ENDIF

C------------ Number of geophysical observations and FORPY initialization.           
      IF (ANY(IODEVICE(:,1).EQ.6)) THEN !Geophysical observations.
         NDEVGEO = COUNT(IODEVICE(:,1).EQ.6)
      ELSE
         NDEVGEO = 0
      END IF

      IF (NDEVGEO.GT.0) THEN
          IERROR = FORPY_INITIALIZE()
          IF (IERROR == NO_NUMPY_ERROR) THEN
              WRITE(*,*) "THIS EXAMPLE NEEDS NUMPY..."
              STOP
          ENDIF
      ENDIF

           CALL SIM_JAC
     &(A_COUPL_DSC,A_COUPL_DSCF,ACTH     ,AFLU        ,AFLUDSC  
     &,AFLUDSCF  ,ALFA     
     &,AREA      ,ATRA      ,ATRADSC     ,ATRADSCF ,BCOUPLED
     &,BETAC     ,BFLU      ,BIBI                  ,BM_ND_FL
     &,BM_ND_TT  ,BM_ZN_FL  ,BM_ZN_TT    ,BTRA
     &,BUOYANCY  ,DBUOYANCY    
     &,CAUDAL    ,CAUX1     ,CAUX2       ,CCALAN   ,CCALIT
     &,CFLU      ,CFPAREL   ,CFPARNP     ,CONCFLOW ,COORD    ,CPREV1
     &,CPREV2    ,CREF      ,DAT_VD      ,DBFLUDFLU,DBFLUDTRA 
     &,DELTAITER   ,DENSITY  ,DENSREF
     &,DERVISC   
     &,DERC      ,DERH      ,DFLU        ,DFLUDFLU ,DFLUDTRA ,DNODALRH 
     &,DPARELDC  ,DPARELDH  ,DQDFLU      ,DQDTRA 
     &,DTMXDS    ,DTPREVINV ,DTRA        ,DTRADFLU ,DTRADTRA ,DVDC    
     &,DVDH      ,DVOBS     ,DWDH        ,FILENAME
     &,FNT       ,GP_COORD  ,GRAVEL      ,GRDFF    ,HAUX1   ,HAUX2  
     &,HBASE     ,HCALAN    ,HCALIT      ,HPREV1  ,HPREV2
     &,IAD_S     ,IADD_S    ,IADN_S
     &,IAD_D     ,IADD_D    ,IADN_D
     &,IAFD_S    ,IAFDD_S   ,IAFDN_S
     &,IAFD_D    ,IAFDD_D   ,IAFDN_D
     &,IA_COUPLED_DSC_COLS  ,IA_COUPLED_DSC_ROWS
     &,IAFLUDSC_COLS,IAFLUDSC_ROWS    ,IATRADSC_COLS,IATRADSC_ROWS,IBCOD
     &                      ,IBTCO      ,IDIMAFLU ,IDIMATRA,IDIMBB
     &,IDIMCFLU  ,IDIMDENS  ,IDIMDERC ,IDIMDERH ,IDIMDFLU  ,IDIMDQ
     &,IDIMDTRA  ,IDIMFNT   ,IDIMGRAVEL
     &,IDIMQ     ,IDIMWORK  ,IFLAGS     
     &,INDEXNOD  ,INORPAR   ,IOBMCMP
     &,IOCONSRC  ,IODENS_INI
     &,IODEVICE  ,IODIM
     &,IODIRECT  ,IOEQT     ,IOFLLI      ,IOFLSAT  ,IOFMLF  ,IOFMLT
     &           ,IOINV     ,IOPINVDT    ,IOPTS
     &,IORTS     ,IOTRLI    ,IOTRS       
     &,IOWRITE   ,IPAR_DIR  ,IPARTNER    ,ISOLEQ   ,ISOT    ,ISOZ
     &,ISPARSE   ,ITERGLMX  ,ITPTVAR 
     &,ITYPAFLU  ,ITYPCFLU  ,ITYPDFLU
     &,ITYPATRA  ,ITYPDTRA  
     &,ITYPTRADSC,ITYPACOUPLDSC          ,ITYPFLUDSC
     &,IVPAR
     &,IXPARNP   ,KINT      ,KXX         ,LDIM
     &,LINMET    ,LMXNDL    ,LNNDEL      ,GRADLOC  ,LTYPE   
     &,LXPAREL   ,IPAR_INV(4) ,MAXNB    ,MAXNBF
     &,MAXNN     ,MAXPG
     &,MIN_STOP  ,NBAND     ,NBAND1
     &,NDEVS                ,NFLAGS      ,NFNL
     &,NFNLPAR   ,NFNLPRG
     &,NFNLTIP   ,NFTPAR    ,NINT        ,NMAXF     ,NMAXT  ,NOOBSIT
     &,NOPTS     ,NPAR
     &,NPARALG   ,NPAREL    ,NPARF       ,NPARNP   ,NPBFL   ,NPBMX
     &,NPBTP     ,NPPEL     ,NPPNP       ,NROW     ,NTDMT   ,NTYPAR
     &,NTYPEL    ,NUMEL     ,NUMITER     ,NUMNP    ,NUMTIT
     &,NUMTNOD   ,NUMTOBS   ,NWRITE      ,NZONE_PAR,NZPAR   ,NZPRG
     &,NZTRA     ,PAR_DIR   ,PARACD      ,PARZ     ,PAREL   ,PARNP
     &,POINTWEIGHT,PRGC     ,QXYZ        ,SOLUTION ,SOURCE  ,TIME
     &,TIT       ,TOBS      ,VAR_REF     ,VD       ,VISCOSITY,VISCREF
     &,VJAC      ,VOBSC     ,WATVOL      ,WORK     ,WTOBSN  ,WTOBST
     &,XNORVD    ,DVDP      ,IOLG_PAR    ,IOCTRA   ,HINI,WSPECHEAT
     &,WTHERMCON
     ;,IDIMWGT   ,WGT_PAR  ,IPNT_PAR,IPOS     ,DERIV
     &,NDEVGEO, PARNAME)

           IF (IFLAGS(4).EQ.1 .AND. NUMITER.EQ.1) CLOSE(69)

           IF (IOFLLI.NE.0.AND.IOWRITE(13).NE.0) THEN

             ITTOTV=ITTOTV+ITACUMF
  
             WRITE(MAINF,*)
     ;       ' GENERAL ACCUMULATIVE RAPSON ITERATIONS=',ITTOTV

           ENDIF

C_______________________Selects the minimization algorithm.
C_______________________Only Marquardt's method operative

           IF (IPAR_INV(4).EQ.1 .AND. IFLAGS(4).EQ.1) RETURN


           IF (IOINV.GT.0) THEN

C______________________________ Computes the contribution of observations to
C______________________________ objective function


*       write(84,*) ' objf'
*       write(84,*) VOBS
             CALL OBJ_VAR
     ;(FNEW     ,FOLD     ,IDIMCOV  ,ISUMFO   ,NDEVS
     ;,NSTAT    ,NUMTOBS  ,OBJCON   ,OBJFLO   ,OBJGEO  ,OBJHED   ,OBJHUM
     ;,COVINV   ,FOBJ_WGT ,VOBS     ,VOBSC    ,MEASTYP)

             CALL MINIMIZATION
     ;(ALF      ,FNEW     ,FOLD     ,GNORM     ,GNORM1    ,IDIMCOV
     ;,IOMIN    ,ISUMFO   ,MAINF    ,MIN_STOP  ,NBANDCOV  ,NDEVS
     ;,NFLAGS   ,NITERF1  ,NITERF2  ,NPAR      ,NPARALG   ,NSTAT
     ;,NUMITER  ,NUMTOBS  ,NWRITE   ,NZPAR     ,OBJCON    ,OBJHED    
     ;,OBJPAR   ,XMAXIM   ,COVINV   ,COVPAR    ,DLT_PAR   ,FOBJ_WGT  
     ;,GRAD     ,HESS     ,HESSAUX  ,IFLAGS    ,IOWRITE   ,IPAR_INV
     ;,PAR_INV  ,PARAUX   ,PARC      ,PARGOOD   ,PARM     ,PARZ
     ;,VJAC     ,VOBS     ,VOBSC     ,WORK      ,WGT_UNK  ,MEASTYP
     ;, PHI     ,OBJGEO)

C______________________No convergence

             IF (MIN_STOP.EQ.0) THEN    ! No convergence yet 
           
C_______________________Reduce the parameters according to the maximum

                CALL MAX_PAR
     ;(ALF     ,MAINF    ,NFLAGS  ,NPAR    ,PAR_INV(9) ,PAR_INV(8)
     ;,XMAXIM  ,DLT_PAR  ,INDPAR  ,IFLAGS  ,PARAUX)

C_______________________ Assign the new computed parameters (updates PARC)

                DO IPAR=1,NPAR
                   PARC(IPAR)=PARC(IPAR)+DLT_PAR(IPAR)
                END DO

C_______________________ Updates the zonal parameters

                CALL UPDATE_PARZ
     ;(IDIMWGT    ,NPAR    ,NZPAR    ,DLT_PAR    ,IPNT_PAR    
     ;,IVPAR      ,PARZ    ,WGT_PAR)

C------------------------- Computes matrix diffusion dimensionless terms
C------------------------- They can change as a consequence of a change
C------------------------- in model parameters (tipically after an inv.
C------------------------- prob. iteration). Dependence on time will be 
C------------------------- treated later, within the transient times loop.

               IF (IOEQT.NE.1.AND.NTDMT.NE.0) 
     ;            CALL UPD_PAR_DMT (PARZ,NZPAR)

             END IF    ! MIN_STOP.EQ.0
           END IF      ! IOINV.GT.0

C------------------------- Rewinds auxiliar files

           IF (IOEQT.NE.2 .AND. IOWRITE(7).NE.0) REWIND(94)
           IF (IOEQT.NE.1 .AND. IOWRITE(8).NE.0) REWIND(93)
           IF (IOEQT.NE.2) REWIND(52)

           IF (IOINV.LE.0) MIN_STOP=1         ! Leaves the loop

         END DO !General convergence loop

C______________________________ Writes the parameters' history

         IF (IOINV.GT.0) CALL WRI_PARAM_HISTORY
     ;(IOWRITE(15) ,ISOT     ,MAINF  ,NTYPAR
     ;,NZPAR       ,INORPAR  ,IVPAR  ,NZONE_PAR    ,PARGOOD)

C______________________________ Last simulation to compute mass balance
C______________________________ and state variable (possibly last iteration
C______________________________ was bad. Lastgood set of parameters was 
C______________________________ recovered in Marquardt's subroutine)

       IF (IOINV.GT.0) THEN
          IOBMCMP=1                     ! Mass balance calculation option 

C_________________________ Rewinds observation file

          IF (IOWRITE(5).EQ.2.OR.IOWRITE(6).EQ.2) REWIND(81)

C_________________________ Initializes array ALFA, containing leakage coeff. 

          CALL ZERO_ARRAY(ALFA,NUMNP)

           CALL SIM_JAC
     &(A_COUPL_DSC,A_COUPL_DSCF,ACTH     ,AFLU        ,AFLUDSC  
     &,AFLUDSCF  ,ALFA     
     &,AREA      ,ATRA      ,ATRADSC     ,ATRADSCF ,BCOUPLED
     &,BETAC     ,BFLU      ,BIBI                  ,BM_ND_FL
     &,BM_ND_TT  ,BM_ZN_FL  ,BM_ZN_TT    ,BTRA
     &,BUOYANCY  ,DBUOYANCY    
     &,CAUDAL    ,CAUX1     ,CAUX2       ,CCALAN   ,CCALIT
     &,CFLU      ,CFPAREL   ,CFPARNP     ,CONCFLOW ,COORD    ,CPREV1
     &,CPREV2    ,CREF      ,DAT_VD      ,DBFLUDFLU,DBFLUDTRA 
     &,DELTAITER   ,DENSITY  ,DENSREF
     &,DERVISC   
     &,DERC      ,DERH      ,DFLU        ,DFLUDFLU ,DFLUDTRA ,DNODALRH 
     &,DPARELDC  ,DPARELDH  ,DQDFLU      ,DQDTRA 
     &,DTMXDS    ,DTPREVINV ,DTRA        ,DTRADFLU ,DTRADTRA ,DVDC    
     &,DVDH      ,DVOBS     ,DWDH        ,FILENAME
     &,FNT       ,GP_COORD  ,GRAVEL      ,GRDFF    ,HAUX1   ,HAUX2  
     &,HBASE     ,HCALAN    ,HCALIT      ,HPREV1  ,HPREV2
     &,IAD_S     ,IADD_S    ,IADN_S
     &,IAD_D     ,IADD_D    ,IADN_D
     &,IAFD_S    ,IAFDD_S   ,IAFDN_S
     &,IAFD_D    ,IAFDD_D   ,IAFDN_D
     &,IA_COUPLED_DSC_COLS  ,IA_COUPLED_DSC_ROWS
     &,IAFLUDSC_COLS,IAFLUDSC_ROWS    ,IATRADSC_COLS,IATRADSC_ROWS,IBCOD
     &                      ,IBTCO      ,IDIMAFLU ,IDIMATRA,IDIMBB
     &,IDIMCFLU  ,IDIMDENS  ,IDIMDERC ,IDIMDERH ,IDIMDFLU  ,IDIMDQ
     &,IDIMDTRA  ,IDIMFNT   ,IDIMGRAVEL
     &,IDIMQ     ,IDIMWORK  ,IFLAGS     
     &,INDEXNOD  ,INORPAR   ,IOBMCMP
     &,IOCONSRC  ,IODENS_INI
     &,IODEVICE  ,IODIM
     &,IODIRECT  ,IOEQT     ,IOFLLI      ,IOFLSAT  ,IOFMLF  ,IOFMLT
     &           ,IOINV     ,IOPINVDT    ,IOPTS
     &,IORTS     ,IOTRLI    ,IOTRS       
     &,IOWRITE   ,IPAR_DIR  ,IPARTNER    ,ISOLEQ   ,ISOT    ,ISOZ
     &,ISPARSE   ,ITERGLMX  ,ITPTVAR 
     &,ITYPAFLU  ,ITYPCFLU  ,ITYPDFLU
     &,ITYPATRA  ,ITYPDTRA  
     &,ITYPTRADSC,ITYPACOUPLDSC          ,ITYPFLUDSC
     &,IVPAR
     &,IXPARNP   ,KINT      ,KXX         ,LDIM
     &,LINMET    ,LMXNDL    ,LNNDEL      ,GRADLOC  ,LTYPE   
     &,LXPAREL   ,IPAR_INV(4) ,MAXNB    ,MAXNBF
     &,MAXNN     ,MAXPG
     &,MIN_STOP  ,NBAND     ,NBAND1
     &,NDEVS                ,NFLAGS      ,NFNL
     &,NFNLPAR   ,NFNLPRG
     &,NFNLTIP   ,NFTPAR    ,NINT        ,NMAXF     ,NMAXT  ,NOOBSIT
     &,NOPTS     ,NPAR
     &,NPARALG   ,NPAREL    ,NPARF       ,NPARNP   ,NPBFL   ,NPBMX
     &,NPBTP     ,NPPEL     ,NPPNP       ,NROW     ,NTDMT   ,NTYPAR
     &,NTYPEL    ,NUMEL     ,NUMITER     ,NUMNP    ,NUMTIT
     &,NUMTNOD   ,NUMTOBS   ,NWRITE      ,NZONE_PAR,NZPAR   ,NZPRG
     &,NZTRA     ,PAR_DIR   ,PARACD      ,PARZ     ,PAREL   ,PARNP
     &,POINTWEIGHT,PRGC     ,QXYZ        ,SOLUTION ,SOURCE  ,TIME
     &,TIT       ,TOBS      ,VAR_REF     ,VD       ,VISCOSITY,VISCREF
     &,VJAC      ,VOBSC     ,WATVOL      ,WORK     ,WTOBSN  ,WTOBST
     &,XNORVD    ,DVDP      ,IOLG_PAR    ,IOCTRA   ,HINI,WSPECHEAT
     &,WTHERMCON
     ;,IDIMWGT   ,WGT_PAR  ,IPNT_PAR,IPOS     ,DERIV
     &,NDEVGEO, PARNAME)
       END IF

      IF (NDEVGEO.GT.0) THEN
          CALL FORPY_FINALIZE
      ENDIF
       IF(IFLAGS(3).EQ.1) CALL IO_SUB('INVER',1)
       RETURN
       END
