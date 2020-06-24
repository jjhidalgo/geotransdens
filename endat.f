       SUBROUTINE ENTDAT 
     ;(EPSFLU   ,EPSTRA   ,IDIMCOV  ,IDIMFNT  ,IDIMQ    ,IERROR   
     ;,IOCNSF   ,IOCNST   ,IODIM    ,IOEQT    ,IOFLLI
     ;,IOFLSAT  ,IOFMLF   ,IOFMLT   ,IOINV    ,IORTS
     ;,IOTRLI   ,IOTRS              ,ISOT     ,LMXNDL   ,MAINF    
     ;,NBAND    ,NDEVS    ,NFNL     ,NFNT     ,NINT     ,NPAR     
     ;,NPAREL   ,NPARF    ,NPARFPRG ,NPARNP             ,NPBFL    
     ;,NPBMX    ,NPBTP    ,NROW     ,NTDMT    ,NTYPAR   ,NUMEL    
     ;,NUMNP    ,NUMTIT   ,NUMTNOD  ,NUMTOBS  ,NWRITE   ,NZDMT    
     ;,NZPAR    ,THETAF   ,THETAT   ,ACTH     ,AREA     
     ;,BTRA     ,BUDAT    ,CAUDAL   ,CCALIT   ,CFPAREL  
     ;,CFPARNP  ,COORD    ,COVINV   ,COVPAR   ,DTMXDS   ,EXTNBU   
     ;,FNT      ,GRAV     ,HCALIT   ,HCALAN   ,IBCOD    ,IBTCO    
     ;,IFLAGS   ,INDEXNOD ,INDPAR   ,INORPAR  ,IOBUTYP  ,IOCALBU  
     ;,IODEVICE ,IOLG_PAR ,IOTINT   ,MEASTYP  ,IOUTYP   ,IOWRITE ,ISOLEQ   
     ;,ISOZ     ,IVPAR    ,IXPARNP  ,KINT     ,KXX      ,LDIM     
     ;,LNNDEL   ,LTYPE    ,LXPAREL  ,NBUW     ,NFNLPAR  ,NFNLPRG  
     ;,NFNLTIP  ,NFTPAR   ,NOBUF    ,NOOBSIT  ,NZONE_PAR,PARACD   
     ;,PARZ     ,PARM     ,QXYZ     ,STPAR    ,TIME     ,TIT      
     ;,TOBS     ,VD       ,VOBS     ,WTOBSBU  ,WTOBSN   ,WTOBST   
     ;,WTOBSU   ,XNORVD   ,DEVNAME  ,FILENAME
     ;,IDMBLCVP ,NBLCVP   ,NFLAGS
     ;,IO_KG_GS         ,IOPT_GS      ,MXDISC_GS  ,MXNZON_GS
     ;,NGROUP_ZN        ,MXGRPZN      ,DIVZN_GS   ,SUPBL_GS
     ;,POSDIS_GS        ,IPOLDRIFT_GS,SEARCH_GS
     ;,TRIM_GS          ,EXDRZN_GS    ,MXMEASPP_GS,IDIMVAR_GS
     ;,MXNPP_GS         ,MXNVAR_GS    ,POSMEAS_GS ,VMEAS_GS
     ;,VSTATS_GS        ,VARIO_GS
     ;,IVARIO_GS        ,IDIMIVARIO_GS,IDIMWGT    ,IPNT_PAR
     ;,WGT_PAR          ,IOINV_GS     ,POSZN_GS   ,COORDGR_GS
     ;,NPARDET          ,LDIM_GS      ,WGT_UNK    ,PARC ,PAR_WGT
     ;,IOSMFL       ,IOSMTP     ,PARC_GS
     ;!NUEVOS
     ;,IODENS_INI,ITPTVAR,    BETAC,CREF,DENSF,TEMPF,VISCREF,WSPECHEAT
     &,WTHERMCON)
*****************************************************************************
*
* PURPOSE
*
*     Reads all input data
*
* DESCRIPTION
*
*     Reads all input data and stores the dimension of every element on
*     array LDIM
*
* EXTERNAL VARIABLES: ARRAYS
*
*  ACTH                   Aquifer thickness of every element. Cross sectional   
*                         area for 1-D elements, thickness for 2-D elements.    
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  BTRA                   Right hand side of transport discretized equation     
*  CAUDAL                 Input/output flow at every node.                      
*  CCAL                   Computed concentration at every node                  
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  CFPARNP                Array containing node coefinient of node j            
*                         corresponding to INpar index parameter zone.          
*  COBS                   Measured concentrations at the observation points.    
*  COORD                  Nodal coordinates                                     
*  COVINV                 Inverse of the covariance matrix                      
*  DEVNAME                Device name                                           
*  DTMXDS                 Maximum allowed time increment for nonlinear          
*                         problems                                              
*  EXTNBU                 Measure (length, area or volume) of basic unit        
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  FNT                    Array containing time functions values                
*  GRAV                   Gravity array direction                               
*  HCAL                   Computed heads at every node                          
*  IBCOD                  Flow boundary condition index                         
*  IBTCO                  Transport boundary condition index                    
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INDEXNOD               Index relating nodes                                  
*  INDPAR                 Array of 0's and 1's. 0 means that the parameter is   
*                         estimated aritmethically, otherwise logarithmically.  
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARZ, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IOCALBU                Calculation method for basic unit                     
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Status for calc. of obs.                    
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*  IOLG_PAR               Array containing all logarithmic options of           
*                         estimated parameters                                  
*  IOWRITE                Array containing all output options                   
*  ISOLEQ                 Array containing the type of head/concentration
*                         solution desired for the user at each obs. time
*  ISOZ                   Anisotropy of every transmissivity zone               
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  IXPARNP                Array containing zone number of each node j,          
*                         corresponding to INpar index parameter zone.          
*  KINT                   Number of solution time increments, between           
*                         successive observation times.                         
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  NBUW                   Number of basic unit weights                          
*  NFNLPAR                Vector containing non-linear function order           
*                         afecting every parameter at each zone.                
*  NFNLPRG                Generic parameter zone number for every nonlinear     
*                         function                                              
*  NFNLTIP                Type of non-linear function                           
*  NFTPAR                 Vector containing time function number at every       
*                         parameter zone                                        
*  NOBUF                  Number of first basic unit                            
*  NOOBSIT                Observation number to which an integration time       
*                         belongs to                                            
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PARACD                 Agreement parameters                                  
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*  PARM                   Vector containing measured values for all             
*                         parameters                                            
*  QXYZ                   Products between the different components of          
*                         Darcy's velocity divided by its norm                  
*  STPAR                  Vector containing standard deviation errors of        
*                         all parameters prioo information                      
*  TIME                   Observation times.                                    
*  TIT                    Integration time                                      
*  TOBS                   Time of observation                                   
*  VD                     Darcy's velocity                                      
*  VOBS                   Observation value                                     
*  WTOBSN                 Weight for node used to calculate observation         
*  WTOBST                 Weight for integration time                           
*  XNORVD                 Euclidean norm of Darcy's velocity                    
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  EPSFLU                 Time weighting parameter for nonlinear flow problems  
*  EPSTRA                 Time weighting parameter for nonlinear transport      
*                         problems                                              
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  IDIMQ                  Used to dimension array QXYZ                          
*  IERROR                 Current number of errors on input data                
*  IOALFT                 If zero, leakage do not vary in time
*  IOCNSF                 Scheme for storage term in flow problem               
*  IOCNST                 Scheme for mass storage term in transport problem     
*  IODIM                  Maximum dimension of any element included             
*                         in the problem                                        
*  IOEQT                  Type of problem to be solved                          
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOFLSAT                Indicates the possibility that one part of the domain 
*                         reaches unsaturated state.                            
*  IOFMLF                 Flow Formulation number                               
*  IOFMLT                 Transport formulation number                          
*  IOINV                  Inverse problem option                                

*  IORTS                  Transport regime
*  IOSMFL                 If 1, flow problems are simultaneous, otherwise 0
*  IOSMTP                 If 1, transport problems are simultaneous, otherwise 0
*  IOTRLI                 Idem to IOFLLI, in the case of transport.             
*  IOTRS                  Flow regime                                           
*  ISOT                   Maximum hydraulic conductivity tensor anisotropy      
*                         degree in the problem.                                
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NBAND                  Half Bandwith (maximum difference between the         
*                         numbers of two nodes belonging to the same element)   
*  NDEVS                  Number of devices
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NFNT                   Number of time functions used for describing time     
*                         dependence of all transient parameters                
*  NINT                   Number of observation times                           
*  NPAR                   Total number of parameters to be estimated            
*  NPAREL                 Number of element parameters in current problem       
*  NPARF                  Number of transient parameters to be estimated        
*  NPARFPRG               Number of uncertain generic parameter zones involved  
*                         in the non-linear flow inverse problem                
*  NPARNP                 Number of nodal parameters in current problem         
*  NPBFL                  Number of simultaneous flow problems                  
*  NPBMX                  MAX(NPBFL,NPBTP)
*  NPBTP                  Number of simultaneous transport problems             
*  NROW                   Current record number                                 
*  NTDMT                  Matrix diffusion option
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NUMTIT                 Total number of integration times                     
*  NUMTNOD                Total number of nodes used for calculating obs.       
*  NUMTOBS                Total number of observations                          
*  NWRITE                 Number of output options (used for dimensioning)      
*  NZDMT                  Number of matrix diffusion zones                      
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  THETAF                 Time weighting parameter for flow problems            
*  THETAT                 Time weighting parameter for transport problems       
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  AREA_ELEM                                                                    
*  ENDATINICOND           Reads initial conditions                              
*  ENDATOBS               Reads all variables related to measurements           
*  ENTDATIX_COOR                                                                
*  ENTDATIX_ZON                                                                 
*  ENTDATLX_ELEM                                                                
*  ENTDATLX_ZON                                                                 
*  ENTDATNZ               Reads zone parameters                                 
*  ENTDAT_DMT                                                                   
*  ENTDAT_FUNNOLI_GRAV    Reads generic parameter zone number, agreement        
*                         parameters and/or gravity direction                   
*  ENTFLOW                Reads flow                                            
*  ENTVEL                 Reads Darcy's velocity if the flow problem is not     
*                         solved                                                
*  IO_SUB                 Writes a message at the end or the beginning of
*                         a subroutine
*  LDIMEN                 Computes the dimension of a given element             
*  LEC_CFE                Reads elements coeficients                            
*  LEC_CFN                Reads node coeficients                                
*  LEC_FT                 Reads time and time functions from TIM file           
*
* HISTORY
*
*     AMS        1988     First coding
*     SCR      4-1997     Revision and verification
*     AMS      1-1998     Common elimination and addition of header
*     AAR      7-2001     Revision and inclusion of observations related var.
*
*****************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       CHARACTER FILENAME(20)*20,DEVNAME(NDEVS)*10

       DIMENSION 
     ;  ACTH(NUMEL), BTRA(NUMNP), CAUDAL(NUMNP), CCALIT(NUMNP),
     ;  CFPAREL(NUMEL,NPAREL), CFPARNP(NUMNP,NPARNP), 
     ;  COORD(NUMNP,3),DTMXDS(NINT) ,FNT(IDIMFNT,NINT), GRAV(3),
     ;  IBCOD(NUMNP,NPBMX), IBTCO(NUMNP,NPBMX), IFLAGS(NFLAGS),
     ;  INDPAR(NPAR),INORPAR(NPAR),IOLG_PAR(NTYPAR,2),IOWRITE(NWRITE),
     ;  IVPAR(NZPAR,4), IXPARNP(NUMNP,NPARNP,NPBMX), 
     ;  KINT(NINT), KXX(LMXNDL,NUMEL), LDIM(NUMEL), LNNDEL(NUMEL), 
     ;  LTYPE(NUMEL),LXPAREL(NUMEL,NPAREL,NPBMX), NFNLPAR(NZPAR), 
     ;  NFNLPRG(8,NFNL),NFNLTIP(8,NFNL), NFTPAR(NZPAR), HCALIT(NUMNP),
     ;  NZONE_PAR(NTYPAR),PARACD(3,NFNL), PARZ(NZPAR), PARM(NZPAR), 
     ;  QXYZ(IDIMQ,NUMEL),STPAR(NZPAR), TIME(NINT),
     ;  VD(IODIM,NUMEL), XNORVD(NUMEL), ISOLEQ(NINT,4),covinv(idimcov),
     ;  COVPAR(NPAR*(NPAR+1)/2),MEASTYP(NUMTOBS),INDEXNOD(NUMTNOD),
     ;  IOPT_GS(MXGRPZN,20),WGT_UNK(NPAR),
     ;  IO_KG_GS(MXGRPZN,16),PARC(NPAR),PARC_GS(MXNPP_GS,NGROUP_ZN),
     ;  LDIM_GS(NUMEL)



C------------------------- FIRST EXECUTABLE STATEMENT.
C------------------------- Output message to know that the program is already
C------------------------- inthis routine

       IF (IFLAGS(3).EQ.1) CALL IO_SUB ('ENTDAT',0)

C------------------------- Reads node numbers,coordinates,IX and boundary conditions
C------------------------- Starts grid definition file

       NROW=0         
       CALL ENTDATIX_COOR
     ; (IERROR   ,IDALF    ,IDALFT   ,IDCHP    ,IDCHPT   ,IDCON
     ; ,IDCONT   ,IDDMT    ,IOWRITE(1),IDQQP   ,IDQQPT   ,IOEQT
     ; ,IOFLLI   ,IOTRLI   ,IOWRITE(2),11      ,MAINF
     ; ,NROW     ,NTDMT    ,NUMNP    ,BTRA     ,CCALIT
     ; ,FILENAME ,HCALAN   ,COORD(1,1),COORD(1,2),COORD(1,3), IDCLK)

C------------------------- Stores the value of IOEQT

       IOEQT_ORIG=IOEQT

       DO IPROB=1,NPBMX                             !  MAX(NPBFL,NPBTP)
          
C------------------------- Location in array IXPARNP

          INCHP =1
          INCHPT=2
          INQQP =3
          INQQPT=4
          INALF =5
          INALFT=6
          INCON =7
          INCONT=8
          INDMT =9
	  INCLK = 10

C------------------------- Next assignations are to avoid reading of flow or 
C------------------------- transport data when the number of problems of this
C------------------------- type is surpassed

          IF (IPROB.GT.NPBFL) THEN
             IOEQT=2
             IPBFL=1
          ELSE
             IPBFL=IPROB
          ENDIF
          IF (IPROB.GT.NPBTP) THEN
             IOEQT=1
             IPBTP=1
          ELSE
             IPBTP=IPROB
          ENDIF

          CALL ENTDATIX_ZON
     ; (IDALF    ,IDALFT   ,IDCHP    ,IDCHPT   ,IDCON    ,IDCONT
     ; ,IDDMT    ,IDQQP    ,IDQQPT   ,IERROR   ,INALF    ,INALFT
     ; ,INCHP    ,INCHPT   ,INCON    ,INCONT   ,INDMT    ,IOWRITE(1)
     ; ,INQQP    ,INQQPT   ,IOEQT    ,IORTS    ,IOTRS
     ; ,IOWRITE(2),IPROB   ,11       ,MAINF    ,NPARNP   ,NROW
     ; ,NTDMT    ,NUMNP    ,NZONE_PAR(6),NZONE_PAR(4),NZONE_PAR(13)
     ; ,NZONE_PAR(16)      ,NZONE_PAR(5)       ,FILENAME
     ; ,IBCOD(1,IPBFL)     ,IBTCO(1,IPBTP)     ,IXPARNP(1,1,IPROB)
     & ,IDCLK,NZONE_PAR(18),INCLK)

       ENDDO

C------------------------- Recovers the value of IOEQT

       IOEQT=IOEQT_ORIG

C------------------------- Reads elements number KXX,LX

         CALL ENTDATLX_ELEM
     ; (IERROR   ,IOWRITE(1),IOEQT   ,IOWRITE(2),11      ,LDARR
     ; ,LDARRT   ,LDCOE    ,LDCRD    ,LDDFM    ,LDDSP    ,LDFOD
     ; ,LDPOR    ,LDSTG    ,LDTRA    ,LMXNDL   ,MAINF    ,NBAND
     ; ,NROW     ,NUMEL    ,NUMNP    ,ACTH     ,FILENAME ,KXX
     ; ,LNNDEL   ,LTYPE    ,COORD(1,1),COORD(1,2) )

       DO IPROB=1,NPBMX                             !  MAX(NPBFL,NPBTP)
          
C------------------------- Next assignations are to avoid reading of flow or 
C------------------------- transport data when the number of problems of this
C------------------------- type is surpassed

          IF (IPROB.GT.NPBFL) IOEQT=2
          IF (IPROB.GT.NPBTP) IOEQT=1

          CALL ENTDATLX_ZON
     ; (IERROR   ,IOWRITE(1),IOEQT   ,IOFLSAT  ,IOTRS    ,IOWRITE(2)
     ; ,IPROB    ,11       ,LDARR    ,LDARRT   ,LDCOE    ,LDCRD
     ; ,LDDFM    ,LDDSP    ,LDFOD    ,LDPOR    ,LDSTG    ,LDTRA
     ; ,MAINF    ,NROW     ,NUMEL    ,NZONE_PAR(3)
     ; ,NZONE_PAR(13) ,NZONE_PAR(12) ,NZONE_PAR(9) ,NZONE_PAR(7) 
     ; ,NZONE_PAR(11) ,NZONE_PAR(10) ,NZONE_PAR(2) ,NZONE_PAR(1)
     ; ,FILENAME ,LDIM     ,LTYPE    ,LXPAREL(1,3,IPROB) 
     ; ,LXPAREL(1,4,IPROB) ,LXPAREL(1,10,IPROB) ,LXPAREL(1,9,IPROB) 
     ; ,LXPAREL(1,6,IPROB) ,LXPAREL(1,5,IPROB)  ,LXPAREL(1,8,IPROB) 
     ; ,LXPAREL(1,7,IPROB) ,LXPAREL(1,2,IPROB)  ,LXPAREL(1,1,IPROB)
     & ,LXPAREL(1,12,IPROB) )

C------------------------- Stores transmissivity zone dimensions

*** OJO, solo funciona con 1 problema!!!!!!!!!!

          IF (IOINV_GS.GT.0) THEN
             DO I=1,NZONE_PAR(1)
                LDIM_GS(I)=LDIM(I)
             END DO
          END IF

       ENDDO !IPROB=1,MAX(NPBFL,NPBTP)

C------------------------- Recovers the value of IOEQT

       IOEQT=IOEQT_ORIG

C------------------------- Reads node coeficients

       NROW=0
       CALL LEC_CFN
     ; (IERROR       ,5            ,6            ,1
     ; ,2            ,7            ,8            ,IOWRITE(1)
     ; ,3            ,4            ,IOEQT        ,IORTS
     ; ,IOTRS        ,IOWRITE(2)   ,12
     ; ,MAINF        ,NPARNP       ,NROW         ,NUMNP
     ; ,NZONE_PAR(6) ,NZONE_PAR(4) ,NZONE_PAR(13),NZONE_PAR(5)
     ; ,CFPARNP      ,FILENAME     ,IBCOD        ,IBTCO
     & ,INCLK        ,NZONE_PAR(18))

C------------------------- Reads elements coeficients

       CALL LEC_CFE
     ; (IERROR       ,3            ,4            ,10
     ; ,9            ,6            ,5            ,8
     ; ,7            ,IOWRITE(1)   ,2            ,1
     ; ,IOEQT        ,IOFLSAT      ,IOTRS        ,IOWRITE(2)
     ; ,12           ,MAINF        ,NPAREL       ,NROW
     ; ,NUMEL        ,NZONE_PAR(3) ,NZONE_PAR(13),NZONE_PAR(12)
     ; ,NZONE_PAR(9) ,NZONE_PAR(7) ,NZONE_PAR(11),NZONE_PAR(10)
     ; ,NZONE_PAR(2) ,NZONE_PAR(1) ,CFPAREL      ,FILENAME
     & ,11           ,NZONE_PAR(19))

C------------------------- Reads zone parameters

       CALL ENTDATNZ
     ; (IDIMWGT   ,IERROR   ,IOWRITE(1),IODIM   ,IOEQT     ,IOFLLI   
     ; ,IOFLSAT   ,IOINV    ,IOTRLI   ,IORTS    ,IOTRS     ,IOWRITE(2)
     ; ,ISOT      ,12       ,MAINF    ,MXGRPZN  ,NFNL      ,NPAR     
     ; ,NPARF     ,NPARFPRG ,NROW     ,NTYPAR   ,NUMEL     ,NZPAR    
     ; ,FILENAME  ,INDPAR   ,INORPAR  ,IOLG_PAR ,IOPT_GS   ,IPNT_PAR  
     ; ,ISOZ      ,IVPAR    ,LDIM     ,NFNLPAR  ,NFTPAR    ,NZONE_PAR 
     ; ,PAR_WGT   ,PARC     ,PARM     ,PARZ     ,STPAR     ,WGT_PAR   
     ; ,NGROUP_ZN ,NPARDET  ,IOINV_GS ,WGT_UNK)    


C------------------------- Reads matrix diffusion zones

       IF (IOEQT.NE.1. AND .NTDMT.NE.0) CALL ENTDAT_DMT
     ; (NZDMT    ,IVPAR    ,PARZ     ,INORPAR   ,NZONE_PAR
     ; ,NTYPAR   ,NZPAR   ,IXPARNP(1,1,1)       ,NUMNP
     ; ,NPAR     ,IOWRITE(1)         ,MAINF     ,IOINV
     ; ,NPBMX    ,NPARNP)

C------------------------- Reads nonlinear parameters and/or gravity direction

       IF(IOFLLI.NE.0.OR.IOTRLI.NE.0 .OR. IODENS_INI.EQ.1)
     ;   CALL ENTDAT_FUNNOLI_GRAV
     ;       (IERROR     ,IOWRITE(1) ,IODENS_INI ,IOFLLI     ,IOTRLI
     ;       ,IOWRITE(2) ,12         ,MAINF      ,NFNL       ,NROW
     &       ,NZONE_PAR(14)          ,FILENAME   ,GRAV       ,NFNLPRG
     7       ,NFNLTIP  ,PARACD)


       IF(IODENS_INI.EQ.1 .OR. ITPTVAR.EQ.1) THEN

          CALL ENTDAT_DENSITY
     &        (BETAC    ,CREF     ,DENSF    ,FILENAME ,IERROR
     &        ,IOWRITE(1)         ,IODENS_INI         ,IOWRITE(2)
     &        ,ITPTVAR  ,12       ,MAINF    ,NROW     ,TEMPF
     &        ,VISCREF  ,WSPECHEAT,WTHERMCON)

      ELSE

          BETAC = 0D0
          CREF  = 0D0
          DENSF = 1D0
          VISCREF = 0D0
          TEMPF = 0D0
          WSPECHEAT = 0D0
          WTHERMCON = 0D0

	  END IF !IODENS_INI.EQ.1 .OR. ITPTVAR.EQ.1

C------------------------- Reads time and time functions if a transient 
C------------------------- problem is solved

       IF (IOTRS+IORTS.NE.0) THEN

C------------------------- Starts record counter 

          NROW=0

          CALL LEC_FT 
     ; (EPSFLU   ,EPSTRA   ,IDIMFNT  ,IERROR   ,IOWRITE(1),IOCNSF
     ; ,IOCNST   ,IOEQT    ,IOFLLI   ,IOFMLF   ,IOFMLT    ,IORTS
     ; ,IOTRLI   ,IOTRS    ,IOWRITE(2),13      ,MAINF     ,NFNT
     ; ,NINT     ,NROW     ,THETAF   ,THETAT   ,DTMXDS    ,FILENAME
     ; ,FNT      ,ISOLEQ   ,KINT     ,TIME)

       ELSE

C------------------------- Assigns ISOLEQ array at time 1 in steady-state 
C------------------------- case

          ISOLEQ(1,1)=1           ! Steady flow
          ISOLEQ(1,2)=1           ! Steady transport 
          ISOLEQ(1,3)=1           ! Number of flow problem
          ISOLEQ(1,4)=1           ! Number of transport problem

       ENDIF

C------------------------- Reads initial conditions
C------------------------- Starts initial condition file

       NROW=0 

       NFL_SIM=MAX(1,IOSMFL*NPBFL)  ! # of simultaneous flow pb
       NTP_SIM=MAX(1,IOSMTP*NPBTP)  ! # of simultaneous tpt pb

       IF (IOTRS.EQ.1.OR.IORTS.EQ.1) THEN
          CALL ENDATINICOND 
     & (IERROR   ,IOWRITE(1),IORTS   ,IOTRS    ,IOWRITE(2),15
     & ,MAINF    ,NFL_SIM   ,NTP_SIM ,NUMNP    ,CCALIT      ,FILENAME
     & ,HCALIT)
       END IF

C------------------------- Sets LDIM array containing elements dimensions

       DO L=1,NUMEL
          ITIP=LTYPE(L)
          LDIM(L)=LDIMEN(ITIP)
       END DO

       IF (IOEQT.EQ.2) THEN
          CALL ENTVEL 
     ; (IDIMQ    ,IERROR   ,IOWRITE(1),IODIM   ,IOWRITE(2),15
     ; ,MAINF    ,NROW     ,NUMEL    ,ACTH     ,FILENAME ,LDIM
     ; ,LTYPE    ,QXYZ     ,VD       ,XNORVD)

          CALL ENTFLOW
     ; (IERROR   ,IOWRITE(1),IOWRITE(2),15  ,MAINF    ,NROW
     ; ,NUMNP    ,CAUDAL   ,FILENAME)

       ENDIF

C------------------------- Reads observation variables. First, computation of
C------------------------- element size is done, as it is needed for spatial
C------------------------- averaging of observation devices

       CALL AREA_ELEM
     ;(IERROR     ,IOWRITE(2) ,LMXNDL   ,MAINF      ,NUMEL    
     ;,NUMNP      ,AREA       ,KXX      ,LTYPE      ,COORD(1,1) 
     ;,COORD(1,2) ,COORD(1,3) ,FILENAME)


       IF (IOINV.GT.0.OR.IOWRITE(5).GT.0.OR.IOWRITE(6).GT.0
     ;               .OR.IOWRITE(11).NE.0.OR.IOWRITE(12).NE.0) THEN

          NROW=0 
          CALL ENDATOBS
     ;(IDIMCOV    ,IERROR   ,IOWRITE(1) ,IOINV      ,IOWRITE(2),14    
     ;,LMXNDL     ,MAINF    ,NDEVS      ,NROW       ,NUMEL     ,NUMNP      
     ;,NUMTIT     ,NUMTNOD  ,NUMTOBS    ,TIME(NINT) ,AREA      ,BUDAT    
     ;,COVINV     ,DEVNAME  ,EXTNBU     ,INDEXNOD   ,IOBUTYP   ,IOCALBU    
     ;,IODEVICE   ,IOTINT   ,MEASTYP    ,IOUTYP     ,KXX       ,LNNDEL
     ;,LTYPE      ,NBUW     ,NOBUF      ,NOOBSIT    ,TIT       ,TOBS
     ;,VOBS       ,WTOBSBU  ,WTOBSN     ,WTOBST     ,WTOBSU  ,COORD(1,1)
     ;,COORD(1,2) ,COORD(1,3),FILENAME)

       ENDIF

C-------------------------  Builds a priori covariance matrix of param. to be 
C-------------------------  estimated DETERMINISTICALLY

       IF (IOINV.GT.0 .AND. NPARDET.GT.0) THEN

           IF (NBLCVP.GT.NPARDET) THEN
             WRITE(MAINF,2000) NPARDET,NBLCVP,NPARDET
 2000        FORMAT(/,' ERROR: NUMBER OF DETERMINISTICALLY ESTIMATED'
     ;                ' PARAMETERS IS: ',I5,'. YOU DEFINED A PRIORI'
     ;                ' COVARIANCE MATRIX OF PARAM. WITH: ',I5,' BLOCKS'
     ;             ,/,' THAT IS NOT COHERENT. MAXIMUM NUMBER OF BLOCKS'
     ;                ' CAN BE: ',I5)
             STOP ' CRITICAL STOP. CHECK RES.OUT'
           
           ELSE

             IDMBLDIM=NPARDET                           ! Diagonal matrix
             IF (NBLCVP.GT.0) IDMBLDIM=NBLCVP           ! Non diagonal matrix
             CALL COVPAR_MATRIX
     ;(IDMBLDIM   ,IERROR   ,IOWRITE(1) ,IOWRITE(2) ,17
     ;,MAINF      ,MXGRPZN  ,NBLCVP     ,NPARDET    ,NZPAR
     ;,COVPAR     ,FILENAME ,IDMBLCVP   ,IVPAR(1,3) ,IOPT_GS
     ;,IVPAR(1,1) ,STPAR)

           END IF ! NBLCVP.GT.NPARDET
       END IF ! IOINV.GT.0 .AND. NPARDET.GT.0

C_______________ Reads all input data related to geoestatistical inverse problem

       IF (IOINV_GS.GT.0) CALL ENTDAT_GROUPS_ZONES
     ;(IDIMIVARIO_GS ,IDIMVAR_GS  ,IDIMWGT      ,IERROR     ,IOWRITE(1)
     ;,IODIM         ,IOINV       ,IOWRITE(2)   ,ISOT       ,16
     ;,LMXNDL        ,MAINF       ,MXDISC_GS    ,MXGRPZN    ,MXMEASPP_GS
     ;,MXNPP_GS      ,MXNVAR_GS   ,MXNZON_GS    ,NGROUP_ZN  ,NPAREL
     ;,NPBMX         ,NTYPAR      ,NUMEL        ,NUMNP      ,NZPAR
     ;,AREA          ,COORD       ,COORDGR_GS   ,DIVZN_GS   ,EXDRZN_GS
     ;,FILENAME      ,INORPAR     ,IO_KG_GS     ,IOPT_GS    ,IPNT_PAR
     ;,IPOLDRIFT_GS  ,IVARIO_GS   ,IVPAR        ,KXX        ,LNNDEL
     ;,LTYPE         ,LXPAREL     ,NZONE_PAR    ,PARZ       ,POSDIS_GS
     ;,POSMEAS_GS    ,POSZN_GS    ,SEARCH_GS    ,SUPBL_GS   ,TRIM_GS
     ;,VARIO_GS      ,VMEAS_GS    ,VSTATS_GS    ,WGT_PAR    ,PARC_GS)

C------------------------- Ckecks number of errors (IERROR)
C------------------------- If IERROR >0 program stops.

       IF (IERROR.NE.0) THEN

          WRITE(MAINF,3000) IERROR

 3000     FORMAT(//5X,
     ;          'ERRORS FOUND ON INPUT DATA',/,
     ;          'NUMBER OF ERRORS = ......',I5)
       
          STOP  'ERRORS ON INPUT DATA'
       ENDIF
       
C------------------------- NPARF is set to one (when zero) because it is used
C------------------------- to dimension some arrays

       IF (NPARF.EQ.0) NPARF=1

C------------------------- NPARF is set to NPAR (when density is variable)
C------------------------- because in that case, flow and tpt. are coupled.

       IF (IODENS_INI.EQ.1) NPARF = NPAR

C------------------------- Output message to know that the program is leaving
C------------------------- this routine

       IF (IFLAGS(3).EQ.1) CALL IO_SUB ('ENTDAT',1)

       RETURN 
       END
