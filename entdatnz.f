       SUBROUTINE ENTDATNZ
     ; (IDIMWGT   ,IERROR   ,INPWR    ,IODIM    ,IOEQT     ,IOFLLI   
     ; ,IOFLSAT   ,IOINV    ,IOTRLI   ,IORTS    ,IOTRS     ,IOWAR    
     ; ,ISOT      ,IUPAR    ,MAINF    ,MXGRPZN  ,NFNL      ,NPAR     
     ; ,NPARF     ,NPARFPRG ,NROW     ,NTYPAR   ,NUMEL     ,NZPAR    
     ; ,FILENAME  ,INDPAR   ,INORPAR  ,IOLG_PAR ,IOPT_GS   ,IPNT_PAR  
     ; ,ISOZ      ,IVPAR    ,LDIM     ,NFNLPAR  ,NFTPAR    ,NZONE_PAR 
     ; ,PAR_WGT   ,PARC     ,PARM     ,PARZ     ,STPAR     ,WGT_PAR   
     ; ,NGROUP_ZN ,IPARDET  ,IOINV_GS ,WGT_UNK)    

*****************************************************************************
*
* PURPOSE
*      Reads zonal information of all parameters
*
* DESCRIPTION
*      Reads zonal related values of all parameters (transmissivity, storage,
*      recharge, etc): initial values, prior estimates, estimation index, 
*      number of nonlinear and time functions and standard deviation
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  INDPAR                 Array of 0's and 1's. 0 means that the parameter is   
*                         estimated aritmethically, otherwise logarithmically.  
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IOLG_PAR               Array containing all logarithmic options of           
*                         estimated parameters                                  
*  ISOZ                   Anisotropy of every transmissivity zone               
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  LDIM                   Vector containing fisical dimension of j-th element   
*  NFNLPAR                Vector containing non-linear function order           
*                         afecting every parameter at each zone.                
*  NFTPAR                 Vector containing time function number for every      
*                         parameter zone                                        
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*  PARM                   Vector containing measured values for all             
*                         parameters                                            
*  STPAR                  Vector containing standard deviation errors of        
*                         all parameters prioo information                      
*  WGT_UNK                Array containing WEIGHTS (for objective function
*                         calculation) of problems unknowns
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IODIM                  Maximum dimension of any element included             
*                         in the problem                                        
*  IOEQT                  Type of problem to be solved                          
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOFLSAT                Indicates the possibility that one part of the domain 
*                         reaches unsaturated state.                            
*  IOINV                  Inverse problem option                                
*  IOTRLI                 Idem to IOFLLI, in the case of transport.             
*  IORTS                  Transport regime                                           
*  IOTRS                  Flow regime                                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  ISOT                   Maximum hydraulic conductivity tensor anisotropy      
*                         degree in the problem.                                
*  IUPAR                  Unit number of file PAR                               
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NPAR                   Total number of parameters to be estimated            
*  NPARF                  Number of transient parameters to be estimated        
*  NPARFPRG               Number of uncertain generic parameter zones involved  
*                         in the non-linear flow inverse problem                
*  NROW                   Current record number                                 
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  NZTRA                  Number of transmissivity zones                        
*
* INTERNAL VARIABLES: SCALARS
*
*  IOLGTRA                Transmissivity log-scaling index                      
*  IPAR                   Current number of parameters to estimate
*  IPOS                   Location of parameters in arrays PARC, PARM, etc
*  IZTRAC                 Number of transmissivity zonal values
*  NPARFDIM               Number of parameters to estimate as given in file DIM
*  NPARFPAR               Number of parameters to estimate as given in file PAR,
*                         it must coincide with NPARFDIM
*  NPFNL                  Counts the number of nonlinear flow functions read
*  NPTNL                  Counts the number of nonlinear transp. functions read
*  NZALF                  Number of leakance zones                              
*  NZARR                  Number of areal recharge zones                        
*  NZCHP                  Number of prescribed head zones                       
*  NZCOE                  Number of external concentration zones                
*  NZCRD                  Number of retardation Coefficient zones               
*  NZDFM                  Number of molecular difusion zones                    
*  NZDSP                  Number of dispersivity zones                          
*  NZFOD                  Number of zones of first order decay                  
*  NZPOR                  Number of porosity zones                              
*  NZPRG                  Total number of generic parameter zones               
*  NZQQP                  Number of prescribed flow zones                       
*  NZSTG                  Number of storage Coefficient zones                   
*  VARNAME                Stores the name of the parameters
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  READ_COV_MAT           Reads covariance matrix
*  READ_PAR               Reads and checks zonal values of a given parameter
*                         type except transmissivities
*  READ_TRA               Reads and checks zonal values of transmissivity
*
* HISTORY
*
*     AMS        1988     First coding
*     SCR      4-1997     Revision and verification
*     AMS      1-1998     Common elimination and modification
*
*****************************************************************************

       IMPLICIT NONE

       INTEGER*4 IDIMWGT,IERROR,IOWAR,IPARDET,NZPAR,NPAR,NTYPAR,IOEQT
     ;          ,IOINV,ISOT,IODIM,INPWR,IOFLLI,IOTRS,IUPAR,MAINF,MXGRPZN
     ;          ,NFNL,NUMEL,IOTRLI,IORTS,IOFLSAT,NPARF,NPARFPRG
     ;          ,NGROUP_ZN,IOINV_GS
     ;          ,ISOZ,IPNT_PAR(NZPAR*IDIMWGT),NZONE_PAR(NTYPAR)
     ;          ,IOLG_PAR(NTYPAR,2),INORPAR(NTYPAR),INDPAR(NPAR)
     ;          ,IOPT_GS(MXGRPZN,20),IVPAR(NZPAR,4),LDIM(NUMEL)
     ;          ,NFTPAR(NZPAR),NFNLPAR(NZPAR)

       REAL*8 WGT_PAR(NZPAR*IDIMWGT),PARZ(NZPAR),PARM(NPAR),STPAR(NZPAR)
     ;       ,WGT_UNK(NPAR),PAR_WGT(NTYPAR),PARC(NPAR)

       INTEGER*4 NPFNL,IOLGTRA,IPOS,IZTRAC,IOLGVAR,IOPBLI
     ,          ,IOTIM,ISTART,NZVAR,NPTNL,NROW,NZONTOT
     ;          ,NPARFDIM,NPARFPAR,NPPFLOW,IGROUP,NPPTOTAL
       REAL*4 ERNUM
       CHARACTER FILENAME(20)*20,VARNAME*20

C------------------------- FIRST EXECUTABLE STATEMENT.
C------------------------- Initializes number of order of FLOW estimated param.
C------------------------- and finds number of pilot points related to flow and
C------------------------- transport parameters

       NPARFDIM=NPARF-NPARFPRG
       NPPTOTAL=0   ! Total number of pilot points
       NPPFLOW=0    ! Flow pilot points

       IF (IOINV_GS.NE.0) THEN
         NPPTOTAL=0
         NPPFLOW=0
         DO IGROUP=1,NGROUP_ZN
           IF (IOPT_GS(IGROUP,2).EQ.1) THEN   ! Pilot points used
             IF (IOPT_GS(IGROUP,1).LE.4) 
     ;          NPPFLOW=NPPFLOW+IOPT_GS(IGROUP,5)
             NPPTOTAL=NPPTOTAL+IOPT_GS(IGROUP,5)
           END IF ! IOPT_GS(IGROUP,2).EQ.1
         END DO ! IGROUP=1,NGROUP_ZN

C------------------------- Echoes partition errors

         IF (NPARF.LT.NPPFLOW) CALL ERROR
     ;(IERROR,IOWAR,MAINF,FILENAME,
     ;'NUMBER OF FLOW ESTIMATED PARAMETERS (NPARF) IS SMALLER THAN '//
     ;'THE NUMBER OF PILOT POINTS FOR FLOW',NROW,6,IUPAR,2,ERNUM+0.01)
 
         IF (NPAR.LT.NPPTOTAL) CALL ERROR
     ;(IERROR,IOWAR,MAINF,FILENAME,
     ;'NUMBER OF ESTIMATED PARAMETERS (NPAR) IS SMALLER THAN '//
     ;'NUMBER OF PILOT POINTS',NROW,6,IUPAR,2,ERNUM+0.01)
     
       END IF ! IOINV_GS.NE.0

C------------------------- Initialization of counter parameters

       IPARDET=0
       ISTART=0
       NZONTOT=0
       NPFNL=0

C------------------------- Reads transmissivity zones

       IF (IOEQT.NE.2) THEN
          IOLGTRA=IOLG_PAR(1,1)
          IZTRAC=MAX0(ISOT,IODIM)*NZONE_PAR(1)    !Aux for TRAC dimensioning
          CALL READ_TRA
     ;(IDIMWGT    ,IERROR     ,INPWR      ,IOINV      ,IOLGTRA   
     ;,IOFLLI     ,IOWAR      ,IPARDET    ,ISOT       ,IOTRS     
     ;,IUPAR      ,IZTRAC     ,MAINF      ,MXGRPZN    ,NFNL       
     ;,NGROUP_ZN  ,NPAR       ,NPFNL      ,NUMEL      ,NZONE_PAR(1) 
     ;,PAR_WGT(1) ,FILENAME   ,INDPAR     ,INORPAR    ,IVPAR(1,4)
     ;,IOPT_GS    ,IVPAR(1,2) ,IPNT_PAR   ,IVPAR(1,1) ,ISOZ       
     ;,IVPAR(1,3) ,LDIM       ,NFNLPAR    ,NFTPAR     ,STPAR      
     ;,PARC       ,PARM       ,PARZ       ,WGT_PAR    ,WGT_UNK,NROW)

C------------------------- Reads storage coeficient zones
          
          NZONTOT=IZTRAC
          NZVAR=NZONE_PAR(2)
          IF (NZVAR.NE.0) THEN
             ISTART=NZONTOT*IDIMWGT+1
             ERNUM=6.20
             IOLGVAR=IOLG_PAR(2,1)
             IOPBLI=IOFLLI
             IOTIM=IOTRS
             VARNAME='STORAGE COEFFICIENT '
             IPOS=INORPAR(7)+1
             CALL READ_PAR
     ;(ERNUM         ,IDIMWGT      ,IERROR       ,INPWR         ,IOINV         
     ;,IOLGVAR       ,IOPBLI       ,IOTIM        ,IPARDET       ,IOWAR    
     ;,ISTART        ,IUPAR        ,MAINF        ,MXGRPZN       ,NFNL     
     ;,NGROUP_ZN     ,NPAR         ,NPFNL        ,NZVAR         ,VARNAME
     ;,PAR_WGT(2)    ,FILENAME     ,INDPAR       ,IVPAR(IPOS,4) ,IOPT_GS  
     ;,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)           ,IVPAR(IPOS,1) 
     ;,IVPAR(IPOS,3) ,NFNLPAR(IPOS),NFTPAR(IPOS) ,STPAR(IPOS)   ,PARC          
     ;,PARM          ,PARZ(IPOS)   ,WGT_UNK      ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(2).NE.0

C------------------------- Reads recharge coeficient zones

          NZONTOT=NZONTOT+NZONE_PAR(2)
          NZVAR=NZONE_PAR(3)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=6.30      
            IOLGVAR=IOLG_PAR(3,1)
            IOPBLI=IOFLLI
            IOTIM=IOTRS
            IPOS=INORPAR(8)+1
            VARNAME='      RECHARGE      '
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR        ,INPWR    ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM         ,IPARDET  ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF         ,MXGRPZN  ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL         ,NZVAR    ,VARNAME  
     ;,PAR_WGT(3)    ,FILENAME      ,INDPAR        ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)        
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS) ,NFTPAR(IPOS)  
     ;,STPAR(IPOS)   ,PARC          ,PARM          ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(3).NE.0

C------------------------- Reads prescribed head

          NZONTOT=NZONTOT+NZONE_PAR(3)
          NZVAR=NZONE_PAR(4)
          IF (NZVAR.NE.0)  THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=6.40      
            IOLGVAR=IOLG_PAR(4,1)
            IOPBLI=IOFLLI
            IOTIM=IOTRS
            IPOS=INORPAR(9)+1
            VARNAME='   PRESCRIBED HEAD  '
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR        ,INPWR    ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM         ,IPARDET  ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF         ,MXGRPZN  ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL         ,NZVAR    ,VARNAME  
     ;,PAR_WGT(4)    ,FILENAME      ,INDPAR        ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)            
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS) ,NFTPAR(IPOS)  
     ;,STPAR(IPOS)   ,PARC          ,PARM          ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(4).NE.0

C------------------------- Reads prescribed flow

          NZONTOT=NZONTOT+NZONE_PAR(4)
          NZVAR=NZONE_PAR(5)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=6.50      
            IOLGVAR=IOLG_PAR(5,1)
            IOPBLI=IOFLLI
            IOTIM=IOTRS
            IPOS=INORPAR(10)+1
            VARNAME='  PRESCRIBED FLOW   '
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR        ,INPWR    ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM         ,IPARDET  ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF         ,MXGRPZN  ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL         ,NZVAR    ,VARNAME  
     ;,PAR_WGT(5)    ,FILENAME      ,INDPAR        ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)      
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS) ,NFTPAR(IPOS)  
     ;,STPAR(IPOS)   ,PARC          ,PARM          ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(5).NE.0

C------------------------- Reads leakage

          NZONTOT=NZONTOT+NZONE_PAR(5)
          NZVAR=NZONE_PAR(6)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=6.60      
            IOLGVAR=IOLG_PAR(6,1)
            IOPBLI=IOFLLI
            IOTIM=IOTRS
            IPOS=INORPAR(11)+1
            VARNAME='       LEAKAGE      '
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR        ,INPWR    ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM         ,IPARDET  ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF         ,MXGRPZN  ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL         ,NZVAR    ,VARNAME  
     ;,PAR_WGT(6)    ,FILENAME      ,INDPAR        ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)      
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS) ,NFTPAR(IPOS)  
     ;,STPAR(IPOS)   ,PARC          ,PARM          ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(6).NE.0

          IF (IOINV.GT.0) NPARFPAR=IPARDET ! Number of determ. estimated param.

       ENDIF !IOEQT.NE.2

C------------------------- Checks if the input number of flow estimated 
C------------------------- parameters (in file DIM) coincides with the
C------------------------- actual number of parameters to estimate (in file PAR)
C------------------------- plus pilot points devoted to flow problems

       NPARFPAR=NPARFPAR+NPPFLOW
       IF (NPARFPAR.NE.NPARFDIM .AND. IOINV.GT.0) THEN

          WRITE(MAINF,3000) NPARFDIM,NPARFPAR

 3000     FORMAT(//,
     ;     ' INPUT NUMBER OF FLOW PARAM. FOR EST. (dim)'
     ;     ' :',I5,/
     ;     ' REAL NUMBER OF FLOW PARAM. FOR EST. (par + pilot points)'
     ;     ' :',I5)

          STOP ' CRITICAL STOP. CHECK RES.OUT'
       ENDIF

C------------------------- Checks if the input number of flow estimated 
C------------------------- parameters (in file DIM) coincides with the
C------------------------- actual number of parameters to estimate (in file PAR)
C------------------------- 

       IF (IOEQT.EQ.1) THEN

          IF (NPARFPAR+NPARFPRG.NE.NPAR .AND. IOINV.GT.0) THEN
             WRITE(MAINF,3100) NPAR,NPARFPAR+NPARFPRG
 3100        FORMAT(//,
     ;     ' INPUT NUMBER OF PARAMETERS FOR ESTIMATION : ',I5,/,
     ;     ' REAL NUMBER OF PARAMETERS FOR ESTIMATION: ',I5)
             STOP ' CRITICAL STOP. CHECK FILE RES.OUT'
          ENDIF
       ENDIF
       

C------------------------- TRANSPORT PARAMETERS

C------------------------- Initializes the counter of transport nonlinear 
C------------------------- functions

       NPTNL=0

       IF (IOEQT.NE.1) THEN

C------------------------- Reads longitudinal dispersivity

          NZONTOT=NZONTOT+NZONE_PAR(6)
          NZVAR=NZONE_PAR(7)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=6.70      
            IOLGVAR=IOLG_PAR(7,1)
            IOPBLI=IOTRLI
            IOTIM=IORTS
            IPOS=INORPAR(12)+1
            VARNAME=' LONG. DISPERSIVITY '
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR        ,INPWR    ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM         ,IPARDET  ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF         ,MXGRPZN  ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL         ,NZVAR    ,VARNAME  
     ;,PAR_WGT(7)    ,FILENAME      ,INDPAR        ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)            
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS) ,NFTPAR(IPOS)  
     ;,STPAR(IPOS)   ,PARC          ,PARM          ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(7).NE.0

C------------------------- Reads transversal dispersivity

          NZONTOT=NZONTOT+NZONE_PAR(7)
          NZVAR=NZONE_PAR(7)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=6.70      
            IOLGVAR=IOLG_PAR(7,1)
            IOPBLI=IOTRLI
            IOTIM=IORTS
            IPOS=INORPAR(13)+1
            VARNAME='TRANSV. DISPERSIVITY'
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR       ,INPWR    ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM        ,IPARDET  ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF        ,MXGRPZN  ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL        ,NZVAR    ,VARNAME  
     ;,PAR_WGT(8)    ,FILENAME      ,INDPAR       ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)            
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS) ,NFTPAR(IPOS)  
     ;,STPAR(IPOS)   ,PARC          ,PARM          ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(7).NE.0

C------------------------- Reads molecular difusion

          NZONTOT=NZONTOT+NZONE_PAR(8)
          NZVAR=NZONE_PAR(9)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=6.80      
            IOLGVAR=IOLG_PAR(9,1)
            IOPBLI=IOTRLI
            IOTIM=IORTS
            IPOS=INORPAR(14)+1
            VARNAME='  MOLECULAR DIFF.   '
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR        ,INPWR    ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM         ,IPARDET  ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF         ,MXGRPZN  ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL         ,NZVAR    ,VARNAME  
     ;,PAR_WGT(9)    ,FILENAME      ,INDPAR        ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)            
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS) ,NFTPAR(IPOS) 
     ;,STPAR(IPOS)   ,PARC          ,PARM          ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(9).NE.0

       ENDIF !IOEQT.NE.1

C------------------------- Reads porosity

       IF (IOEQT.NE.1 .OR. IOFLSAT.NE.0) THEN

          NZONTOT=NZONTOT+NZONE_PAR(9)
          NZVAR=NZONE_PAR(10)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=6.90      
            IOLGVAR=IOLG_PAR(10,1)
            IOPBLI=IOTRLI
            IOTIM=IORTS
            IPOS=INORPAR(15)+1
            VARNAME='      POROSITY      '
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR        ,INPWR    ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM         ,IPARDET  ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF         ,MXGRPZN  ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL         ,NZVAR    ,VARNAME
     ;,PAR_WGT(10)   ,FILENAME      ,INDPAR        ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)            
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS) ,NFTPAR(IPOS)  
     ;,STPAR(IPOS)   ,PARC          ,PARM          ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          ELSE

            CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;      'TRANSPORT OR UNSATURATED FLOW PROBLEM '//
     ;      'AND POROSITY IS NOT DEFINED'
     ;      ,NROW,0,IUPAR,2,6.70)

          END IF ! NZONE_PAR(10).NE.0

       ENDIF ! IOEQT.NE.1 .OR. IOFLSAT.NE.0
 
       IF (IOEQT.NE.1) THEN

C------------------------- Reads first order decay coeficient

          NZONTOT=NZONTOT+NZONE_PAR(10)
          NZVAR=NZONE_PAR(11)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=7.00      
            IOLGVAR=IOLG_PAR(11,1)
            IOPBLI=IOTRLI
            IOTIM=IORTS
            IPOS=INORPAR(16)+1
            VARNAME='  FIRST ORDER DECAY '
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR        ,INPWR   ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM         ,IPARDET ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF         ,MXGRPZN ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL         ,NZVAR   ,VARNAME  
     ;,PAR_WGT(11)   ,FILENAME      ,INDPAR        ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)            
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS) ,NFTPAR(IPOS)  
     ;,STPAR(IPOS)   ,PARC          ,PARM          ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)
 
          END IF ! NZONE_PAR(11) 

C------------------------- Reads retardation coeficient

          NZONTOT=NZONTOT+NZONE_PAR(11)
          NZVAR=NZONE_PAR(12)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=7.10      
            IOLGVAR=IOLG_PAR(12,1)
            IOPBLI=IOTRLI
            IOTIM=IORTS
            IPOS=INORPAR(17)+1
            VARNAME=' RETARDATION COEFF. '
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR      ,INPWR   ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM       ,IPARDET ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF       ,MXGRPZN ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL       ,NZVAR   ,VARNAME
     ;,PAR_WGT(12)   ,FILENAME      ,INDPAR      ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)             
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS)        ,NFTPAR(IPOS)
     ;,STPAR(IPOS)   ,PARC          ,PARM        ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(12)

C------------------------- Reads external concentration

          NZONTOT=NZONTOT+NZONE_PAR(12)
          NZVAR=NZONE_PAR(13)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=7.20      
            IOLGVAR=IOLG_PAR(13,1)
            IOPBLI=IOTRLI
            IOTIM=IORTS
            IPOS=INORPAR(18)+1
            VARNAME='  EXTERNAL CONCENT. '
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR        ,INPWR   ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM         ,IPARDET ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF         ,MXGRPZN ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL         ,NZVAR   ,VARNAME  
     ;,PAR_WGT(13)   ,FILENAME      ,INDPAR        ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)            
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS) ,NFTPAR(IPOS)  
     ;,STPAR(IPOS)   ,PARC          ,PARM          ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(13).NE.0

C------------------------- Reads concentration leakage coeficient

          NZONTOT=NZONTOT+NZONE_PAR(17)
          NZVAR=NZONE_PAR(18)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=7.30      
            IOLGVAR=IOLG_PAR(16,1)
            IOPBLI=IOTRLI
            IOTIM=IORTS
            IPOS=INORPAR(21)+1
            VARNAME='CONC. LEAKAGE COEF. '
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR      ,INPWR   ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM       ,IPARDET ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF       ,MXGRPZN ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL       ,NZVAR   ,VARNAME
     ;,PAR_WGT(18)   ,FILENAME      ,INDPAR      ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)             
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS)        ,NFTPAR(IPOS)
     ;,STPAR(IPOS)   ,PARC          ,PARM        ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(18).NE.0

       ENDIF    !IOEQT.NE.1

C------------------------- Reads group parameters for non linear problems

       IF (IOFLLI.NE.0 .OR. IOTRLI.NE.0) THEN

          NZONTOT=NZONTOT+NZONE_PAR(13)
          NZVAR=NZONE_PAR(14)
          IF (NZVAR.NE.0) THEN
            ISTART=NZONTOT*IDIMWGT+1
            ERNUM=7.20      
            IOLGVAR=IOLG_PAR(14,1)
            IOPBLI=IOFLLI
            IOTIM=IOTRS
            IPOS=INORPAR(19)+1
            VARNAME='  GENERIC PARAMETERS'
            CALL READ_PAR
     ;(ERNUM         ,IDIMWGT       ,IERROR        ,INPWR   ,IOINV         
     ;,IOLGVAR       ,IOPBLI        ,IOTIM         ,IPARDET ,IOWAR    
     ;,ISTART        ,IUPAR         ,MAINF         ,MXGRPZN ,NFNL     
     ;,NGROUP_ZN     ,NPAR          ,NPFNL         ,NZVAR   ,VARNAME  
     ;,PAR_WGT(14)   ,FILENAME      ,INDPAR        ,IVPAR(IPOS,4)
     ;,IOPT_GS       ,IVPAR(IPOS,2) ,IPNT_PAR(ISTART)            
     ;,IVPAR(IPOS,1) ,IVPAR(IPOS,3) ,NFNLPAR(IPOS) ,NFTPAR(IPOS)  
     ;,STPAR(IPOS)   ,PARC          ,PARM          ,PARZ(IPOS)    
     ;,WGT_UNK       ,WGT_PAR(ISTART) ,NROW)

          END IF ! NZONE_PAR(14).NE.0

       ENDIF ! IOFLLI.NE.0 .OR. IOTRLI.NE.0

C------------------------- Checks the total number of parameters with a 
C------------------------- nonlinear dependence

       IF (IOFLLI.NE.0. AND .NPFNL.EQ.0) THEN

          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;              'NON LINEAR FLOW BUT NON LINEAR FUNCTIONS ARE NOT'//
     ;              ' DEFINED',NROW,5,IUPAR,0,0.00)

       ELSE IF (IOTRLI.NE.0. AND .NPTNL.EQ.0) THEN
          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;              'NON LINEAR TRANSPORT BUT NON LINEAR FUNCTIONS'
     ;          //' ARE NOT DEFINED',0,0,IUPAR,0,0.00)
       ENDIF   

C------------------------- Checks NPAR

       IF (IOINV.GT.0) THEN
         NPARFPAR=IPARDET+NPPTOTAL
         IF (NPARFPAR+NPARFPRG.NE.NPAR .AND. IOINV.GT.0) THEN
           WRITE(MAINF,3100) NPAR,NPARFPAR+NPARFPRG
           STOP ' CRITICAL STOP. CHECK FILE RES.OUT'
          ENDIF
       END IF

       RETURN
       END
 
