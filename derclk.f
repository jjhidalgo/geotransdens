      SUBROUTINE DERCLK
     &        (CAUX1    ,CCALAN   ,CCALIT   ,CFPARNP  ,DERC
     &        ,DTIM     ,EPSTRA   ,FNT      ,IDIMDERC ,IDIMFNT
     &        ,IFLAGS   ,INDSSTR  ,INEW     ,INORPAR  ,INTI
     &        ,IOCAP    ,IOFMLT   ,IOTRLI   ,IXPARNP  ,IVPAR
     &        ,KXX      ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLPAR
     &        ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR
     &        ,NPARNP   ,NPPNP    ,NTYPAR   ,NUMEL    ,NUMNP
     &        ,NZONE_PAR,NZPAR    ,PARACD   ,PARC     ,PARNP
     ;        ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV
     &        ,ITPTVAR  ,WSPECHEAT)

******************************************************************************
*
* PURPOSE
*
*         Computes the RHS of inverse problem for concentration leak
*
* DESCRIPTION
*
*         Computes the derivatives w.r.t. concentration leak
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,
*                         volume for 3-D)
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  CCALIT                 Computed concentration at every node                  
*  CCALAN                 Computed concentrations in the previous time step.    
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
*  FNT                    Array containing time functions values                
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LNNDEL                 Number of nodes at every element                      
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
*  WATVOL                 Array containing the water content of every element   
*
* INTERNAL VARIABLES: ARRAYS
*
*  CFPARAM                Used to store element coefficient
*  DERIV                  Stores the derivative of parameter w.r.t. zonal and/
*                         or generic coefficient
*  INDEX                  Location in array *PAR*
*  IPOS                   Location in array DERC of derivative
*  XPARAM                 Auxiliar array used in nonlinear functions
*
* EXTERNAL VARIABLES: SCALARS
*
*  DTIM                   Time used in time function
*  EPSTRA                 Time weighting parameter for nonlinear transport      
*                         problems                                              
*  IDIMFNT                First dimension of array FNT, it coincides with       
*  INDSSTR                Current problem state. If 1, transient state,         
*                         if 0, steady-state                                    
*  INEW                   Location in array DERC
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOFMLT                 Transport formulation number                          
*  IOTRLI                 Idem to IOFLLI, in the case of transport.             
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFLAGS                 Maximum number of allowed flags                       
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NINT                   Number of observation times                           
*  NPAR                   Total number of parameters to be estimated            
*  NPAREL                 Number of element parameters in current problem       
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*
* INTERNAL VARIABLES: SCALARS
*
*  L                      Element number
*  LD                     Element dimension
*  NNUD                   Number of nodes of the current element                
*  NPTOT                  Total number of parameters to be estimated (usually 1,
*                         but it can be larger if generic parameters are to be 
*                         estimated)
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  DER_PARAM              Computes the derivatives of a given parameter w.r.t.
*                         its zonal parameter and its generic parameters
*
* HISTORY
*
*
*
*****************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW     ,INTI
     &          ,IOCAP    ,IOFMLT   ,IOTRLI   ,LMXNDL   ,NFLAGS
     &          ,NFNL     ,NINT     ,NPAR     ,NPARNP   ,NPPNP
     &          ,NTYPAR   ,NUMEL    ,NUMNP    ,NZPAR

      INTEGER*4::IDIMWGT  ,IPNT_PAR ,ITPTVAR

      REAL*8::DTIM   ,EPSTRA ,WSPECHEAT



      INTEGER*4::IFLAGS(NFLAGS)   ,INORPAR(NTYPAR),IXPARNP(NUMNP,NPARNP)
     &          ,IVPAR(NZPAR)     ,KXX(LMXNDL,NUMEL) ,NFNLPAR(NZPAR)
     &          ,NFNLPRG(8,NFNL)  ,NFNLTIP(NFNL)     ,NFTPAR(NZPAR)
     &          ,NZONE_PAR(NTYPAR)


      REAL*8::CAUX1(NUMNP)               ,CCALAN(NUMNP)
     &       ,CCALIT(NUMNP)              ,CFPARNP(NUMNP,NPARNP)
     &       ,DERC(NUMNP,NPAR,IDIMDERC)  ,FNT(IDIMFNT,NINT)
     &       ,PARACD(3,NFNL)             ,PARC(NZPAR)
     &       ,PARNP(NUMNP,NPPNP)         ,WGT_PAR(IDIMWGT)

C------------------------- Internal

      INTEGER*4::I    ,IC   ,IP   ,JJ   ,NCNF ,NPTOT,NZCLK

      REAL*8::CEXT_CNOD,CLK

      INTEGER*4::IPOS(NPAR) ,INDEX(12)

      REAL*8::CFPARAM(12) ,DERIV(NPAR) ,XPARAM(8)


C------------------------- First executable statement


      DO I=1,NUMNP

          NZCLK = IXPARNP(I,10)
          JJ = INORPAR(21) + NZCLK
          IP = IVPAR(JJ)

          IF (IP.NE.0 .OR. IOTRLI.NE.0) THEN


              IF (IOTRLI.NE.0) THEN

                  NCNF = NFNLPAR(JJ)

              ELSE

                  NCNF = 0

              END IF !IOTRLI.NE.0

              INDEX(1) = JJ

              CFPARAM(1) = CFPARNP(I,10)

              CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLT   ,IP       ,I
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,1        ,NPTOT    ,1        ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),CCALIT   ,CCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

              IF (NPTOT.GT.0) THEN

                  CEXT_CNOD = PARNP(I,4) - CAUX1(I)

                  IF (ITPTVAR.EQ.1) DERIV(:) = DERIV(:)/WSPECHEAT

                  DO IC=1,NPTOT

                      CLK = DERIV(IC)*CEXT_CNOD

                      DERC(I,IPOS(IC),INEW) = DERC(I,IPOS(IC),INEW)
     &                          + CLK

                  END DO !IC=1,NPTOT

              END IF !NPTOT.GT.0

          END IF !IP.NE.0 .OR. IOTRLI.NE.0

      END DO !I=1,NUMNP

      END SUBROUTINE DERCLK
