      SUBROUTINE DERDFM
     &        (BIBI     ,CAUX1    ,CCALAN   ,CCALIT   ,CFPAREL
     &        ,DENSITY  ,DERC     ,DTIM     ,EPSTRA   ,FNT
     &        ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,INDSSTR
     &        ,INEW     ,INORPAR  ,INTI     ,IOCAP    ,IODENS
     &        ,IOFMLT   ,IOTRLI   ,IOVRWC   ,ITPTVAR  ,IVPAR
     &        ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LXPAREL
     &        ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &        ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPPEL
     &        ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &        ,PARACD   ,PARC     ,PAREL    ,WATVOL   ,WSPECHEAT
     &        ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR)

*******************************************************************************
*
* PURPOSE
*
*   Manages the computation of the RHS of inverse problem for molecular diff.
*
* DESCRIPTION
*
*   Manages the computation of the RHS of inverse problem for molecular diff.
*   by calling the corresponding rutines by NODES or by ELEMENTS as needed 
*   depending on the way WATVOL is computed
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  DERDFM_I           Derivatives with WATVOL by nodes
*  DERDFM_L           Derivatives with WATVOL by elements
*
* HISTORY
*
*     AMS      6-2002     First coding
*
*******************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMBB   ,IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW
     &          ,INTI     ,IOCAP    ,IODENS   ,IOFMLT   ,IOTRLI
     &          ,IOVRWC,ITPTVAR  ,LMXNDL   ,NFLAGS   ,NFNL   ,NINT
     &          ,NPAR     ,NPAREL   ,NPPEL    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZPAR    ,IDIMWGT  ,IPNT_PAR

      REAL*8::DTIM     ,EPSTRA   ,WSPECHEAT    ,WGT_PAR(IDIMWGT)

      

      INTEGER*4::IFLAGS(NFLAGS)       ,INORPAR(NTYPAR),IVPAR(NZPAR)
     &          ,KXX(LMXNDL,NUMEL)    ,LDIM(NUMEL)    ,LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)        ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)

      REAL*8::BIBI(IDIMBB,NUMEL)        ,CAUX1(NUMNP)
     &       ,CCALIT(NUMNP)             ,CCALAN(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)     ,DENSITY(NUMEL)
     &       ,DERC(NUMNP,NPAR,IDIMDERC) ,FNT(IDIMFNT,NINT)
     &       ,PARACD(3,NFNL)            ,PARC(NZPAR)
     &       ,PAREL(NUMEL,NPPEL)        ,WATVOL(LMXNDL,NUMEL,3)


C------------------------- First executable statement

      IF (IOVRWC.LE.1) THEN

          CALL DERDFM_L
     &        (BIBI     ,CAUX1    ,CCALAN   ,CCALIT   ,CFPAREL
     &        ,DENSITY  ,DERC     ,DTIM     ,EPSTRA   ,FNT
     &        ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IDIMWGT  ,IFLAGS
     &        ,INDSSTR  ,INEW     ,INORPAR  ,INTI     ,IOCAP
     &        ,IODENS   ,IOFMLT   ,IOTRLI   ,IPNT_PAR ,ITPTVAR
     &        ,IVPAR    ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL
     &        ,LXPAREL  ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &        ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR     ,NPAREL
     &        ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &        ,NZPAR    ,PARACD   ,PARC     ,PAREL    ,WATVOL
     &        ,WGT_PAR  ,WSPECHEAT)



      ELSE

          CALL DERDFM_I
     &        (BIBI     ,CAUX1    ,CCALAN   ,CCALIT   ,CFPAREL
     &        ,DENSITY  ,DERC     ,DTIM     ,EPSTRA   ,FNT
     &        ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IDIMWGT  ,IFLAGS
     &        ,INDSSTR  ,INEW     ,INORPAR  ,INTI     ,IOCAP
     &        ,IODENS   ,IOFMLT   ,IOTRLI   ,IPNT_PAR ,ITPTVAR
     &        ,IVPAR    ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL
     &        ,LXPAREL  ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &        ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR     ,NPAREL
     &        ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &        ,NZPAR    ,PARACD   ,PARC     ,PAREL    ,WATVOL
     &        ,WGT_PAR  ,WSPECHEAT)

      END IF !IOVRWC.LE.1

      END SUBROUTINE DERDFM

*****************************************************************************
*****************************************************************************

      SUBROUTINE DERDFM_L
     &          (BIBI     ,CAUX1    ,CCALAN   ,CCALIT   ,CFPAREL
     &          ,DENSITY  ,DERC     ,DTIM     ,EPSTRA   ,FNT
     &          ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IDIMWGT  ,IFLAGS
     &          ,INDSSTR  ,INEW     ,INORPAR  ,INTI     ,IOCAP
     &          ,IODENS   ,IOFMLT   ,IOTRLI   ,IPNT_PAR ,ITPTVAR
     &          ,IVPAR    ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL
     &          ,LXPAREL  ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &          ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR     ,NPAREL
     &          ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &          ,NZPAR    ,PARACD   ,PARC     ,PAREL    ,WATVOL
     &          ,WGT_PAR  ,WSPECHEAT)
******************************************************************************
*
* PURPOSE
*
*         Computes the RHS of inverse problem for molecular diffusion
*
* DESCRIPTION
*
*         Computes the derivatives w.r.t. molecular diffusion. WATVOL is 
*         defined by elements
*
* EXTERNAL VARIABLES: ARRAYS
*
*  BIBI                   Array containing the product of interpolation         
*                         functions gradient, for a given element               
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  CCAL                   Computed concentration at every node                  
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
*  IDIMBB                 Used to dimension array BIBI. Is equal to IDIMQ times 
*                         the maximum possible anisotropy of the problem        
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
*  NPBTP                  Number of simultaneous transport problems             
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
*  KNT                    Counter to situate in array BIBI
*  L                      Element number
*  LANI                   Maximum possible anisotropy of a given element. Equal 
*                         to LDIM*(LDIM+1)/2                                    
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
*     AMS      3-2002     First coding (starting from TRANSIN-II)
*
*****************************************************************************


      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMBB   ,IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW
     &          ,INTI     ,IOCAP    ,IODENS   ,IOFMLT   ,IOTRLI
     &          ,ITPTVAR  ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT
     &          ,NPAR     ,NPAREL   ,NPPEL    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZPAR    ,IDIMWGT  ,IPNT_PAR

      REAL*8::DTIM     ,EPSTRA   ,WSPECHEAT

      

      INTEGER*4::IFLAGS(NFLAGS)       ,INORPAR(NTYPAR),IVPAR(NZPAR)
     &          ,KXX(LMXNDL,NUMEL)    ,LDIM(NUMEL)    ,LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)        ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)

      REAL*8::BIBI(IDIMBB,NUMEL)        ,CAUX1(NUMNP)
     &       ,CCALIT(NUMNP)             ,CCALAN(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)     ,DENSITY(NUMEL)
     &       ,DERC(NUMNP,NPAR,IDIMDERC) ,FNT(IDIMFNT,NINT)
     &       ,PARACD(3,NFNL)            ,PARC(NZPAR)
     &       ,PAREL(NUMEL,NPPEL) ,WATVOL(1,NUMEL,3)    ,WGT_PAR(IDIMWGT)


C------------------------- Internal


      INTEGER*4::I    ,I1   ,I2   ,IC   ,IP   ,JJ   ,K1   ,K2   ,KNT  ,L
     &          ,LANI ,LD   ,NCNF ,NNUD ,NPTOT,NZON

      REAL*8::C1       ,C2       ,FACTOR   ,POR1
     &       ,S   ,WTV

      INTEGER*4::IPOS(NPAR),INDEX(12)
      REAL*8::XPARAM(8),CFPARAM(12),DERIV(NPAR)

C------------------------- First executable statement


      DO L=1,NUMEL

          NZON = LXPAREL(L,6)
          JJ = INORPAR(14) + NZON
          IP = IVPAR(JJ)

          IF (IP.NE.0 .OR. IOTRLI.NE.0) THEN

              NNUD = LNNDEL(L)
              LD = LDIM(L)

              IF (IOTRLI.NE.0) THEN

                  NCNF = NFNLPAR(JJ)

             ELSE

                  NCNF = 0

             END IF !IOTRLI.NE.0

             INDEX(1) = JJ
             CFPARAM(1) = CFPAREL(L,6)

             CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLT   ,IP       ,L
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,1        ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),CCALIT   ,CCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

              IF (NPTOT.GT.0) THEN

                  LANI = (LD*(LD+1))/2
                  KNT = 0

                  IF (ITPTVAR.EQ.0) THEN

                      IF (IODENS.EQ.1) THEN

                          WTV = DENSITY(L)*WATVOL(1,L,2)

                      ELSE

                          WTV = WATVOL(1,L,2)

                      END IF !IODENS.EQ.1

                  
                  ELSE

                      POR1 = 1D0 - PAREL(L,12)
                      FACTOR = POR1*WSPECHEAT

                  END IF !ITPTVAR.EQ.0

                  DO K1=1,NNUD-1

                      I1 = KXX(K1,L)
                      C1 = CAUX1(I1)

                      DO K2=K1+1,NNUD

                          I2 = KXX(K2,L)
                          C2 = CAUX1(I2)
                          S = 0D0

                          DO I=1,LD

                              IF (ITPTVAR.EQ.0) THEN

                                  S = S + BIBI(KNT+I,L)*WTV*(C2-C1)

                              ELSE

                                  S = S + BIBI(KNT+I,L)*FACTOR*(C2-C1)

                              END IF !ITPTVAR.EQ.0

                          END DO !I=1,LD


                          DO IC=1,NPTOT

                              DERC(I1,IPOS(IC),INEW)=
     &                            DERC(I1,IPOS(IC),INEW) - S*DERIV(IC)
                              DERC(I2,IPOS(IC),INEW)=
     &                            DERC(I2,IPOS(IC),INEW) + S*DERIV(IC)

                          END DO !IC=1,NPTOT

                          KNT = KNT + LANI

                      END DO !K2=K1+1,NNUD
                  END DO ! K1=1,NNUD-1
              ENDIF ! NPTOT.GT.0
          END IF ! IP.NE.0 .OR. IOTRLI.NE.0
      END DO !L=!,NUMEL

      END SUBROUTINE DERDFM_L

*****************************************************************************
*****************************************************************************

      SUBROUTINE DERDFM_I
     &          (BIBI     ,CAUX1    ,CCALAN   ,CCALIT   ,CFPAREL
     &          ,DENSITY  ,DERC     ,DTIM     ,EPSTRA   ,FNT
     &          ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IDIMWGT  ,IFLAGS
     &          ,INDSSTR  ,INEW     ,INORPAR  ,INTI     ,IOCAP
     &          ,IODENS   ,IOFMLT   ,IOTRLI   ,IPNT_PAR ,ITPTVAR
     &          ,IVPAR    ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL
     &          ,LXPAREL  ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &          ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR     ,NPAREL
     &          ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &          ,NZPAR    ,PARACD   ,PARC     ,PAREL    ,WATVOL
     &          ,WGT_PAR  ,WSPECHEAT)

******************************************************************************
*
* PURPOSE
*
*         Computes the RHS of inverse problem for molecular diffusion
*
* DESCRIPTION
*
*         Computes the derivatives w.r.t. molecular diffusion. WATVOL is 
*         defined by nodes
*
* EXTERNAL VARIABLES: ARRAYS
*
*  BIBI                   Array containing the product of interpolation         
*                         functions gradient, for a given element               
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  CCAL                   Computed concentration at every node                  
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
*  IDIMBB                 Used to dimension array BIBI. Is equal to IDIMQ times 
*                         the maximum possible anisotropy of the problem        
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
*  NPBTP                  Number of simultaneous transport problems             
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
*  KNT                    Counter to situate in array BIBI
*  L                      Element number
*  LANI                   Maximum possible anisotropy of a given element. Equal 
*                         to LDIM*(LDIM+1)/2                                    
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
*     AMS      3-2002     First coding (starting from TRANSIN-II)
*
*****************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMBB   ,IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW
     &          ,INTI     ,IOCAP    ,IODENS   ,IOFMLT   ,IOTRLI
     &          ,ITPTVAR  ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT
     &          ,NPAR     ,NPAREL   ,NPPEL    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZPAR    ,IDIMWGT  ,IPNT_PAR

      REAL*8::DTIM     ,EPSTRA   ,WSPECHEAT    ,WGT_PAR(IDIMWGT)

      

      INTEGER*4::IFLAGS(NFLAGS)       ,INORPAR(NTYPAR),IVPAR(NZPAR)
     &          ,KXX(LMXNDL,NUMEL)    ,LDIM(NUMEL)    ,LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)        ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)

      REAL*8::BIBI(IDIMBB,NUMEL)        ,CAUX1(NUMNP)
     &       ,CCALIT(NUMNP)             ,CCALAN(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)     ,DENSITY(NUMEL)
     &       ,DERC(NUMNP,NPAR,IDIMDERC) ,FNT(IDIMFNT,NINT)
     &       ,PARACD(3,NFNL)            ,PARC(NZPAR)
     &       ,PAREL(NUMEL,NPPEL)        ,WATVOL(LMXNDL,NUMEL,3)


C------------------------- Internal


      INTEGER*4::I    ,I1   ,I2   ,IC   ,IP   ,JJ   ,K1   ,K2   ,KNT  ,L
     &          ,LANI ,LD   ,NCNF ,NNUD ,NPTOT,NZON

      REAL*8::C1       ,C2       ,FACTOR   ,POR1     ,S        ,WTV_AVG

      INTEGER*4::IPOS(NPAR),INDEX(12)
      REAL*8::XPARAM(8),CFPARAM(12),DERIV(NPAR)

C------------------------- First executable statement


      DO L=1,NUMEL

          NZON = LXPAREL(L,6)
          JJ = INORPAR(14)+NZON
          IP = IVPAR(JJ)

          IF (IP.NE.0 .OR. IOTRLI.NE.0) THEN
              NNUD = LNNDEL(L)
              LD = LDIM(L)

              IF (IOTRLI.NE.0) THEN

                  NCNF = NFNLPAR(JJ)

              ELSE

                  NCNF=0

              END IF !IOTRLI.NE.0

              INDEX(1) = JJ
              CFPARAM(1) = CFPAREL(L,6)

              CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLT   ,IP       ,L
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,1        ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),CCALIT   ,CCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

              IF (NPTOT.GT.0) THEN

C------------------------- WTV_AVG is the average WATVOL in the element

                  WTV_AVG = 0D0

                  IF (ITPTVAR.EQ.0) THEN

                      DO K1=1,NNUD
                          WTV_AVG = WTV_AVG + WATVOL(K1,L,2)
                      END DO !K1=1,NNUD

                      WTV_AVG = WTV_AVG/NNUD

                      IF (IODENS.EQ.1) THEN

                          WTV_AVG = WTV_AVG*DENSITY(L)

                      END IF !IODENS.EQ.1

                  ELSE

                      POR1 = 1D0 - PAREL(L,12)
                      FACTOR = POR1*WSPECHEAT

                  END IF !ITPTVAR.EQ.0

                  LANI = (LD*(LD+1))/2
                  KNT = 0

                  DO K1=1,NNUD-1

                      I1 = KXX(K1,L)
                      C1 = CAUX1(I1)
              
                      DO K2=K1+1,NNUD

                          I2 = KXX(K2,L)
                          C2 = CAUX1(I2)
                          S = 0D0

                          DO I=1,LD
                              
                              IF (ITPTVAR.EQ.0) THEN

                                  S = S + BIBI(KNT+I,L)*WTV_AVG*(C2-C1)

                              ELSE

                                  S = S + BIBI(KNT+I,L)*FACTOR*(C2-C1)

                              END IF !ITPTVAR.EQ.1

                          END DO !I=1,LD

                          DO IC=1,NPTOT

                              DERC(I1,IPOS(IC),INEW)=
     &                              DERC(I1,IPOS(IC),INEW) - DERIV(IC)*S
                              DERC(I2,IPOS(IC),INEW)=
     &                              DERC(I2,IPOS(IC),INEW) + DERIV(IC)*S

                          END DO !IC=1,NPTOT

                          KNT = KNT + LANI

                      END DO !K2=K1+1,NNUD

                  END DO ! K1=1,NNUD-1

              END IF ! NPTOT.GT.0

          END IF ! IP.NE.0 .OR. IOTRLI.NE.0

      END DO !L=!,NUMEL

      END SUBROUTINE DERDFM_I
