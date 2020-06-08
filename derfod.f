      SUBROUTINE DERFOD
     &          (ACTH     ,AREA     ,BETAC    ,CAUX1    ,CCALAN
     &          ,CCALIT   ,CFPAREL  ,CREF     ,DENSREF  ,DERC
     &          ,DERCS    ,DTIM     ,EPSTRA   ,FNT      ,ICMPDERS
     &          ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,INDSSTR  ,INEW
     &          ,INORPAR  ,INTI     ,IOCAP    ,IODENS   ,IOFMLT
     &          ,IOTRLI   ,IOVRWC   ,IVPAR    ,KXX      ,LDIM
     &          ,LMXNDL   ,LNNDEL   ,LXPAREL  ,NFLAGS   ,NFNL
     &          ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT
     &          ,NPAR     ,NPAREL   ,NPPEL    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &          ,PAREL    ,WATVOL   ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR
     &          ,IPOS     ,DERIV)

******************************************************************************
*
* PURPOSE
*
*         Computes the RHS of inverse problem for retardation
*
* DESCRIPTION
*
*         Computes the derivatives w.r.t. retardation
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
*     AMS      3-2002     First coding (starting from TRANSIN-II)
*
*****************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::ICMPDERS  ,IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW
     &          ,INTI      ,IOCAP    ,IODENS   ,IOFMLT   ,IOTRLI
     &          ,IOVRWC    ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT
     &          ,NPAR     ,NPAREL   ,NPPEL    ,NTYPAR    ,NUMEL
     &          ,NUMNP    ,NZPAR    ,IDIMWGT  ,IPNT_PAR

      REAL*8::BETAC    ,CREF     ,DENSREF  ,DTIM     ,EPSTRA         

      REAL*8::DENS

      INTEGER*4::IFLAGS(NFLAGS)       ,INORPAR(NTYPAR),IVPAR(NZPAR) 
     &          ,KXX(LMXNDL,NUMEL)    ,LDIM(NUMEL)    ,LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)        ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)


      REAL*8::ACTH(NUMEL)        ,AREA(NUMEL)    ,WGT_PAR(IDIMWGT)
     &       ,CAUX1(NUMNP)               ,CCALAN(NUMNP)
     &       ,CCALIT(NUMNP)              ,CFPAREL(NUMEL,NPAREL)
     &       ,DERC(NUMNP,NPAR,IDIMDERC)  ,DERCS(NUMNP,NPAR,2)
     &       ,FNT(IDIMFNT,NINT)          ,PARACD(3,NFNL)
     &       ,PARC(NZPAR)                ,PAREL(NUMEL,NPPEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)


C------------------------- Internal

      INTEGER*4::I    ,IC   ,IP    ,JJ   ,K1   ,L    ,LD   ,NCNF ,NNUD
     &          ,NPTOT,NZON

      REAL*8::AREALN   ,C1       ,DENSNODE ,FOD      ,RETARD   ,WTV_I

      INTEGER*4 IPOS(NPAR) ,INDEX(12)

      REAL*8 CFPARAM(12) ,DERIV(NPAR) ,XPARAM(8)


C------------------------- First executable statement


      DO L=1,NUMEL

          NZON = LXPAREL(L,8)
          JJ = INORPAR(16) + NZON
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

              CFPARAM(1) = CFPAREL(L,8)

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

                  RETARD = PAREL(L,14)*ACTH(L)
                  
                  AREALN = AREA(L)/NNUD

                  DO K1=1,NNUD

					I = KXX(K1,L)
                      C1 = CAUX1(I)

					IF (IODENS.EQ.1) THEN

						DENSNODE = DENS(DENSREF,BETAC,C1,CREF)

					ELSE

						DENSNODE = 1D0

					END IF !IODENS.EQ.1

                      IF (IOVRWC.LT.2) THEN

                          WTV_I = WATVOL(1,L,2)

	                ELSE

                          WTV_I = WATVOL(K1,L,2)

	                END IF !IOVRWC.LT.2

	                FOD = DENSNODE*AREALN*(RETARD + WTV_I)*C1

                      DO IC=1,NPTOT

                          DERC(I,IPOS(IC),INEW)=DERC(I,IPOS(IC),INEW)-
     &                          DERIV(IC)*FOD

C------------------------- In a radioactive chain, part of the derivatives
C------------------------- of the RHS (simulation) of the "son" w.r.t the
C------------------------- parameters of the "father" is added to the RHS
C------------------------- (inversion) of the "son"

                          IF (ICMPDERS.EQ.1) THEN

                              DERCS(I,IPOS(IC),INEW)=
     &                                DERCS(I,IPOS(IC),INEW)
     &                                + DERIV(IC)*FOD

                          END IF !ICMPDERS.EQ.1

                      END DO !IC=1,NPTOT

                  END DO !K1=1,NNUD

              END IF !NPTOT.GT.0

          END IF !IP.NE.0 .OR. IOTRLI.NE.0

      END DO !L=1,NUMEL

      END SUBROUTINE DERFOD
