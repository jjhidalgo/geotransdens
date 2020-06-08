      SUBROUTINE DERCRD
     &          (ACTH     ,AREA     ,BETAC    ,CAUX1    ,CAUX2
     &          ,CCALAN   ,CCALIT   ,CFPAREL  ,CREF     ,DENSREF
     &          ,DERC     ,DERCS    ,DTIM     ,EPSTRA   ,FNT
     &          ,ICMPDERS ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,INDSSTR
     &          ,INEW     ,INORPAR  ,INTI     ,IOCAP    ,IODENS
     &          ,IOFMLT   ,IOTRLI   ,IOVRWC   ,ITPTVAR  ,IVPAR
     &          ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LXPAREL
     &          ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPPEL
     &          ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &          ,PARACD   ,PARC     ,PAREL    ,THETAT   ,WSPECHEAT
     &          ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)

******************************************************************************
*     
*  PURPOSE
*     
*   Computes the RHS of inverse problem for retardation
*     
*  DESCRIPTION
*     
*     Computes the derivatives w.r.t. retardation
*     
*     EXTERNAL VARIABLES: ARRAYS
*     
*  AREA                   Element size (length for 1-D elem, area for 2-D,
*                         volume for 3-D)
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  CAUX2                  Array containing diference of concentrations in two
*                         consecutives times, related to time increment
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
*                      PARACD                 Agreement parameters                                  
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*  PAREL                  Parameter values at every element and current time
*                         for all nodal parameters (each value is computed as
*                         the product of up to four terms:
*                         elem. coeff*zonal value*time funct.*nonl. funct. )
*     
*  INTERNAL VARIABLES: ARRAYS
*     
*  CFPARAM                Used to store element coefficient
*  DERIV                  Stores the derivative of parameter w.r.t. zonal and/
*                         or generic coefficient
*  INDEX                  Location in array *PAR*
*  IPOS                   Location in array DERC of derivative
*  XPARAM                 Auxiliar array used in nonlinear functions
*     
*  EXTERNAL VARIABLES: SCALARS
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
*  NPBTP                  Number of simultaneous transport problems             
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*     
*  INTERNAL VARIABLES: SCALARS
*     
*  L                      Element number
*  LD                     Element dimension
*  NNUD                   Number of nodes of the current element                
*  NPTOT                  Total number of parameters to be estimated (usually 1,
*                         but it can be larger if generic parameters are to be 
*                         estimated)
*     
*  FUNCTIONS AND SUBROUTINES REFERENCED
*     
*  DER_PARAM              Computes the derivatives of a given parameter w.r.t.
*                         its zonal parameter and its generic parameters
*     
*  HISTORY
*     
*  AMS      3-2002     First coding (starting from TRANSIN-II)
*     
*****************************************************************************


      IMPLICIT NONE


C-------------------------External

      INTEGER*4::ICMPDERS ,IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW
     &          ,INTI     ,IOCAP    ,IODENS   ,IOFMLT   ,IOTRLI
     &          ,IOVRWC   ,ITPTVAR  ,LMXNDL   ,NFLAGS   ,NFNL
     &          ,NINT     ,NPAR     ,NPAREL   ,NPPEL    ,NTYPAR
     &          ,NUMEL    ,NUMNP    ,NZPAR    ,IDIMWGT  ,IPNT_PAR


      REAL*8::BETAC    ,CREF     ,DENSREF  ,DTIM     ,EPSTRA
     &       ,THETAT   ,WSPECHEAT

      REAL*8 DENS    ,WGT_PAR(IDIMWGT)

      INTEGER*4::IFLAGS(NFLAGS)       ,INORPAR(NTYPAR),IVPAR(NZPAR)
     &          ,KXX(LMXNDL,NUMEL)    ,LDIM(NUMEL)    ,LNNDEL(NUMEL)        
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)        ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)

      REAL*8::ACTH(NUMEL)                 ,AREA(NUMEL)
     &       ,CAUX1(NUMNP)                ,CAUX2(NUMNP)
     &       ,CCALAN(NUMNP)               ,CCALIT(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)       ,DERC(NUMNP,NPAR,IDIMDERC)
     &       ,DERCS(NUMNP,NPAR,2)         ,FNT(IDIMFNT,NINT)
     &       ,PARC(NZPAR)                 ,PAREL(NUMEL,NPPEL)
     &       ,PARACD(3,NFNL)

      
C-------------------------Internal

      INTEGER*4::I    ,IC   ,INODE,IP   ,JJ   ,K    ,KNODE,L    ,LD
     &          ,NCNF,NNUD ,NPTOT,NZON

      REAL*8::AREALN   ,CTH1     ,CAVG     ,CNODE    ,DENSCRD  ,DENSFOD
     &       ,DENSL    ,FACTOR   ,FOD      ,TH1      ,VAUX

      INTEGER*4::INDEX(12)  ,IPOS(NPAR)

      REAL*8::CFPARAM(12)   ,DERIV(NPAR)   ,XPARAM(8)

C-------------------------First executable statement

      FOD = 0D0
      FACTOR = 0D0

      TH1 = 1D0 - THETAT
      DENSCRD = 1D0
      DENSFOD = 1D0
      DENSL = 1D0

      DO L=1,NUMEL

          NZON = LXPAREL(L,9)
          JJ = INORPAR(17) + NZON
          IP = IVPAR(JJ)

          IF (IP.NE.0 .OR. IOTRLI.NE.0) THEN

              NNUD = LNNDEL(L)
              LD = LDIM(L)

              IF (IOTRLI.NE.0) THEN

                  NCNF = NFNLPAR(JJ)

              ELSE

                  NCNF = 0

              END IF

              INDEX(1) = JJ
              CFPARAM(1) = CFPAREL(L,9)

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

C------------------------- Computes density at k+1-thetat
C------------------------- elementwise, if needed.

                  IF (IODENS.GT.0 .AND. IOVRWC.LE.1) THEN

                      CAVG = 0D0

                      DO I=1,NNUD

                          INODE = KXX(I,L)
                          CAVG = CAVG + TH1*CCALIT(INODE)
     &                                + THETAT*CCALAN(INODE)

                      END DO !I=1,NNUD

                      CAVG = CAVG/NNUD
                      DENSL = DENS(DENSREF,BETAC,CAVG,CREF)

                  END IF !IODENS.GT.0 .AND. IOVRWC.LE.1


                  AREALN = ACTH(L)*AREA(L)/NNUD


                  IF (ITPTVAR.EQ.0) THEN

                      FOD = PAREL(L,13)
                      FACTOR = 1D0

                  ELSE

                      FOD = 0D0
                      FACTOR = 1D0/WSPECHEAT

                  END IF !ITPTVAR.EQ.0


                  DO K=1,NNUD

                      KNODE = KXX(K,L)

                      IF (IODENS.EQ.1) THEN

                          CNODE = CAUX1(KNODE)
                          DENSFOD =  DENS(DENSREF,BETAC,CNODE,CREF)

                          IF (IOVRWC.LE.1) THEN

                              DENSCRD = DENSL

                          ELSE IF (IOVRWC.EQ.2) THEN

                              CTH1 = TH1*CCALIT(KNODE)
     &                             + THETAT*CCALAN(KNODE)
                              DENSCRD = DENS(DENSREF,BETAC,CTH1,CREF)

                          END IF !IOVRWC

                      ELSE

                          DENSFOD = 1D0
                          DENSCRD = 1D0

                      END IF !IODENS.EQ.1

                      VAUX = (DENSFOD*FOD*CAUX1(KNODE)
     &                       + DENSCRD*FACTOR*CAUX2(KNODE))*AREALN

                      DO IC=1,NPTOT

                          DERC(KNODE,IPOS(IC),INEW) =
     &                        DERC(KNODE,IPOS(IC),INEW) - DERIV(IC)*VAUX

C------------------------- In a radioactive chain, part of the derivatives
C------------------------- of the RHS (simulation) of the "son" w.r.t the    
C------------------------- parameters of the "father" is added to the RHS
C------------------------- (inversion) of the "son"

                          IF (ICMPDERS.EQ.1) THEN

                              DERCS(KNODE,IPOS(IC),INEW) =
     &                        DERCS(KNODE,IPOS(IC),INEW) + DERIV(IC)*
     &                        FOD*CAUX1(KNODE)

                          END IF !ICMPDERS.EQ.1

                      END DO !IC=1,NPTOT

                  END DO !K=1,NNUD

              END IF !NPTOT.GT.0

          END IF !IP.NE.0 .OR. IOTRLI.NE.0

      END DO !L=1,NUMEL

      END SUBROUTINE DERCRD
