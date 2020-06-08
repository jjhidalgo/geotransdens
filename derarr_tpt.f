      SUBROUTINE DERARR_TPT
     &          (DTIM     ,EPSFLU   ,IDIMDERC ,IDIMFNT  ,INARR
     &          ,INDSSTR  ,INEW     ,INTI     ,IOCAP    ,IOFLLI
     &          ,IOFMLF   ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT
     &          ,NPAR     ,NPAREL   ,NPPEL    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZPAR    ,AREA     ,BETAC    ,CAUX1
     &          ,CFPAREL  ,DERC     ,FNT      ,HCALIT   ,HCALAN
     &          ,IFLAGS   ,INORPAR  ,IVPAR    ,KXX      ,LDIM
     &          ,LNNDEL   ,LXPAREL  ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NZONE_PAR,PARACD   ,PARC     ,PAREL
     &          ,IODENS   ,DENSREF  ,CREF     ,IDIMWGT  ,WGT_PAR
     &          ,IPNT_PAR ,IPOS     ,DERIV)

******************************************************************************
*
* PURPOSE
*
*      Computes the RHS of inverse problem for areal recharge (explicit dep.)
*
* DESCRIPTION
*
*      Computes the RHS of inverse problem for areal recharge (explicit dep.)
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,
*                         volume for 3-D)
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
*  FNT                    Array containing time functions values                
*  HCAL                   Computed heads at every node                          
*  HCALAN                 Head level at previous time                           
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

      INTEGER*4::IDIMDERC ,IDIMFNT  ,INARR    ,INDSSTR  ,INEW
     &          ,INTI     ,IOCAP    ,IOFLLI   ,IOFMLF   ,LMXNDL
     &          ,NFLAGS   ,NFNL     ,NINT     ,NPAR     ,NPAREL
     &          ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZPAR
     &          ,IODENS   ,IDIMWGT  ,IPNT_PAR

      REAL*8::BETAC,DTIM,EPSFLU,CREF,DENSREF

      REAL*8::DENS

C------------------------- External

      INTEGER*4::IFLAGS(NFLAGS),INORPAR(NTYPAR) ,IVPAR(NZPAR)
     &          ,KXX(LMXNDL,NUMEL) ,LDIM(NUMEL),LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL) ,NFTPAR(NZPAR),NZONE_PAR(NTYPAR)

      REAL*8::AREA(NUMEL) ,PAREL(NUMEL,NPPEL),CAUX1(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL) ,DERC(NUMNP,NPAR,IDIMDERC)
     &       ,FNT(IDIMFNT,NINT),PARACD(3,NFNL),PARC(NZPAR)
     &       ,HCALIT(NUMNP),HCALAN(NUMNP)    ,WGT_PAR(IDIMWGT)

C------------------------- Internal


      INTEGER*4::K1,I1,IC,L,NZON,JJ,IP,NNUD,LD,NCNF,NPTOT
      REAL*8::AREALN,CEXT,CINT,DENSEXT,FACTOR,RECHARGE,CEXT_CINT

      INTEGER*4::IPOS(NPAR),INDEX(12)
      REAL*8::XPARAM(8),CFPARAM(12),DERIV(NPAR)

C------------------------- First executable statement

      DO L=1,NUMEL

          NZON = LXPAREL(L,INARR)
          JJ = INORPAR(8) + NZON
          IP = IVPAR(JJ)
          RECHARGE = PAREL(L,8)

          IF (IP.NE.0 .OR. IOFLLI.NE.0) THEN

              NNUD = LNNDEL(L)
              LD = LDIM(L)

              IF (IOFLLI.NE.0) THEN
                  NCNF = NFNLPAR(JJ)
              ELSE
                  NCNF = 0
              END IF !IOFLLI.NE.0

              INDEX(1) = JJ
              CFPARAM(1) = CFPAREL(L,INARR)

              CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,L
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,1        ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV  ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

              IF (NPTOT.GT.0) THEN

                  IF (RECHARGE.GT.0) THEN

                      CEXT = PAREL(L,15)

                  ELSE

                      CEXT = 0D0

                  END IF !RECHARGE.GT.0   

                  IF (IODENS.EQ.1) THEN

                      DENSEXT = DENS(DENSREF,BETAC,CEXT,CREF)

                  END IF !IODENS.EQ.1 .AND. PAREL(L,8).GT.0

                  AREALN = AREA(L)/NNUD

                  DO K1=1,NNUD

                      I1 = KXX(K1,L)

                      CINT = CAUX1(I1)

                      IF (IODENS.EQ.1) THEN

                          CEXT_CINT = DENSEXT*(CEXT - CINT)

                      ELSE

                          CEXT_CINT = CEXT - CINT

                      END IF !IODENS.EQ.1


                        FACTOR = AREALN*CEXT_CINT

                      DO IC=1,NPTOT

                          DERC(I1,IPOS(IC),INEW) =DERC(I1,IPOS(IC),INEW)
     &                                           +DERIV(IC)*FACTOR

                      END DO !IC=1,NPTOT

                  END DO !K1=1,NNUD

              END IF !NPTOT.GT.0

          END IF ! IP.NE.0 .OR. IOFLLI.NE.0

      END DO !L=1,NUMEL

      END SUBROUTINE DERARR_TPT
