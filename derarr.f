      SUBROUTINE DERARR
     &          (AREA     ,BETAC    ,CFPAREL  ,CREF     ,DENSREF
     &          ,DERH     ,DTIM     ,EPSFLU   ,FNT      ,HCALIT
     &          ,HCALAN   ,IDIMDERH ,IDIMFNT  ,IFLAGS   ,INARR
     &          ,INDSSTR  ,INEW     ,INORPAR  ,INTI     ,IODENS
     &          ,IOFLLI   ,IOFMLF   ,IVPAR    ,KXX      ,LMXNDL
     &          ,LNNDEL   ,LXPAREL  ,NFLAGS   ,NFNL     ,NFNLPAR
     &          ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAREL
     &          ,NPARF    ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP
     &          ,NZONE_PAR,NZPAR    ,PARACD   ,PARC     ,PAREL
     &          ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,NPAR     ,IPOS
     &          ,DERIV)

*****************************************************************************
*
* PURPOSE
*
*     Compute the contribution of the derivatives of BFLU array with respect
*     to areal recharge to the RHS of the derivatives of h  wrt parameters
*
* DESCRIPTION
*
*     Compute the contribution of the derivatives of AFLU matrix with respect
*     to areal recharge to the RHS of the derivatives of h  wrt parameters
*     More precisely, adds to the RHS the expression
*
*                              d BFLU     
*                         ----------------
*                         d areal recharge
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  DERH                   Nodal head derivatives with respect to estimated      
*                         flow parameters.                                      
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
*  LNNDEL                 Number of nodes at every element                      
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
*  DTIM                   Value to compute the time function value at the
*                         current time (k+theta)
*  EPSFLU                 Time weighting parameter for nonlinear flow problems  
*  FNT                    Array containing time functions values                
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  INARR                  Index for areal recharge                              
*                         in array variables (LXPAREL and CFPAREL)              
*  INDSSTR                Current problem state. If 1, transient state,         
*                         if 0, steady-state                                    
*  INEW                   Index to locate the RHS of the derivatives of head
*                         at current time with respect to parameters inside
*                         array DERH
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOFMLF                 Flow Formulation number                               
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFLAGS                 Maximum number of allowed flags                       
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NINT                   Number of observation times                           
*  NPAREL                 Number of element parameters in current problem       
*  NPARF                  Number of transient parameters to be estimated        
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
*  IOPAR                  
*  IZON                   Areal recharge zone number of the current element
*  JJ                     Location of areal recharge values of the current
*                         areal recharge zone in arrays PARC, NFNLPAR, etc
*  L                      Current element
*  NCNF                   Nonlinear function number of the current
*                         areal recharge zone
*  NNUD                   Number of nodes of the current element
*  NPTOT                  Total number of parameters to be estimated in the
*                         current areal recharge zone (includes zonal plus
*                         generic parameters)
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  DER_PARAM              Computes the derivatives of the parameter with
*                         respect to its zonal value plus its generic parameters
*
* HISTORY
*
*     SCR      ?-1997     First coding
*     AMS     12-1998     Full modification
*
*****************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMDERH,IDIMFNT,INARR,INDSSTR,INEW,INTI,IOCAP
     &          ,IODENS,IOFLLI,IOFMLF,LMXNDL,NCNF,NFLAGS,NFNL,NINT
     &          ,NPAREL,NPARF,NPPEL,NPTOT,NTYPAR,NUMEL,NUMNP,NZPAR
     ;          ,IDIMWGT  ,IPNT_PAR ,NPAR

      REAL*8::BETAC,CREF,DENS,DENSREF,DTIM,EPSFLU

      INTEGER*4::IFLAGS(NFLAGS),INDEX(12),INORPAR(NTYPAR)
     &          ,IVPAR(NZPAR),KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR),NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL),NFTPAR(NZPAR),NZONE_PAR(NTYPAR)

      REAL*8::AREA(NUMEL),CFPAREL(NUMEL,NPAREL)
     &       ,DERH(NUMNP,NPARF,IDIMDERH),FNT(IDIMFNT,NINT)
     &       ,HCALIT(NUMNP),HCALAN(NUMNP),PARACD(3,NFNL),PARC(NZPAR)
     &       ,PAREL(NUMEL,NPPEL),WGT_PAR(IDIMWGT)

C------------------------- Internal


      INTEGER*4::I,IC,IP,IZON,JJ,K,L,NNUD,NZONECEXT

      

      REAL*8::AREALN,CEXT,DENSREC,RECHARGE,RHOAREALN


      INTEGER*4::IPOS(NPAR)

      REAL*8::CFPARAM(12),XPARAM(8),DERIV(NPAR)

C------------------------- FIRST EXECUTABLE STATEMENT.

C------------------------- Cross over all elements

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          IZON = LXPAREL(L,INARR)

          IF (IZON.NE.0) THEN

              IP = IVPAR(INORPAR(8)+IZON)

              IF (IP.GT.0 .OR. IOFLLI.NE.0) THEN

C------------------------- Computes derivatives

                  JJ = INORPAR(8) + IZON            ! Auxiliar variables

                  IF (IOFLLI.NE.0) THEN

                      NCNF = NFNLPAR(JJ)

                  ELSE

                      NCNF=0

                  END IF !IOFLLI.NE.0

                  INDEX(1) = JJ

                  CFPARAM(1) = CFPAREL(L,INARR)

                  CALL DER_PARAM
     & (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     & ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,L
     & ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     & ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,1        ,NUMEL
     & ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     & ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     & ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ; ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

                  IF (NPTOT.NE.0) THEN

                      AREALN = AREA(L)/NNUD
                      
                      NZONECEXT = LXPAREL(L,10) 
                      RECHARGE = PAREL(L,8)

C------------------------- Density dependent contribution.

                      IF (NZONECEXT.GT.0 .AND. RECHARGE.GT.0) THEN
                          CEXT = PAREL(L,15)
                      ELSE
                          CEXT = 0D0
                      END IF !NZONECEXT.GT.0 .AND. RECHARGE.GT.0
       
                      IF (IODENS.EQ.1) THEN
                          DENSREC =  DENS(DENSREF,BETAC,CEXT,CREF)
                      ELSE
                          DENSREC = 1D0
                      END IF !IODENS.EQ.1

                      RHOAREALN = DENSREC*AREALN

C------------------------- Cross all element nodes
 
                      DO K=1,NNUD

                          I=KXX(K,L)
                      

C------------------------- Computes element contribution to derivate

                           DO IC=1,NPTOT
                              DERH(I,IPOS(IC),INEW) =
     &                             DERH(I,IPOS(IC),INEW)
     &                           + RHOAREALN*DERIV(IC)
                          END DO !IC=1,NPTOT

                      END DO !K=1,NNUD

                  END IF ! NPTOT.NE.0

              END IF ! IP.GT.0 .OR. IOFLLI.NE.0

          END IF ! IZON.NE.0

      END DO ! L=!,NUMEL

      END SUBROUTINE DERARR
