      SUBROUTINE DERPOR
     &          (ACTH     ,AREA     ,BETAC    ,BIBI     ,CAUX1
     &          ,CAUX2    ,CCALAN   ,CCALIT   ,CFPAREL  ,CREF
     &          ,DENSITY  ,DENSREF  ,DERC     ,DERCS    ,DTIM
     &          ,EPSTRA   ,FNT      ,ICMPDERS ,IDIMBB   ,IDIMDERC
     &          ,IDIMFNT  ,IFLAGS   ,INDSSTR  ,INEW     ,INORPAR
     &          ,INTI     ,IOCAP    ,IODENS   ,IOFMLT   ,IOTRLI
     &          ,IOVRWC   ,ITPTVAR  ,IVPAR    ,KXX      ,LDIM
     &          ,LMXNDL   ,LNNDEL   ,LXPAREL  ,NFLAGS   ,NFNL
     &          ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT
     &          ,NPAR     ,NPAREL   ,NPPEL    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &          ,PAREL    ,THETAT   ,WSPECHEAT,WTHERMCON,IDIMWGT
     &          ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)

*******************************************************************************
*
* PURPOSE
*
*   Computes generic derivatives of WATVOL w.r.t. porosity and assembles it.
*
* DESCRIPTION
*
*   Computes generic derivatives of WATVOL (either by nodes or elements)
*   w.r.t. porosity and assembles it in its RHS.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  ACTH                   Aquifer thickness of every element. Cross sectional   
*                         area for 1-D elements, thickness for 2-D elements.    
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  BIBI                   Array containing the product of interpolation         
*                         functions gradient, for a given element               
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  CAUX2                  Array containing diference of concentrations in two   
*                         consecutives times, related to time increment         
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
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
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMBB                 Used to dimension array BIBI. Is equal to IDIMQ times 
*                         the maximum possible anisotropy of the problem        
*  INDSSTR                Current problem state. If 1, transient state,         
*                         if 0, steady-state                                    
*  INEW                   Index for the RHS of derivatives of concentrations
*  LMXNDL                 Maximum number of nodes per element                   
*  NFLAGS                 Maximum number of allowed flags                       
*  NINT                   Number of observation times                           
*  NPAR                   Total number of parameters to be estimated            
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
*
* INTERNAL VARIABLES: SCALARS
*
*  DER_WTV                Derivative of WATVOL in the current element
*  KNT                    Counter to locate elements in BIBI array
*  L                      Current element
*  LANI                   Maximum possible anisotropy of a given element. Equal
*                         to LDIM*(LDIM+1)/2
*  LD                     Dimension of the current element
*  NNUD                   Number of nodes of the current element
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY
*
*     AMS      6-2002     First coding
*
*******************************************************************************

      IMPLICIT NONE


C------------------------- External

      INTEGER*4::ICMPDERS ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IDIMWGT
     &          ,INDSSTR  ,INEW     ,INTI     ,IOCAP    ,IODENS
     &          ,IOFMLT   ,IOTRLI   ,IOVRWC   ,IPNT_PAR ,ITPTVAR
     &          ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT     ,NPAR
     &          ,NPAREL   ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP
     &          ,NZPAR

      REAL*8::BETAC    ,CREF     ,DENSREF  ,DTIM     ,EPSTRA
     &       ,THETAT   ,WSPECHEAT,WTHERMCON

      REAL*8::DENS

      INTEGER*4::IFLAGS(NFLAGS)       ,INORPAR(NTYPAR),IVPAR(NZPAR)
     &          ,KXX(LMXNDL,NUMEL)    ,LDIM(NUMEL)    ,LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)        ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)


      REAL*8::ACTH(NUMEL)              ,AREA(NUMEL)
     &       ,BIBI(IDIMBB,NUMEL)       ,CAUX1(NUMNP)
     &       ,CAUX2(NUMNP)             ,CCALAN(NUMNP)
     &       ,CCALIT(NUMNP)            ,CFPAREL(NUMEL,NPAREL)
     &       ,DENSITY(NUMEL)           ,DERC(NUMNP,NPAR,IDIMDERC)
     &       ,DERCS(NUMNP,NPAR,2)      ,PARC(NZPAR)
     &       ,FNT(IDIMFNT,NINT)        ,PAREL(NUMEL,NPPEL)
     &       ,PARACD(3,NFNL)           ,WGT_PAR(IDIMWGT)

C------------------------- Internal

      INTEGER*4::I    ,I1   ,I2   ,IC   ,IP   ,JJ   ,K1   ,K2   ,KNT
     &          ,L    ,LANI ,LD   ,NCNF ,NNUD ,NPTOT,NZ

      REAL*8::ACTHL    ,AREALN   ,C1       ,C2       ,C2_C1
     &       ,CAVG     ,CNODE    ,CONCTH1  ,DENSL    ,DENSLK1TH
     &       ,DENSNODE ,DENSNK1TH,DER_WTV  ,FOD      ,S
     &       ,SDER     ,STG      ,SWCOND   ,TH1      ,WSPHEATINV

      INTEGER*4::INDEX(12),IPOS(NPAR)

      REAL*8::CFPARAM(12),DERIV(NPAR),XPARAM(8)


C------------------------- First executable statement

      TH1 = 1D0 - THETAT

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          NZ = LXPAREL(L,7)
          JJ = INORPAR(15) + NZ
          IP = IVPAR(JJ)

          IF (IP.GT.0 .OR. IOTRLI.NE.0) THEN

              IF (IOTRLI.NE.0) THEN
                  NCNF = NFNLPAR(JJ)
              ELSE
                  NCNF = 0
              END IF !IOTRLI.NE.0

              INDEX(1) = JJ

              CFPARAM(1) = CFPAREL(L,7)

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

                  LD = LDIM(L)
                  LANI = LD*(LD+1)/2
                  AREALN = AREA(L)/NNUD
                  ACTHL = ACTH(L)

C------------------------- Computes element average density at k+thetat
C------------------------- for diffusion-dispersion terms and at
C------------------------- k + 1 - thetat for storage term if transient
C------------------------- transport is solved.

                  IF (IODENS.EQ.1) THEN

                      DENSL = DENSITY(L)

                      IF (INDSSTR.NE.0) THEN

                          CAVG = 0D0

                          DO I=1,NNUD
                              I1 = KXX(I,L)
                              CAVG = CAVG + TH1*CCALIT(I1)
     &                                    + THETAT*CCALAN(I1)
                          END DO !I=1,NNUD

                          CAVG = CAVG/NNUD
                          DENSLK1TH = DENS(DENSREF,BETAC,CAVG,CREF)

                      END IF !INDSSTR.NE.0

                  ELSE

                      DENSL = 1D0
                      DENSLK1TH  = 1D0

                  END IF !IODENS.EQ.1

                  IF (ITPTVAR.EQ.1) THEN

                      WSPHEATINV = 1D0/(WSPECHEAT)

                  END IF !ITPTVAR.EQ.1

                  !DO K=1,NNUD

                  !    DER_WTV(K) = PAREL(L,12)/PARC(INORPAR(15)+NZ)*ACTH(L)
                  !    DER_WTV_AVG = DER_WTV_AVG + DER_WTV(K)

                  !END DO !K=1,NNUD

                  !DER_WTV_AVG = DER_WTV_AVG/NNUD

                  KNT = 0

C------------------------- Diffusion

                  DO K1=1,NNUD-1

                      I1 = KXX(K1,L)
                      C1 = CAUX1(I1)

                      DO K2=K1+1,NNUD

                          I2 = KXX(K2,L)
                          C2 = CAUX1(I2)
                          S = 0D0
                          SWCOND = 0D0
                          C2_C1 = (C2-C1)

                          DO I=1,LD


                              IF (ITPTVAR.EQ.0) THEN !Solute

                                  S = S+PAREL(L,11)*BIBI(KNT+I,L)*C2_C1

                              ELSE !energy

                                  S = S+PAREL(L,11)*BIBI(KNT+I,L)*C2_C1

                                  SWCOND = SWCOND
     &                                   + WTHERMCON*BIBI(KNT+I,L)*C2_C1

                              END IF !ITPTVAR.EQ.0

                          END DO !I=1,LD


                          DO IC=1,NPTOT

                              DER_WTV = DERIV(IC)*ACTHL*DENSL

                              IF (ITPTVAR.EQ.0) THEN

                                  SDER = S*DER_WTV

                              ELSE

                                  SDER = WSPHEATINV
     &                                   *(SWCOND*DER_WTV - DERIV(IC)*S)
                             
                              END IF !ITPTVAR.EQ.0


                              DERC(I1,IPOS(IC),INEW) =
     &                                     DERC(I1,IPOS(IC),INEW) - SDER

                              DERC(I2,IPOS(IC),INEW) =
     &                                     DERC(I2,IPOS(IC),INEW) + SDER

                          END DO !IC=1,NPTOT

                          KNT =KNT + LANI

                      END DO !K2=K1+1,NNUD

                  END DO ! K1=1,NNUD-1

C------------------------- First order decay (always nodewise) and storage

                  DO K1=1,NNUD

                      I1 = KXX(K1,L)

C------------------------- Computes node density at k+thetat for FOD
C------------------------- terms and at k + 1 - thetat for storage term.

                      IF (IODENS.EQ.1) THEN

                          CNODE = CAUX1(I1)
                          DENSNODE = DENS(DENSREF,BETAC,CNODE,CREF)

                          CONCTH1 = TH1*CCALIT(I1) + THETAT*CCALAN(I1)
                          DENSNK1TH = DENS(DENSREF,BETAC,CONCTH1,CREF)

                      ELSE

                          DENSNODE = 1D0
                          DENSNK1TH = 1D0

                      END IF !IODENS.EQ.1
                      DO IC=1,NPTOT

                          DER_WTV = DERIV(IC)*ACTHL*DENSNODE

                          FOD = PAREL(L,13)*DER_WTV*AREALN*CAUX1(I1)

                          DERC(I1,IPOS(IC),INEW) =
     &                                      DERC(I1,IPOS(IC),INEW) - FOD

C------------------------- Adds storativity term for transient transport

                          IF (INDSSTR.NE.0) THEN

                              IF (IOVRWC.LT.2) THEN

                                  DER_WTV = DERIV(IC)*ACTHL*DENSLK1TH

                              ELSE

                                  DER_WTV = DERIV(IC)*ACTHL*DENSNK1TH

                  	        END IF !IOVRWC.LT.2

                              STG = DER_WTV*AREALN*CAUX2(I1)

                              DERC(I1,IPOS(IC),INEW)=
     &                                      DERC(I1,IPOS(IC),INEW) - STG

                          END IF !INDSSTR.NE.0


C------------------------- In a radioactive chain, part of the derivatives
C------------------------- of the RHS (simulation) of the "son" w.r.t the
C------------------------- parameters of the "father" is added to the RHS
C------------------------- (inversion) of the "son"
 
                          IF (ICMPDERS.EQ.1) THEN

                              DERCS(I1,IPOS(IC),INEW) =
     &                                     DERCS(I1,IPOS(IC),INEW) + FOD

	                    END IF !ICMPDERS.EQ.1
                      END DO !IC=1,NPTOT
                  END DO !K1=1,NNUD

              END IF !NPTOT.GT.0

          END IF ! IP > 0, ie., porosity of this element is estimated

      END DO !L=1,NUMEL

      END SUBROUTINE DERPOR
