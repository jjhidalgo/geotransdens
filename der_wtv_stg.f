       SUBROUTINE DER_WTV_STG
     &           (AREA     ,BETAC    ,BIBI     ,CAUX1    ,CAUX2
     &           ,CFPAREL  ,CREF     ,DENSITY  ,DENSREF  ,DERC
     &           ,DTIM     ,EPSFLU   ,HCALAN   ,HCALIT   ,HINI
     &           ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW
     &           ,INORPAR  ,INTI     ,IOCAP    ,IODENS   ,IOFLLI
     &           ,IOFMLF   ,IOVRWC   ,IVPAR    ,KXX      ,LDIM
     &           ,LMXNDL   ,LNNDEL   ,LXPAREL  ,NFLAGS   ,NFNL
     &           ,NINT     ,NPAR     ,NPAREL   ,NPPEL    ,NPZON
     &           ,NTYPAR   ,NUMEL    ,NUMNP    ,NZPAR    ,PARC
     &           ,PAREL    ,THETAT   ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR)

*******************************************************************************
*
* PURPOSE
*
*   Manages the derivatives of WATVOL w.r.t. storativity
*
* DESCRIPTION
*
*   Manages the derivatives of WATVOL w.r.t. storativity by calling the 
*   necessary rutines by NODES or by ELEMENTS
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  DER_WTV_STG_I       Derivatives of WATVOL by nodes
*  DER_WTV_STG_L       Derivatives of WATVOL by elements
*
* HISTORY
*
*     AMS      6-2002     First coding
*
*******************************************************************************

      IMPLICIT NONE


C------------------------- External
      
      INTEGER*4::IDIMBB   ,IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW
     &          ,INTI     ,IOCAP    ,IODENS   ,IOFLLI   ,IOFMLF
     &          ,IOVRWC   ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT,NPAR
     &          ,NPAREL   ,NPPEL    ,NPZON    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZPAR    ,IDIMWGT  ,IPNT_PAR

      REAL*8::BETAC,CREF,DENSREF,DTIM,EPSFLU,THETAT


      INTEGER*4::IFLAGS(NFLAGS)       ,INORPAR(NTYPAR),IVPAR(NZPAR)
     &          ,KXX(LMXNDL,NUMEL)    ,LDIM(NUMEL)    ,LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)        ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)



      REAL*8::AREA(NUMEL)               ,BIBI(IDIMBB,NUMEL)
     &       ,CAUX1(NUMNP)              ,CAUX2(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)     ,DENSITY(NUMEL)
     &       ,DERC(NUMNP,NPAR,IDIMDERC) ,FNT(IDIMFNT,NINT)
     &       ,HCALAN(NUMNP)             ,HCALIT(NUMNP)
     &       ,HINI(NUMNP)               ,PARACD(3,NFNL)
     &       ,PARC(NZPAR)               ,PAREL(NUMEL,NPPEL)
     &       ,WGT_PAR(IDIMWGT)


       IF (IOVRWC.EQ.1) THEN ! Elementwise

          CALL DER_WTV_STG_L
     &        (AREA    ,BETAC    ,BIBI     ,CAUX1    ,CAUX2
     &        ,CFPAREL  ,CREF     ,DENSITY  ,DENSREF  ,DERC
     &        ,DTIM     ,EPSFLU   ,FNT      ,HCALAN   ,HCALIT
     &        ,HINI     ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IDIMWGT  
     &        ,IFLAGS   ,INDSSTR  ,INEW     ,INORPAR  ,INTI
     &        ,IOCAP    ,IODENS   ,IOFLLI   ,IOFMLF   ,IPNT_PAR
     &        ,IVPAR    ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL
     &        ,LXPAREL  ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &        ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   
     &        ,NPPEL    ,NPZON    ,NTYPAR   ,NUMEL    ,NUMNP    
     &        ,NZONE_PAR,NZPAR    ,PARACD   ,PARC     ,PAREL
     &        ,THETAT   ,WGT_PAR)


       ELSE IF (IOVRWC.EQ.2) THEN         ! Nodewise

          CALL DER_WTV_STG_I
     &        (AREA     ,BETAC    ,BIBI     ,CAUX1    ,CAUX2
     &        ,CFPAREL  ,CREF     ,DENSREF  ,DERC     ,DTIM
     &        ,EPSFLU   ,FNT      ,HCALAN   ,HCALIT   ,HINI
     &        ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IDIMWGT  ,IFLAGS
     &        ,INDSSTR  ,INEW     ,INORPAR  ,INTI     ,IOCAP
     &        ,IODENS   ,IOFLLI   ,IOFMLF   ,IPNT_PAR ,IVPAR
     &        ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LXPAREL
     &        ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &        ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPPEL
     &        ,NPZON    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &        ,NZPAR    ,PARACD   ,PARC     ,PAREL    ,THETAT
     &        ,WGT_PAR)

       END IF !IOVRWC.EQ.1

      
       END SUBROUTINE DER_WTV_STG

************************************************************************
************************************************************************

      SUBROUTINE DER_WTV_STG_L
     &          (AREA    ,BETAC    ,BIBI     ,CAUX1    ,CAUX2
     &         ,CFPAREL  ,CREF     ,DENSITY  ,DENSREF  ,DERC
     &         ,DTIM     ,EPSFLU   ,FNT      ,HCALAN   ,HCALIT
     &         ,HINI     ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IDIMWGT  
     &         ,IFLAGS   ,INDSSTR  ,INEW     ,INORPAR  ,INTI
     &         ,IOCAP    ,IODENS   ,IOFLLI   ,IOFMLF   ,IPNT_PAR
     &         ,IVPAR    ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL
     &         ,LXPAREL  ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &         ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   
     &         ,NPPEL    ,NPZON    ,NTYPAR   ,NUMEL    ,NUMNP    
     &         ,NZONE_PAR,NZPAR    ,PARACD   ,PARC     ,PAREL
     &         ,THETAT   ,WGT_PAR)

*******************************************************************************
*
* PURPOSE
*
*   Computes derivatives of WATVOL w.r.t. storativity and assembles it.
*
* DESCRIPTION
*
*   Computes derivatives of WATVOL (BY ELEMENTS) w.r.t. storativity and 
*   assembles it in its RHS. 
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
*  HAUX1                  Array containing heads, ponderated by THETAF time
*                         weight. Is equal to THETAF*HCAL+(1-THETAF)*HCALAN
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
      
      INTEGER*4::IDIMBB   ,IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW
     &          ,INTI     ,IOCAP    ,IODENS   ,IOFLLI   ,IOFMLF
     &          ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT     ,NPAR
     &          ,NPAREL   ,NPPEL    ,NPZON    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZPAR    ,IDIMWGT  ,IPNT_PAR

      REAL*8::BETAC,CREF,DENSREF,DTIM,EPSFLU,THETAT

      REAL*8::DENS

      INTEGER*4::IFLAGS(NFLAGS)       ,INORPAR(NTYPAR),IVPAR(NZPAR)
     &          ,KXX(LMXNDL,NUMEL)    ,LDIM(NUMEL)    ,LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)        ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)



      REAL*8::AREA(NUMEL)               ,BIBI(IDIMBB,NUMEL)
     &       ,CAUX1(NUMNP)              ,CAUX2(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)     ,DENSITY(NUMEL)
     &       ,DERC(NUMNP,NPAR,IDIMDERC) ,FNT(IDIMFNT,NINT)
     &       ,HCALAN(NUMNP)             ,HCALIT(NUMNP)
     &       ,HINI(NUMNP)               ,PARACD(3,NFNL)
     &       ,PARC(NZPAR)               ,PAREL(NUMEL,NPPEL)
     &       ,WGT_PAR(IDIMWGT)

     

C------------------------- Internal

      INTEGER*4::I    ,I1   ,I2   ,IP   ,IPAR ,IPS  ,JJ   ,K1   ,K2
     &          ,KNT  ,L   ,LANI  ,LD   ,NCNF ,NNUD ,NPTOT,NZ

      REAL*8::AREALN ,C1     ,C2     ,CNODE    ,DENSL  ,DENSNODE
     &       ,DER_WTV,DER_WTVK1TH    ,FOD      ,HK1TH  ,HKTH
     &       ,S      ,SK1TH  ,STG      ,THT1

      INTEGER*4::INDEX(12)  ,IPOS(NPAR)

      REAL*8::CFPARAM(12)   ,DERIV_STG(NPAR) ,XPARAM(8)

C------------------------- First Executable statement

      THT1 = 1D0 - THETAT

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          NZ = LXPAREL(L,2)
          JJ = INORPAR(7)+NZ
          IP = IVPAR(JJ)
          DENSL =DENSITY(L)

          IF (IP.GT.0 .OR. IOFLLI.NE.0) THEN

              LD = LDIM(L)
              LANI = LD*(LD+1)/2
              
              IF (IOFLLI.NE.0) THEN

                  NCNF = NFNLPAR(JJ)

              ELSE

                  NCNF = 0

              END IF !IOFLLI.NE.0

              INDEX(1) = JJ

              CFPARAM(1) = CFPAREL(L,2)

              CALL DER_PARAM
     & (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)  
     & ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,L
     & ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     & ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,NPZON   ,NUMEL    
     & ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV_STG  ,FNT      
     & ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG  
     & ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ; ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

              IF (NPTOT.GT.0) THEN

                  DO IPAR=1,NPTOT

                      KNT = 0
                      S = 0D0
	                SK1TH = 0D0
                      IPS = IPOS(IPAR)

                      DO K1=1,NNUD

                          I = KXX(K1,L)

	                    HKTH = THETAT*HCALIT(I) + THT1*HCALAN(I)
		                HK1TH = THT1*HCALIT(I) + THETAT*HCALAN(I)

                          S = S + HKTH - HINI(I)
                          SK1TH = SK1TH + HK1TH - HINI(I)

                      END DO !K1=1,NNUD

                      DER_WTV = DERIV_STG(IPAR)*S/NNUD
                      DER_WTVK1TH = DERIV_STG(IPAR)*SK1TH/NNUD

C------------------------- Diffusion

                      DO K1=1,NNUD-1

                          I1 = KXX(K1,L)
                          C1 = CAUX1(I1)

                          DO K2=K1+1,NNUD

                              I2 = KXX(K2,L)
                              C2 = CAUX1(I2)
                              S = 0D0

                              DO I=1,LD

                                  S = S + PAREL(L,11)*BIBI(KNT+I,L)
     &                                    *DER_WTV*(C2-C1)
                              END DO !I=1,LD

                              S = DENSL*S

                              DERC(I1,IPS,INEW) = DERC(I1,IPS,INEW) - S
                              DERC(I2,IPS,INEW) = DERC(I2,IPS,INEW) + S
                              KNT = KNT + LANI

                          END DO !K2=K1+1,NNUD

                      END DO ! K1=1,NNUD-1

C------------------------- First order decay and storage

                      AREALN = AREA(L)/NNUD
                      FOD = PAREL(L,13)*DER_WTV*AREALN
                      STG = DER_WTVK1TH*AREALN

                      DO K1=1,NNUD

                          I1 = KXX(K1,L)

	                    IF (IODENS.EQ.1) THEN

                              CNODE = CAUX1(I1)
	                        DENSNODE = DENS(DENSREF,BETAC,CNODE,CREF)

	                    END IF !IODENS.EQ.1

                          DERC(I1,IPS,INEW) = DERC(I1,IPS,INEW)
     &                                      - DENSNODE*FOD*CAUX1(I1)

                          IF (INDSSTR.NE.0) THEN

                              DERC(I1,IPS,INEW) = DERC(I1,IPS,INEW)
     &                                          - DENSL*STG*CAUX2(I1)

                          END IF !INDSSTR.NE.0

                      END DO !K1=1,NNUD

                  END DO !IPAR=1,NPTOT

              END IF !NPTOT.GT.0

          END IF ! IP > 0, ie., storativity of this element is estimated

      END DO ! L

      END SUBROUTINE DER_WTV_STG_L

************************************************************************
************************************************************************

      SUBROUTINE DER_WTV_STG_I
     &          (AREA     ,BETAC    ,BIBI     ,CAUX1    ,CAUX2
     &          ,CFPAREL  ,CREF     ,DENSREF  ,DERC     ,DTIM
     &          ,EPSFLU   ,FNT      ,HCALAN   ,HCALIT   ,HINI
     &          ,IDIMBB   ,IDIMDERC ,IDIMFNT  ,IDIMWGT  ,IFLAGS
     &          ,INDSSTR  ,INEW     ,INORPAR  ,INTI     ,IOCAP
     &          ,IODENS   ,IOFLLI   ,IOFMLF   ,IPNT_PAR ,IVPAR
     &          ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LXPAREL
     &          ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPPEL
     &          ,NPZON    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &          ,NZPAR    ,PARACD   ,PARC     ,PAREL    ,THETAT
     &          ,WGT_PAR)
*******************************************************************************
*
* PURPOSE
*
*   Computes derivatives of WATVOL w.r.t. storativity and assembles it.
*
* DESCRIPTION
*
*   Computes derivatives of WATVOL (BY NODES) w.r.t. storativity and 
*   assembles it in its RHS. 
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
*  HAUX1                  Array containing heads, ponderated by THETAF time
*                         weight. Is equal to THETAF*HCAL+(1-THETAF)*HCALAN
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


      INTEGER*4::IDIMBB,IDIMDERC,IDIMFNT,INDSSTR,INEW,INTI,LMXNDL,NFLAGS
     &          ,NFNL,NINT,NPAR,NPAREL,NPPEL,NTYPAR,NUMEL,NUMNP,NZPAR
     &          ,IOCAP,NPZON,IOFLLI,IOFMLF,IODENS,IDIMWGT  ,IPNT_PAR

      REAL*8::BETAC,CREF,DENSREF,DTIM,EPSFLU,THETAT

      REAL*8::DENS

      INTEGER*4::IFLAGS(NFLAGS)       ,INORPAR(NTYPAR),IVPAR(NZPAR)
     &          ,KXX(LMXNDL,NUMEL)    ,LDIM(NUMEL)    ,LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)        ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)



      REAL*8::AREA(NUMEL)            ,BIBI(IDIMBB,NUMEL)
     &       ,CAUX1(NUMNP)           ,CAUX2(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)  ,DERC(NUMNP,NPAR,IDIMDERC)
     &       ,FNT(IDIMFNT,NINT)      ,HCALAN(NUMNP)
     &       ,HCALIT(NUMNP)          ,HINI(NUMNP)
     &       ,PARACD(3,NFNL)         ,PARC(NZPAR)
     &       ,PAREL(NUMEL,NPPEL)     ,WGT_PAR(IDIMWGT)

C------------------------- Internal

      INTEGER*4::I    ,I1   ,I2   ,IP   ,IPAR ,IPS  ,JJ   ,K1   ,K2
     &          ,KNT  ,L   ,LANI  ,LD   ,NCNF ,NNUD ,NPTOT,NZ

      REAL*8::AREALN   ,C1       ,C2       ,CNODE    ,DENSNODE
     &       ,DER_WTV_AVG        ,FOD      ,HK1TH    ,HKTH
     &       ,S       ,STG       ,THT1

      INTEGER*4::INDEX(12)  ,IPOS(NPAR)

      REAL*8::CFPARAM(12) ,DERIV_STG(NPAR),DER_WTV(10),DER_WTVK1TH(10)
     &       ,XPARAM(8)


C------------------------- First Executable statement

      DENSNODE = 1D0
      THT1 = 1D0 - THETAT

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          NZ = LXPAREL(L,2)
          JJ = INORPAR(7)+NZ
          IP = IVPAR(JJ)

          IF (IP.GT.0 .OR. IOFLLI.NE.0) THEN

              LD = LDIM(L)
              LANI = LD*(LD+1)/2

              IF (IOFLLI.NE.0) THEN

                  NCNF = NFNLPAR(JJ)

              ELSE

                  NCNF = 0

              END IF !IOFLLI.NE.0

              INDEX(1) = JJ

              CFPARAM(1) = CFPAREL(L,2)


              CALL DER_PARAM
     & (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)  
     & ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,L
     & ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF   ,NFNLTIP
     & ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,NPZON   ,NUMEL    
     & ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV_STG  ,FNT      
     & ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG  
     & ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ; ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

              IF (NPTOT.GT.0) THEN

                  DO IPAR=1,NPTOT

                      DER_WTV_AVG = 0D0
                      KNT = 0
                      IPS = IPOS (IPAR)
                      AREALN = AREA(L)/NNUD

                      DO K1=1,NNUD

                          I = KXX(K1,L)

                         HKTH = THETAT*HCALIT(I) + THT1*HCALAN(I)
                         HK1TH = THT1*HCALIT(I) + THETAT*HCALAN(I)

                          DER_WTV(K1) = DERIV_STG(IPAR)*(HKTH - HINI(I))

                         DER_WTVK1TH(K1) = 
     &                                 DERIV_STG(IPAR)*(HK1TH - HINI(I))


                          DER_WTV_AVG = DER_WTV_AVG + DER_WTV(K1)

                      END DO

                      DER_WTV_AVG = DER_WTV_AVG/NNUD

C------------------------- Diffusion

                      DO K1=1,NNUD-1

                          I1 = KXX(K1,L)
                          C1 = CAUX1(I1)

                            DO K2=K1+1,NNUD

                              I2 = KXX(K2,L)
                              C2 = CAUX1(I2)
                              S = 0D0

                              DO I=1,LD
                                  S = S + PAREL(L,11)*BIBI(KNT+I,L)
     &                                    *DER_WTV_AVG*(C2-C1)
                          END DO !I=1,LD

                              DERC(I1,IPS,INEW) = DERC(I1,IPS,INEW) - S
                              DERC(I2,IPS,INEW) = DERC(I2,IPS,INEW) + S
                              KNT = KNT + LANI

                          END DO !K2=K1+1,NNUD

                      END DO ! K1=1,NNUD-1

C------------------------- First order decay and storage
                  
                      
                      DO K1=1,NNUD

                         I1 = KXX(K1,L)

                         IF (IODENS.EQ.1) THEN

                              CNODE = CAUX1(I1)
                             DENSNODE = DENS(DENSREF,BETAC,CNODE,CREF)

                         END IF !IODENS.EQ.1

                          STG = DER_WTVK1TH(K1)*AREALN
                          
                          FOD = PAREL(L,13)*DER_WTV(K1)*AREALN

                          DERC(I1,IPS,INEW) = DERC(I1,IPS,INEW)
     &                                      - DENSNODE*FOD*CAUX1(I1)

                          IF (INDSSTR.NE.0) THEN

                              DERC(I1,IPS,INEW) = DERC(I1,IPS,INEW)
     &                                          - DENSNODE*STG*CAUX2(I1)

                          END IF !IF INDSSTR.NE.0

                      END DO !K1=1,NNUD

                  END DO !IPAR=1,NPTOT

              END IF !NPTOT.GT.0

          END IF ! IP > 0, ie., storativity of this element is estimated

      END DO ! L=!,NUMEL


      END SUBROUTINE DER_WTV_STG_I
