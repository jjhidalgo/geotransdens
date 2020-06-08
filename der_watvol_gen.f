      SUBROUTINE DER_WATVOL_GEN
     &          (IDIMBB   ,IDIMDERC ,IDIMDERH ,INEW     ,IOVRWC   ,IPAR
     &          ,LMXNDL   ,NPAR     ,NPARF    ,NPPEL    ,NUMEL    ,NUMNP
     &          ,THETAT   ,AREA     ,BIBI     ,CAUX1    ,CAUX2    ,DERC
     &          ,DERH     ,KXX      ,LDIM     ,LNNDEL   ,PAREL    ,INEWH
     &          ,IOLDH)

*******************************************************************************
*
* PURPOSE
*
*   Manages the derivatives of WATVOL w.r.t. flow parameters except POR & STG
*
* DESCRIPTION
*
*   Manages the derivatives of WATVOL w.r.t. flow parameters except POR & STG, 
*   by calling as necessary the corresponding rutines by NODES or by ELEMENTS
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  DER_WATVOL_GEN_I       Derivatives of WATVOL by nodes
*  DER_WATVOL_GEN_L       Derivatives of WATVOL by elements
*
* HISTORY
*
*     AMS      6-2002     First coding
*
*******************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMBB   ,IDIMDERC ,IDIMDERH ,INEW     ,INEWH
     &          ,IOLDH    ,IOVRWC   ,IPAR     ,LMXNDL   ,NPAR
     &          ,NPARF    ,NPPEL    ,NUMEL    ,NUMNP

      REAL*8::DTH,THETAT

      INTEGER*4::LDIM(NUMEL),LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)

      REAL*8::AREA(NUMEL)               ,BIBI(IDIMBB,NUMEL)
     &       ,CAUX1(NUMNP)              ,CAUX2(NUMNP)
     &       ,DERC(NUMNP,NPAR,IDIMDERC) ,DERH(NUMNP,NPARF,IDIMDERH)
     &       ,PAREL(NUMEL,NPPEL)

      DTH=1.D0-THETAT

      IF (IOVRWC.EQ.1) THEN

          CALL DER_WATVOL_GEN_L
     &        (DTH      ,IDIMBB   ,IDIMDERC ,IDIMDERH ,INEW
     &        ,IPAR     ,LMXNDL   ,NPAR     ,NPARF    ,NPPEL
     &        ,NUMEL    ,NUMNP    ,THETAT   ,AREA     ,BIBI
     &        ,CAUX1    ,CAUX2    ,DERC     ,DERH     ,KXX
     &        ,LDIM     ,LNNDEL   ,PAREL    ,INEWH    ,IOLDH)

      ELSE

          CALL DER_WATVOL_GEN_I
     &        (DTH      ,IDIMBB   ,IDIMDERC ,IDIMDERH ,INEW
     &        ,IPAR     ,LMXNDL   ,NPAR     ,NPARF    ,NPPEL
     &        ,NUMEL    ,NUMNP    ,THETAT   ,AREA     ,BIBI
     &        ,CAUX1    ,CAUX2    ,DERC     ,DERH     ,KXX
     &        ,LDIM     ,LNNDEL   ,PAREL    ,INEWH    ,IOLDH)

      END IF !IOVRWC.EQ.1,2


      END SUBROUTINE DER_WATVOL_GEN

************************************************************************
************************************************************************

       SUBROUTINE DER_WATVOL_GEN_L
     &           (DTH      ,IDIMBB   ,IDIMDERC ,IDIMDERH ,INEW
     &           ,IPAR     ,LMXNDL   ,NPAR     ,NPARF    ,NPPEL
     &           ,NUMEL    ,NUMNP    ,THETAT   ,AREA     ,BIBI
     &           ,CAUX1    ,CAUX2    ,DERC     ,DERH     ,KXX
     &           ,LDIM     ,LNNDEL   ,PAREL    ,INEWH    ,IOLDH)

*******************************************************************************
*
* PURPOSE
*
*   Computes generic derivatives of WATVOL and assembles it in its RHS
*
* DESCRIPTION
*
*   Computes generic derivatives of WATVOL (BY ELEMENTS) and assembles it in 
*   its RHS. It does not include derivatives of WATVOL w.r.t. porosity and 
*   storativity
*
* EXTERNAL VARIABLES: ARRAYS
*
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
*  DERH                   Nodal head derivatives with respect to estimated      
*                         flow parameters.                                      
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LNNDEL                 Number of nodes at every element                      
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
*  DTH                    =THETAT for INTI=1 and =1-THETAT for INTI>1
*  IDIMBB                 Used to dimension array BIBI. Is equal to IDIMQ times 
*                         the maximum possible anisotropy of the problem        
*  INEW                   Index for the RHS of derivatives of concentrations
*  IPAR                   Current flow parameter
*  LMXNDL                 Maximum number of nodes per element                   
*  NPAR                   Total number of parameters to be estimated            
*  NPARF                  Number of transient parameters to be estimated        
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  THETAT                 Time weighting parameter for transport problems       
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

      INTEGER*4::IDIMBB   ,IDIMDERC ,IDIMDERH ,INEW     ,INEWH
     &          ,IOLDH    ,IPAR     ,LMXNDL   ,NPAR     ,NPARF
     &          ,NPPEL    ,NUMEL    ,NUMNP

      REAL*8::DTH,THETAT

      INTEGER*4::LDIM(NUMEL),LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)

      REAL*8::AREA(NUMEL)               ,BIBI(IDIMBB,NUMEL)
     &       ,CAUX1(NUMNP)              ,CAUX2(NUMNP)
     &       ,DERC(NUMNP,NPAR,IDIMDERC) ,DERH(NUMNP,NPARF,IDIMDERH)
     &       ,PAREL(NUMEL,NPPEL)

C------------------------- Internal

      INTEGER*4::I   ,I1  ,I2  ,K   ,K1  ,K2  ,KNT ,L   ,LANI,LD  ,NNUD

      REAL*8::C1     ,C2     ,DER_WTV,DER_WTVK1TH    ,DH_AVG ,DH_AVGK1TH
     &       ,S      ,VAUX1  ,VAUX2

C------------------------- First executable statement

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          LD = LDIM(L)
          LANI = LD*(LD+1)/2
          DH_AVG = 0D0
	  DH_AVGK1TH = 0D0

C------------------------- Derivative of water volume of the element 

          DO K=1,NNUD

              I = KXX(K,L)
              DH_AVG = DH_AVG + (THETAT*DERH(I,IPAR,IOLDH)
     &                                    + DTH*DERH(I,IPAR,INEWH))

              DH_AVGK1TH = DH_AVGK1TH + (DTH*DERH(I,IPAR,IOLDH)
     &                                    + THETAT*DERH(I,IPAR,INEWH))
          END DO !K=1,NNUD


          DER_WTV = PAREL(L,7)*DH_AVG/NNUD
	    DER_WTVK1TH = PAREL(L,7)*DH_AVGK1TH/NNUD

C------------------------- Assembles the RHS

          KNT = 0

C------------------------- Diffusion

          DO K1=1,NNUD-1

              I1 = KXX(K1,L)
              C1 = CAUX1(I1)

              DO K2=K1+1,NNUD

                  I2 = KXX(K2,L)
                  C2 = CAUX1(I2)
                  S = 0D0

                  DO I=1,LD
                      S = S + BIBI(KNT+I,L)
                  END DO !I=1,LD

                  S = S*PAREL(L,11)*DER_WTV
                  DERC(I1,IPAR,INEW) = DERC(I1,IPAR,INEW) - S*(C2-C1)
                  DERC(I2,IPAR,INEW) = DERC(I2,IPAR,INEW) + S*(C2-C1)
                  KNT = KNT + LANI

              END DO !K2=K1+1,NNUD
          END DO ! K1=1,NNUD-1

C------------------------- First order decay and storage

          VAUX1 = DER_WTVK1TH*AREA(L)/NNUD    ! Storage
          VAUX2 = PAREL(L,13)*VAUX1           ! FOD

          DO K1=1,NNUD

              I1 = KXX(K1,L)
              DERC(I1,IPAR,INEW) = DERC(I1,IPAR,INEW) - VAUX2*CAUX1(I1)
     &                           - VAUX1*CAUX2(I1)
          END DO !K1=1,NNUD

      END DO !L=1,NUMEL

      END SUBROUTINE DER_WATVOL_GEN_L

************************************************************************
************************************************************************

      SUBROUTINE DER_WATVOL_GEN_I
     &          (DTH      ,IDIMBB   ,IDIMDERC ,IDIMDERH ,INEW
     &          ,IPAR     ,LMXNDL   ,NPAR     ,NPARF    ,NPPEL
     &          ,NUMEL    ,NUMNP    ,THETAT   ,AREA     ,BIBI
     &          ,CAUX1    ,CAUX2    ,DERC     ,DERH     ,KXX
     &          ,LDIM     ,LNNDEL   ,PAREL    ,INEWH    ,IOLDH)

*******************************************************************************
*
* PURPOSE
*
*   Computes generic derivatives of WATVOL and assembles it in its RHS
*
* DESCRIPTION
*
*   Computes generic derivatives of WATVOL (BY NODES) and assembles it in 
*   its RHS. It does not include derivatives of WATVOL w.r.t. porosity and 
*   storativity
*
* EXTERNAL VARIABLES: ARRAYS
*
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
*  DERH                   Nodal head derivatives with respect to estimated      
*                         flow parameters.                                      
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LNNDEL                 Number of nodes at every element                      
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
*  DTH                    =THETAT for INTI=1 and =1-THETAT for INTI>1
*  IDIMBB                 Used to dimension array BIBI. Is equal to IDIMQ times 
*                         the maximum possible anisotropy of the problem        
*  INEW                   Index for the RHS of derivatives of concentrations
*  IPAR                   Current flow parameter
*  LMXNDL                 Maximum number of nodes per element                   
*  NPAR                   Total number of parameters to be estimated            
*  NPARF                  Number of transient parameters to be estimated        
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  THETAT                 Time weighting parameter for transport problems       
*
* INTERNAL VARIABLES: SCALARS
*
*  DER_WTV_AVG            Average of WATVOL in the current element
*  DER_WTV                Derivative of WATVOL in the nodes of the current 
*                         element
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

      INTEGER*4::IDIMBB   ,IDIMDERC ,IDIMDERH ,INEW     ,INEWH
     &          ,IOLDH    ,IPAR     ,LMXNDL   ,NPAR     ,NPARF
     &          ,NPPEL    ,NUMEL    ,NUMNP

      REAL*8::DTH,THETAT

      INTEGER*4::LDIM(NUMEL),LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)

      REAL*8::AREA(NUMEL)               ,BIBI(IDIMBB,NUMEL)
     &       ,CAUX1(NUMNP)              ,CAUX2(NUMNP)
     &       ,DERC(NUMNP,NPAR,IDIMDERC) ,DERH(NUMNP,NPARF,IDIMDERH)
     &       ,PAREL(NUMEL,NPPEL)

C------------------------- Internal

      INTEGER*4::I   ,I1  ,I2  ,K   ,K1  ,K2  ,KNT ,L   ,LANI,LD  ,NNUD

	REAL*8::C1     ,C2     ,DER_WTV_AVG    ,S      ,VAUX1    ,VAUX2

	REAL*8::DER_WTV(10),DER_WTVK1TH(10)

C------------------------- First executable statement

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          LD = LDIM(L)
          LANI = LD*(LD+1)/2
          DER_WTV_AVG = 0D0

          DO K=1,NNUD

              I = KXX(K,L)
              DER_WTV(K) = PAREL(L,7)*
     &              (THETAT*DERH(I,IPAR,IOLDH) + DTH*DERH(I,IPAR,INEWH))
              DER_WTV_AVG = DER_WTV_AVG+DER_WTV(K)

              DER_WTVK1TH(K) = PAREL(L,7)*
     &              (DTH*DERH(I,IPAR,IOLDH) + THETAT*DERH(I,IPAR,INEWH))

          END DO !K=1,NNUD

          DER_WTV_AVG = DER_WTV_AVG/NNUD

C------------------------- Derivative of water volume of the element 


C------------------------- Assembles the RHS

          KNT = 0

C------------------------- Diffusion

          DO K1=1,NNUD-1

              I1 = KXX(K1,L)
              C1 = CAUX1(I1)

              DO K2=K1+1,NNUD

                  I2 = KXX(K2,L)
                  C2 = CAUX1(I2)
                  S = 0D0

                  DO I=1,LD

                      S =S+PAREL(L,11)*BIBI(KNT+I,L)*DER_WTV_AVG*(C2-C1)

                  END DO !I=1,LD

                  DERC(I1,IPAR,INEW) = DERC(I1,IPAR,INEW) - S
                  DERC(I2,IPAR,INEW) = DERC(I2,IPAR,INEW) + S
                  KNT = KNT + LANI

              END DO !K2=K1+1,NNUD

          END DO !K1=1,NNUD-1

C------------------------- First order decay and storage

          VAUX1 = AREA(L)/NNUD      ! Storage
          VAUX2 = PAREL(L,13)*VAUX1 ! FOD

          DO K1=1,NNUD

              I1 = KXX(K1,L)
              DERC(I1,IPAR,INEW) = DERC(I1,IPAR,INEW)
     &                           - VAUX2*DER_WTV(K1)*CAUX1(I1)
     &                            -VAUX1*DER_WTVK1TH(K1)*CAUX2(I1)
          END DO !K1=1,NNUD

      END DO !L=1,NUMEL

      END SUBROUTINE DER_WATVOL_GEN_I
