      SUBROUTINE COMP_ATRA
     & (AREA     ,ATRA     ,BETAC    ,BIBI     ,CAUX1    ,DAT_VD
     & ,DENSITY  ,DENSREF  ,DTRADFLU ,DTRADTRA ,DVDH     ,DVDC
     & ,DWDH     ,EPSFLU   ,EPSTRA   ,GRDFF    ,IDIMBB   ,IDIMDENS
     & ,IDIMDQ   ,IDIMQ    ,IODENS   ,IODIM    ,IOINV    ,IOVRWC
     & ,ITPTVAR  ,LINMET   ,LMXNDL   ,NPPEL    ,NTYPAR   ,NUMEL
     & ,NUMNP    ,KXX      ,LDIM     ,LNNDEL   ,LTYPE    ,NZONE_PAR
     & ,PAREL    ,QXYZ     ,THETAT   ,VD       ,WATVOL   ,WSPECHEAT
     & ,WTHERMCON,XNORVD)


********************************************************************************
*
* PURPOSE
*
*  Computes of ATRA matrix without boundary conditions
*
*
* DESCRIPTION
*
*  Computes the transport coefficient matrix (ATRA) that includes 
*  dispersion, advection, first order decay term, sink/sources and age term 
*  whitout boundary conditions under any circumstances (linear/nonlinear, 
*  steady/transient)
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  ATRA                   Matrix of finite elements equations for transport     
*                         problem. No boundary conditions are included.         
*  BIBI                   Array containing the product of interpolation         
*                         functions gradient, for a given element               
*
*  DENSITY                Array containing the density of every element.
*
*  GRDFF                  Array containing the product between interpolation    
*                         functions integrals and interp. functions gradient    
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PAREL                  
*  QXYZ                   Products between the different components of          
*                         Darcy's velocity divided by its norm                  
*  VD                     Darcy's velocity                                      
*  WATVOL                                                                       
*  XNORVD                 Euclidean norm of Darcy's velocity                    
*
* INTERNAL VARIABLES: ARRAYS
*
*  DAUX                                                                         
*  SUM                                                                          
*  VV                                                                           
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMBB                 Used to dimension array BIBI                          
*  IDIMQ                  Used to dimension array QXYZ                          
*  IOCNST                 Scheme for mass storage term in transport problem     
*                             1.Lumped.
*  IODIM                  Maximum dimension of any element included             
*                         in the problem                                        
*  IOINV                  Inverse problem option                                
*  LMXNDL                 Maximum number of nodes per element
*  ITPTVAR                Indicates the type of state variable of transport.
*                             0. Mass fraction
*                             1. Temperature.
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NPPEL                                                                        
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*
* INTERNAL VARIABLES: SCALARS
*
*  DISPIJ                                                                       
*  LANI                                                                         
*  LTY                                                                          
*  NNUD                   Number of nodes of the current element                
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  COMDADQ                                                                      
*  ZERO_ARRAY             Initializes to zero all components of a given array
*
* HISTORY
*
*     AMS        1988     First coding
*     JCA      4-1998     Revision: Inclusion of header and some small modific.
*     AMS     12-1998     Revision: Inclusion of nonlinear calls and addition 
*                         of comments
*     JHG      5-2003     Inclusion of density dependence and elementwise
*                         computation of ATRA and its derivatives.
*
******************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMBB   ,IDIMDENS ,IDIMDQ   ,IDIMQ    ,IOCALCDEVF
     &          ,IODENS   ,IODIM    ,IOINV    ,IOVRWC   ,ITPTVAR
     &          ,LMXNDL   ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP

      REAL*8::BETAC,DENSREF,EPSFLU,EPSTRA,THETAT,WSPECHEAT,WTHERMCON


      INTEGER*4::KXX(LMXNDL, NUMEL)  ,LDIM(NUMEL)   ,LINMET(3,2)
     &          ,LNNDEL(NUMEL)       ,LTYPE(NUMEL)  ,NZONE_PAR(NTYPAR)

      REAL*8::AREA(NUMEL)                  ,ATRA(NUMEL,LMXNDL*LMXNDL)
     &       ,BIBI(IDIMBB, NUMEL)          ,CAUX1(NUMNP)
     &       ,DAT_VD(IODIM,IDIMDQ,NUMEL)   ,DENSITY(IDIMDENS)
     &       ,DTRADFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DTRADTRA(NUMEL,LMXNDL*LMXNDL),DVDC(LMXNDL,IODIM,NUMEL)
     &       ,DVDH(LMXNDL,IODIM,NUMEL)     ,DWDH(1,NUMEL)
     &       ,GRDFF(IODIM, LMXNDL, NUMEL)  ,PAREL(NUMEL,NPPEL)
     &       ,QXYZ(IDIMQ, NUMEL)           ,VD(IODIM, NUMEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
     &       ,XNORVD(NUMEL)

C------------------------- Internal

      INTEGER*4::I,I1,I2,IDIM,IND,J,K,L,LANI,LTY,NNUD

      REAL*8::AREAL    ,AREALN
     &       ,DFML     ,DISPIJ   ,EFTHRMCON,POR1     ,S
     &       ,WATVAVG

      REAL*8::DAUX(6),SUM(6),VV(6),WATV(LMXNDL)


C------------------------- FIRST EXECUTABLE STATEMENT.

C------------------------- Initializes ATRA matrix to zero

      ATRA = 0D0


C------------------------- Cross over elements

      DO L=1,NUMEL

          LTY = LTYPE(L)
          NNUD = LNNDEL(L)
          AREAL = AREA(L)
          AREALN = AREAL/NNUD   
          IDIM = LDIM(L)
          LANI = IDIM*(IDIM+1)/2    ! Maximum anisotropy of the current element
          WATV(:) = 0D0
          WATVAVG = 0D0

C------------------------- Computes the averge of water content.

          IF (IOVRWC.LE.1) THEN !Elementwise

              WATV(:) = WATVOL(1,L,2)
              WATVAVG = WATVOL(1,L,2)

          ELSE !Nodewise


              DO K=1,NNUD

                  WATV(K) = WATVOL(K,L,2)
                  WATVAVG = WATVAVG + WATV(K)

              END DO !K=1,NNUD

              WATVAVG = WATVAVG/NNUD

          END IF !IOVRWC.LE.1

C------------------------- Computes the hydrodynamic dispersion tensor

C------------------------- Moecular Diffusion.

C------------------------- If solving solute transport...

          IF (ITPTVAR.EQ.0) THEN

C------------------------- ... molecular diffusion is computed
C------------------------- acording to water content and molecular
C------------------------- diffusion coeficient.

              DFML = WATVAVG*PAREL(L,11)

C------------------------- If solving energy transport...

          ELSE

C------------------------- ...equivalent term to molecular diffusion
C------------------------- is computed according to effective thermal
C------------------------- conductivity

              POR1 = 1D0 - PAREL(L,12)

              EFTHRMCON = WATVAVG*WTHERMCON + POR1*PAREL(L,11)
              IF (IODENS.EQ.1) THEN

                  DFML = EFTHRMCON/(DENSITY(L)*WSPECHEAT)
              ELSE
C------------------------- If density is constant, water density appears
C------------------------- in this term (see TRANSDENS Guia rapida, 3.7).
                  DFML = EFTHRMCON/(DENSREF*WSPECHEAT)

              END IF!IODENS.EQ.1

          END IF !ITPTVAR.EQ.1

C------------------------- Hydrodynamic Dispersion.

C------------------------- PAREL(L,9)  <=> Long. dispersion coef.
C------------------------- PAREL(L,10) <=> Transv. dispersion coef.

C------------------------- Diagonal terms.

          DO I=1,IDIM 

              DAUX(I) = DFML + PAREL(L,9)*QXYZ(I,L)

              DO J=0,IDIM-2

                  DAUX(I) = DAUX(I)+PAREL(L,10)*QXYZ( MOD(I+J,IDIM)+1,L)

              END DO !J=0,IDIM-2

          END DO !I=1,IDIM

C------------------------- Non-diagonal terms.

          DO I=IDIM+1,LANI

              DAUX(I) = (PAREL(L,9)-PAREL(L,10))*QXYZ(I,L)

          END DO !I=IDIM+1,LANI


C------------------------- Contribution of velocity at current element 

          DO I=1,NNUD

              S = 0D0

              DO J=1,IDIM

                  S = S + VD(J,L)*GRDFF(J,I,L)

              END DO !J=1,IDIM

              VV(I) = S*AREALN


C------------------------- SUM is not related to contribution of
C------------------------- velocity, but current loop is used to
C------------------------- initialize it to zero for optimization
C------------------------- purposes.

              SUM(I) = 0D0
                        
          END DO !I=1,NNUD


C------------------------- Update of ATRA matrix

          IND = 0

          DO I=1,NNUD-1

              DO J=I+1,NNUD

                  DISPIJ = 0D0

                  DO K=1,LANI

                      DISPIJ = DISPIJ + DAUX(K)*BIBI(IND+K,L)

                  END DO !K=1,LANI

C------------------------- Diagonal contribution

                  SUM(I) = SUM(I) + DISPIJ
                  SUM(J) = SUM(J) + DISPIJ


C------------------------- Non diagonal terms 

                  I1 = (I-1)*NNUD + J
                  I2 = (J-1)*NNUD + I

                  ATRA(L,I1) = ATRA(L,I1) + DISPIJ + VV(J)
                  ATRA(L,I2) = ATRA(L,I2) + DISPIJ + VV(I)

                  IND = IND + LANI

              END DO ! J=I+1,NNUD

C------------------------- Diagonal terms

              I1 = (I-1)*NNUD + I

              ATRA(L,I1) = ATRA(L,I1) - SUM(I) + VV(I)

          END DO ! I=1,NNUD-1

C------------------------- Last diagonal term

          I2 = (NNUD-1)*NNUD + NNUD

          ATRA(L,I2) = ATRA(L,I2) - SUM(NNUD) + VV(NNUD)

C------------------------- Finally, contribution of density is added

          IF (IODENS.GT.0) THEN

              ATRA(L,:) = DENSITY(L)*ATRA(L,:)

          END IF !IODENS.GT.0

C------------------------- Derivatives of dispersivity part of ATRA matrix 
C------------------------- with respect to velocity components (qx,qy,qz).
C------------------------- Only if flow parameters are estimated with tpt inv.
C------------------------- or transport is solved with Newton-Raphson method
C------------------------- if there is dispersivity, and velocity is not null.

          IF (IOINV.EQ.3 .OR. LINMET(2,2).EQ.2.OR.LINMET(3,2).EQ.2) THEN

              IF  (NZONE_PAR(7).NE.0 .AND. XNORVD(L).GE.1.D-25) THEN

                  CALL DER_ATRA_VD
     &            (DFML  ,PAREL(L,9),PAREL(L,10),IDIMBB   ,IDIMDQ
     &            ,IDIMQ ,IODIM    ,L        ,LANI     ,IDIM     ,NNUD
     &            ,NUMEL ,XNORVD(L),BIBI     ,DAT_VD   ,DAUX     ,QXYZ
     &            ,VD)

              END IF
          END IF !IOINV.EQ.3 .OR. LINMET(2,2).EQ.2.OR.LINMET(3,2).EQ.2)


C------------------------- Derivatives of ATRA w.r.t. state variables.
C------------------------- Only if Newton's method is used.
       
          IF (LINMET(2,2).EQ.2.OR.LINMET(3,2).EQ.2) THEN

              IOCALCDEVF = 0

              IF(LINMET(3,2).EQ.2) THEN
                  IOCALCDEVF=1  !Calculate dtradflu also
              END IF

              CALL COMP_DER_ATRA
     &           (AREA     ,ATRA     ,BETAC    ,CAUX1    ,DAT_VD
     &           ,DENSITY  ,DENSREF  ,DTRADFLU ,DTRADTRA ,DVDH     ,DVDC
     &           ,DWDH     ,EPSFLU   ,EPSTRA   ,IDIMDQ   ,IOCALCDEVF
     &           ,IODENS   ,IODIM    ,IOVRWC   ,ITPTVAR  ,GRDFF
     &           ,KXX      ,L        ,LDIM     ,LMXNDL   ,LNNDEL
     &           ,LTYPE    ,NPPEL   ,NUMEL     ,NUMNP    ,PAREL
     &           ,THETAT   ,WATVOL  ,WSPECHEAT ,WTHERMCON)


          END IF !LINMET(2,2).EQ.2.OR.LINMET(3,2).EQ.2

      END DO ! L=1,NUMEL

      END SUBROUTINE COMP_ATRA
