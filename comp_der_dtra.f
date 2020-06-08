      SUBROUTINE COMP_DER_DTRA
     &          (AREA     ,BETAC    ,CAUX1    ,CAUX2    ,CREF
     &          ,DENSITY  ,DENSREF  ,DTRA     ,DTRADFLU ,DTRADTRA
     &          ,DWDH     ,IDIMDTRA ,IOVRWC   ,ITPTVAR  ,KXX
     &          ,LMXNDL   ,LNNDEL   ,NUMEL    ,NUMNP    ,THETAT
     &          ,WATVOL)

********************************************************************************
*
* PURPOSE
*
*  Manages the computation of the derivatives of DTRA matrix w.r.t. state variables.
*
*
* DESCRIPTION
*
*  Manages the computation of the derivatives of DTRA matrix w.r.t. state variables.
*  The derivative is stored elementwise.
*
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CAUX2                  Array containing diference of concentrations in two   
*                         consecutives times, related to time increment         
*  CNST                   Interpolation functions gradient for a given element  
*                         nodes.
*  DTRADFLU               Derivative of transport matrix w.r.t. fresh water head.
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element
*  PAREL
*
* INTERNAL VARIABLES: ARRAYS
*
*  DER_DTRA               Derivative of DTRA (storage matrix) w.r.t. fresh water head.
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMDTRA               Used to dimension array DTRA (second dimension).
*                             IDIMDTRA=LMXNDL if consistent scheme.
*                             IDIMDTRA=LMXNDL * (LMXNDL+1)/2 if lumped scheme.
*  IOCNST                 Scheme for mass storage term in transport problem
*                             IOCNST = 1 Lumped scheme.
*                             IOCNST = 0 Consistent scheme (diagonal matrix)
*  LMXNDL                 Maximum number of nodes per element
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  TINC                   Current time increment  
*                              
*
* INTERNAL VARIABLES: SCALARS
*
*
*  NNUD                   Number of nodes of the current element                
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY
*
*     JHG      8-2003     First coding
********************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER::IDIMDTRA ,IOVRWC   ,ITPTVAR  ,LMXNDL   ,NNUD     ,NUMEL
     &        ,NUMNP         

      REAL*8::BETAC,CREF,DENSREF,THETAT

      REAL*8::DENS

      INTEGER*4::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

      REAL*8::AREA(NUMEL),CAUX1(NUMNP),CAUX2(NUMNP),DENSITY(NUMEL)
     &       ,DTRA(NUMEL,IDIMDTRA)
     &       ,DTRADFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DTRADTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DWDH(MAX(1,(IOVRWC-1)*2*LMXNDL),NUMEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
    
C------------------------- Internal

      INTEGER::I      ,II_POS ,IJ_POS ,INODE  ,IPOS   ,IPOS1  ,J      ,L

      REAL*8::AREAL    ,AREALN   ,DENSNODE ,DDTRADC  ,DDTRADH
     &  ,FACTOR,DER_DTRA,WTV,THT1


      REAL*8::DER_DFLU(2)

C------------------------- First executable statement

      THT1 = 1D0 - THETAT

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          AREAL = AREA(L)
          AREALN = AREAL/NNUD

          DER_DTRA = 0D0
          DER_DFLU(1:2) = 0D0

          DO I=1,NNUD

                  INODE = KXX(I,L)

              SELECT CASE (IOVRWC)

          CASE (0,1) ! Constant, Elementwise

                  IF (ITPTVAR.EQ.0) THEN ! Mass fraction

                      DER_DTRA = THT1*BETAC*DTRA(L,I)/NNUD

                  ELSE !temperature


                      FACTOR = THT1*BETAC*DENSITY(L)*AREALN/NNUD

                      DER_DTRA =  WATVOL(1,L,3)*FACTOR

                  END IF !ITPTVAR.EQ.0


                  DER_DFLU(1:2) = DENSITY(L)*THT1*DWDH(1,L)*AREALN



              CASE (2)

                  IF (ITPTVAR.EQ.0) THEN ! Mass fraction

                      DER_DTRA = THT1*BETAC*DTRA(L,I)

                  ELSE !Temperature

                      WTV = WATVOL(I,L,3)
                      DENSNODE = DENS(DENSREF,BETAC,CAUX1(INODE),CREF)
                      FACTOR = THT1*BETAC*DENSNODE*AREALN

                      DER_DTRA = WTV*FACTOR

                  END IF !ITPTVAR.EQ.0

C------------------------- Derivatives of watvol when stored nodewise
C------------------------- are dependent of the node.

                      IPOS = 2*(I-1) + 1 !@Wi/@Hi
                      IPOS1 = IPOS + 1   !!@Wi/@Hj j.ne.i

                      DER_DFLU(1) = DENSNODE*THT1*DWDH(IPOS,L)*AREALN
                      DER_DFLU(2) = DENSNODE*THT1*DWDH(IPOS1,L)*AREALN

              END SELECT !IOVRWC


              II_POS = (I-1)*NNUD + I

              DDTRADC = DER_DTRA*CAUX2(INODE)
              DTRADTRA(L,II_POS) = DTRADTRA(L,II_POS) + DDTRADC

              DO J=1,NNUD

                  IF (I.EQ.J .OR. IOVRWC.NE.2) THEN

                      DDTRADH = DER_DFLU(1)*CAUX2(INODE)

                  ELSE

                      DDTRADH = DER_DFLU(2)*CAUX2(INODE)

                  END IF !I.EQ.J


                  IJ_POS = (I-1)*NNUD + J
                  DTRADFLU(L,IJ_POS) = DTRADFLU(L,IJ_POS) + DDTRADH

C------------------------- When density is node wise there are no cross
C------------------------- derivatives.

                  IF (IOVRWC.NE.2 .AND. J.NE.I) THEN

                      DTRADTRA(L,IJ_POS) = DTRADTRA(L,IJ_POS) + DDTRADC

                  END IF !IOVRWC.NE.2 .AND. J.NE.I

              END DO !J=1,NNUD

          END DO !I=1,NNUD
      
      END DO !L=1,NUMEL
      
      END SUBROUTINE COMP_DER_DTRA
