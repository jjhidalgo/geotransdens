      SUBROUTINE COMP_AFLU
     &(IDIMAFLU   ,IDIMBB     ,IOCALCDEVF ,IOCALCDEVT ,LMXNDL      
     &,NPAREL     ,NPPEL      ,NUMEL      ,NUMNP      ,NZTRA
     &,AFLU       ,BETAC      ,BIBI       ,DENSITY    ,DFLUDTRA
     &,DFLUDFLU   ,DPARELDH   ,DPARELDC   ,EPSFLU     ,EPSTRA
     &,HAUX1      ,ISOZ       ,KXX        ,LDIM       ,LNNDEL
     &,LXPAREL    ,PAREL      ,THETAT)
  
***************************************************************************
*
*  PURPOSE
*  To calculate the matrix AFLU and its derivatives w.r.t head and conc. if
*  necesary.
*
* 
* 
*  DESCRIPTION
*  
*     Step 1:Initializes a few vectors and matricesInitialize some arrays
*            and vectors
*
*     Step 2: Computation of AFLU and its derivatives. Following are some
*             hints about the algorithm.
*      
*             Hydraulic conductivity has to be computed takin into account
*             the degree of anisotroy.
*
*             AFLU is stored elementwise. Due to its simetry only the
*             upper part of the matrix has to be stored. Diagonal
*             coefficients are equal to the sum of the row with opposite
*             sign.
*
*             DFLUDFLU does not depends on the node with respect the
*             derivation is made then only one value per node is stored.
*             The same applies to DFLUDTRA.
*
*             AFLUi,j = Integral [ rho_average * GradNi K gradNj ] =
*                     = rho_avg * Integral [GradNi K gradNj] =
*                     = rho_avg * BIBI*TRACT
*
*             DFLUDFLU(L,I) = SUMj[dAFLUi,j/dH *Hj] =
*             = SUMj[ Integral [ rho_average * GradNi dK/dH gradNj ] *Hj ] =
*             = rho_avg * BIBI*DPARELDH*Hj
*
*             DFLUDFLU(L,I) = SUMj[dAFLUi,j/dC *Hj] =
*             = SUMj[ Integral [ rho_average * GradNi dK/dC gradNj ] *Hj ] +
*             + (Beta_c/N)*SUMj [Ai,j * Hj]
*             = rho_avg * (BIBI*DPARELDC + (BETAC/NNUD)*BIBI*TRACT)*Hj
*
*             Since AFLU is part of DFLUDLTRA, the loop takes advantage of
*             the computaions of AFLU to add the cotribution to DFLUDTRA.
*             This is made through the variable GNKGN, which stands 
*             for Grad_Ni�Kij�Grad_Nj.
*
*             The product Matrix*StateVariable take advantage of the
*             symetric storage. This loop can also be found in 
*             PROD_MAT_VEC for matrix type number 6.
*
*         
*
*  NEW  VARIABLES
*
*     IODENS       indicates wether density is constant or variable
*     IOCALCDEVF   flag indicating wether DFLUDFLU should be calculated
*     IOCALCDEVT   flag indicating wether DFLUDTRA should be calculated
*
****************************************************************************

      IMPLICIT NONE

C-------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4 LMXNDL,NUMEL,NUMNP,IDIMAFLU,NPPEL,IDIMBB 
     &         ,NZTRA,NPAREL,IOCALCDEVF ,IOCALCDEVT
      
      REAL*8::EPSFLU,EPSTRA,THETAT

C-------------------- EXTERNAL VARIABLES: ARRAYS

      INTEGER*4  LDIM(NUMEL),LNNDEL(NUMEL),ISOZ(NZTRA)
     &           ,KXX(LMXNDL,NUMEL),LXPAREL(NUMEL,NPAREL)

      REAL*8::AFLU(NUMEL,IDIMAFLU),DFLUDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DENSITY(NUMEL)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL),HAUX1(NUMNP)
     &       ,PAREL(NUMEL,NPPEL),BIBI(IDIMBB,NUMEL)
     &       ,DPARELDH(NPPEL,NUMEL)
     &       ,DPARELDC(NPPEL,NUMEL) 


C-------------------- INTERNAL VARIABLES: SCALARS

      INTEGER*4 NZONE,ISZ,ISMAX,NNUD,ISOTTR,L,M,IS
     &         ,I,J,K,IJ_POS,IK_POS,JK_POS

      REAL*8 BETAC,H1,H2,GNKGN,SUMGNKGN,DERWRTH,DERWRTW
                      

C-------------------- INTERNAL VARIABLES: ARRAYS

      REAL*8 TRACT(9)

C-------------------- Step 1:Initializes a few vectors and matrices
      
      AFLU = 0D0

      IF(IOCALCDEVF .NE.0) DFLUDFLU = 0D0
      IF(IOCALCDEVT .NE.0) DFLUDTRA = 0D0


C-------------------- Step 2: Computation of AFLU and its derivatives.

C--------------------------- For each element...

      DO L=1, NUMEL

C--------------------------- ...its properties are set: hid. cond. zone,
C--------------------------- anisotropy degree, number of hid. cond.
C--------------------------- components and number of nodes.

          NZONE=LXPAREL(L,1)
          ISZ=ISOZ(NZONE)            
          ISMAX=MAX(ISZ,LDIM(L))
          NNUD=LNNDEL(L)
          ISOTTR=(LDIM(L)*(LDIM(L)+1))/2
                
C--------------------------- Then the value of hid. cond. components and
C--------------------------- their derivatives w.r.t. stae variables
C--------------------------- (if needed) are computed according to the
C--------------------------- anisotropy degree.

          DO IS=1,ISMAX
              TRACT(IS)=PAREL(L,IS)
              IF (IS.GT.1. .AND. (ISZ.LT.2.OR.IS.NE.3) ) THEN      
                  IF (IOCALCDEVF.EQ.1) DPARELDH(IS,L)=DPARELDH(IS-1,L)
                  IF (IOCALCDEVT.EQ.1) DPARELDC(IS,L)=DPARELDC(IS-1,L)
              END IF 
          END DO     ! IS=1,ISMAX (Loop over T components)

C--------------------------- Then AFLU and its derivatives are computed
C--------------------------- according to the description in the
C--------------------------- subroutine header.

          M=0        ! Pointer to BIBI array


          DO I=1,NNUD-1

C--------------------------- Node 'I' head stored

              IF (MAX(IOCALCDEVT,IOCALCDEVF).EQ.1) THEN
                  H1 = HAUX1(KXX(I,L))
              ELSE
                  H1 = 0D0
              END IF
                

              DO J=I+1,NNUD

C--------------------------- Node 'J' head stored

                  IF (MAX(IOCALCDEVT,IOCALCDEVF) .EQ. 1) THEN
                      H2 = HAUX1(KXX(J,L))
                  ELSE
                      H2 = 0D0
                  END IF
       
C--------------------------- Cumulative variables initialized

                  SUMGNKGN = 0D0
                  DERWRTH = 0D0
                  DERWRTW = 0D0

C--------------------------- Loop over anisotropy components

                  DO IS=1,ISMAX

C--------------------------- Contribution to AFLU.
C--------------------------- GNKGN stands for Grad_Ni�Kij�Grad_Nj.

                      GNKGN = TRACT(IS)*BIBI(M+IS,L)*DENSITY(L)

                      SUMGNKGN = SUMGNKGN + GNKGN

C--------------------------- Contribution to the derivative w.r.t. hk+1


                      IF (IOCALCDEVF.EQ.1) THEN
                          DERWRTH = DERWRTH
     &                    +EPSFLU*DPARELDH(IS,L)*BIBI(M+IS,L)*DENSITY(L)
                      END IF

C--------------------------- Contribution to the derivative w.r.t. ck+1.

                      IF (IOCALCDEVT.EQ.1) THEN
                          DERWRTW = DERWRTW
     &                   +(EPSTRA*DPARELDC(IS,L)*BIBI(M+IS,L)*DENSITY(L)
     &                            + THETAT*BETAC*GNKGN/NNUD)
                      END IF

                  END DO ! IS=1,ISMAX

C--------------------------- Derivatives.

                  DO K=1,NNUD

                      IK_POS = (I - 1)*NNUD + K
                      JK_POS = (J - 1)*NNUD + K

C--------------------------- Derivatives w.r.t. flow

                      IF (IOCALCDEVF.EQ.1) THEN

                          DFLUDFLU(L,IK_POS) = DFLUDFLU(L,IK_POS)
     &                                       + DERWRTH*(H2-H1)

                          DFLUDFLU(L,JK_POS) = DFLUDFLU(L,JK_POS)
     &                                       + DERWRTH*(H1-H2)
                      END IF !IOCALCDEVF.EQ.1

C--------------------------- Derivatives w.r.t. transport

                      IF (IOCALCDEVT.EQ.1) THEN 
                          DFLUDTRA(L,IK_POS) = DFLUDTRA(L,IK_POS) 
     &                                       + DERWRTW*(H2-H1) 
              
                          DFLUDTRA(L,JK_POS) = DFLUDTRA(L,JK_POS)
     &                                       + DERWRTW*(H1-H2)

                      END IF !IOCALCDEVT.EQ.1


                  END DO !K=1,NNUD
                  
                  
C--------------------------- A_IJ (Symetric storage without diagonal)

                 IJ_POS =(I - 1)*NNUD + J - I*(I+1)/2
                 AFLU(L,IJ_POS) = AFLU(L,IJ_POS) + SUMGNKGN
                     
                 M=M+ISOTTR  !Update pointer to BIBI

              END DO !J=I+1,NNUD

          END DO  !I=1,NNUD-1               
                  
      END DO !L=1,NUMEL

      END SUBROUTINE COMP_AFLU
