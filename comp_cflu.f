      SUBROUTINE COMP_CFLU
     &          (AREA     ,BETAC    ,CAUX1    ,CAUX2    ,CFLU
     &          ,CREF     ,DENSITY  ,DENSREF  ,DFLUDFLU ,DFLUDTRA
     &          ,DWDH     ,IDIMCFLU ,IOCALCDEV,IOVRWC   ,KXX
     &          ,LMXNDL   ,LNNDEL   ,NUMEL    ,NUMNP    ,THETAT
     &          ,WATVOL)
*********************************************************************
*
* PURPOSE
* to calculate the matrix CFLU. CFLU is a symetric matrix
* The derivative of Cij of an element is equal to beta_w*Dij of the 
* element. 
*
* WATVOL stored elementwise or nodewise
*
*  NEW VARIABLES
*
*  IOCALCDEV: indivcates of derivatives of some sort need to be calculated
*
********************************************************************************

      IMPLICIT NONE


* EXTERNAL VARIABLES: SCALARS
      REAL*8::BETAC,CREF,DENSREF,THETAT

      INTEGER*4::NUMEL,IOCALCDEV,IOVRWC,LMXNDL,NUMNP,IDIMCFLU


* EXTERNAL VARIABLES: ARRAYS

      REAL*8::AREA(NUMEL)
     &       ,CFLU(NUMEL, IDIMCFLU)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL),CAUX1(NUMNP),CAUX2(NUMNP)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL),DENSITY(NUMEL)
     &       ,DWDH(MAX(1,(IOVRWC-1)*2*LMXNDL),NUMEL)

      INTEGER*4::LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)

      IF (IOVRWC.LE.1) THEN

          CALL COMP_CFLU_L
     &        (AREA     ,BETAC    ,CAUX2    ,CFLU     ,DENSITY
     &        ,DFLUDFLU ,DFLUDTRA ,DWDH     ,IDIMCFLU ,IOCALCDEV
     &        ,IOVRWC   ,KXX      ,LMXNDL   ,LNNDEL   ,NUMEL
     &        ,NUMNP    ,THETAT   ,WATVOL)

      ELSE

          CALL COMP_CFLU_I
     &        (AREA     ,BETAC    ,CAUX1    ,CAUX2    ,CFLU
     &        ,CREF     ,DENSREF  ,DFLUDFLU ,DFLUDTRA ,DWDH
     &        ,IDIMCFLU ,IOCALCDEV,IOVRWC   ,KXX      ,LMXNDL
     &        ,LNNDEL   ,NUMEL    ,NUMNP    ,THETAT   ,WATVOL)

      END IF !IOVRWC.LE.1

      END SUBROUTINE COMP_CFLU

*****************************************************************************
*****************************************************************************

      SUBROUTINE COMP_CFLU_L
     &          (AREA     ,BETAC    ,CAUX2    ,CFLU     ,DENSITY
     &          ,DFLUDFLU ,DFLUDTRA ,DWDH     ,IDIMCFLU ,IOCALCDEV
     &          ,IOVRWC   ,KXX      ,LMXNDL   ,LNNDEL   ,NUMEL
     &          ,NUMNP    ,THETAT,WATVOL)

*********************************************************************
*
* PURPOSE
* to calculate the matrix CFLU. CFLU is a symetric matrix
* The derivative of Cij of an element is equal to beta_w*Dij of the 
* element. 
*
* WATVOL stored elementwise.
*
*  NEW VARIABLES
*
*  IOCALCDEV: indivcates of derivatives of some sort need to be calculated
*
********************************************************************************

      IMPLICIT NONE


C-------------------- External

      REAL*8::BETAC,THETAT

      INTEGER*4::NUMEL,IOCALCDEV,IOVRWC,LMXNDL,NUMNP,IDIMCFLU


      REAL*8::AREA(NUMEL)
     &       ,CFLU(NUMEL, IDIMCFLU)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL),CAUX2(NUMNP)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL),DENSITY(NUMEL)
     &       ,DWDH(1,NUMEL)

      INTEGER*4::LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)

C-------------------- Internal

      REAL*8::AREAL,CFLUI,DCFLUDC,DCFLUDH,WATV
      INTEGER*4::I,IJ_POS,J,L,NNUD,INODE

C-------------------- First executable statement

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          AREAL = AREA(L) 

          WATV = WATVOL(1,L,2)

c---------------------------------------------------------------
c______________ Deshabilitado temporalmente hasta que se entre en
c______________ STATE_VARIABLE_INIT en flujo.
c
c         !calculating mean head in element
c         DO INODE =1,NNUD
c            NODEI=KXX(INODE,L)
c            H=(1D0/NNUD)*HCALIT(NODEI)
c            HAN=(1D0/NNUD)*HCALAN(NODEI)  
c         ENDDO 
c         WATV=WATV+PAREL(L,7)*(H-HAN)*EPSFLU
c---------------------------------------------------------------

          CFLUI = WATV*BETAC*AREAL*DENSITY(L)/NNUD

          CFLU(L,1:IDIMCFLU) = CFLUI


C-------------------- Calculate derivatives.

          IF (IOCALCDEV.EQ.1) THEN

              DO I=1,NNUD

                  INODE = KXX(I,L)

                  DCFLUDH = BETAC*DENSITY(L)*DWDH(1,L)*CAUX2(INODE)
     &                      *AREAL/NNUD

                  DCFLUDC = THETAT*CFLUI*BETAC*CAUX2(INODE)/NNUD

                  DO J=1,NNUD

                      IJ_POS = (I-1)*NNUD + J

                      DFLUDTRA(L,IJ_POS) = DFLUDTRA(L,IJ_POS) + DCFLUDC
                      DFLUDFLU(L,IJ_POS) = DFLUDFLU(L,IJ_POS) + DCFLUDH

                  END DO !J=!,NNUD

              END DO !I=1,NNUD


          END IF !IOCALCDEV.EQ.1
                                                   
      END DO !L=1,NUMEL


      END SUBROUTINE COMP_CFLU_L

*****************************************************************************
*****************************************************************************

      SUBROUTINE COMP_CFLU_I
     &          (AREA     ,BETAC    ,CAUX1    ,CAUX2    ,CFLU
     &          ,CREF     ,DENSREF  ,DFLUDFLU ,DFLUDTRA ,DWDH
     &          ,IDIMCFLU ,IOCALCDEV,IOVRWC   ,KXX      ,LMXNDL
     &          ,LNNDEL   ,NUMEL    ,NUMNP    ,THETAT   ,WATVOL)

*********************************************************************
*
* PURPOSE
* to calculate the matrix CFLU. CFLU is a symetric matrix
* The derivative of Cij of an element is equal to beta_w*Dij of the 
* element. 
*
* WATVOL stored nodewise.
*
*  NEW VARIABLES
*
*  IOCALCDEV: indivcates of derivatives of some sort need to be calculated
*
********************************************************************************

      IMPLICIT NONE

C--------------------------- External

      REAL*8::BETAC    ,CREF     ,DENSREF  ,THETAT

      REAL*8::DENS

      INTEGER*4::IDIMCFLU ,IOCALCDEV,IOVRWC   ,LMXNDL
     &          ,NUMEL    ,NUMNP

      REAL*8::AREA(NUMEL)
     &       ,CFLU(NUMEL, IDIMCFLU)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL),CAUX1(NUMNP),CAUX2(NUMNP)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DWDH(2*LMXNDL,NUMEL)

      INTEGER*4::LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)

C---------------------------Internal

      REAL*8::AREAL    ,AREALN   ,CNODE    ,DCFLUDC  ,DCFLUDH
     &       ,DENSNODE ,FACTOR

      INTEGER*4::I        ,IJ_POS   ,IPOS     ,IPOS1    ,J        ,L
     &          ,NNUD     ,INODE


      REAL*8::WATV(LMXNDL)


      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          AREAL = AREA(L)
          AREALN = AREAL/NNUD

          WATV(1:NNUD) = WATVOL(1:NNUD,L,2)

c---------------------------------------------------------------
c______________ Deshabilitado temporalmente hasta que se entre en
c______________ STATE_VARIABLE_INIT en flujo.
c
c         WATV=WATV+PAREL(L,7)*(H-HAN)*EPSFLU
c---------------------------------------------------------------

          DO I=1,NNUD

              INODE = KXX(I,L)
              CNODE = CAUX1(INODE)
              DENSNODE = DENS(DENSREF,BETAC,CNODE,CREF)
              FACTOR = BETAC*DENSNODE*AREALN

              CFLU(L,I) = WATV(I)*FACTOR

C-------------------- Calculate derivatives. 

              IF (IOCALCDEV.EQ.1) THEN

                  INODE = KXX(I,L)
                  IPOS = 2*(I-1) + 1
                  IPOS1 = IPOS + 1


                  DO J=1,NNUD

                      IF (I.EQ.J) THEN

                          DCFLUDH = FACTOR*DWDH(IPOS,L)*CAUX2(INODE)
                          DCFLUDC = THETAT*CFLU(L,I)*BETAC*CAUX2(INODE)

                      ELSE

                          DCFLUDH = FACTOR*DWDH(IPOS1,L)*CAUX2(INODE)
                          DCFLUDC = 0D0

                      END IF ! I.EQ.J

                      IJ_POS = (I-1)*NNUD + J

                      DFLUDFLU(L,IJ_POS) = DFLUDFLU(L,IJ_POS) + DCFLUDH

                      DFLUDTRA(L,IJ_POS) = DFLUDTRA(L,IJ_POS) + DCFLUDC

                  END DO !J=!,NNUD

              END IF !IOCALCDEV.EQ.1

          END DO !I=1,NNUD

      END DO !L=1,NUMEL

      END SUBROUTINE COMP_CFLU_I
