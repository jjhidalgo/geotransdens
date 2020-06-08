      SUBROUTINE COMP_DER_ATRA
     &           (AREA     ,ATRA     ,BETAC    ,CAUX1    ,DAT_VD
     &           ,DENSITY  ,DENSREF  ,DTRADFLU ,DTRADTRA ,DVDH     ,DVDC
     &           ,DWDH     ,EPSFLU   ,EPSTRA   ,IDIMDQ   ,IOCALCDEVF
     &           ,IODENS   ,IODIM   ,IOVRWC    ,ITPTVAR  ,GRDFF
     &           ,KXX      ,L       ,LDIM      ,LMXNDL   ,LNNDEL
     &           ,LTYPE    ,NPPEL   ,NUMEL     ,NUMNP    ,PAREL
     &           ,THETAT   ,WATVOL   ,WSPECHEAT,WTHERMCON)

********************************************************************************
*
* PURPOSE
*
*      Computes derivatives of ATRA matrix w. r. t. state variables.
*
* DESCRIPTION
*
*      Computes 
*
* EXTERNAL VARIABLES: ARRAYS
*
*                        
*
********************************************************************************

      IMPLICIT NONE


C------------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4::IDIMDQ   ,IOCALCDEVF         ,IODENS   ,IODIM
     &          ,IOVRWC   ,ITPTVAR  ,LMXNDL   ,NPPEL    ,NUMEL
     &          ,NUMNP

      REAL*8::BETAC    ,DENSREF,EPSFLU   ,EPSTRA   ,THETAT   ,WSPECHEAT
     &       ,WTHERMCON

C------------------------- EXTERNAL ARRAYS

      INTEGER*4::LDIM(NUMEL)   ,KXX(LMXNDL,NUMEL) ,LNNDEL(NUMEL)
     &          ,LTYPE(NUMEL)
      REAL*8::AREA(NUMEL)   ,CAUX1(NUMNP),DAT_VD(IODIM,IDIMDQ,NUMEL)
     &       ,DENSITY(NUMEL),DTRADFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DTRADTRA(NUMEL,LMXNDL*LMXNDL),DVDH(LMXNDL,IODIM,NUMEL)
     &       ,DVDC(LMXNDL,IODIM,NUMEL)  ,PAREL(NUMEL,NPPEL)
     &       ,GRDFF(IODIM,LMXNDL,NUMEL) ,ATRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DWDH(MAX(1,(IOVRWC-1)*2*LMXNDL),NUMEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)


C------------------------- INTERNAL SCALARS

      INTEGER*4::I          ,IDIM       ,II_DER_POS ,II_POS
     &           ,IJ_DER_POS,IJ_POS     ,IM_POS     ,INODE      ,IPOS
     &           ,J         ,JM_POS     ,JNODE      ,K          ,L
     &           ,LMXNDLSQR ,LTY        ,M          ,NDIM
     &           ,NNUD      ,NODEJ

      REAL*8::AMJ_C       ,AMJ_H    ,AREAL    ,AREALN   ,DENSL
     &       ,DERWATVDH   ,EFTHERMCON         ,FACTOR   ,POR1
     &       ,RHOSPHEATINV,SUMH     ,SUMW     ,TOTH     ,TOTW
     &       ,WATVAVG

C------------------------- Internal Arrays

      REAL*8::VEC(LMXNDL)



      LTY = LTYPE(L)
      NNUD = LNNDEL(L)
      IDIM = LDIM(L)
      AREAL = AREA(L)
      AREALN = AREAL/NNUD
      SUMH = 0D0
      SUMW = 0D0

      IF (IODENS.GT.0) THEN

          DENSL = DENSITY(L)

      ELSE

          IF (ITPTVAR.EQ.1) THEN

              DENSL = DENSREF

          ELSE

              DENSL =1D0

          END IF !ITPTVAR.EQ.1

      END IF!IODENS.GT.0

      LMXNDLSQR = LMXNDL*LMXNDL

      DTRADTRA(L,1:LMXNDLSQR) = 0D0

      IF (IOCALCDEVF.EQ.1) THEN

          DTRADFLU(L,1:LMXNDLSQR) = 0D0

      END IF

C------------------------- Contribution of velocity. 


      DO M=1,NNUD               !node to which we derive
        
        AMJ_C = 0D0
        AMJ_H = 0D0
        
         DO NDIM = 1,IDIM       !Loop over dimensions

            DO J=1,NNUD         !Loop over caux1
               
               JNODE = KXX(J,L)

               AMJ_C = AMJ_C +  DVDC(M,NDIM,L)*GRDFF(NDIM,J,L)
     &        *CAUX1(JNODE)*AREALN

               AMJ_H = AMJ_H +  DVDH(M,NDIM,L)*GRDFF(NDIM,J,L)
     &        *CAUX1(JNODE)*AREALN
                 
              
            END DO !J=1,NNUD

         END DO !NDIM=1,IDIM

         DO I=1,NNUD

            IPOS =  (I-1)*NNUD + M
            DTRADTRA(L,IPOS) = DTRADTRA(L,IPOS) + EPSTRA*AMJ_C

            IF(IOCALCDEVF.EQ.1) THEN
                DTRADFLU(L,IPOS) = DTRADFLU(L,IPOS) + EPSFLU*AMJ_H
            END IF

         END DO !I=1,NNUD

      ENDDO !M=1,NNUD

C------------------------- Contribution of dispersivity tensor.


      DO I=1,NNUD-1

         INODE = KXX(I,L)

         DO J=I+1,NNUD

            JNODE = KXX(J,L)
 
            DO M=1,NNUD

               TOTH = 0D0
               TOTW = 0D0

               IJ_POS = (I-1)*NNUD + J - I*(I+1)/2  
               
               DO NDIM=1,IDIM 
                  TOTH = TOTH + DAT_VD(NDIM,IJ_POS,L)*DVDH(M,NDIM,L)
                  TOTW = TOTW + DAT_VD(NDIM,IJ_POS,L)*DVDC(M,NDIM,L)
               END DO !NDIM=1,IDIM

               IM_POS = (I-1)*NNUD + M
               JM_POS = (J-1)*NNUD + M

               DTRADTRA(L,IM_POS)=DTRADTRA(L,IM_POS) 
     &         +  TOTW*EPSTRA*(CAUX1(JNODE)-CAUX1(INODE))

               DTRADTRA(L,JM_POS)=DTRADTRA(L,JM_POS)
     &         +  TOTW*EPSTRA*(CAUX1(INODE)-CAUX1(JNODE))

               IF (IOCALCDEVF.EQ.1) THEN
                  DTRADFLU(L,IM_POS)=DTRADFLU(L,IM_POS) 
     &           +  TOTH*EPSFLU*(CAUX1(JNODE)-CAUX1(INODE))

                  DTRADFLU(L,JM_POS)=DTRADFLU(L,JM_POS)
     &           +  TOTH*EPSFLU*(CAUX1(INODE)-CAUX1(JNODE))

               END IF !IOCALCDEVF.EQ.1

            END DO !M=1,NNUD
        END DO !J=I+1,NNUD
      END DO !I=1,NNUD-1


C------------------------- Contribution of molecular difusion.
C------------------------- PAREL(L,7)  --> Storage
C------------------------- PAREL(L,11) --> Diff. coef. / solid thermal cond.
C------------------------- WATVOL stored either elementwise or nodewise.

      IF (ITPTVAR.EQ.1) THEN

C------------------------- Computes the averge of water content.

          IF (IOVRWC.LE.1) THEN !Elementwise

              WATVAVG = WATVOL(1,L,2)

          ELSE !Nodewise

              WATVAVG = 0D0

              DO K=1,NNUD

                  WATVAVG = WATVAVG + WATVOL(K,L,2)

              END DO !K=1,NNUD

              WATVAVG = WATVAVG/NNUD

          END IF !IOVRWC.LE.1

          POR1 = 1D0 - PAREL(L,12)
          EFTHERMCON = WATVAVG*WTHERMCON + POR1*PAREL(L,11)
          RHOSPHEATINV = 1D0/(DENSL*WSPECHEAT)

      END IF !ITPTVAR.EQ.1

      DO I=1,NNUD

          INODE = KXX(I,L)
          II_POS = (I-1)*NNUD + I

C------------------------- Since molecular diffusion term is computed with
C------------------------- the average of WATVOL /when stored nodewise)
C------------------------- the derivative of that average has to be computed.

          IF (IOVRWC.LE.1) THEN

              DERWATVDH = DWDH(1,L)
              II_DER_POS = 1
              IJ_DER_POS = 1

          ELSE IF (IOVRWC.EQ.2) THEN

              II_DER_POS = 2*(I-1) + 1
              IJ_DER_POS = II_DER_POS + 1
              DERWATVDH = (DWDH(II_DER_POS,L) + 
     &                     (NNUD-1)*DWDH(IJ_DER_POS,L))/NNUD

          END IF !IOVRWC.EQ.1,2

          DERWATVDH = THETAT*DERWATVDH

C------------------------- If solving solute transport...

          IF (ITPTVAR.LT.1) THEN


              IF (IOCALCDEVF.EQ.1) THEN

                  DTRADFLU(L,II_POS) = DTRADFLU(L,II_POS)
     &            + CAUX1(INODE)*DERWATVDH*PAREL(L,11)*AREALN

              END IF !IOCALCDEVF.EQ.1

C------------------------- If solving energy transport...
          ELSE

              IF (IOCALCDEVF.EQ.1) THEN

                  DTRADFLU(L,II_POS) = DTRADFLU(L,II_POS)
     &            + CAUX1(INODE)*DERWATVDH*WTHERMCON*AREALN

              END IF !IOCALCDEVF.EQ.1

                  DTRADTRA(L,II_POS) = DTRADTRA(L,II_POS)
     &            - CAUX1(INODE)*BETAC*EFTHERMCON*RHOSPHEATINV*AREALN


          END IF ! ITPTVAR.EQ.1


          DO J=1,NNUD

              IF (J.NE.I) THEN

                  IJ_POS = (I-1)*NNUD + J

C------------------------- Contribution of Molecular Diffusion

                  IF (ITPTVAR.EQ.0) THEN

                      DTRADFLU(L,IJ_POS) = DTRADFLU(L,IJ_POS) +
     &                   PAREL(L,11)*DERWATVDH*CAUX1(INODE)*AREALN

                  ELSE

                      IF (IOCALCDEVF.EQ.1) THEN

                          DTRADFLU(L,IJ_POS) = DTRADFLU(L,IJ_POS)
     &                    + CAUX1(INODE)*DERWATVDH*WTHERMCON*AREALN

                      END IF !IOCALCDEVF.EQ.1

                          DTRADTRA(L,IJ_POS) = DTRADTRA(L,IJ_POS)
     &        - CAUX1(INODE)*THETAT*BETAC*EFTHERMCON*RHOSPHEATINV*AREALN

                  END IF !ITPTVAR.EQ.0


              END IF !J.NE.I

          END DO !J=1,NNUD

      END DO !I=1,NNUD

C------------------------- All the terms must be multiplicated by the density

      DTRADFLU(L,1:LMXNDLSQR) = DENSL*DTRADFLU(L,1:LMXNDLSQR)
      DTRADTRA(L,1:LMXNDLSQR) = DENSL*DTRADTRA(L,1:LMXNDLSQR)


C------------------------- Derivatives w.r.t. concentration.
C------------------------- BETAC/NNUD * ATRA *CAUX! is added to DTRADTRA
C------------------------- (Only if density is not constant).

      IF (IODENS.GT.0) THEN

          VEC = 0D0
          FACTOR = THETAT*BETAC/NNUD
          DO I=1,NNUD
              DO J=1,NNUD
                  NODEJ = KXX(J,L)
                  IJ_POS = (I-1)*NNUD + J
                  VEC(I) = VEC(I)+ATRA(L,IJ_POS)*CAUX1(NODEJ)
              END DO !J=1,NNUD
          END DO !I=1,NNUD

          DO I=1,NNUD
              DO J=1,NNUD
                  IJ_POS = (I-1)*NNUD + J
                  DTRADTRA(L,IJ_POS) = DTRADTRA(L,IJ_POS)+VEC(I)*FACTOR
              END DO !J=1,NNUD
          END DO !I=1,NNUD

      END IF !IODENS.GT.0

      END SUBROUTINE COMP_DER_ATRA
