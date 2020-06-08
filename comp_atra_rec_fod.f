      SUBROUTINE COMP_ATRA_REC_FOD
     &          (ACTH     ,AREA     ,ATRA     ,BETAC    ,CAUX1
     &          ,CREF     ,DENSREF  ,DTRADFLU ,DTRADTRA ,DWDH
     &          ,IODENS   ,IONEWT   ,IOVRWC   ,ITPTVAR  ,LINMET
     &          ,LMXNDL   ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP
     &          ,KXX      ,LNNDEL   ,NZONE_PAR,PAREL    ,THETAT
     &          ,WATVOL)


********************************************************************************
*
* PURPOSE
*
*  Adds to ATRA matrix recharge and first order decay contributions.
*
*
*
*
******************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IODENS   ,IONEWT   ,IOVRWC   ,ITPTVAR  ,LMXNDL
     &          ,NPPEL    ,NTYPAR   ,NUMEL    ,NUMNP

      REAL*8::BETAC,CREF,DENSREF,THETAT

      REAL*8::DENS

      INTEGER*4::KXX(LMXNDL, NUMEL)  ,LINMET(3,2)
     &          ,LNNDEL(NUMEL)       ,NZONE_PAR(NTYPAR)

      REAL*8::ACTH(NUMEL)                             ,AREA(NUMEL)
     &       ,ATRA(NUMEL,LMXNDL*LMXNDL)               ,CAUX1(NUMNP)
     &       ,DTRADFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DTRADTRA(NUMEL,LMXNDL*LMXNDL)           ,DWDH(2,NUMEL)
     &       ,PAREL(NUMEL,NPPEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)

C------------------------- Internal

      INTEGER*4::I        ,II_DER_POS         ,IJ_DER_POS
     &          ,II_POS   ,IJ_POS   ,INODE    ,J        ,K
     &          ,L        ,NNUD

      REAL*8::AREAL    ,AREALN   ,CNODE    ,CREC     ,DENSNODE ,DENSREC
     &       ,DFODDH   ,DFODDW   ,FOD      ,FODNODE  ,RECNODE  ,WATVAVG
     &       ,RECHARGE ,RETARD

      REAL*8::WATV(LMXNDL)


C------------------------- FIRST EXECUTABLE STATEMENT.


      DENSNODE = 1D0
      DENSREC = 1D0

C------------------------- Cross over elements

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          AREAL = AREA(L)
          AREALN = AREAL/NNUD   
          FOD = 0D0
          RECHARGE = 0D0
          WATV(:) = 0D0
          WATVAVG = 0D0
          RETARD = PAREL(L,14)*ACTH(L)
          RECNODE = PAREL(L,8)*AREALN
          FODNODE = PAREL(L,13)*AREALN

          IF (RECNODE.GT.0) THEN

              CREC = PAREL(L,15)

          ELSE

              CREC = 0D0

          END IF !RECNODE.GT.0

C------------------------- Computes the averge of water content.

          IF (IOVRWC.LE.1) THEN !Elementwise

              WATV(:) = WATVOL(1,L,2)

          ELSE !Nodewise


              DO K=1,NNUD

                  WATV(K) = WATVOL(K,L,2)

              END DO !K=1,NNUD


          END IF !IOVRWC.LE.1


          DO I=1,NNUD

              INODE = KXX(I,L)
              CNODE = CAUX1(INODE)
              II_POS = (I-1)*NNUD + I

C------------------------- Areal recharge (Only Picard).

              IF (NZONE_PAR(3).NE.0 .AND. IONEWT.EQ.0) THEN
                
                  IF (IODENS.EQ.1) THEN

                      DENSREC =  DENS(DENSREF,BETAC,CREC,CREF)

                  END IF !IODENS.EQ.1

                  RECHARGE = RECNODE*DENSREC

              END IF !NZONE_PAR(3).NE.0 .AND. IONEWT.EQ.0


C------------------------- First order decay coefficient and retardation.
C------------------------- FOD*CRD*B*AREA/NNUD
C------------------------- Lambda * (WATVOL + alfa_s) * NjNi

              IF (NZONE_PAR(11).NE.0 .AND. ITPTVAR.NE.1) THEN

                  IF (IODENS.EQ.1) THEN

                      DENSNODE =  DENS(DENSREF,BETAC,CNODE,CREF)

                  END IF !IODENS.EQ.1

                  FOD = FODNODE*DENSNODE*(WATV(I) + RETARD)


C------------------------- Derivatives of FOD w.r.t. concentration.
C------------------------- Only if Newton's method is used.

                  IF (LINMET(2,2).EQ.2 .OR. LINMET(3,2).EQ.2) THEN

                      DFODDW = THETAT*BETAC*FOD*CAUX1(INODE)

                      DTRADTRA(L,II_POS) = DTRADTRA(L,II_POS) + DFODDW

                  END IF !LINMET(2,2).EQ.2 .OR. LINMET(3,2).EQ.2

C------------------------- Derivatives of FOD w.r.t. head.
C------------------------- Only coupled Newton's method.

                  IF (LINMET(3,2).EQ.2) THEN

                      IF (IOVRWC.LE.1) THEN

                          II_DER_POS = 1
                          IJ_DER_POS = 1

                      ELSE IF (IOVRWC.EQ.2) THEN

                          II_DER_POS = 2*(I-1) + 1
                          IJ_DER_POS = II_DER_POS + 1

                      END IF

                      DFODDH = FODNODE*DENSNODE*CAUX1(INODE)
     &                        *THETAT*DWDH(II_DER_POS,L)

                      DTRADFLU(L,II_POS) = DTRADFLU(L,II_POS) + DFODDH

                      DO J=1,NNUD

                          IF (J.NE.I) THEN

                              IJ_POS = (I-1)*NNUD + J

                              DFODDH = FODNODE*DENSNODE*CAUX1(INODE)
     &                        *DWDH(IJ_DER_POS,L)

                              DTRADFLU(L,IJ_POS) = DTRADFLU(L,IJ_POS)
     &                                            + DFODDH

                          END IF !J.NE.I
                      END DO !J=1,NNUD
                  END IF !LINMET(3,2).EQ.2
              END IF !NZONE_PAR(11).NE.0 .AND. ITPTVAR.NE.1

C------------------------- Update of ATRA matrix

              ATRA(L,II_POS) = ATRA(L,II_POS) + FOD + RECHARGE


          END DO ! I=1,NNUD

      END DO ! L=1,NUMEL

      END SUBROUTINE COMP_ATRA_REC_FOD
