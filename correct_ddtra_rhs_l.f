      SUBROUTINE CORRECT_DDTRA_RHS_L
     &          (AREA     ,BETAC    ,CAUX2    ,CCALAN   ,CCALIT
     &          ,CREF     ,DENSREF  ,DERC     ,DERH     ,DPARELDH
     &          ,DTRA     ,EPSTRA   ,HCALAN   ,HINI     ,IDIMDTRA
     &          ,INEWT    ,IODENS   ,IOFLLI   ,IOLD     ,IOLDT
     &          ,IOVRWC   ,ITPTVAR  ,KXX      ,LMXNDL   ,LNNDEL
     &          ,NPAR     ,NPPEL    ,NUMEL    ,NUMNP    ,PAREL
     &          ,THETAT   ,WATVOL)
***********************************************************************
* PURPOSE
*
* Corrects the contribution of DTRADTRA and DTRADFLU to the inverse
* problem right hand side due to the reverse integration scheme used in
* DTRA.
*
*
* DESCRIPTION
*
* In the direct problem DTRA is computed at K+1-EPS. DTRADFLU and
* DTRADTRA contain the derivatives of DTRA w. r. t. hk+1 and ck+1,
* respectively.
*
* The call in JAC_C to RHS_NOLI_IN adds to the inverseproble RHS a wrong
* term since the conversien from derivatives w. r. t. hk+1 to hk cannot
* be done with the factor (EPS-1)/EPS in th case of DTRA.
*
* This subroutine removes the wrong term added in JAC_C and adds the
* contribution of the derivatives of DTRA w. r. t. time k to the inverse
* problem right hand side.
*
* INTERNAL VARIABLES:
*
*
* EXTERNAL VARIABLES:
*
*
* HISTORY
*
*
***********************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER::IDIMDTRA ,INEWT    ,IODENS   ,IOFLLI   ,IOLD     ,IOLDT
     &        ,IOVRWC   ,ITPTVAR  ,LMXNDL   ,NPAR     ,NPPEL    ,NUMEL
     &        ,NUMNP         

      REAL*8::BETAC,CREF,DENSREF,EPSTRA,THETAT

      REAL*8::DENS

      INTEGER*4::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

      REAL*8::AREA(NUMEL)           ,CAUX2(NUMNP)
     &       ,CCALAN(NUMNP)         ,CCALIT(NUMNP)
     &       ,DERC(NUMNP,NPAR,2)    ,DERH(NUMNP,NPAR,2)
     &       ,DTRA(NUMEL,IDIMDTRA)
     &       ,HCALAN(NUMNP)         ,HINI(NUMNP)
     &       ,DPARELDH(NPPEL,NUMEL) ,PAREL(NUMEL,NPPEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
    
C------------------------- Internal

      INTEGER::I      ,INODE  ,J      ,JNODE  ,K      ,KNODE  ,L
     &        ,NNUD

      REAL*8::AREAL    ,AREALN   ,CAVG     ,DELTAHAVG,DDTRADC  ,DDTRADH
     &       ,DER_DFLU ,DERWTVDH ,DIFFEPS  ,EPS1     ,EPSFACT  ,FACTOR
     &       ,DER_DTRA ,RHO      ,THT1

C---------------------- First executable statment.

      EPSFACT = (1D0 - 2D0*EPSTRA)/EPSTRA
      EPS1 = 1D0 - EPSTRA
	DIFFEPS = DABS(EPS1-EPSTRA)
      THT1 = 1D0 - THETAT
 
      DO L=1,NUMEL


          DER_DTRA = 0D0

          NNUD = LNNDEL(L)
          AREAL = AREA(L)
          AREALN = AREAL/NNUD

          IF (IODENS.GT.0) THEN
              CAVG = 0D0
              DO K=1,NNUD

                  KNODE = KXX(K,L)

                  CAVG = CAVG + THT1*CCALIT(KNODE)+THETAT*CCALAN(KNODE)

              END DO !K=1,NNUD

              CAVG = CAVG/NNUD
              RHO = DENS(DENSREF,BETAC,CAVG,CREF)

          ELSE

              RHO = 1D0

          END IF !IODENS.GT.0

C------------------------- Correction of derivatives w. r. t. concentration

          IF (IODENS.GT.0 .AND. DIFFEPS.GT.1D-5) THEN

              IF (ITPTVAR.EQ.0) THEN ! Mass fraction

                  DER_DTRA = EPSFACT*BETAC*DTRA(L,1)/NNUD

              ELSE !temperature

                  FACTOR = BETAC*RHO*AREALN/NNUD

                  DER_DTRA = EPSFACT*WATVOL(1,L,3)*FACTOR

              END IF !ITPTVAR.EQ.0


              DO I=1,NNUD

                  INODE = KXX(I,L)
                  DDTRADC = DER_DTRA*CAUX2(INODE)
 
                  DERC(INODE,1:NPAR,INEWT) = DERC(INODE,1:NPAR,INEWT)
     &                                + DDTRADC*DERC(INODE,1:NPAR,IOLDT)

C------------------------- Cross derivatives.

                  DO J=1,NNUD

                      IF (J.NE.I) THEN

                          JNODE = KXX(J,L)

                          DERC(INODE,1:NPAR,INEWT) =
     &                                          DERC(INODE,1:NPAR,INEWT)
     &                                + DDTRADC*DERC(JNODE,1:NPAR,IOLDT)

	                END IF !J.NE.I

                  END DO !J=1,NNUD

              END DO !I=1,NNUD

          END IF !IODENS.GT.0 .AND. DIFFEPS.GT.1D-5

C------------------------- Correction of derivatives w. r. t. head

          IF (IOVRWC.GT.0 .AND. (DIFFEPS.GT.1D-5 .OR. IOFLLI.GT.0)) THEN
          
              DER_DFLU = 0D0

              DELTAHAVG = 0D0

              DO K=1,NNUD

                  KNODE = KXX(K,L)
                  DELTAHAVG = DELTAHAVG + HCALAN(KNODE) - HINI(KNODE)

              END DO !K=1,NNUD

              DELTAHAVG = DELTAHAVG/NNUD

              DERWTVDH = EPSFACT *PAREL(L,7)/NNUD
     &                  - EPS1*DPARELDH(7,L)*DELTAHAVG

              DER_DFLU = RHO*THT1*DERWTVDH*AREALN

              DO I=1,NNUD

                  INODE = KXX(I,L)
                  DDTRADH = DER_DFLU*CAUX2(INODE)
                  DERC(INODE,1:NPAR,INEWT) = DERC(INODE,1:NPAR,INEWT)
     &                                 + DDTRADH*DERH(INODE,1:NPAR,IOLD)

C------------------------- Cross derivatives.

                  DO J=1,NNUD

	                IF (J.NE.I) THEN

                          JNODE = KXX(J,L)

                          DERC(INODE,1:NPAR,INEWT) = 
     &                                          DERC(INODE,1:NPAR,INEWT)
     &                                 + DDTRADH*DERH(JNODE,1:NPAR,IOLD)
                      END IF !J.NE.I

                  END DO !J=1,NNUD

              END DO !I=1,NNUD

         END IF !IOVRWC.GT.0 .AND. (DIFFEPS.GT.1D-5 ...

      END DO !L=1,NUMEL

      END SUBROUTINE CORRECT_DDTRA_RHS_L
