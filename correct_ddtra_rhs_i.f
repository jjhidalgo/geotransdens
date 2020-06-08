      SUBROUTINE CORRECT_DDTRA_RHS_I
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

      REAL*8::BETAC    ,CREF     ,DENSREF  ,EPSTRA   ,THETAT

      REAL*8::DENS

      INTEGER*4::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

      REAL*8::AREA(NUMEL)            ,CAUX2(NUMNP)
     &       ,CCALAN(NUMNP)          ,CCALIT(NUMNP)
     &       ,DERC(NUMNP,NPAR,2)     ,DERH(NUMNP,NPAR,2)
     &       ,DTRA(NUMEL,IDIMDTRA)   ,HCALAN(NUMNP)
     &       ,HINI(NUMNP)            ,DPARELDH(NPPEL,NUMEL)
     &       ,PAREL(NUMEL,NPPEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
    
C------------------------- Internal

      INTEGER::I      ,INODE  ,J      ,JNODE  ,L      ,NNUD

      REAL*8::AREAL    ,AREALN   ,CONC     ,DELTAH   ,DDTRADC  ,DDTRADH
     &       ,EPS1     ,EPSFACT  ,DER_DTRA ,DIFFEPS  ,RHO      ,THT1

      REAL*8::DER_DFLU(2),DERWTVDH(2)

C---------------------- First executable statment.

      EPSFACT = (1D0 - 2D0*EPSTRA)/EPSTRA
      EPS1 = 1D0 - EPSTRA
      THT1 = 1D0 - THETAT
	DIFFEPS = DABS(EPS1-EPSTRA)
 
      DO L=1,NUMEL


          NNUD = LNNDEL(L)
          AREAL = AREA(L)
          AREALN = AREAL/NNUD

          IF (IODENS.GT.0) THEN

              CONC = THT1*CCALIT(INODE) + THETAT*CCALAN(INODE)
              RHO = DENS(DENSREF,BETAC,CONC,CREF)
          ELSE

              RHO = 1D0

          END IF !IODENS.GT.0

C------------------------- Correction of derivatives w. r. t. concentration
C------------------------- (when density is nodewise there are no cross
C------------------------- derivatives).

          IF (IODENS.GT.0 .AND. DIFFEPS.GT.1D-5) THEN

              DO I=1,NNUD

                  INODE = KXX(I,L)
                  DER_DTRA = 0D0

                  IF (ITPTVAR.EQ.0) THEN ! Mass fraction

                      DER_DTRA = EPSFACT*BETAC*DTRA(L,I)

                  ELSE !Temperature

                      DER_DTRA = EPSFACT*BETAC*RHO* WATVOL(I,L,3)*AREALN

                  END IF !ITPTVAR.EQ.0

                  DDTRADC = DER_DTRA*CAUX2(INODE)
 
                  DERC(INODE,1:NPAR,INEWT) = DERC(INODE,1:NPAR,INEWT)
     &                                + DDTRADC*DERC(INODE,1:NPAR,IOLDT)

              END DO !I=1,NNUD

          END IF !IODENS.GT.0 .AND. DIFFEPS.GT.1D-5

C------------------------- Correction of derivatives w. r. t. head

          IF (DIFFEPS.GT.1D-5 .OR. IOFLLI.GT.0) THEN
          
              DER_DFLU(1:2) = 0D0

              DO I=1,NNUD

                  INODE = KXX(I,L)

C------------------------- Derivatives of watvol when stored nodewise
C------------------------- dependen on the node.

                  DELTAH = (HCALAN(INODE) - HINI(INODE))

                  DERWTVDH(1) = EPSFACT*PAREL(L,7)
     &                           - EPS1*DPARELDH(7,L)*DELTAH

                  DERWTVDH(2) = -EPS1*DPARELDH(7,L)*DELTAH

                  DER_DFLU(1)=RHO*EPSFACT*DERWTVDH(1)*AREALN

                  DER_DFLU(2)=RHO*EPSFACT*DERWTVDH(2)*AREALN

                  DO J=1,NNUD

                      JNODE = KXX(J,L)

                      IF (I.EQ.J) THEN

                          DDTRADH = DER_DFLU(1)*CAUX2(INODE)

                      ELSE

                          DDTRADH = DER_DFLU(2)*CAUX2(INODE)

                      END IF !I.EQ.J

                      DERC(INODE,1:NPAR,INEWT) =DERC(INODE,1:NPAR,INEWT)
     &                                 + DDTRADH*DERH(INODE,1:NPAR,IOLD)

                  END DO !J=1,NNUD

              END DO !I=1,NNUD

         END IF !DIFFEPS.GT.1D-5 .OR. IOFLLI.GT.0

      END DO !L=1,NUMEL

      END SUBROUTINE CORRECT_DDTRA_RHS_I
