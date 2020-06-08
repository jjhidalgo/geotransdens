      SUBROUTINE CORRECT_DDTRA_RHS
     &          (AREA     ,BETAC    ,CAUX2    ,CCALAN    ,CCALIT
     &          ,CREF     ,DENSREF  ,DERC     ,DERH      ,DPARELDH
     &          ,DTRA     ,EPSTRA   ,HCALAN   ,HINI      ,IDIMDTRA
     &          ,INEWT    ,IODENS   ,IOFLLI   ,IOLD      ,IOLDT
     &          ,IOVRWC   ,ITPTVAR  ,KXX      ,LMXNDL    ,LNNDEL
     &          ,NPAR     ,NPPEL    ,NUMEL     ,NUMNP    ,PAREL
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
* The call in JAC_C to RHS_NOLI_IN adds to the inverse problem RHS a wrong
* term since the conversion from derivatives w. r. t. hk+1 to hk cannot
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

      INTEGER*4::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

      REAL*8::AREA(NUMEL)           ,CAUX2(NUMNP)
     &       ,CCALAN(NUMNP)         ,CCALIT(NUMNP)
     &       ,DERC(NUMNP,NPAR,2)    ,DERH(NUMNP,NPAR,2)
     &       ,DTRA(NUMEL,IDIMDTRA)  ,HCALAN(NUMNP)
     &       ,HINI(NUMNP)           ,DPARELDH(NPPEL,NUMEL)
     &       ,PAREL(NUMEL,NPPEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
    

C---------------------- First executable statment.


      IF (IOVRWC.LT.2) THEN

          CALL CORRECT_DDTRA_RHS_L
     &        (AREA     ,BETAC    ,CAUX2    ,CCALAN   ,CCALIT
     &        ,CREF     ,DENSREF  ,DERC     ,DERH     ,DPARELDH
     &        ,DTRA     ,EPSTRA   ,HCALAN   ,HINI     ,IDIMDTRA
     &        ,INEWT    ,IODENS   ,IOFLLI   ,IOLD     ,IOLDT
     &        ,IOVRWC   ,ITPTVAR  ,KXX      ,LMXNDL   ,LNNDEL
     &        ,NPAR     ,NPPEL    ,NUMEL    ,NUMNP    ,PAREL
     &        ,THETAT   ,WATVOL)

      ELSE IF (IOVRWC.EQ.2) THEN

          CALL CORRECT_DDTRA_RHS_I
     &        (AREA     ,BETAC    ,CAUX2    ,CCALAN   ,CCALIT
     &        ,CREF     ,DENSREF  ,DERC     ,DERH     ,DPARELDH
     &        ,DTRA     ,EPSTRA   ,HCALAN   ,HINI     ,IDIMDTRA
     &        ,INEWT    ,IODENS   ,IOFLLI   ,IOLD     ,IOLDT
     &        ,IOVRWC   ,ITPTVAR  ,KXX      ,LMXNDL   ,LNNDEL
     &        ,NPAR     ,NPPEL    ,NUMEL    ,NUMNP    ,PAREL
     &        ,THETAT   ,WATVOL)

      END IF !IOVRWC.LT.2 ...

      END SUBROUTINE CORRECT_DDTRA_RHS
