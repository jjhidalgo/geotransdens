      SUBROUTINE DER_WATVOL_H
     &          (DPARELDH ,DWDH     ,EPSTRA   ,HCALAN   ,HCALIT
     &          ,IOVRWC   ,KXX      ,LMXNDL   ,LNNDEL   ,NPPEL
     &          ,NUMEL    ,NUMNP    ,PAREL)

     
********************************************************************************
*
* PURPOSE
*
*  Manages the computation of the derivatives of WATVOL evaluated at k+th
*  w.r.t. head k + 1, DIVIDED BY THETAT.
*
*
* DESCRIPTION
*
*  Manages the computation of the derivatives of WATVOL evaluated at k+alfa
*  w.r.t. head at k + 1, DIVIDED BY THETHAT
*
*  IOVRWC = 0
*
*    WATVOL is constant. Therefore its derivatives are null.
*
*  IOVRWC = 1
*
*    WATVOL is variable elementwise.
*
*        W(k+th) = W(k) + th*PAREL(L,7)*dltHavg
*
*    The derivative is computed as
*    (hk1m =h^(k+1)_m; hkm =h^(k)_m; hem = h^(k+e)_m)
*
*
*      dW(k+th)/dhk1m = dW(k)/dh1m + th*dPAREL/dhk1m *dltHavg + th*PAREL(L,7)*DdltHavg/Dhk1m =
*                  =      0    + th*e*dPAREL/dhem *dltHavg + th*PAREL(L,7)*1/(N) =
*                  = th*(e*dPAREL/dhm *dltHavg + PAREL(L,7)/N)
*
*      dW(k+th)/dhkm = dW(k)/dhkm + th*dPAREL/dhkm *dltHavg + th*PAREL(L,7)*DdltHavg/Dhkm =
*                  =   dW(k)/dhkm + th*(1-e)*dPAREL/dhem *dltHavg - th*PAREL(L,7)*1/(N) =
*                  = th*((1-e)*dPAREL/dhm *dltHavg - PAREL(L,7)/N)
*
*
*    The derivative does not depends on m, therefore, it is stored in a vector of
*    size NUMEL*2, where position 1 stores dW(k+th)/dhk1 or dW(k+th)/dhk and position
*    2 stores dW(k)/dhk.
*    THE STORED VALUE IS DIVIDED BY TH IN ORDER TO ALLOW THE DIRECT AND REVERSE INTEGRATION IN TIME
*
*  IOVRWC = 2
*
*    WATVOL is variable nodewise.
*
*        W(k+e) = W(k) + e*PAREL(L,7)*dltHavg_i
*
*    The derivative is computed as (hm = h^(k+e)_m)
*
*      dW(k+e)/dhk1m = dW(k)/dhk1m + th*e*dPAREL/dhem *dltHavg_i + th*PAREL(L,7)*DdltHavg_i/Dhk1m =
*                  =      0    + th*e*dPAREL/dhem *dltHavg + th*PAREL(L,7)*Dim =
*                  = th*(e*dPAREL/dhem *dltHavg + PAREL(L,7)*Dim)
*
*      dW(k+e)/dhkm = dW(k)/dhkm + th*e*dPAREL/dhem *dltHavg_i + th*PAREL(L,7)*DdltHavg_i/Dhkm =
*                  =  dW(k)/dhkm + th*(1-e)*dPAREL/dhem *dltHavg - th*PAREL(L,7)*Dim =
*                  = th*((1-e)*dPAREL/dhem *dltHavg - PAREL(L,7)*Dim)
*
*      Where Dim is equal to 1 if i=m, and zero otherwise.
*
*    There are two different values for the derivative one if m=i and m/=i 
*    (two values for each node). Therefore, it is stored in a vector of
*    size NUMEL*2*N 
*    THE STORED VALUE IS DIVIDED BY TH IN ORDER TO ALLOW THE DIRECT AND REVERSE INTEGRATION IN TIME
*
*    In this vector, the derivative of node I w. r. t. m=I is in position
*    4*(I-1)+1 and the derivatives w. r. t. m/=I are in position 4*(I-1)+2
*    Positions 4*(I-1)+3 and 4*(I-1)+4 store dW(k)/dhkm (m=I and m/=I respectively).
*
*  STORAGE
*
*    DWDH(MAX(1,(IOVRWC-1)*2*LMXNDL),NUMEL)
*
*      Elementwise (IOVRWC = 1): DWDH(1,NUMEL)
*                   DWDH(1,NUMEL) -> @W/hk+1
*
*      Nodewise (IOVRWC = 2): DWDH(2*LMXNDL,NUMEL)
*                   DWDH(1,NUMEL) -> @Wi/hik+1 
*                   DWDH(2,NUMEL) -> @Wi/hjk+1
*                Positions are not 1 to 2, but
*                IPOS = 2*(I-1) + (0--1)    I=1,NNUD
*
*
* HISTORY
*
*     JHG      7-2005     First coding
********************************************************************************

      IMPLICIT NONE


C------------------------- External
      INTEGER*4::IOVRWC,LMXNDL,NPPEL,NUMEL,NUMNP
      
      REAL*8::EPSTRA

      REAL*8::DWDH(MAX(1,(IOVRWC-1)*2*LMXNDL),NUMEL),HCALAN(NUMNP)
     &       ,HCALIT(NUMNP),PAREL(NUMEL,NPPEL),DPARELDH(NPPEL,NUMEL)

      INTEGER*4::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)


C------------------------- Internal

      REAL*8::HAVG,HI
      INTEGER*4::I,INODE,IPOS,L,NNUD



C------------------------- First executable statement.

      DO L=1,NUMEL

          NNUD = LNNDEL(L)

          SELECT CASE(IOVRWC)

              !CASE(0) !Constant, elementwise

                  !DWDH(1,L) = PAREL(L,7)

              CASE(1) !Variable, elementwise

                  HAVG = 0D0

                  DO I=1,NNUD

                      INODE = KXX(I,L)

                      HAVG = HAVG + (HCALIT(INODE)-HCALAN(INODE))

                  END DO !I=1,NNUD

                  HAVG = HAVG/NNUD

                  DWDH(1,L) = EPSTRA*DPARELDH(7,L)*HAVG+PAREL(L,7)/NNUD

              CASE(2)

                  DO I=1,NNUD !Variable, nodewise

                      INODE = KXX(I,L)

                      HI = (HCALIT(INODE)-HCALAN(INODE))

                      IPOS = 2*(I-1) + 1

                      DWDH(IPOS,L) = EPSTRA*DPARELDH(7,L)*HI+PAREL(L,7)

                      DWDH(IPOS+1,L) = EPSTRA*DPARELDH(7,L)*HI

                  END DO !I=1,NNUD


          END SELECT !IOVRWC
      END DO !L=1,NUMEL

      END SUBROUTINE DER_WATVOL_H
