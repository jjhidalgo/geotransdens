      SUBROUTINE DER_PRESC_HEAD
     &          (CFPARNP  ,DERH     ,DTIM     ,EPSFLU   ,FNT
     &          ,HCALAN   ,HCALIT   ,IDIMDERH ,IDIMFNT  ,IFLAGS
     &          ,INCHP    ,INDSSTR  ,INEW     ,INORPAR  ,INTI
     &          ,IOCAP    ,IOFLLI   ,IOFMLF   ,IVPAR    ,IXPARNP
     &          ,KXX      ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLPAR
     &          ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT     ,NODE
     &          ,NPAR     ,NPARNP   ,NPZON    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &          ,THETAF   ,TINC     ,TINTERVOBS         ,IDIMWGT
     &          ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)

*******************************************************************************
*
* PURPOSE
*
*   Computes flow RHS of derivatives with respect to prescribed head for a given
*   node.
*
* DESCRIPTION
*
*   Computes flow RHS of derivatives with respect to prescribed head for a given
*   node.
*
* HISTORY
*
*     JHG     12-2005     First coding
*
*******************************************************************************

      IMPLICIT NONE
     
C------------------------- External

      INTEGER*4::IDIMDERH ,IDIMFNT  ,INCHP    ,INDSSTR  ,INEW
     &          ,INTI     ,IOCAP    ,IOFLLI   ,LMXNDL   ,NFLAGS
     &          ,NFNL     ,NINT     ,NODE     ,NPAR     ,NPARNP
     &          ,NPZON    ,NTYPAR   ,NUMEL    ,IOFMLF   ,NUMNP
     &          ,NZPAR    ,IDIMWGT  ,IPNT_PAR


      REAL*8::DTIM     ,EPSFLU   ,THETAF   ,TINC     ,TINTERVOBS


      INTEGER*4::IFLAGS(NFLAGS)        ,INORPAR(NTYPAR)  ,IVPAR(NZPAR)
     &          ,IXPARNP(NUMNP,NPARNP) ,KXX(LMXNDL,NUMEL),NFNLPAR(NZPAR)
     &          ,NFNLPRG(8,NFNL)       ,NFNLTIP(NFNL)    ,NFTPAR(NZPAR)
     &          ,NZONE_PAR(NTYPAR)

      REAL*8::CFPARNP(NUMNP,NPARNP)  ,DERH(NUMNP,NPAR,IDIMDERH)
     &       ,FNT(IDIMFNT,NINT)      ,HCALAN(NUMNP)
     &       ,HCALIT(NUMNP)          ,PARC(NZPAR)
     &       ,PARACD(3,NFNL)         ,WGT_PAR(IDIMWGT)
       
C------------------------- Internal

      INTEGER*4::IC     ,IP     ,IZON   ,JJ     ,NCNF   ,NPTOT

      REAL*8::DTIMAUX

      INTEGER*4::INDEX(12)  ,IPOS(NPAR)

      REAL*8::CFPARAM(12)   ,DERIV(NPAR)   ,XPARAM(8)

C------------------------- First executable statement




      IZON = IXPARNP(NODE,INCHP)
      JJ = INORPAR(9) + IZON
      IP = IVPAR(JJ)
        
C------------------------- Derivatives have to be computed if IP>0 or
C------------------------- even when the zonal parameter is not estimated, it
C------------------------- may depend on some generic parameters that may have
C------------------------- to be estimated

      IF (IP.GT.0 .OR. IOFLLI.NE.0) THEN

          IF (IOFLLI.NE.0) THEN 
              NCNF=NFNLPAR(JJ)
          ELSE
              NCNF=0
          END IF !IOFLLI.NE.0


          INDEX(1) = JJ

          CFPARAM(1) = CFPARNP(NODE,INCHP)

C------------------------- Derivative must be compute at time k+1, not at
C------------------------- k+theta. So, DTIM = 1D0.

          DTIMAUX = DTIM + (1D0-THETAF)*TINC/TINTERVOBS

          CALL DER_PARAM
     &    (CFPARAM  ,DTIMAUX  ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,NODE
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,1        ,NPTOT    ,NPZON    ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)



          IF (NPTOT.GT.0) THEN

              DO IC=1,NPTOT

                  DERH(NODE,IPOS(IC),INEW) = DERIV(IC)

              END DO !IC=1,NPTOT

          END IF !NPTOT.GT.0

      END IF !IP.GT.0 .OR. IOFLLI.NE.0 <==> Estimation?

	END SUBROUTINE DER_PRESC_HEAD
