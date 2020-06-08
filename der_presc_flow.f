      SUBROUTINE DER_PRESC_FLOW
     &          (BETAC    ,CAUX1    ,CFPARNP  ,CREF     ,DENSREF
     &         ,DERH      ,DTIM     ,EPSFLU   ,FNT      ,HCALAN
     &         ,HCALIT    ,IDIMDERH ,IDIMFNT  ,IFLAGS   ,INDSSTR
     &         ,INEW      ,INORPAR  ,INQQP    ,INTI     ,IOCAP
     &         ,IODENS    ,IOFLLI   ,IOFMLF   ,IVPAR    ,IXPARNP
     &         ,KXX       ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLPAR
     &         ,NFNLPRG   ,NFNLTIP  ,NFTPAR   ,NINT     ,NODE
     &         ,NPAR      ,NPARNP   ,NPPNP    ,NPZON    ,NTYPAR
     &         ,NUMEL     ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD
     &         ,PARC      ,PARNP    ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR
     &         ,IPOS     ,DERIV)

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

      INTEGER*4::IDIMDERH ,IDIMFNT  ,INDSSTR  ,INEW     ,INQQP
     &          ,INTI     ,IOCAP    ,IODENS   ,IOFLLI   ,IOFMLF
     &          ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT     ,NODE
     &          ,NPAR     ,NPARNP   ,NPZON    ,NTYPAR   ,NUMEL
     &          ,NPPNP   ,NUMNP     ,NZPAR

      REAL*8::BETAC,CREF,DENSREF,DTIM     ,EPSFLU

	REAL*8::DENS

      INTEGER*4::IFLAGS(NFLAGS)        ,INORPAR(NTYPAR)  ,IVPAR(NZPAR)
     &          ,IXPARNP(NUMNP,NPARNP) ,KXX(LMXNDL,NUMEL),NFNLPAR(NZPAR)
     &          ,NFNLPRG(8,NFNL)       ,NFNLTIP(NFNL)    ,NFTPAR(NZPAR)
     &          ,NZONE_PAR(NTYPAR)     ,IDIMWGT  ,IPNT_PAR

      REAL*8::CAUX1(NUMNP)              ,CFPARNP(NUMNP,NPARNP)
     &       ,DERH(NUMNP,NPAR,IDIMDERH) ,FNT(IDIMFNT,NINT)
     &       ,HCALAN(NUMNP)             ,HCALIT(NUMNP)
     &       ,PARC(NZPAR)               ,PARACD(3,NFNL)
     &       ,PARNP(NUMNP,NPPNP)        ,WGT_PAR(IDIMWGT)
       
C------------------------- Internal

      INTEGER*4::IC     ,IP     ,IZON   ,JJ     ,NCNF   ,NPTOT

      REAL*8::DENSNODE

      INTEGER*4::INDEX(12)  ,IPOS(NPAR)

      REAL*8::CFPARAM(12)   ,DERIV(NPAR)   ,XPARAM(8)

C------------------------- First executable statement



C-------------------------- PARNP(NODE,4) ==> conc. ext
C-------------------------- PARNP(NODE,2) ==> caudal

      IF(IODENS.EQ.1) THEN

          IF (PARNP(NODE,2).GT.0D0) THEN

              DENSNODE = DENS(DENSREF,BETAC,PARNP(NODE,4),CREF)

          ELSE IF (PARNP(NODE,2).LT.0D0) THEN

	        DENSNODE = DENS(DENSREF,BETAC,CAUX1(NODE),CREF)

          ELSE

              DENSNODE = DENS(DENSREF,BETAC,PARNP(NODE,4),CREF)
     &                 + DENS(DENSREF,BETAC,CAUX1(NODE),CREF)

              DENSNODE = DENSNODE/2D0

	    END IF !PARNP(I,2).GT.0D0

      ELSE
      
          DENSNODE = 1D0

      END IF !IODENS.EQ.1

      IZON = IXPARNP(NODE,INQQP)
      JJ = INORPAR(10)+IZON
      IP = IVPAR(JJ)

      IF (IP.NE.0.OR.IOFLLI.NE.0) THEN

          IF (IOFLLI.NE.0) THEN
              NCNF=NFNLPAR(JJ)
          ELSE
              NCNF=0
          END IF !IOFLLI.NE.0

          INDEX(1) = JJ

          CFPARAM(1) = CFPARNP(NODE,INQQP)

          CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,NODE
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,1        ,NPTOT    ,NPZON    ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

          IF (NPTOT.GT.0) THEN 

              DO IC=1,NPTOT
                  DERH(NODE,IPOS(IC),INEW) = DERH(NODE,IPOS(IC),INEW)
     &                                     + DERIV(IC)*DENSNODE
              END DO !NPTOT.GT.0

          END IF !NPTOT.GT.0

      END IF !IP.NE.0.OR.IOFLLI.NE.0 <==> Estimation?

	END SUBROUTINE DER_PRESC_FLOW
