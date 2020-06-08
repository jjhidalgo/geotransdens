      SUBROUTINE DER_LEAK
     &          (BETAC    ,CAUX1    ,CFPARNP  ,CREF     ,DENSREF
     &          ,DERH     ,DTIM     ,EPSFLU   ,FNT      ,HAUX1
     &          ,HCALAN   ,HCALIT   ,IDIMDERH ,IDIMFNT  ,IFLAGS
     &          ,INALF    ,INCHP    ,INDSSTR  ,INEW     ,INORPAR
     &          ,INTI     ,IOCAP    ,IODENS   ,IOFLLI   ,IOFMLF
     &          ,IVPAR    ,IXPARNP  ,KXX      ,LMXNDL   ,NFLAGS
     &          ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR
     &          ,NINT     ,NODE     ,NPAR     ,NPARNP   ,NPPNP
     &          ,NPZON    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &          ,NZPAR    ,PARACD   ,PARC     ,PARNP    ,IDIMWGT
     &          ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)

*******************************************************************************
*
* PURPOSE
*
*   Computes flow RHS of derivatives with respect to leakage for a given node.
*
* DESCRIPTION
*
*   Computes flow RHS of derivatives with respect to leakage  for a given node.
*
* HISTORY
*
*     JHG     12-2005     First coding
*
*******************************************************************************

      IMPLICIT NONE
     
C------------------------- External

      INTEGER*4::IDIMDERH ,IDIMFNT  ,INALF    ,INCHP    ,INDSSTR
     &          ,INEW     ,INTI     ,IOCAP    ,IODENS   ,IOFLLI
     &          ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT     ,NODE
     &          ,NPAR     ,NPARNP   ,NPPNP    ,NPZON    ,NTYPAR
     &          ,NUMEL    ,IOFMLF   ,NUMNP    ,NZPAR


      REAL*8::BETAC    ,CREF      ,DENSREF  ,DTIM     ,EPSFLU

      REAL*8::DENS

      INTEGER*4::IFLAGS(NFLAGS)        ,INORPAR(NTYPAR)  ,IVPAR(NZPAR)
     &          ,IXPARNP(NUMNP,NPARNP) ,KXX(LMXNDL,NUMEL),NFNLPAR(NZPAR)
     &          ,NFNLPRG(8,NFNL)       ,NFNLTIP(NFNL)    ,NFTPAR(NZPAR)
     &          ,NZONE_PAR(NTYPAR)    ,IDIMWGT  ,IPNT_PAR

      REAL*8::CAUX1(NUMNP)              ,CFPARNP(NUMNP,NPARNP)
     &       ,DERH(NUMNP,NPAR,IDIMDERH) ,FNT(IDIMFNT,NINT)
     &       ,HCALAN(NUMNP)             ,HCALIT(NUMNP)
     &       ,HAUX1(NUMNP)              ,PARNP(NUMNP,NPPNP)
     &       ,PARC(NZPAR)  ,PARACD(3,NFNL)    ,WGT_PAR(IDIMWGT) 
       
C------------------------- Internal

      INTEGER*4::IC     ,IP     ,IZH    ,IZON   ,JH     ,JJ
     &          ,NCNF   ,NPTOT

      REAL*8::ALF      ,ALFAX    ,CAUDLEAK ,CEXT     ,CHP
     &       ,DENSNODE ,HEAD

      INTEGER*4::INDEX(12)  ,IPOS(NPAR)

      REAL*8::CFPARAM(12)   ,DERIV(NPAR)   ,XPARAM(8)

C------------------------- First executable statement


C-------------------------- PARNP(NODE,4) ==> conc. ext

      ALFAX = PARNP(NODE,3)
      HEAD = PARNP(NODE,1)
      CAUDLEAK = ALFAX*(HEAD - HAUX1(NODE))
	CEXT = PARNP(NODE,4)

      IF (IODENS.EQ.1) THEN

          IF (CAUDLEAK.GT.0D0) THEN

              DENSNODE = DENS(DENSREF,BETAC,CEXT,CREF)

          ELSE IF (CAUDLEAK.LT.0D0) THEN

	        DENSNODE = DENS(DENSREF,BETAC,CAUX1(NODE),CREF)

          ELSE

	        DENSNODE = DENS(DENSREF,BETAC,PARNP(NODE,4),CREF)
     &                 + DENS(DENSREF,BETAC,CAUX1(NODE),CREF)

              DENSNODE = DENSNODE/2D0

	    END IF !CAUDLEAK.GT.0D0

      ELSE

          DENSNODE = 1D0

      END IF !IODENS.EQ.1


C------------------------- Derivative with respect to leakage coefficient

      IZON = IXPARNP(NODE,INALF)
      JJ = INORPAR(11) + IZON
      IP = IVPAR(JJ)

      IF (IP.NE.0.OR.IOFLLI.NE.0) THEN


C------------------------- Stores in CHP the prescribed head value at node I

          CHP = PARNP(NODE,1)


          IF (IOFLLI.NE.0) THEN
              NCNF=NFNLPAR(JJ)
          ELSE
              NCNF=0
          END IF !IOFLLI.NE.0

          INDEX(1) = JJ


          CFPARAM(1) = CFPARNP(NODE,INALF)

          CALL DER_PARAM
     & (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     & ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,NODE
     & ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     & ,NFTPAR   ,NINT     ,1        ,NPTOT    ,NPZON    ,NUMEL
     & ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     & ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     & ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ; ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

          IF (NPTOT.GT.0) THEN 

              DO IC=1,NPTOT
                  DERH(NODE,IPOS(IC),INEW) = DERH(NODE,IPOS(IC),INEW)
     &                            + DERIV(IC)*(CHP-HAUX1(NODE))*DENSNODE

              END DO !IC=1,NPTOT

          END IF !NPTOT.GT.0

      END IF !IP.NE.0.OR.IOFLLI.NE.0 <==> Estimated?

C------------------------- Derivative with respect to leakage prescribed head
           
	IZH = IXPARNP(NODE,INCHP)
	JH = INORPAR(9) + IZH
      IP = IVPAR(JH)

      IF (IP.NE.0.OR.IOFLLI.NE.0) THEN

C------------------------- Stores in ALF the leakage coefficient value at node I

          ALF = PARNP(NODE,3)

          IF (IOFLLI.NE.0) THEN
              NCNF=NFNLPAR(JH)
          ELSE
              NCNF=0
          END IF !IOFLLI.NE.0

          INDEX(1) = JH

          CFPARAM(1)=CFPARNP(NODE,INCHP)

          CALL DER_PARAM
     & (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     & ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,NODE
     & ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     & ,NFTPAR   ,NINT     ,1        ,NPTOT    ,1        ,NUMEL
     & ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JH),DERIV ,FNT
     & ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     & ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ; ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

          IF (NPTOT.GT.0) THEN

              DO IC=1,NPTOT

                  DERH(NODE,IPOS(IC),INEW) = DERH(NODE,IPOS(IC),INEW)
     &                      + DERIV(IC)*ALF*DENSNODE

              END DO !IC=1,NPTOT

          END IF !NPTOT.GT.0

      END IF !IP.NE.0.OR.IOFLLI.NE.0 <==> Estimated?

      END SUBROUTINE DER_LEAK
