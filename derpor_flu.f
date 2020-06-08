      SUBROUTINE DERPOR_FLU
     &          (ACTH     ,AREA     ,BETAC    ,CAUX1    ,CAUX2
     &          ,CCALAN   ,CCALIT   ,CFPAREL  ,CREF     ,DENSITY
     &          ,DENSREF  ,DERH     ,DTIM     ,EPSTRA   ,FNT
     &          ,IDIMDERH ,IDIMFNT  ,IFLAGS   ,INDSSTR  ,INEW
     &          ,INORPAR  ,INTI     ,IOCAP    ,IOFMLT   ,IOTRLI
     &          ,IOVRWC   ,IVPAR    ,KXX      ,LMXNDL   ,LNNDEL
     &          ,LXPAREL  ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &          ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR     ,NPAREL
     &          ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &          ,PARACD   ,PARC     ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR
     &          ,IPOS     ,DERIV)

*********************************************************************
*
* PURPOSE
*
*    Computes the derivatives of flow equation w. r. t. porosity
*
*********************************************************************
  
      IMPLICIT NONE


C--------------------  External

      INTEGER*4::IDIMDERH ,IDIMFNT  ,INDSSTR  ,INEW     ,INTI
     &          ,IOCAP    ,IOFMLT   ,IOTRLI   ,IOVRWC   ,LMXNDL
     &          ,NFLAGS   ,NFNL     ,NINT     ,NPAR     ,NPAREL
     7          ,NTYPAR   ,NUMEL    ,NUMNP    ,NZPAR


      REAL*8::BETAC    ,CREF     ,DENSREF  ,DTIM      ,EPSTRA

      REAL*8::DENS

      INTEGER*4::IFLAGS(NFLAGS)   ,INORPAR(NTYPAR)  ,IVPAR(NZPAR)
     &          ,KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)  ,LXPAREL(NUMEL,NPAREL)
     &          ,NFNLPAR(NZPAR)   ,NFNLPRG(8,NFNL)  ,NFNLTIP(NFNL)
     &          ,NFTPAR(NZPAR)    ,NZONE_PAR(NTYPAR) ,IDIMWGT  


      REAL*8::ACTH(NUMEL)               ,AREA(NUMEL)   ,CAUX1(NUMNP)
     &       ,CAUX2(NUMNP)              ,CCALAN(NUMNP) ,CCALIT(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)     ,DENSITY(NUMEL)
     &       ,DERH(NUMNP,NPAR,IDIMDERH) ,PARC(NZPAR)  ,FNT(IDIMFNT,NINT)
     &       ,PARACD(3,NFNL)            ,WGT_PAR(IDIMWGT)


C--------------------  Internal

      INTEGER*4::INDEX(12),IPOS(NPAR)    ,IPNT_PAR
      
      REAL*8::AREAL,CNODE,DENSNODE,FACTOR,VAUX

      INTEGER*4::I,IC,INODE,IP,JJ,L,NCNF,NNUD,NPTOT,NZ

      REAL*8::CFPARAM(12),DERIV(NPAR),XPARAM(8)

C--------------------  First executable statement.

      DO L=1,NUMEL

          NZ = LXPAREL(L,7)
          JJ = INORPAR(15) + NZ
          IP = IVPAR(JJ)


          IF (IP.GT.0 .OR. IOTRLI.NE.0) THEN

              NNUD = LNNDEL(L)

              IF (IOTRLI.NE.0) THEN
                  NCNF = NFNLPAR(JJ)
              ELSE
                  NCNF = 0
              END IF !IOTRLI.NE.0

              INDEX(1) = JJ

              CFPARAM(1) = CFPAREL(L,7)

              CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLT   ,IP       ,L
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,1        ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),CCALIT   ,CCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

              IF (NPTOT.GT.0) THEN

              
                  NNUD = LNNDEL(L)
                  AREAL = AREA(L)

                  FACTOR = ACTH(L)*BETAC*AREAL/NNUD

                  DO I=1,NNUD

                      INODE = KXX(I,L)

                      IF (IOVRWC.LT.2) THEN

                          DENSNODE = DENSITY(L)

                      ELSE

                          CNODE = CAUX1(INODE)
                      
                          DENSNODE = DENS(DENSREF,BETAC,CNODE,CREF)

                      END IF !IOVRWC.LT.2


                      DO IC=1,NPTOT

                          VAUX = FACTOR*DENSNODE*CAUX2(INODE)

                          DERH(INODE,IPOS(IC),INEW) = 
     &                                         DERH(INODE,IPOS(IC),INEW)
     &                                         - DERIV(IC)*VAUX
                      END DO !IC=1,NPAR
                  END DO !I=1,NNUD

              END IF !NPTOT.GT.0

          END IF !IP.GT.0 .OR. IOTRLI.NE.0

      END DO !L=1,NUMEL
  
      
      END SUBROUTINE DERPOR_FLU
