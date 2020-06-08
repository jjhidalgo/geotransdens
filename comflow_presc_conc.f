      SUBROUTINE COMFLOW_PRESC_CONC
     &          (AFLU     ,AREA     ,ARRC     ,ATRA     ,BETAC
     &          ,BUOYANCY ,CAUDAL   ,CAUX1    ,CAUX2    ,CFLU
     &          ,COECEL   ,COORD    ,CREF     ,DFLU     ,DENSITY
     &          ,DENSREF  ,DTRA     ,FLUX     ,GP_COORD ,GRADLOC
     &          ,HAUX1    ,HAUX2    ,IAD_S    ,IADN_S   ,IBCOD
     &          ,IBTCO    ,IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU
     &          ,IDIMDTRA ,IODENS   ,IODIM    ,IORECATRA,ISOLFL
     &          ,ISOLTR   ,ISOZ     ,ITYPAFLU ,ITYPATRA ,ITYPCFLU
     &          ,ITYPDFLU ,ITYPDTRA ,KXX      ,LDIM     ,LMXNDL
     &          ,LNNDEL   ,LTYPE    ,LXARR    ,LXPAREL  ,MAXNB
     &          ,MAXPG    ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL
     &          ,NUMEL    ,NUMNP    ,NZARR    ,NZTRA    ,PAREL
     &          ,PARNP    ,POINTWEIGHT        ,THETAT)

********************************************************************************
*
*     Computes nodal fluxes at prescribed concentration nodes
*
*     M = [(Div q) + DFLU*DeltaH + CFLU*DeltaC + ATRA]*Ck+theta + DTRA*DeltaC +
*         - DESNREC*Recharge*CREC - DESNEXT*Q*CEXT
*
*
*  ISOLFL                 If 0, no flow has been solved at current time.
*                         If 1, steady flow has been solved at current time.
*                         If 2, transient flow has been solved at current time
*  ISOLTR                 If 0, no transport has been solved at current time.
*                         If 1, steady transport has been solved at current time
*                         If 2, trans. transport has been solved at current time
*
********************************************************************************

      IMPLICIT NONE


C------------------------- External

      INTEGER*4::IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU ,IDIMDTRA
     &          ,IODENS   ,IODIM    ,IORECATRA,ISOLFL   ,ISOLTR
     &          ,ITYPAFLU ,ITYPATRA ,ITYPCFLU ,ITYPDFLU ,ITYPDTRA
     &          ,LMXNDL   ,MAXNB    ,MAXPG    ,NPAREL   ,NPPEL
     &          ,NPPNP    ,NTYPEL   ,NUMEL    ,NUMNP    ,NZARR
     &          ,NZTRA

      REAL*8::BETAC    ,CREF     ,DENSREF  ,THETAT

      REAL*8::DENS

      INTEGER*4::IAD_S(MAXNB,NUMNP)  ,IADN_S(NUMNP)
     &          ,IBCOD(NUMNP)        ,IBTCO(NUMNP)
     &          ,ISOZ(NZTRA)         ,KXX(LMXNDL,NUMEL)
     &          ,LDIM(NUMEL)         ,LNNDEL(NUMEL)
     &          ,LTYPE(NUMEL)        ,LXARR(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL)

      REAL*8::AFLU(NUMEL,IDIMAFLU)          ,ATRA(NUMEL,IDIMATRA)
     &       ,AREA(NUMEL)                   ,ARRC(NUMEL)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL)  ,CAUDAL(NUMNP)
     &       ,CAUX1(NUMNP)
     &       ,CAUX2(NUMNP)                  ,CFLU(NUMEL,IDIMCFLU)
     &       ,COECEL(NUMEL)                 ,COORD(NUMNP,3)
     &       ,DENSITY(NUMEL)                ,DTRA(NUMEL,IDIMDTRA)
     &       ,FLUX(NUMNP)                   ,GP_COORD(6,8,IODIM)
     &       ,GRADLOC(IODIM,LMXNDL,MAXPG)   ,HAUX1(NUMNP)
     &       ,PAREL(NUMEL,NPPEL)            ,POINTWEIGHT(MAXPG,NTYPEL)
     &       ,dflu(numel,idimdflu),haux2(numnp),parnp(NUMNP,NPPNP)


C------------------------- Internal

      INTEGER*4::I        ,INODE    ,L        ,NNUD     ,NZARRL

      REAL*8::AREALN   ,CAUD     ,CNODE    ,CONC     ,COE
     &       ,DENSCAUD ,DENSNODE ,DENSREC  ,REC      ,RECFLUX

C--------------------  First executable statement.

C-------------------- First, computes  (Div q).

      CALL COMP_DIV_Q
     &        (AFLU     ,AREA     ,BUOYANCY ,COORD    ,DENSITY
     &        ,FLUX     ,GP_COORD ,GRADLOC  ,HAUX1    ,IAD_S
     &        ,IADN_S   ,IDIMAFLU ,IODENS   ,IODIM    ,ISOZ
     &        ,ITYPAFLU ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL
     &        ,LTYPE    ,LXPAREL  ,MAXNB    ,MAXPG    ,NPAREL
     &        ,NPPEL    ,NTYPEL   ,NUMEL    ,NUMNP    ,NZTRA
     &        ,PAREL    ,POINTWEIGHT        ,THETAT)

C-------------------- Contributions only for transient problems.

      IF (ISOLFL .NE. 1) THEN

C-------------------- Second, computes DFLU*DeltaH
          CALL PROD_MAT_VEC
     &        (1D0      ,IAD_S    ,IADN_S   ,IDIMDFLU ,NUMEL    ,NUMNP
     &        ,1        ,ITYPDFLU ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &        ,KXX      ,LNNDEL   ,DFLU     ,FLUX     ,HAUX2)

C-------------------- Third, computes CFLU*DeltaC, if density dependent flow.

          IF (IODENS.EQ.1) THEN

              CALL PROD_MAT_VEC
     &          (1D0      ,IAD_S    ,IADN_S   ,IDIMCFLU ,NUMEL    ,NUMNP
     &          ,1        ,ITYPCFLU ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &          ,KXX      ,LNNDEL   ,CFLU     ,FLUX     ,CAUX2)

         END IF !IODENS.EQ.1

      END IF !ISOLFL .NE. 1

C-------------------- Computes [(Div q) + DFLU*DeltaH + CFLU*DeltaC]*Ck+theta

      FLUX = FLUX*CAUX1

C--------------------  Fourth, computes  (ATRA*ck+theta)

      CALL PROD_MAT_VEC
     &    (1D0      ,IAD_S    ,IADN_S   ,IDIMATRA ,NUMEL    ,NUMNP
     &    ,1        ,ITYPATRA ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &    ,KXX      ,LNNDEL   ,ATRA     ,FLUX     ,CAUX1)     


C-------------------- Fifth, computes  (DTRA*DeltaC), if transient

      IF (ISOLTR .NE. 1) THEN

          CALL PROD_MAT_VEC
     &        (1D0      ,IAD_S    ,IADN_S   ,IDIMDTRA ,NUMEL    ,NUMNP
     &        ,1        ,ITYPDTRA ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &        ,KXX      ,LNNDEL   ,DTRA     ,FLUX     ,CAUX2)

      END IF !ISOLTR .NE. 1

C-------------------- Nodal flow contribution except for nodes with
C-------------------- prescribed head where it is null. At these nodes
C-------------------- caudal stores the in/out flow for transport
C-------------------- purposes.

      DO I=1,NUMNP

          IF (IBCOD(I).NE.1) THEN

              CAUD = CAUDAL(I)

          ELSE
               CAUD = 0D0

          END IF !IBCOD(I).NE.1

          IF (CAUD.GT.0) THEN

              COE = PARNP(I,4)

          ELSE

              COE = CAUX1(I)

          END IF !CAUD.GT.0

          IF (IODENS.EQ.1) THEN

              DENSCAUD = DENS(DENSREF,BETAC,COE,CREF)

          ELSE

              DENSCAUD = 1D0

          END IF !IODENS.EQ.1

          FLUX(I) = FLUX(I) - DENSCAUD*CAUD*COE

      END DO !I=1,NUMNP

C-------------------- Then recharge contribution (Only prescribed conc. nodes).

      IF (NZARR.GT.0) THEN

          DO L=1,NUMEL

              NNUD = LNNDEL(L)
              NZARRL = LXARR(L)

              IF (NZARRL.GT.0 .AND. ANY(IBTCO(KXX(1:NNUD,L)).EQ.1)) THEN

                  AREALN = AREA(L)/NNUD

                  REC = ARRC(L)*AREALN

                  IF (IODENS.EQ.1 .AND. REC.GT.0) THEN
                  
                      CONC = COECEL(L)

                      DENSREC = DENS(DENSREF,BETAC,CONC,CREF)

                  ELSE

                      DENSREC = 1D0
                      DENSNODE = 1D0
                      CONC = 0D0
 
                  END IF !IODENS.EQ.1 .AND. REC.GT.0


                  DO I=1,NNUD

                      INODE = KXX(I,L)

                      IF (IBTCO(INODE).EQ.1) THEN

                          RECFLUX = - REC*DENSREC*CONC

                          IF(IORECATRA.EQ.0) THEN

                              CNODE = CAUX1(INODE)

                              IF (IODENS.EQ.1) THEN

                                  DENSNODE =
     &                                   DENS(DENSREF,BETAC,CNODE,CREF)

                              END IF !IODENS.EQ.1

                              RECFLUX = RECFLUX - REC*DENSNODE*CNODE

                          END IF !IORECATRA.EQ.0

                          FLUX(INODE) = FLUX(INODE) + RECFLUX

                      END IF !IBTCO(INODE):EQ.1

                  END DO !I=1,NNUD

              END IF !NZARRL.GT.0 .AND. ANY presc. conc. node

          END DO !L=1,NUMEL

      END IF !NZARR.GT.0

C--------------------  Set FLUX to zero in nodes where there is not prescribed
C--------------------  concentration.

      DO I=1,NUMNP

          IF (IBTCO(I).NE.1) THEN

              FLUX(I) = 0D0

          END IF !IBTCO(1.).NE.1

      END DO !I=1,NUMNP

      END SUBROUTINE COMFLOW_PRESC_CONC
