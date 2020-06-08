      SUBROUTINE COMFLOW_PRESC_HEAD
     &          (AFLU     ,AREA     ,BETAC    ,BUOYANCY ,CAUDAL
     &          ,CAUX2    ,CFLU     ,COORD    ,CREF     ,DBUOYANCY
     &          ,DENSITY  ,DENSREF  ,DFLU     ,DFLUDFLU ,DFLUDTRA
     &          ,DPARELDH ,GP_COORD ,GRADLOC  ,GRAVEL   ,GRDFF
     &          ,HAUX1    ,HAUX2    ,IBCOD    ,IDIMAFLU ,IDIMDFLU
     &          ,IODENS   ,IODIM    ,IOFLLI   ,ISOLFL   ,ISOZ
     &          ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &          ,LXARR    ,LXPAREL  ,MAXPG    ,NPAREL   ,NPPEL
     &          ,NTYPEL   ,NUMEL    ,NUMNP    ,NZTRA    ,PAREL
     &          ,POINTWEIGHT        ,THETAT)

********************************************************************************
*
*     Computes nodal fluxes at prescribed head nodes
*     DENSITY CONTRIBUTION NOT REMOVED
*
*
*
********************************************************************************

      IMPLICIT NONE

      

C------------------------- External

      INTEGER*4::IDIMAFLU ,IDIMDFLU ,IODENS   ,IODIM    ,IOFLLI
     &          ,ISOLFL   ,LMXNDL   ,MAXPG    ,NPAREL   ,NPPEL
     &          ,NTYPEL   ,NUMEL    ,NUMNP    ,NZTRA

      REAL*8::BETAC    ,CREF     ,DENSREF  ,THETAT

      REAL*8::DENS

      INTEGER*4::IBCOD(NUMNP) ,ISOZ(NZTRA)           ,KXX(LMXNDL,NUMEL)
     &          ,LDIM(NUMEL)  ,LNNDEL(NUMEL)         ,LTYPE(NUMEL)
     &          ,LXARR(NUMEL) ,LXPAREL(NUMEL,NPAREL)

      REAL*8::AFLU(NUMEL,IDIMAFLU)                 ,AREA(NUMEL)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL)
     &       ,CAUDAL(NUMNP)                        ,CAUX2(NUMNP)
     &       ,CFLU(NUMEL,IDIMDFLU)                 ,COORD(NUMNP,3)
     &       ,DBUOYANCY(IODIM,LMXNDL*LMXNDL,NUMEL)
     &       ,DENSITY(NUMEL)                       ,DFLU(NUMEL,IDIMDFLU)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DPARELDH(NPPEL,NUMEL)                ,GP_COORD(6,8,IODIM)
     &       ,GRADLOC(IODIM,LMXNDL,MAXPG)          ,GRAVEL(NUMEL,3)
     &       ,GRDFF(IODIM,LMXNDL,NUMEL)            ,HAUX1(NUMNP)
     &       ,HAUX2(NUMNP)                         ,PAREL(NUMEL,NPPEL)
     &       ,POINTWEIGHT(MAXPG,NTYPEL)

C------------------------- Internal

      INTEGER*4::I      ,INOD_IJ,INODE  ,J      ,JNODE  ,L      ,NNUD
     &          ,NZ
      REAL*8::CEXT     ,DENSEXT

      REAL*8,ALLOCATABLE::CAUDBYNCY(:)


C------------------------- First executable statement
      DO L=1,NUMEL

          NNUD=LNNDEL(L) 
          NZ=LXARR(L)

          DO I=1,NNUD

              INODE=KXX(I,L)
              
              IF (I.NE.NNUD) THEN 

                  DO J=I+1,NNUD

                      JNODE = KXX(J,L)


                      INOD_IJ = (I-1)*NNUD + J- I*(I+1)/2

                      IF (IBCOD(INODE).EQ.1) THEN   
                          CAUDAL(INODE) = CAUDAL(INODE)
     &                    + AFLU(L,INOD_IJ)*(HAUX1(JNODE)-HAUX1(INODE))

                      END IF

                      IF (IBCOD(JNODE).EQ.1) THEN   
                          CAUDAL(JNODE) = CAUDAL(JNODE)
     &                    + AFLU(L,INOD_IJ)*(HAUX1(INODE)-HAUX1(JNODE))

                      END IF

                   END DO !J=I+1,NNUD

               END IF !(I.NE.NNUD)
               
C------------------------Presc. head. Transient state.

              IF (IBCOD(INODE).EQ.1) THEN

                  IF (ISOLFL.GT.1) THEN

                      CAUDAL(INODE) = CAUDAL(INODE)
     &                                + DFLU(L,I)*HAUX2(INODE)


C------------------------Density-dependent flow contribution.

                      IF (IODENS.GT.0) THEN
                          
                          CAUDAL(INODE) = CAUDAL(INODE)
     &                                    + CFLU(L,I)*CAUX2(INODE)

                      END IF !IODENS.GT.0

                  END IF !(ISOLFL.LT.1)

C------------------------- Correction for areal recharge
C------------------------- (includes evaporation)

                  IF (NZ.NE.0) THEN

                      IF (IODENS.EQ.1 .AND. PAREL(L,8).GT.0) THEN

                          CEXT = PAREL(L,15)
                          DENSEXT = DENS(DENSREF,BETAC,CEXT,CREF)

                      ELSE IF (IODENS.EQ.1 .AND. PAREL(L,8).LT.0) THEN

                          DENSEXT = DENS(DENSREF,BETAC,0D0,CREF)

	                ELSE

                          DENSEXT = 1D0

                      END IF !IODENS.EQ.1 .AND. PAREL(L,8).GT.0

                      CAUDAL(INODE) = CAUDAL(INODE)
     &                                - DENSEXT*PAREL(L,8)*AREA(L)/NNUD

                  END IF !(NZ.NE.0)

              END IF !IBCOD(INODE).EQ.1

          END DO !I=1,NNUD

      END DO !L=1,NUMEL

C------------------------ Buoyancy contribution in density dependent flow.

      IF (IODENS.EQ.1 .OR. IOFLLI.NE.0) THEN

          ALLOCATE(CAUDBYNCY(NUMNP))
          CAUDBYNCY = 0D0

	END IF !IODENS.EQ.1 .OR. IOFFLI.NE.0

      IF (IODENS.EQ.1) THEN

          CALL COMP_BFLU_BUOYANCY
     &        (AREA     ,BETAC    ,CAUDBYNCY,BUOYANCY ,COORD
     &        ,DBUOYANCY,DENSITY  ,DFLUDTRA ,GP_COORD ,GRADLOC
     &        ,0        ,IODIM    ,ISOZ     ,KXX      ,LDIM
     &        ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPAREL  ,MAXPG
     &        ,NPAREL   ,NPPEL    ,NTYPEL   ,NUMEL    ,NUMNP
     &        ,NZTRA    ,PAREL    ,POINTWEIGHT        ,THETAT)

      
      ELSE IF (IODENS.EQ.0 .AND. IOFLLI.NE.0) THEN

          CALL COMP_FLOWGRAV
     &          (AREA       ,CAUDBYNCY  ,DFLUDFLU   ,DPARELDH   ,GRAVEL
     &          ,GRDFF      ,0          ,IODIM      ,KXX        ,LDIM
     &          ,LMXNDL     ,LNNDEL     ,NUMNP      ,NPPEL      ,NUMEL
     &          ,PAREL)

      END IF !IODENS.EQ.0 .AND. IOFLLI.NE.0

      IF (IODENS.EQ.1 .OR. IOFLLI.NE.0) THEN

          DO I=1,NUMNP

              IF (IBCOD(I).EQ.1) THEN

                  CAUDAL(I) = CAUDAL(I) - CAUDBYNCY(I)

              END IF !IBCOD(I).EQ.1

          END DO !I=1,NUMNP

          DEALLOCATE(CAUDBYNCY)

      END IF !IODENS.EQ.1 .OR. IOFFLI.NE.0

	END SUBROUTINE COMFLOW_PRESC_HEAD
