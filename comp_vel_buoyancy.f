      SUBROUTINE COMP_VEL_BUOYANCY
     &          (AREA     ,BUOYANCY ,COORD    ,DBUOYANCY,DVDC_L
     &          ,GP_COORD ,GRADLOC  ,IODIM    ,IOCALCDEV,ISOZ
     &          ,KXX      ,L        ,LDIM     ,LMXNDL   ,LNNDEL
     &          ,LTYPE    ,LXPAREL  ,MAXPG    ,NPAREL   ,NPPEL
     &          ,NUMEL    ,NUMNP    ,NZTRA    ,PAREL    ,POINTWEIGHT
     &          ,VEL)

********************************************************************************
*
* Computes contribution of buoyancy term to Darcy's velocity for a given element.
*
********************************************************************************


      IMPLICIT NONE

C--------------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4::IOCALCDEV,IODIM,LMXNDL,MAXPG,NPAREL,NPPEL,NUMEL,NUMNP
     &          ,NZTRA

C---------------------------  EXTERNAL VARIABLES: ARRAYS
      REAL*8::AREA(NUMEL)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL) ,COORD(NUMNP,IODIM)
     &       ,DBUOYANCY(IODIM,LMXNDL*LMXNDL,NUMEL)
     &       ,DVDC_L(LMXNDL,IODIM) 
     &       ,GP_COORD(6,8,IODIM)          ,GRADLOC(IODIM,LMXNDL,MAXPG)
     &       ,PAREL(NUMEL,NPPEL)           ,POINTWEIGHT(6,8)

      INTEGER*4::ISOZ(NZTRA)   ,KXX(LMXNDL,NUMEL) ,LDIM(NUMEL)
     &          ,LNNDEL(NUMEL) ,LTYPE(NUMEL)      ,LXPAREL(NUMEL,NPAREL)


C---------------------------  INTERNAL VARIABLES: SCALARS
      REAL*8::DETERMNTJ

      INTEGER*4::I        ,IDM      ,IDIM     ,IGAUSS 
     &          ,ISMAX    ,ISZ      ,J        ,K,KJ_POS        ,KDIM
     &          ,L        ,LD       ,LTYP     ,NG       ,NNUD
     &          ,NZONE


C---------------------------  INTERNAL VARIABLES: ARRAYS
      INTEGER*4::IND(3,3,3),NUMGP(6)

      REAL*8::DGRAV(IODIM,LMXNDL)      ,DSUMA(IODIM,LMXNDL)
     &       ,DSUMATOT(IODIM,LMXNDL)   ,GRAV(IODIM)
     &       ,SUMA(IODIM)              ,SUMATOT(IODIM)
     &       ,TRANSDGRAV(IODIM,LMXNDL) ,TRANSGRAV(IODIM)
     &       ,TRACT(9)                 ,TRANSINV(IODIM,IODIM)
     &       ,VEL(IODIM)               ,XGAUS(MAXPG)
     &       ,XYZNUD(LMXNDL,IODIM)     ,YGAUS(MAXPG)
     &       ,ZGAUS(MAXPG)



      DATA ((IND(I,J,3),I=1,3),J=1,3)/1,4,5,4,2,6,5,6,3/
      DATA ((IND(I,J,2),I=1,2),J=1,2)/1,3,3,2/,IND(1,1,1)/1/

      DATA (NUMGP(I),I=1,2)/1,1/,(NUMGP(I),I=3,4)/4,4/
      DATA  NUMGP(5)/1/,NUMGP(6)/6/

C--------------------------- For each element...

      SUMATOT= 0D0
      DSUMATOT = 0D0
      TRACT = 0D0

C--------------------------- sets element properties

      LTYP = LTYPE(L)
      LD = LDIM(L)
      NNUD = LNNDEL(L)
      NZONE = LXPAREL(L,1)      ! T zone of current element
      ISZ = ISOZ(NZONE)         ! Anisotropy degree of the zone
      ISMAX = MAX(ISZ,LD)       ! Number of T-tensor components

      XYZNUD(1:NNUD,1:IODIM) = COORD(KXX(1:NNUD,L),1:IODIM)

C--------------------------- Identifies the value of the T components

      TRACT(1:ISMAX) = PAREL(L,1:ISMAX) 

C--------------------------- and the number of Gauss points.
          
      NG = NUMGP(LTYP)

C--------------------------- Computes local gradident at Gauss points.

      XGAUS = 0D0
      YGAUS = 0D0
      ZGAUS = 0D0

      IF (IODIM.GE.1) THEN
          XGAUS(1:NG) = GP_COORD(LTYP,1:NG,1)
      END IF

      IF (IODIM.GE.2) THEN
          YGAUS(1:NG) = GP_COORD(LTYP,1:NG,2)
      END IF
      IF (IODIM.EQ.3) THEN
          ZGAUS(1:NG) = GP_COORD(LTYP,1:NG,3)
      END IF

      CALL COMP_GRAD_LOC
     &    (GRADLOC  ,IODIM    ,LMXNDL   ,LTYP     ,NG       ,XGAUS
     &    ,YGAUS    ,ZGAUS)


C--------------------------- Computes buoyancy contribution 
C--------------------------- to BFLU and its derivatives.

      DO IGAUSS=1,NG

C--------------------------- Computes "local gravity" at Gauss point
C--------------------------- as a sumatory of bouyancy coefficients by
C--------------------------- shape funtions (evaluated at Gauss point).
          SUMA=0D0
          DSUMA = 0D0
          GRAV = 0D0
          DGRAV = 0D0
          TRANSGRAV =0D0
          TRANSDGRAV =0D0
              
          DO K=1,NNUD
              DO IDM=1,LD
                  GRAV(IDM) = GRAV(IDM)
     &                      +BUOYANCY(IDM,K,L)*GRADLOC(IDM,K,IGAUSS)

                  IF (IOCALCDEV.EQ.1) THEN
                              
                      DO J=1,NNUD
                          KJ_POS = (K - 1)*NNUD + J
                          DGRAV(IDM,J) = DGRAV(IDM,J)
     &                                 + DBUOYANCY(IDM,KJ_POS,L)
     &                                   *GRADLOC(IDM,K,IGAUSS)
                      END DO !J=1,NNUD

                  END IF !IOCALCDEV.EQ.1

              END DO !IDM=1,LD
          END DO !K=1,NNUD

C--------------------------- Computes inverse of the Jacobian
C--------------------------- at Gauss point.
                  
          CALL COMP_TRANS_MAT
     &        (AREA(L)  ,DETERMNTJ,XYZNUD   ,GRADLOC
     &        ,LD       ,IGAUSS   ,LTYP     ,MAXPG    ,NNUD
     &        ,TRANSINV)


          DETERMNTJ = DABS(DETERMNTJ)


C--------------------------- Transformation of the gravity 
C--------------------------- and transformation of its derivative.


          DO J=1,LD
              DO K=1,LD

                  TRANSGRAV(J) = TRANSGRAV(J) + TRANSINV(J,K)*GRAV(K)

                  IF (IOCALCDEV.EQ.1) THEN

                      TRANSDGRAV(J,1:NNUD) = TRANSDGRAV(J,1:NNUD) 
     &                                   + TRANSINV(J,K)*DGRAV(K,1:NNUD)
                  END IF !IOCALCDEV.EQ.1

              END DO !DO J=1,LD
          END DO !DO K=1,LD

C--------------------------- Finally the product of TRANSGRAV*K

          DO IDIM=1,LD

              DO KDIM=1,LD 

                  SUMA(IDIM) = SUMA(IDIM)
     &                       + TRACT(IND(IDIM,KDIM,LD))*TRANSGRAV(IDIM)

                  IF (IOCALCDEV.EQ.1) THEN
                      DSUMA(IDIM,1:NNUD) = DSUMA(IDIM,1:NNUD)
     &                                   + TRACT(IND(IDIM,KDIM,LD))
     &                                     *TRANSDGRAV(IDIM,1:NNUD)
                  END IF !IOCALCDEV.EQ.1

              END DO !KDIM=1,LD
          END DO !IDIM=1,LD

C--------------------------- Contribution of the Gauss point is weighted
C--------------------------- and the determinant of the transformation.

          SUMATOT = SUMATOT + SUMA*POINTWEIGHT(LTYP,IGAUSS)*DETERMNTJ

          IF (IOCALCDEV.EQ.1) THEN

              DSUMATOT=DSUMATOT+DSUMA*POINTWEIGHT(LTYP,IGAUSS)*DETERMNTJ

          END IF !IOCALCDEV.EQ.1

      END DO !IGAUS=1,NG

C--------------------------- The computed values are stored.
C--------------------------- The minus sign comes from Darcy's law.

      VEL(1:LD) = VEL(1:LD) - SUMATOT(1:LD)/AREA(L)

      IF (IOCALCDEV.EQ.1) THEN
          DO J=1,NNUD
              DVDC_L(J,1:LD) = DVDC_L(J,1:LD) - DSUMATOT(1:LD,J)/AREA(L)
          END DO
      END IF !IOCALCDEV.EQ.1

      END SUBROUTINE COMP_VEL_BUOYANCY

