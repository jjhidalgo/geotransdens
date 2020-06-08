      SUBROUTINE COMP_BFLU_BUOYANCY
     &          (AREA     ,BETAC    ,BFLU     ,BUOYANCY ,COORD
     &          ,DBUOYANCY,DENSITY  ,DFLUDTRA ,GP_COORD ,GRADLOC
     &          ,IOCALCDEV,IODIM    ,ISOZ     ,KXX      ,LDIM
     &          ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPAREL  ,MAXPG
     &          ,NPAREL   ,NPPEL    ,NTYPEL   ,NUMEL    ,NUMNP
     &          ,NZTRA    ,PAREL    ,POINTWEIGHT        ,THETAT)

********************************************************************************
*
*  GRADLOC: GRADLOC(j,i,k) contains the gradient in local coordinates 
*           of shape function j in direction i evaluted in node k.
*           The notation in inline comments is: @Ni/@Xj (k)
*           Where Xj = X, Y, Z (J=1,2,3, respectively)
*           i.e. Derivative of shape function "i", with respect
*           direction "Xj", evaluates in node "k".
*
*  BUOYANCY  Buoyancy coefficients. Its dimensions are (IODIM,LMXNDL,NUMEL)
*            For each element (NUMEL), contains the coeffcient for each node
*            (LMXNDL) in each dimension (IODIM).
*
*  DBUOYANCY  Derivatives of buoyancy coefficients w. r. t. concentration.
*             Contains for each element (NUMEL), de derivatives of de
*             coefficient in each dimension (IODIM) w. r. t. the concentration
*             of the nodes of the element. There are NxN derivatives (N=#nodes).
*             The position (I-1)*N + J contains the derivative of coefficient I
*             w. r. t. concentration of node J.
*
*
********************************************************************************          
       
      IMPLICIT NONE

C--------------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4::IOCALCDEV,IODIM,LMXNDL,MAXPG,NPAREL,NPPEL,NTYPEL,NUMEL
     &          ,NUMNP,NZTRA

      REAL*8::BETAC,THETAT

C---------------------------  EXTERNAL VARIABLES: ARRAYS

      REAL*8::AREA(NUMEL)                  ,BFLU(NUMNP)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL) ,COORD(NUMNP,3)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DBUOYANCY(IODIM,LMXNDL*LMXNDL,NUMEL)
     &       ,DENSITY(NUMEL)
     &       ,GP_COORD(6,8,IODIM)          ,GRADLOC(IODIM,LMXNDL,MAXPG)
     &       ,PAREL(NUMEL,NPPEL)           ,POINTWEIGHT(MAXPG,NTYPEL)
            
      INTEGER*4::ISOZ(NZTRA)   ,KXX(LMXNDL,NUMEL) ,LDIM(NUMEL)
     &          ,LNNDEL(NUMEL) ,LTYPE(NUMEL)      ,LXPAREL(NUMEL,NPAREL)


C---------------------------  INTERNAL VARIABLES: SCALARS

      REAL*8::DETERMNTJ,SUM,SUMTOT

      INTEGER*4::I        ,IDM      ,IGAUSS   ,IJPOS    ,INODE
     &          ,ISMAX    ,ISZ      ,J        ,K        ,KJ_POS     ,L
     &          ,LD       ,LTYP     ,NG       ,NNUD     ,NZONE


C---------------------------  INTERNAL VARIABLES: ARRAYS

      REAL*8::DGRAV(IODIM,LMXNDL)     ,DSUM(LMXNDL)     ,DSUMTOT(LMXNDL)
     &       ,GRAV(IODIM)             ,TRACT(6)
     &       ,TRANSDGRAV(IODIM,LMXNDL),TRANSGRAV(IODIM)
     &       ,TRANSGRD(IODIM)         ,TRANSINV(IODIM,IODIM)
     &       ,XGAUS(MAXPG)            ,XYZNUD(LMXNDL,IODIM)
     &       ,YGAUS(MAXPG)            ,ZGAUS(MAXPG)


      INTEGER*4::IND(3,3,3), NUMGP(6)

      DATA ((IND(I,J,3),I=1,3),J=1,3)/1,4,5,4,2,6,5,6,3/
      DATA ((IND(I,J,2),I=1,2),J=1,2)/1,3,3,2/,IND(1,1,1)/1/

      DATA (NUMGP(I),I=1,2)/1,1/,(NUMGP(I),I=3,4)/4,4/
      DATA  NUMGP(5)/1/,NUMGP(6)/6/

      TRACT = 0D0
      XGAUS = 0D0
      YGAUS = 0D0
      ZGAUS = 0D0
C--------------------------- For each element...

      DO L=1,NUMEL

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

          XGAUS(1:NG) = GP_COORD(LTYP,1:NG,1)
          
          IF (IODIM.GT.1) THEN
              YGAUS(1:NG) = GP_COORD(LTYP,1:NG,2)
          END IF

          IF (IODIM.GT.2) THEN
              ZGAUS(1:NG) = GP_COORD(LTYP,1:NG,3)
          END IF

          CALL COMP_GRAD_LOC
     &        (GRADLOC  ,IODIM    ,LMXNDL   ,LTYP     ,NG       ,XGAUS
     &        ,YGAUS    ,ZGAUS)


C--------------------------- Computes buoyancy contribution 
C--------------------------- to BFLU and its derivatives.

          DO I=1,NNUD

              INODE = KXX(I,L)
              SUMTOT = 0D0
              DSUMTOT = 0D0
              
              DO IGAUSS=1,NG

C--------------------------- Computes "local gravity" (and its derivative)
C--------------------------- at Gauss point as a sumatory of bouyancy
C--------------------------- coficients by shape funtions
C--------------------------- (evaluated at Gauss point).

                  SUM = 0D0
                  DSUM = 0D0
                  GRAV = 0D0
                  TRANSGRD = 0D0
                  DGRAV =0D0
                  TRANSGRAV = 0D0
                  TRANSDGRAV = 0D0

                  DO K=1,NNUD
                      DO IDM=1,LD
                          GRAV(IDM) = GRAV(IDM)
     &                     +BUOYANCY(IDM,K,L)*GRADLOC(IDM,K,IGAUSS)

                          IF (IOCALCDEV.EQ.1) THEN
                              
                              DO J=1,NNUD
                                  KJ_POS = (K - 1)*NNUD + J
                                  DGRAV(IDM,J) = DGRAV(IDM,J)
     &                                 + DBUOYANCY(IDM,KJ_POS,L)
     &                                 *GRADLOC(IDM,K,IGAUSS)
                              END DO !J=1,NNUD
                          END IF
                      END DO !IDM=1,LD
                  END DO !K=1,NNUD

C--------------------------- Computes inverse of the Jacobian
C--------------------------- at Gauss point ad its determinant.
                  
                  CALL COMP_TRANS_MAT
     &                (AREA(L)  ,DETERMNTJ,XYZNUD   ,GRADLOC  ,LD
     &                ,IGAUSS   ,LTYP     ,MAXPG    ,NNUD     ,TRANSINV)

                  DETERMNTJ = DABS(DETERMNTJ)

C--------------------------- Transformation of the gradient of 
C--------------------------- shape function I at Gauss point,
C--------------------------- transformation of the gravity 
C--------------------------- and transformation of its derivative.


                  DO J=1,LD
                      DO K=1,LD

                          TRANSGRD(J)=TRANSGRD(J) 
     &                    + TRANSINV(J,K)*GRADLOC(K,I,IGAUSS)

                          TRANSGRAV(J)=TRANSGRAV(J) 
     &                    + TRANSINV(J,K)*GRAV(K)

                          IF (IOCALCDEV.EQ.1) THEN

                              TRANSDGRAV(J,1:NNUD)=TRANSDGRAV(J,1:NNUD)
     &                        + TRANSINV(J,K)*DGRAV(K,1:NNUD)

                          END IF !IOCALCDEV.EQ.1

                      END DO !DO J=1,LD
                  END DO !DO K=1,LD


C--------------------------- Finally the product of TRANSGRAV*K*TRANSGRD

                  DO J=1,LD
                      DO K=1,LD

                          SUM = SUM + TRACT(IND(J,K,LD))*TRANSGRAV(K)
     &                                *TRANSGRD(J)

                          IF (IOCALCDEV.EQ.1) THEN
                              DSUM(1:NNUD) = DSUM(1:NNUD)
     &                        + TRACT(IND(J,K,LD))*TRANSDGRAV(K,1:NNUD)
     &                        *TRANSGRD(J)
                          END IF

                      END DO !DO J2=1,LD
                  END DO !DO J1=1,LD

C--------------------------- Contribution of the Gauss point is weighted
C--------------------------- and multiplied by the determinant of the
C--------------------------- transformation (Jacobian matrix).

                  SUMTOT= SUMTOT +SUM*POINTWEIGHT(LTYP,IGAUSS)*DETERMNTJ

                  IF (IOCALCDEV.EQ.1) THEN
                      DSUMTOT(1:NNUD) = DSUMTOT(1:NNUD) + DSUM(1:NNUD)
     &                              *POINTWEIGHT(LTYP,IGAUSS)*DETERMNTJ
                  END IF !IOCALCDEV.EQ.1

              END DO !IGAUS=1,NG

              BFLU(INODE) = BFLU(INODE) - SUMTOT*DENSITY(L)

              IF (IOCALCDEV.EQ.1) THEN

                  DO J=1,NNUD

                      IJPOS = (I-1)*NNUD + J

                      DFLUDTRA(L,IJPOS) = DFLUDTRA(L,IJPOS) 
     &                                + (DSUMTOT(J) + BETAC*SUMTOT/NNUD)
     &                                  *DENSITY(L)*THETAT
                  END DO !J=1,NNUD

              END IF !IOCALCDEV.EQ.1

          END DO !I=1,NNUD
      END DO     !DO L=1,NUMEL


      END SUBROUTINE COMP_BFLU_BUOYANCY
