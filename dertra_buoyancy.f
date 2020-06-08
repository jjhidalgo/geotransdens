      SUBROUTINE DERTRA_BUOYANCY
     &          (AREA     ,BUOYANCY ,COORD    ,DENSITY  ,DERH
     &          ,DERIV    ,DERTRA   ,GP_COORD ,GRADLOC  ,IDIMDERH
     &          ,INEW     ,IPOS     ,IODIM    ,ISOZ     ,KXX
     &          ,L        ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &          ,LXPAREL  ,MAXPG    ,NPAREL   ,NPARF    ,NPTOT
     &          ,NTYPEL   ,NUMEL    ,NUMNP    ,NZTRA    ,POINTWEIGHT)

********************************************************************************
*
*  PURPOSE: Computes the derivative  of bouyancy term w.r.t. transmisivity
*           parameter for a given element and a given direction of anisotropy.
*
*
*  EXTERNAL VARIABLES
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
*  DERTRA    Vector containing the derivatives of transmissivity tensor w. r. t.
*            a given anisotropy direction.
*
********************************************************************************          
       
      IMPLICIT NONE

C--------------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4::IDIMDERH,INEW,IODIM,LMXNDL,MAXPG,NPAREL,NPARF,NPTOT
     &          ,NTYPEL,NUMEL ,NUMNP,NZTRA


C---------------------------  EXTERNAL VARIABLES: ARRAYS

      REAL*8::AREA(NUMEL)                  ,BUOYANCY(IODIM,LMXNDL,NUMEL)
     &       ,COORD(NUMNP,3)               ,DENSITY(NUMEL)
     &       ,DERH(NUMNP,NPARF,IDIMDERH)   ,DERIV(12)
     &       ,DERTRA(6)
     &       ,GP_COORD(6,8,IODIM)          ,GRADLOC(IODIM,LMXNDL,MAXPG)
     &       ,POINTWEIGHT(MAXPG,NTYPEL)
            
      INTEGER*4::IPOS(12)                  ,ISOZ(NZTRA)
     &          ,KXX(LMXNDL,NUMEL)         ,LDIM(NUMEL)
     &          ,LNNDEL(NUMEL)             ,LTYPE(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL)


C---------------------------  INTERNAL VARIABLES: SCALARS

      REAL*8::DETERMNTJ,SUM,SUMTOT

      INTEGER*4::I        ,IDM      ,IGAUSS   ,INODE    ,IP
     &          ,ISMAX    ,ISZ      ,J        ,K        ,L
     &          ,LD       ,LTYP     ,NG       ,NNUD     ,NZONE


C---------------------------  INTERNAL VARIABLES: ARRAYS

      REAL*8::GRAV(IODIM)             ,TRANSGRAV(IODIM)
     &       ,TRANSGRD(IODIM)         ,TRANSINV(IODIM,IODIM)
     &       ,XGAUS(MAXPG)            ,XYZNUD(LMXNDL,IODIM)
     &       ,YGAUS(MAXPG)            ,ZGAUS(MAXPG)


      INTEGER*4::IND(3,3,3) ,NUMGP(6)

      DATA ((IND(I,J,3),I=1,3),J=1,3)/1,4,5,4,2,6,5,6,3/
      DATA ((IND(I,J,2),I=1,2),J=1,2)/1,3,3,2/,IND(1,1,1)/1/

      DATA (NUMGP(I),I=1,2)/1,1/,(NUMGP(I),I=3,4)/4,4/
      DATA  NUMGP(5)/1/,NUMGP(6)/6/

      
C--------------------------- For the given element...


C--------------------------- sets element properties

      LTYP = LTYPE(L)
      LD = LDIM(L)
      NNUD = LNNDEL(L)
      NZONE = LXPAREL(L,1)      ! T zone of current element
      ISZ = ISOZ(NZONE)         ! Anisotropy degree of the zone
      ISMAX = MAX(ISZ,LD)       ! Number of T-tensor components


      XYZNUD(1:NNUD,1:IODIM) = COORD(KXX(1:NNUD,L),1:IODIM)


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
     &    (GRADLOC  ,IODIM    ,LMXNDL   ,LTYP     ,NG       ,XGAUS
     &    ,YGAUS    ,ZGAUS)


C--------------------------- Computes buoyancy contribution 
C--------------------------- to inverse problem RHS

      DO I=1,NNUD

          INODE = KXX(I,L)
          SUMTOT = 0D0
          
          DO IGAUSS=1,NG

C--------------------------- Computes "local gravity" at Gauss
C--------------------------- point as a sumatory of bouyancy
C--------------------------- coeficients by shape funtions
C--------------------------- (evaluated at Gauss point).

              SUM = 0D0
              GRAV = 0D0
	        TRANSGRD = 0D0
	        TRANSGRAV = 0D0

              DO K=1,NNUD
                  DO IDM=1,LD
                      GRAV(IDM) = GRAV(IDM)
     &                 +BUOYANCY(IDM,K,L)*GRADLOC(IDM,K,IGAUSS)
                  END DO !IDM=1,LD
              END DO !K=1,NNUD

C--------------------------- Computes inverse of the Jacobian
C--------------------------- at Gauss point ad its determinant.
                  
              CALL COMP_TRANS_MAT
     &            (AREA(L)  ,DETERMNTJ,XYZNUD   ,GRADLOC  ,LD
     &            ,IGAUSS   ,LTYP     ,MAXPG    ,NNUD     ,TRANSINV)

              DETERMNTJ = DABS(DETERMNTJ)

C--------------------------- Transformation of the gradient of 
C--------------------------- shape function I at Gauss point,
C--------------------------- transformation of the gravity 
C--------------------------- and transformation of its derivative.


              DO J=1,LD
                  DO K=1,LD

                      TRANSGRD(J)=TRANSGRD(J) 
     &                + TRANSINV(J,K)*GRADLOC(K,I,IGAUSS)

	                TRANSGRAV(J)=TRANSGRAV(J) 
     &                + TRANSINV(J,K)*GRAV(K)


                  END DO !DO J=1,LD
              END DO !DO K=1,LD


C--------------------------- Finally the product of TRANSGRAV*K*TRANSGRD

              DO J=1,LD
                  DO K=1,LD

                      SUM = SUM + DERTRA(IND(J,K,LD))*TRANSGRAV(K)
     &                            *TRANSGRD(J)

                  
                  END DO !DO J2=1,LD
              END DO !DO J1=1,LD

C--------------------------- Contribution of the Gauss point is weighted
C--------------------------- and multiplied by the determinant of the
C--------------------------- transformation (Jacobian matrix).

              SUMTOT = SUMTOT + SUM*POINTWEIGHT(LTYP,IGAUSS)*DETERMNTJ

          END DO !IGAUS=1,NG

          SUMTOT = SUMTOT*DENSITY(L)

C--------------------------- Contribution added to inverse problem RHS.

          DO IP=1,NPTOT
              
              DERH(INODE,IPOS(IP),INEW) = DERH(INODE,IPOS(IP),INEW)
     &                                  - SUMTOT*DERIV(IP)

          END DO !IP=1,NPTOT

      END DO !I=1,NNUD


      END SUBROUTINE DERTRA_BUOYANCY
