      SUBROUTINE DER_VDTRA_BUOYANCY
     &          (AREA     ,BUOYANCY ,COORD    ,DERTRA   ,DER_VD
     &          ,GP_COORD ,GRADLOC  ,IODIM    ,KXX      ,L
     &          ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE    ,MAXPG
     &          ,NUMEL    ,NUMNP    ,POINTWEIGHT)

********************************************************************************
*
* Computes contribution of buoyancy term to the derivative of Darcy's velocity
* w. r. t. a transmissivity parameter for a given element (direct dependence)
*
********************************************************************************
          
       
      IMPLICIT NONE

C--------------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4::IODIM,LMXNDL,MAXPG,NUMEL,NUMNP


C---------------------------  EXTERNAL VARIABLES: ARRAYS
      REAL*8::AREA(NUMEL)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL) ,COORD(NUMNP,IODIM)
     &       ,DERTRA(6)                    ,DER_VD(IODIM)
     &       ,GP_COORD(6,8,IODIM)          ,GRADLOC(IODIM,LMXNDL,MAXPG)
     &       ,POINTWEIGHT(6,8)
            
      INTEGER*4::KXX(LMXNDL,NUMEL) ,LDIM(NUMEL)
     &          ,LNNDEL(NUMEL)     ,LTYPE(NUMEL)


C---------------------------  INTERNAL VARIABLES: SCALARS
      REAL*8::DETERMNTJ

      INTEGER*4::I        ,IDM      ,IDIM     ,IGAUSS 
     &          ,J        ,K        ,KDIM
     &          ,L        ,LD       ,LTYP     ,NG       ,NNUD


C---------------------------  INTERNAL VARIABLES: ARRAYS
      INTEGER*4::IND(3,3,3)   ,NUMGP(6)

      REAL*8::GRAV(IODIM)              ,SUMA(IODIM)
     &       ,SUMATOT(IODIM)           ,TRANSGRAV(IODIM)
     &       ,TRANSINV(IODIM,IODIM)    ,XGAUS(MAXPG)
     &       ,XYZNUD(LMXNDL,IODIM)     ,YGAUS(MAXPG)
     &       ,ZGAUS(MAXPG)



      DATA ((IND(I,J,3),I=1,3),J=1,3)/1,4,5,4,2,6,5,6,3/
      DATA ((IND(I,J,2),I=1,2),J=1,2)/1,3,3,2/,IND(1,1,1)/1/

      DATA (NUMGP(I),I=1,2)/1,1/,(NUMGP(I),I=3,4)/4,4/
      DATA  NUMGP(5)/1/,NUMGP(6)/6/

    
C--------------------------- For each element...

      SUMATOT= 0D0


C--------------------------- sets element properties

      LTYP = LTYPE(L)
      LD = LDIM(L)
      NNUD = LNNDEL(L)
      
      XYZNUD(1:NNUD,1:IODIM) = COORD(KXX(1:NNUD,L),1:IODIM)

C--------------------------- Identifies the value of the T components

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
          GRAV = 0D0
          TRANSGRAV =0D0
              
          DO K=1,NNUD
              DO IDM=1,LD
                  GRAV(IDM) = GRAV(IDM)
     &                      +BUOYANCY(IDM,K,L)*GRADLOC(IDM,K,IGAUSS)

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

              END DO !DO J=1,LD
          END DO !DO K=1,LD

C--------------------------- Finally the product of TRANSGRAV*K

          DO IDIM=1,LD

              DO KDIM=1,LD 

                  SUMA(IDIM) = SUMA(IDIM)
     &                       + DERTRA(IND(IDIM,KDIM,LD))*TRANSGRAV(IDIM)

              END DO !KDIM=1,LD
          END DO !IDIM=1,LD

C--------------------------- Contribution of the Gauss point is weighted
C--------------------------- and the determinant of the transformation.

          SUMATOT = SUMATOT + SUMA*POINTWEIGHT(LTYP,IGAUSS)*DETERMNTJ

      END DO !IGAUS=1,NG

C--------------------------- The computed values are stored.
C--------------------------- The minus sign comes from Darcy's law.

      DER_VD(1:LD) = DER_VD(1:LD) - SUMATOT(1:LD)/AREA(L)

      END SUBROUTINE DER_VDTRA_BUOYANCY