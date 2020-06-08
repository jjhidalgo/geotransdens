      SUBROUTINE COMP_TRANS_MAT
     &          (AREAL    ,DETJ     ,GCOORD   ,GRADLOC
     &          ,IDIM     ,IPOINT   ,LTYP     ,MAXPG    ,NNUD
     &          ,TRANSINV)

****************************************************************
*
*  PURPOSE
*  To calculate the inverse of the transformation matrix J. 
*  If a vector v in local coordinates is written as v(e,n) and
*  in globbal coordinates as v(x,y) then 
*  v(x,y)=inv(J)*v(e,n)
*  It also computes the determinant of J
*
C EXTERNALS
C     AREA                Element size (length for 1-D elem, area for 2-D, 
C                         volume for 3-D)                  
C     DETJ                D eterminant of Jacobian matrix
C     TRANSINV            TRANSINVation matrix of the element (inverse of Jacobian matrix)
C     IDIMTRANS           Dimension of TRANSINVation matrix; equal to dimension
C                         of highest dimensional element squared.
C     ELDIM               dimension of the element
C     LTYP                type of element
C     GCOORD              Global nodal coordinates
C     XGP, YGP, ZGP       local coordinates of point where the TRANSINVation 
C                         matrix has to be evaluated.
C         
C INTERNALS
C     AREAFCTR            factor multiplying elements of TRANSINVation matrix
C     NODEX               array containing node coordinate in x direction
C     NODEY               array containing node coordinate in y direction
C     NODEZ               array containing node coordinate in z direction
C     NODENR              node number
C     INODE               dummy for do loop 
C     
C
C    HISTORY 
C     LJS   MARCH 03 first coding
*
****************************************************************
      IMPLICIT NONE

C-------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4::IDIM     ,IPOINT   ,LTYP     ,MAXPG
     &          ,NNUD

      REAL*8::AREAL       ,DETJ

C-------------------- EXTERNAL VARIABLES: ARRAYS

      REAL*8::GRADLOC(IDIM,NNUD,MAXPG)  ,GCOORD(NNUD,IDIM)
     &       ,TRANSINV(IDIM,IDIM)
       
C-------------------- INTERNAL VARIABLES: SCALARS

      INTEGER*4::I  ,J  ,K

      REAL*8:: DFACTOR1 ,DFACTOR2
     &        ,GRAD1     ,GRAD2    ,GRAD3    ,FACT1
     &        ,FACT2     ,FACT3    ,TRM1     ,TRM2
     &        ,TRM3


C-------------------- INTERNAL VARIABLES: ARRAYS
      REAL*8 XDIF(NNUD),YDIF(NNUD),ZDIF(NNUD)

      TRANSINV=0D0
      DETJ = 0D0
C========================================================================
C=======================if we have 1d linear elements====================
C========================================================================
      
      SELECT CASE(LTYP)

          CASE(1)

            TRANSINV(1,1) = 1D0/AREAL
            DETJ = AREAL
      
C========================================================================
C=======================if we have 2d linear triangle====================
C========================================================================
          CASE(2,5)

              DETJ = (2D0*AREAL)
                  
              TRANSINV(1,1) = (GCOORD(3,2)-GCOORD(1,2))/DETJ
              TRANSINV(1,2) = (GCOORD(1,2)-GCOORD(2,2))/DETJ
              TRANSINV(2,1) = (GCOORD(1,1)-GCOORD(3,1))/DETJ
              TRANSINV(2,2) = (GCOORD(2,1)-GCOORD(1,1))/DETJ


C========================================================================
C=======================if we have 3d linear tetrahedron=================
C========================================================================
          CASE(4)
         
              XDIF(1:3) = GCOORD(2:4,1) - GCOORD(1,1)
              YDIF(1:3) = GCOORD(2:4,2) - GCOORD(1,2)
              ZDIF(1:3) = GCOORD(2:4,3) - GCOORD(1,3)
         
C-------------------- calculating determinant

              DETJ =   XDIF(1)*( YDIF(2)*ZDIF(3)-YDIF(3)*ZDIF(2) )
     &                -YDIF(1)*( XDIF(2)*ZDIF(3)-XDIF(3)*ZDIF(2) )
     &                +ZDIF(1)*( XDIF(2)*YDIF(3)-XDIF(3)*YDIF(2) )

              TRANSINV(1,1)=( YDIF(2)*ZDIF(3) - YDIF(3)*ZDIF(2))/DETJ
              TRANSINV(1,2)=(-YDIF(1)*ZDIF(3) + YDIF(3)*ZDIF(1))/DETJ
              TRANSINV(1,3)=( YDIF(1)*ZDIF(2) - YDIF(2)*ZDIF(1))/DETJ

              TRANSINV(2,1)=(-XDIF(2)*ZDIF(3) + XDIF(3)*ZDIF(2))/DETJ
              TRANSINV(2,2)=( XDIF(1)*ZDIF(3) - XDIF(3)*ZDIF(1))/DETJ
              TRANSINV(2,3)=(-XDIF(1)*ZDIF(2) + XDIF(2)*ZDIF(1))/DETJ
         
              TRANSINV(3,1)=( XDIF(2)*YDIF(3) - XDIF(3)*YDIF(2))/DETJ
              TRANSINV(3,2)=(-XDIF(1)*YDIF(3) + XDIF(3)*YDIF(1))/DETJ
              TRANSINV(3,3)=( XDIF(1)*YDIF(2) - XDIF(2)*YDIF(1))/DETJ

C========================================================================
C=======================if we have 2d quadrilaterals=====================
C========================================================================

          CASE(3)
         
              DO I=1,4

                  GRAD1 = GRADLOC(1,I,IPOINT)

                  DO J=1,4

                      GRAD2 = GRADLOC(2,J,IPOINT)

                      DFACTOR1 = GRAD1*GRAD2

                      TRM1 = GCOORD(I,1)*GCOORD(J,2)
     &                      - GCOORD(I,2)*GCOORD(J,1)

                      DETJ = DETJ + DFACTOR1*TRM1

                      IF (I.EQ.1) THEN

                          TRANSINV(1,1)= TRANSINV(1,1)
     &                                 + GRADLOC(2,J,IPOINT)*GCOORD(J,2)

                          TRANSINV(1,2)= TRANSINV(1,2)
     &                                 - GRADLOC(1,J,IPOINT)*GCOORD(J,2)

                          TRANSINV(2,1)= TRANSINV(2,1)
     &                                 - GRADLOC(2,J,IPOINT)*GCOORD(J,1)

                          TRANSINV(2,2)= TRANSINV(2,2)
     &                                 + GRADLOC(1,J,IPOINT)*GCOORD(J,1)

                      END IF

                  END DO !J

              END DO !I
          
              TRANSINV = TRANSINV/DETJ
              

C========================================================================
C=======================if we have 3d prisms=============================
C========================================================================
          CASE(6)

              DO I=1,6

                  GRAD1 = GRADLOC(1,I,IPOINT)
       
                  DO J=1,6

                      GRAD2 = GRADLOC(2,J,IPOINT)
               
                      DO K=1,6

                          GRAD3 = GRADLOC(3,K,IPOINT)
                  
                          TRM1 = GCOORD(J,2)*GCOORD(K,3)
     &                          - GCOORD(J,3)*GCOORD(K,2)

                          TRM2 = GCOORD(J,1)*GCOORD(K,3)
     &                          - GCOORD(J,3)*GCOORD(K,1)

                          TRM3 = GCOORD(J,1)*GCOORD(K,2)
     &                          - GCOORD(J,2)*GCOORD(K,1)  
                  
C-------------------- Factors multiplying the determinant

                          DFACTOR1 = GRAD1*GRAD2*GRAD3

                          DFACTOR2 = GCOORD(I,1)*TRM1
     &                             + GCOORD(I,2)*TRM2
     &                             + GCOORD(I,3)*TRM3

                          DETJ = DETJ + DFACTOR1*DFACTOR2


C-------------------- Factors multiplying the matrix elements
C-------------------- This needs to be done only for varying j and k, not i

                          IF(I.EQ.1) THEN

                              FACT1 = GRADLOC(2,J,IPOINT)
     &                                *GRADLOC(3,K,IPOINT)

                              FACT2 = GRADLOC(1,J,IPOINT)
     &                                *GRADLOC(3,K,IPOINT)

                              FACT3 = GRADLOC(1,J,IPOINT)
     &                                *GRADLOC(2,K,IPOINT)

                              TRANSINV(1,1) = TRANSINV(1,1) + FACT1*TRM1
                              TRANSINV(1,2) = TRANSINV(1,2) - FACT2*TRM1
                              TRANSINV(1,3) = TRANSINV(1,3) + FACT3*TRM1
                              TRANSINV(2,1) = TRANSINV(2,1) - FACT1*TRM2
                              TRANSINV(2,2) = TRANSINV(2,2) + FACT2*TRM2
                              TRANSINV(2,3) = TRANSINV(2,3) - FACT3*TRM2
                              TRANSINV(3,1) = TRANSINV(3,1) + FACT1*TRM2
                              TRANSINV(3,2) = TRANSINV(3,2) - FACT2*TRM3
                              TRANSINV(3,3) = TRANSINV(3,3) + FACT3*TRM3

                          ENDIF !I.EQ.1

                      END DO !K=1,6
                  END DO !J=1,6
              END DO !I=1,6

              TRANSINV = TRANSINV/DETJ  

      END SELECT
      
      END SUBROUTINE COMP_TRANS_MAT
