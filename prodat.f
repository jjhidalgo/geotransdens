       SUBROUTINE PRODAT 
     ;(IDIMBB   ,IERROR   ,IODIM    ,IOEQT    ,IOFLLI   
     ;,IOTRLI   ,LDIM     ,LMXNDL   ,MAINF    ,NTDMT    ,NUMEL
     ;,NUMNP    ,ACTH     ,AREA     ,BIBI     ,BTRA     ,CBASE
     ;,CCAL     ,CNST     ,COORD    ,GRAV     ,GRAVEL   ,GRDFF
     ;,HBASE    ,HCAL     ,KXX      ,LXTRA    ,LNNDEL   ,LTYPE    
     ;,VOLNOD   ,X        ,Y        ,Z       
     ;!NUEVOS
     ;,IODENS   ,IPARTNER ,LCOORD   ,POINTWEIGHT,GP_COORD
     &,IDIMGRAVEL)

*****************************************************************************
***  Perform various computations at the element level
***   BIBI(I,L)   Integral grad(FIi).grad(FIj) for element number L. Only the
***               terms above the diagonal are stored (symmetric mode).
***   GRDFF(I,J,L) Derivative of the Jth basis function of the Lth element
***               with respect to the Ith global coordinate (I=1 for X, etc).
***   CNST(i,j,k) Integral FIi.FIj in consistent scheme and local coordinates
***               for element type k.
***   AREA(L)     Size (area in 2-D, length in 1-D or volume in 3-D) of elt.L
***  
***  It also initializas arrays that are needed for the consistent velocity 
***   approximation:
***     IPARTNER
***     LCOORD
***     POINTWEIGHT   
***     GP_COORD
***  
*****************************************************************************

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION KXX(LMXNDL,NUMEL),X(NUMNP),Y(NUMNP),Z(NUMNP),
     .      BIBI(IDIMBB,NUMEL),AREA(NUMEL),ACTH(NUMEL),BTRA(NUMNP),
     .      LNNDEL(NUMEL),LTYPE(NUMEL),CNST(6,6,6),
     .      B(4),C(4),D(4),HCAL(NUMNP),CCAL(NUMNP),HBASE(NUMEL),
     .      CBASE(NUMEL),GRDFF(IODIM,LMXNDL,NUMEL),
     .      GRADN(24,6),GNGN(15,15),PESO(6)
     .      ,LCOORD(6,3,6),IPARTNER(6,3,6),POINTWEIGHT(6,6)
     .      ,GP_COORD(6,8,3)
     &      ,GRAVEL(IDIMGRAVEL,3*MIN(3*IDIMGRAVEL,3))


C------------------------- MODIFY LTYPE ARRAY FOR INVERSE PROBLEM COMPUTATIONS

       DO L=1,NUMEL
          IF (LTYPE(L).EQ.3) THEN              !TRIANGULAR ELEMENT (2-D)
             LTYPE(L)=2
          ELSE IF (LTYPE(L).EQ.5) THEN         !CUADRANGULAR ELEMENT (2-D)
             LTYPE(L)=3
          ELSE IF (LTYPE(L).EQ.9) THEN         !TETRAHEDRON ELEMENT (3-D)
             LTYPE(L)=4
          ELSE IF (LTYPE(L).EQ.10) THEN        !TRIANGULAR ELEMENT (3-D)
             LTYPE(L)=5
          ELSE IF (LTYPE(L).EQ.11)THEN         !TOBLERONE (3-D)
             LTYPE(L)=6
          END IF
       END DO




C------------------------- Define some arrays which are used
C------------------------- for the consistent velocity approximation 

       IF (IODENS.EQ.1) THEN

          IPARTNER(:,:,:) = 0

C------------------------- IPARTNER(I,J,K) The node used for interpolationg 
C------------------------- mass fraction for node I in direction J in element type K


C------------------------- One dimensional elements

          IPARTNER(1,1,1)=2
          IPARTNER(2,1,1)=1

C------------------------- Triangular elements

          IPARTNER(2,1,2)=1
          IPARTNER(3,2,2)=1

C------------------------- Triangular elements, in 3D medium.
C------------------------- Nodes not present are not used.

          IPARTNER(2,1,5)=1
          IPARTNER(3,2,5)=1



C------------------------- Quadrilateral elements.
C------------------------- Nodes not present are not used.

          IPARTNER(2,1,3)=1
          IPARTNER(3,1,3)=4
          IPARTNER(3,2,3)=1
          IPARTNER(4,2,3)=2

          
C------------------------- Tetrahedrical elements.
C------------------------- Nodes not present are not used.

          IPARTNER(2,1,4)=1
          IPARTNER(3,2,4)=1
          IPARTNER(4,3,4)=1

C------------------------- Prismatic elements.

          IPARTNER(2,2,6)=1
          IPARTNER(3,3,6)=1
          IPARTNER(4,1,6)=1
          IPARTNER(5,1,6)=2
          IPARTNER(5,2,6)=4
          IPARTNER(6,1,6)=3
          IPARTNER(6,3,6)=4


C------------------------- lcoord: the local node coordinates. 
C------------------------- lcoord(i,j,k) is the coordinate in the j direction of node i
C------------------------- in element type k

C-------------------------  1d Elements

          LCOORD(1,1,1)=0
          LCOORD(1,2,1)=1
          LCOORD(2,1,1)=0
          LCOORD(2,2,1)=0

C------------------------- Triangular elements

          LCOORD(1,1,2)=0
          LCOORD(1,2,2)=0
          LCOORD(2,1,2)=1
          LCOORD(2,2,2)=0
          LCOORD(3,1,2)=0
          LCOORD(3,2,2)=1
          LCOORD(1,1,5)=0
          LCOORD(1,2,5)=0
          LCOORD(2,1,5)=1
          LCOORD(2,2,5)=0
          LCOORD(3,1,5)=0
          LCOORD(3,2,5)=1

C------------------------- Cuadrilateral elements

          LCOORD(1,1,3)=0
          LCOORD(1,2,3)=0
          LCOORD(2,1,3)=1
          LCOORD(2,2,3)=0
          LCOORD(3,1,3)=1
          LCOORD(3,2,3)=1
          LCOORD(4,1,3)=0
          LCOORD(4,2,3)=1

C------------------------- Tetrahedical  elements.

          LCOORD(1,1,4)=0
          LCOORD(1,2,4)=0
          LCOORD(1,3,4)=0

          LCOORD(2,1,4)=1
          LCOORD(2,2,4)=0
          LCOORD(2,3,4)=0

          LCOORD(3,1,4)=0
          LCOORD(3,2,4)=1
          LCOORD(3,3,4)=0

          LCOORD(4,1,4)=0
          LCOORD(4,2,4)=0
          LCOORD(4,3,4)=1

C------------------------- Prismatical  elements.

          LCOORD(1,1,6)=0
          LCOORD(1,2,6)=0
          LCOORD(1,3,6)=0

          LCOORD(2,1,6)=0
          LCOORD(2,2,6)=1
          LCOORD(2,3,6)=0

          LCOORD(3,1,6)=0
          LCOORD(3,2,6)=0
          LCOORD(3,3,6)=1

          LCOORD(4,1,6)=1
          LCOORD(4,2,6)=0
          LCOORD(4,3,6)=0

          LCOORD(5,1,6)=1
          LCOORD(5,2,6)=1
          LCOORD(5,3,6)=0

          LCOORD(6,1,6)=1
          LCOORD(6,2,6)=0
          LCOORD(6,3,6)=1



C------------------------- Gauss points weights.
C------------------------- POINTWEIGHT(NTYPEL,MAXPG)

C------------------------- One dimensional elements

          POINTWEIGHT(1,1)=1D0

C------------------------- Triangular elements.
C------------------------- In theory, the weight must be 1.0
C------------------------- but since the whole formula is divided
C------------------------- by a factor of 2.0, this factor is 
C------------------------- included in the weigth.

          POINTWEIGHT(2,1)=0.5D0
          POINTWEIGHT(5,1)=0.5D0

C------------------------- Quadrilateral elements.

          POINTWEIGHT(3,1)=0.25D0
          POINTWEIGHT(3,2)=0.25D0
          POINTWEIGHT(3,3)=0.25D0
          POINTWEIGHT(3,4)=0.25D0
                  
C------------------------- Tetrahedrons.
C------------------------- In theory, the weight must be 0.25
C------------------------- but since the whole formula is divided
C------------------------- by a factor of 6.0, this factor is 
C------------------------- included in the weigth.

          POINTWEIGHT(4,1)=0.25D0/6D0
          POINTWEIGHT(4,2)=0.25D0/6D0
          POINTWEIGHT(4,3)=0.25D0/6D0
          POINTWEIGHT(4,4)=0.25D0/6D0

C------------------------- Prismatic elements (toblerones).

          POINTWEIGHT(6,1)=1D0/12D0
          POINTWEIGHT(6,2)=1D0/12D0
          POINTWEIGHT(6,3)=1D0/12D0
          POINTWEIGHT(6,4)=1D0/12D0
          POINTWEIGHT(6,5)=1D0/12D0
          POINTWEIGHT(6,6)=1D0/12D0


C------------------------- Gauss points coordinates.
C------------------------- GP_COORD(6,8,3)


C------------------------- One dimensional elements

          GP_COORD(1,1,1)=0.5D0

C------------------------- Triangular elements

          GP_COORD(2,1,1)=1D0/3D0  !no importan estos valores, porque
          GP_COORD(2,1,2)=1D0/3D0  !en estos elementos el valor no depende
          GP_COORD(2,1,3)=0d0      !de la posición del punto


C------------------------- Triangles in space

          GP_COORD(5,1,1)=1D0/3D0  !no importan estos valores, porque
          GP_COORD(5,1,2)=1D0/3D0  !en estos elementos el valor no depende
          GP_COORD(5,1,3)=1D0/3D0  !de la posición del punto


C------------------------- Quadrilaterals

          GP_COORD(3,1,1)=(1d0+1d0/dsqrt(3d0))/2d0 !0.78886751346
          GP_COORD(3,1,2)=(1d0+1d0/dsqrt(3d0))/2d0 !0.78886751346
          GP_COORD(3,1,3)=0.D0

          GP_COORD(3,2,1)=(1d0-1d0/dsqrt(3d0))/2d0 !0.2113248655
          GP_COORD(3,2,2)=(1d0-1d0/dsqrt(3d0))/2d0 !0.2113248655
          GP_COORD(3,2,3)=0.D0

          GP_COORD(3,3,1)=(1d0+1d0/dsqrt(3d0))/2d0 !0.78886751346
          GP_COORD(3,3,2)=(1d0-1d0/dsqrt(3d0))/2d0 !0.2113248655
          GP_COORD(3,3,3)=0.D0

          GP_COORD(3,4,1)=(1d0-1d0/dsqrt(3d0))/2d0 !0.2113248655
          GP_COORD(3,4,2)=(1d0+1d0/dsqrt(3d0))/2d0 !0.78886751346
          GP_COORD(3,4,3)=0.D0

C------------------------- Tetrahedrons
c          GP_COORD(4,1,1)=
c          GP_COORD(4,1,2)=
c          GP_COORD(4,1,3)=

c          GP_COORD(4,2,1)=
c          GP_COORD(4,2,2)=
c          GP_COORD(4,2,3)=

c          GP_COORD(4,3,1)=
c          GP_COORD(4,3,2)=
c          GP_COORD(4,3,3)=

c          GP_COORD(4,4,1)=
c          GP_COORD(4,4,2)=
c          GP_COORD(4,4,3)=

C------------------------- Prismatic elements (toblerones).

          GP_COORD(6,1,1)=2D0/3D0
          GP_COORD(6,1,2)=1D0/6D0
          GP_COORD(6,1,3)=1D0/2D0 - DSQRT(3D0)/6D0

          GP_COORD(6,2,1)=1D0/6D0
          GP_COORD(6,2,2)=2D0/3D0
          GP_COORD(6,2,3)=1D0/2D0 - DSQRT(3D0)/6D0

          GP_COORD(6,3,1)=1D0/6D0
          GP_COORD(6,3,2)=1D0/6D0
          GP_COORD(6,3,3)=1D0/2D0 - DSQRT(3D0)/6D0

          GP_COORD(6,4,1)=2D0/3D0
          GP_COORD(6,4,2)=1D0/6D0
          GP_COORD(6,4,3)=1D0/2D0 + DSQRT(3D0)/6D0

          GP_COORD(6,5,1)=1D0/6D0
          GP_COORD(6,5,2)=2D0/3D0
          GP_COORD(6,5,3)=1D0/2D0 + DSQRT(3D0)/6D0

          GP_COORD(6,6,1)=1D0/6D0
          GP_COORD(6,6,2)=1D0/6D0
          GP_COORD(6,6,3)=1D0/2D0 + DSQRT(3D0)/6D0


       ENDIF      !if iodens = 1


C------------------------- Provisional for aplying the
C------------------------- OCD method instead GALERKIN
       if (iodim.eq.3)then
         ind=0
         do l=1,numel
           if (ltype(l).eq.4)then
             ind=1
             goto 8
           endif
         end do
 8       if (ind.ne.0) then
           open (unit=116,file='distorsion.dat',status='old',err=9)
           read (116,*)iopocd
           goto 11
 9         write(mainf,333)
 333       format(
     ;     ' When using thetraedron elements, you tave to create the',/,
     ;' file DISTORSION.DAT. It indicates the methood for computing',/, 
     ;' the conductance matrix (OCD or GALERKIN) and the distorsion ',/,
     ;' tensor. I cannot find it.')
           stop
         endif
       endif
 11    continue
        

       DO 100 L=1,NUMEL
          LTYPE1=LTYPE(L)
          I1=KXX(1,L)
          I2=KXX(2,L)
          X1=X(I1)
          X2=X(I2)
          Y1=Y(I1)
          Y2=Y(I2)
          Z1=Z(I1)
          Z2=Z(I2)

C------------------------- COMPUTE THE AVERAGE OF THE THICKNESS ON THE ELEMENT

          IF (LTYPE1.EQ.1) THEN
             BIBI(1,L)=-1/AREA(L)
             GRDFF(1,1,L)=-1/AREA(L)
             GRDFF(1,2,L)=-GRDFF(1,1,L)
             IF (IOFLLI.NE.0)THEN
               HBASE(L)=(HCAL(I1)+HCAL(I2))/2D0
             ENDIF
             IF (IOTRLI.NE.0)THEN
               CBASE(L)=(CCAL(I1)+CCAL(I2))/2D0
             ENDIF

          ELSE IF (LTYPE1.EQ.2) THEN
             I3=KXX(3,L)
             X3=X(I3)
             Y3=Y(I3)
             B1=Y2-Y3
             B2=Y3-Y1
             C1=X3-X2
             C2=X1-X3
             B3=Y1-Y2
             C3=X2-X1
             IF (IOEQT.NE.1) 
     .       ACTH(L)=(BTRA(I1)+BTRA(I2)+BTRA(I3))/3.D+00
             IF (IOFLLI.NE.0)THEN
               HBASE(L)=(HCAL(I1)+HCAL(I2)+HCAL(I3))/3D0
             ENDIF
             IF (IOTRLI.NE.0)THEN
               CBASE(L)=(CCAL(I1)+CCAL(I2)+CCAL(I3))/3D0
             ENDIF

             AREA2=AREA(L)*4.D+00

C------------------------- PRODUCTS OF COEF. BEFORE COMPUTING, AND THEY STORAGE

             BIBI(1,L)=B1*B2/ AREA2      
             BIBI(2,L)=C1*C2/ AREA2
             BIBI(3,L)=(B2*C1+B1*C2) / AREA2
             BIBI(4,L)=B3*B1/ AREA2
             BIBI(5,L)=C3*C1/ AREA2
             BIBI(6,L)=(B1*C3+B3*C1) / AREA2
             BIBI(7,L)=B2*B3/ AREA2
             BIBI(8,L)=C2*C3/ AREA2      
             BIBI(9,L)=(B3*C2+B2*C3) / AREA2      

             GRDFF(1,1,L)=B1/AREA(L)/2
             GRDFF(1,2,L)=B2/AREA(L)/2
             GRDFF(1,3,L)=-(B1+B2)/AREA(L)/2
             GRDFF(2,1,L)=C1/AREA(L)/2
             GRDFF(2,2,L)=C2/AREA(L)/2
             GRDFF(2,3,L)=-(C1+C2)/AREA(L)/2
          ELSE IF (LTYPE1.EQ.3) THEN
             I3=KXX(3,L)
             I4=KXX(4,L)
             X3=X(I3)
             Y3=Y(I3)
             X4=X(I4)
             Y4=Y(I4)
             HL=DABS(X4-X2)
             HH=DABS(Y4-Y2)
             WX=HH/HL/6D0
             WY=HL/HH/6D0
             BIBI(1,L)=-2*WX
             BIBI(2,L)=WY
             BIBI(3,L)=0.D0
             BIBI(4,L)=-WX
             BIBI(5,L)=-WY
             BIBI(6,L)=-1.D0/2
             BIBI(7,L)=WX
             BIBI(8,L)=-2*WY
             BIBI(9,L)=0.D0
             BIBI(10,L)=WX
             BIBI(11,L)=-2*WY
             BIBI(12,L)=0.D0
             BIBI(13,L)=-WX
             BIBI(14,L)=-WY
             BIBI(15,L)=1.D0/2
             BIBI(16,L)=-2*WX
             BIBI(17,L)=WY
             BIBI(18,L)=0.D0

             IF (IOEQT.NE.1) 
     .       ACTH(L)=(BTRA(I1)+BTRA(I2)+BTRA(I3)+BTRA(I4))/4.D+00
             IF (IOFLLI.NE.0)THEN
               HBASE(L)=(HCAL(I1)+HCAL(I2)+HCAL(I3)+HCAL(I4))/4D0
             ENDIF
             IF (IOTRLI.NE.0)THEN
               CBASE(L)=(CCAL(I1)+CCAL(I2)+CCAL(I3)+CCAL(I4))/4D0
             ENDIF

             GRDFF(1,1,L)=1/(X1-X2)/2
             GRDFF(1,2,L)=1/(X2-X1)/2
             GRDFF(1,3,L)=1/(X2-X1)/2
             GRDFF(1,4,L)=1/(X1-X2)/2
             GRDFF(2,1,L)=1/(Y1-Y4)/2
             GRDFF(2,2,L)=1/(Y1-Y4)/2
             GRDFF(2,3,L)=1/(Y4-Y1)/2
             GRDFF(2,4,L)=1/(Y4-Y1)/2

          ELSE IF(LTYPE1.EQ.4)THEN !TETRAHEDRON'S VOLUME

       I3=KXX(3,L)
       I4=KXX(4,L)
       x1=x(i1)
       y1=y(i1)
       z1=z(i1)
       x2=x(i2)
       y2=y(i2)
       z2=z(i2)
       x3=x(i3)
       y3=y(i3)
       z3=z(i3)
       x4=x(i4)
       y4=y(i4)
       z4=z(i4)

       b(1)=-(y3*z4-y4*z3)+(y2*z4-y4*z2)-(y2*z3-y3*z2)
       b(2)=+(y3*z4-y4*z3)-(y1*z4-y4*z1)+(y1*z3-y3*z1)
       b(3)=-(y2*z4-y4*z2)+(y1*z4-y4*z1)-(y1*z2-y2*z1)
       b(4)=+(y2*z3-y3*z2)-(y1*z3-y3*z1)+(y1*z2-y2*z1)

       c(1)=+(x3*z4-x4*z3)-(x2*z4-x4*z2)+(x2*z3-x3*z2)
       c(2)=-(x3*z4-x4*z3)+(x1*z4-x4*z1)-(x1*z3-x3*z1)
       c(3)=+(x2*z4-x4*z2)-(x1*z4-x4*z1)+(x1*z2-x2*z1)
       c(4)=-(x2*z3-x3*z2)+(x1*z3-x3*z1)-(x1*z2-x2*z1)

       d(1)=-(x3*y4-x4*y3)+(x2*y4-x4*y2)-(x2*y3-x3*y2)
       d(2)=+(x3*y4-x4*y3)-(x1*y4-x4*y1)+(x1*y3-x3*y1)
       d(3)=-(x2*y4-x4*y2)+(x1*y4-x4*y1)-(x1*y2-x2*y1)
       d(4)=+(x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1)


       xb=6.0D0*AREA(L)

       gradn(1,1)=b(1)/xb    ! dN1/dx
       gradn(2,1)=c(1)/xb    ! dN1/dy
       gradn(3,1)=d(1)/xb    ! dN1/dz
       gradn(4,1)=b(2)/xb    ! ...
       gradn(5,1)=c(2)/xb 
       gradn(6,1)=d(2)/xb
       gradn(7,1)=b(3)/xb
       gradn(8,1)=c(3)/xb
       gradn(9,1)=d(3)/xb
       gradn(10,1)=b(4)/xb
       gradn(11,1)=c(4)/xb
       gradn(12,1)=d(4)/xb





*..............................................tetrahedron (3-D)
       g1=gradn(1,1)*AREA(L)
       g2=gradn(2,1)*AREA(L)
       g3=gradn(3,1)*AREA(L)

       gngn(1,1)=g1*gradn(4,1)  !xx 12
       gngn(2,1)=g1*gradn(5,1)  !xy
       gngn(3,1)=g1*gradn(6,1)  !xz
       gngn(4,1)=g2*gradn(4,1)  !yx
       gngn(5,1)=g2*gradn(5,1)  !yy
       gngn(6,1)=g2*gradn(6,1)  !yz
       gngn(7,1)=g3*gradn(4,1)  !zx
       gngn(8,1)=g3*gradn(5,1)  !zy
       gngn(9,1)=g3*gradn(6,1)  !zz

       gngn(1,2)=g1*gradn(7,1)  !13
       gngn(2,2)=g1*gradn(8,1)
       gngn(3,2)=g1*gradn(9,1)
       gngn(4,2)=g2*gradn(7,1)
       gngn(5,2)=g2*gradn(8,1)
       gngn(6,2)=g2*gradn(9,1)
       gngn(7,2)=g3*gradn(7,1)
       gngn(8,2)=g3*gradn(8,1)
       gngn(9,2)=g3*gradn(9,1)

       gngn(1,3)=g1*gradn(10,1)  !14
       gngn(2,3)=g1*gradn(11,1)
       gngn(3,3)=g1*gradn(12,1)
       gngn(4,3)=g2*gradn(10,1)
       gngn(5,3)=g2*gradn(11,1)
       gngn(6,3)=g2*gradn(12,1)
       gngn(7,3)=g3*gradn(10,1)
       gngn(8,3)=g3*gradn(11,1)
       gngn(9,3)=g3*gradn(12,1)

       g4=gradn(4,1)*AREA(L)
       g5=gradn(5,1)*AREA(L)
       g6=gradn(6,1)*AREA(L)

       gngn(1,4)=g4*gradn(7,1)    !23
       gngn(2,4)=g4*gradn(8,1)
       gngn(3,4)=g4*gradn(9,1)
       gngn(4,4)=g5*gradn(7,1)
       gngn(5,4)=g5*gradn(8,1)
       gngn(6,4)=g5*gradn(9,1)
       gngn(7,4)=g6*gradn(7,1)
       gngn(8,4)=g6*gradn(8,1)
       gngn(9,4)=g6*gradn(9,1)

       gngn(1,5)=g4*gradn(10,1)  !24
       gngn(2,5)=g4*gradn(11,1)
       gngn(3,5)=g4*gradn(12,1)
       gngn(4,5)=g5*gradn(10,1)
       gngn(5,5)=g5*gradn(11,1)
       gngn(6,5)=g5*gradn(12,1)
       gngn(7,5)=g6*gradn(10,1)
       gngn(8,5)=g6*gradn(11,1)
       gngn(9,5)=g6*gradn(12,1)

       g7=gradn(7,1)*AREA(L)
       g8=gradn(8,1)*AREA(L)
       g9=gradn(9,1)*AREA(L)
       gngn(1,6)=g7*gradn(10,1)  !34
       gngn(2,6)=g7*gradn(11,1)
       gngn(3,6)=g7*gradn(12,1)
       gngn(4,6)=g8*gradn(10,1)
       gngn(5,6)=g8*gradn(11,1)
       gngn(6,6)=g8*gradn(12,1)
       gngn(7,6)=g9*gradn(10,1)
       gngn(8,6)=g9*gradn(11,1)
       gngn(9,6)=g9*gradn(12,1)


       if (iopocd.eq.0)then             !Galerkin
         ibibi=-6
         do icon=1,6
           ibibi=ibibi+6
           bibi(ibibi+1,l)=gngn(1,icon)                 !xx
           bibi(ibibi+2,l)=gngn(5,icon)                 !yy
           bibi(ibibi+3,l)=gngn(9,icon)                 !zz
           bibi(ibibi+4,l)=gngn(2,icon)+gngn(4,icon)    !xy
           bibi(ibibi+5,l)=gngn(3,icon)+gngn(7,icon)    !xz
           bibi(ibibi+6,l)=gngn(6,icon)+gngn(8,icon)    !yz
         enddo
       else                                      !ocd
         iocd=116
         call tetrah_conductance_ocd_matrix 
     ;   (bibi,x1,x2,x3,x4,y1,y2,y3,y4,
     ;   z1,z2,z3,z4,l,idimbb,numel,lxtra,iocd)
       endif



             grdff(1,1,l)=gradn(1,1)  !  b(1)/xb =dN1/dx
             grdff(2,1,l)=gradn(2,1)  !  c(1)/xb =dN1/dy
             grdff(3,1,l)=gradn(3,1)  !  d(1)/xb =dN1/dz
             grdff(1,2,l)=gradn(4,1)  !  b(2)/xb ...
             grdff(2,2,l)=gradn(5,1)  !  c(2)/xb
             grdff(3,2,l)=gradn(6,1)  !  d(2)/xb
             grdff(1,3,l)=gradn(7,1)  !  b(3)/xb
             grdff(2,3,l)=gradn(8,1)  !  c(3)/xb
             grdff(3,3,l)=gradn(9,1)  !  d(3)/xb
             grdff(1,4,l)=gradn(10,1) !  b(4)/xb
             grdff(2,4,l)=gradn(11,1) !  c(4)/xb
             grdff(3,4,l)=gradn(12,1) !  d(4)/xb

          ELSE IF (LTYPE1.EQ.5) THEN                      !TRIANGULOS EN 3-D
             K2=0   !provi mientras se arregla la opcion de no todo en el plano
             IF (K2.EQ.0)THEN                                  !provi
               CALL BIBITRANG3D (X,Y,Z,KXX,L,NUMEL,NUMNP,LMXNDL,
     ;         BIBI,BTRA,IOEQT,HCAL,GRDFF,AREA,HBASE,
     ;         ACTH,IOFLLI,IDIMBB,IODIM)
               GOTO 100
             ENDIF
c             I3=KXX(3,L)  !below operations (tensor of fracture in the space)
c             X3=X(I3)     !are no operatives (yet)
c             Y3=Y(I3)
c             Z1=Z(I1)
c             Z2=Z(I2)
c             Z3=Z(I3)
c
c* .................. Components of jacobian matrix (global to local coordinates)
c             B(1)= X2 - X1
c             C(1)= X3 - X1
c             B(2)= Y2 - Y1
c             C(2)= Y3 - Y1
c             B(3)= Z2 - Z1
c             C(3)= Z3 - Z1
c
c* ........................ Covariant components and determinant of metric tensor
c             D(2)= B(1)*B(1) + B(2)*B(2) + B(3)*B(3)   ! G11
c             D(1)= C(1)*C(1) + C(2)*C(2) + C(3)*C(3)   ! G22 
c             D(3)= B(1)*C(1) + B(2)*C(2) + B(3)*C(3)   ! G12 
c             DETG= D(1)*D(2) - D(3)*D(3)
c
c* .................................... Contravariant components of metric tensor
c             D(1)= D(1)/DETG       ! G11
c             D(2)= D(2)/DETG       ! G22
c             D(3)= -D(3)/DETG      ! G12
c
c* ....................................................... Compute element's area
c             IF (IOEQT.NE.1) 
c     .       ACTH(L)=(BTRA(I1)+BTRA(I2)+BTRA(I3))/3.D+00
c             IF (IOFLLI.NE.0)THEN
c               HBASE(L)=(HCAL(I1)+HCAL(I2)+HCAL(I3))/3D0
c             ENDIF
c             IF (IOTRLI.NE.0)THEN
c               CBASE(L)=(CCAL(I1)+CCAL(I2)+CCAL(I3))/3D0
c             ENDIF
c
c             IF (DETG.LT.0) THEN
c                 CALL ERROR(' NEGATIVE DETG IN ELEMENT@',L,0,0,3)
c                 GOTO 100
c             END IF
c             AREA1=DSQRT(DETG)/2.D+00
c             IF (DABS(AREA1) .LT. 1.D-20) 
c     .         CALL ERROR(' AREA ELEMEN LOWER THAT 1.D-20 IN ELEMENT@',
c     .              L,0,0,3)
c             IF (AREA1 .LT. 0.D+00) 
c     .         CALL ERROR(' AREA ELEMEN NEGATIVE IN ELEMENT (INCORRECT 
c     .              NODES  NUMERATION ) IN ELEMENT@',L,0,0,3)
c             AREA(L)=AREA1
c
c* .......................................... Compute gradient of basis functions
c             DO I=1,3
c                X4 = D(1)*B(I) + D(3)*C(I)
c                Y4 = D(3)*B(I) + D(2)*C(I)
c                GRDFF(I,1,L) = -X4 -Y4
c                GRDFF(I,2,L) = X4
c                GRDFF(I,3,L) = Y4
c             END DO
c
c* ....................................................... Compute element matrix
c             BIBI(1,L)=  AREA1*GRDFF(1,1,L)*GRDFF(1,2,L)  !ojo los bibi(n,l) 
c             BIBI(2,L)=  AREA1*GRDFF(3,1,L)*GRDFF(3,2,L)  !n=2 8 y 14
c             BIBI(3,L)=  AREA1*GRDFF(2,1,L)*GRDFF(2,2,L)  !deben ser de Tzz
c             BIBI(4,L)= AREA1*(GRDFF(1,1,L)*GRDFF(2,2,L) +
c     .                         GRDFF(2,1,L)*GRDFF(1,2,L))
c             BIBI(5,L)= AREA1*(GRDFF(1,1,L)*GRDFF(3,2,L) +
c     .                         GRDFF(3,1,L)*GRDFF(1,2,L))
c             BIBI(6,L)= AREA1*(GRDFF(2,1,L)*GRDFF(3,2,L) +
c     .                         GRDFF(3,1,L)*GRDFF(2,2,L))
c             BIBI(7,L)=  AREA1*GRDFF(1,1,L)*GRDFF(1,3,L)
c             BIBI(8,L)=  AREA1*GRDFF(3,1,L)*GRDFF(3,3,L)
c             BIBI(9,L)=  AREA1*GRDFF(2,1,L)*GRDFF(2,3,L)
c             BIBI(10,L)= AREA1*(GRDFF(1,1,L)*GRDFF(2,3,L) +
c     .                         GRDFF(2,1,L)*GRDFF(1,3,L))
c             BIBI(11,L)= AREA1*(GRDFF(1,1,L)*GRDFF(3,3,L) +
c     .                         GRDFF(3,1,L)*GRDFF(1,3,L))
c             BIBI(12,L)= AREA1*(GRDFF(2,1,L)*GRDFF(3,3,L) +
c     .                         GRDFF(3,1,L)*GRDFF(2,3,L))
c             BIBI(13,L)=  AREA1*GRDFF(1,3,L)*GRDFF(1,2,L)
c             BIBI(14,L)=  AREA1*GRDFF(3,3,L)*GRDFF(3,2,L)
c             BIBI(15,L)=  AREA1*GRDFF(2,3,L)*GRDFF(2,2,L)
c             BIBI(16,L)= AREA1*(GRDFF(1,3,L)*GRDFF(2,2,L) +
c     .                         GRDFF(2,3,L)*GRDFF(1,2,L))
c             BIBI(17,L)= AREA1*(GRDFF(1,3,L)*GRDFF(3,2,L) +
c     .                         GRDFF(3,3,L)*GRDFF(1,2,L))
c             BIBI(18,L)= AREA1*(GRDFF(2,3,L)*GRDFF(3,2,L) +
c     .                         GRDFF(3,3,L)*GRDFF(2,2,L))
          ELSE IF (LTYPE1.EQ.6)THEN                         !TOBLERONE
            NVAL=6                              !NUMBER OF INTEGRATION POINTS
            I3=KXX(3,L)
            I4=KXX(4,L)
            I5=KXX(5,L)
            I6=KXX(6,L)
 
            X1=X(I1)
            Y1=Y(I1)
            Z1=Z(I1)
            X2=X(I2)
            Y2=Y(I2)
            Z2=Z(I2)
            X3=X(I3)
            Y3=Y(I3)
            Z3=Z(I3)
            X4=X(I4)
            Y4=Y(I4)
            Z4=Z(I4)
            X5=X(I5)
            Y5=Y(I5)
            Z5=Z(I5)
            X6=X(I6)
            Y6=Y(I6)
            Z6=Z(I6)
            CALL GRADN_PT_3D (X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,
     .      X4,Y4,Z4,X5,Y5,Z5,X6,Y6,Z6, 
     .      GRADN,PESO,NVAL)

        DO I=1,9
          DO J=1,15
            GNGN(I,J)=0D0
          END DO
        END DO

        N=6
        DO IVAL=1,NVAL
          K=1
          DO I=1,N-1
            DO J=I+1,N
             GNGN(1,K)=GNGN(1,K)+
     ;           GRADN(1+(I-1)*3,IVAL)*GRADN(1+(J-1)*3,IVAL)*PESO(IVAL)
             GNGN(2,K)=GNGN(2,K)+
     ;           GRADN(1+(I-1)*3,IVAL)*GRADN(2+(J-1)*3,IVAL)*PESO(IVAL)
             GNGN(3,K)=GNGN(3,K)+
     ;           GRADN(1+(I-1)*3,IVAL)*GRADN(3+(J-1)*3,IVAL)*PESO(IVAL)
             GNGN(4,K)=GNGN(4,K)+
     ;           GRADN(2+(I-1)*3,IVAL)*GRADN(1+(J-1)*3,IVAL)*PESO(IVAL)
             GNGN(5,K)=GNGN(5,K)+
     ;           GRADN(2+(I-1)*3,IVAL)*GRADN(2+(J-1)*3,IVAL)*PESO(IVAL)
             GNGN(6,K)=GNGN(6,K)+
     ;           GRADN(2+(I-1)*3,IVAL)*GRADN(3+(J-1)*3,IVAL)*PESO(IVAL)
             GNGN(7,K)=GNGN(7,K)+
     ;           GRADN(3+(I-1)*3,IVAL)*GRADN(1+(J-1)*3,IVAL)*PESO(IVAL)
             GNGN(8,K)=GNGN(8,K)+
     ;           GRADN(3+(I-1)*3,IVAL)*GRADN(2+(J-1)*3,IVAL)*PESO(IVAL)
             GNGN(9,K)=GNGN(9,K)+
     ;           GRADN(3+(I-1)*3,IVAL)*GRADN(3+(J-1)*3,IVAL)*PESO(IVAL)
             K=K+1
            END DO
          END DO

*         CREATES BIBI ARRAY
          IBIBI=-6
          DO ICON=1,15
            IBIBI=IBIBI+6
            BIBI(IBIBI+1,L)=GNGN(1,ICON)                 !XX
            BIBI(IBIBI+2,L)=GNGN(5,ICON)                 !YY
            BIBI(IBIBI+3,L)=GNGN(9,ICON)                 !ZZ
            BIBI(IBIBI+4,L)=GNGN(2,ICON)+GNGN(4,ICON)    !XY
            BIBI(IBIBI+5,L)=GNGN(3,ICON)+GNGN(7,ICON)    !XZ
            BIBI(IBIBI+6,L)=GNGN(6,ICON)+GNGN(8,ICON)    !YZ
          ENDDO
        END DO

*............................................... COMPUTES AVERAGE GRADIENTS
        DO I=1,3
         DO J=1,6
           GRDFF(I,J,L)=0.0
         END DO
        END DO

        DO IVAL=1,NVAL
            GRDFF(1,1,L)= GRDFF(1,1,L)+ GRADN(1,ival)               ! DN1/DX
            GRDFF(2,1,L)= GRDFF(2,1,L)+ GRADN(2,ival)               ! DN1/DY  
            GRDFF(3,1,L)= GRDFF(3,1,L)+ GRADN(3,ival)               ! DN1/DZ    
            GRDFF(1,2,L)= GRDFF(1,2,L)+ GRADN(4,ival)               ! DN2/DX    
            GRDFF(2,2,L)= GRDFF(2,2,L)+ GRADN(5,ival)               ! DN2/DY    
            GRDFF(3,2,L)= GRDFF(3,2,L)+ GRADN(6,ival)               ! ...  
            GRDFF(1,3,L)= GRDFF(1,3,L)+ GRADN(7,ival)
            GRDFF(2,3,L)= GRDFF(2,3,L)+ GRADN(8,ival)
            GRDFF(3,3,L)= GRDFF(3,3,L)+ GRADN(9,ival)
            GRDFF(1,4,L)= GRDFF(1,4,L)+ GRADN(10,ival)
            GRDFF(2,4,L)= GRDFF(2,4,L)+ GRADN(11,ival)
            GRDFF(3,4,L)= GRDFF(3,4,L)+ GRADN(12,ival)
            GRDFF(1,5,L)= GRDFF(1,5,L)+ GRADN(13,ival)
            GRDFF(2,5,L)= GRDFF(2,5,L)+ GRADN(14,ival)
            GRDFF(3,5,L)= GRDFF(3,5,L)+ GRADN(15,ival)
            GRDFF(1,6,L)= GRDFF(1,6,L)+ GRADN(16,ival)
            GRDFF(2,6,L)= GRDFF(2,6,L)+ GRADN(17,ival)
            GRDFF(3,6,L)= GRDFF(3,6,L)+ GRADN(18,ival)
        END DO

        DO I=1,3
         DO J=1,6
          GRDFF(I,J,L)=GRDFF(I,J,L)/DFLOAT(NVAL)
         END DO
        END DO


        ENDIF
 100   CONTINUE
       
* ............ Asigna coeficientes de la integral FIi.FIj en esquema consistente

***  Initialize to zero

       DO L=1,5
          DO J=1,4
             DO I=1,4
                CNST(I,J,L)=0D0
             END DO
          END DO
       END DO

***  LINEAR ELEMENT

       CNST(1,1,1)=1D0/3D0
       CNST(1,2,1)=1D0/6D0
       CNST(2,1,1)=1D0/6D0
       CNST(2,2,1)=1D0/3D0

***  TRIANGULAR ELEMENT IN 2-D

       CNST(1,1,2)=1D0/6D0
       CNST(1,2,2)=1D0/12D0
       CNST(1,3,2)=1D0/12D0
       CNST(2,1,2)=1D0/12D0
       CNST(2,2,2)=1D0/6D0
       CNST(2,3,2)=1D0/12D0
       CNST(3,1,2)=1D0/12D0
       CNST(3,2,2)=1D0/12D0
       CNST(3,3,2)=1D0/6D0

***  RECTANGULAR ELEMENT

       CNST(1,1,3)=1D0/9D0
       CNST(1,2,3)=1D0/18D0
       CNST(1,3,3)=1D0/36D0
       CNST(1,4,3)=1D0/18D0
       CNST(2,1,3)=1D0/18D0
       CNST(2,2,3)=1D0/9D0
       CNST(2,3,3)=1D0/18D0
       CNST(2,4,3)=1D0/36D0
       CNST(3,1,3)=1D0/36D0
       CNST(3,2,3)=1D0/18D0
       CNST(3,3,3)=1D0/9D0
       CNST(3,4,3)=1D0/18D0
       CNST(4,1,3)=1D0/18D0
       CNST(4,2,3)=1D0/36D0
       CNST(4,3,3)=1D0/18D0
       CNST(4,4,3)=1D0/9D0

***  TETRAHEDRON ELEMENT

       CNST(1,1,4)=1D0/10D0
       CNST(1,2,4)=1D0/20D0
       CNST(1,3,4)=1D0/20D0
       CNST(1,4,4)=1D0/20D0
       CNST(2,1,4)=1D0/20D0
       CNST(2,2,4)=1D0/10D0
       CNST(2,3,4)=1D0/20D0
       CNST(2,4,4)=1D0/20D0
       CNST(3,1,4)=1D0/20D0
       CNST(3,2,4)=1D0/20D0
       CNST(3,3,4)=1D0/10D0
       CNST(3,4,4)=1D0/20D0
       CNST(4,1,4)=1D0/20D0
       CNST(4,2,4)=1D0/20D0
       CNST(4,3,4)=1D0/20D0
       CNST(4,4,4)=1D0/10D0

***  TRIANGULAR ELEMENT IN 3-D

       CNST(1,1,5)=1D0/6D0
       CNST(1,2,5)=1D0/12D0
       CNST(1,3,5)=1D0/12D0
       CNST(2,1,5)=1D0/12D0
       CNST(2,2,5)=1D0/6D0
       CNST(2,3,5)=1D0/12D0
       CNST(3,1,5)=1D0/12D0
       CNST(3,2,5)=1D0/12D0
       CNST(3,3,5)=1D0/6D0

***  TOBLERONE ELEMENT 3-D

       CNST(1,1,6)=1D0/18
       CNST(1,2,6)=1D0/36
       CNST(1,3,6)=1D0/36
       CNST(1,4,6)=1D0/36
       CNST(1,5,6)=1D0/72
       CNST(1,6,6)=1D0/72
       CNST(2,1,6)=1D0/36
       CNST(2,2,6)=1D0/18
       CNST(2,3,6)=1D0/36
       CNST(2,4,6)=1D0/72
       CNST(2,5,6)=1D0/36
       CNST(2,6,6)=1D0/72
       CNST(3,1,6)=1D0/36
       CNST(3,2,6)=1D0/36
       CNST(3,3,6)=1D0/18
       CNST(3,4,6)=1D0/72
       CNST(3,5,6)=1D0/72
       CNST(3,6,6)=1D0/36
       CNST(4,1,6)=1D0/36
       CNST(4,2,6)=1D0/72
       CNST(4,3,6)=1D0/72
       CNST(4,4,6)=1D0/18
       CNST(4,5,6)=1D0/36
       CNST(4,6,6)=1D0/36
       CNST(5,1,6)=1D0/72
       CNST(5,2,6)=1D0/36
       CNST(5,3,6)=1D0/72
       CNST(5,4,6)=1D0/36
       CNST(5,5,6)=1D0/18
       CNST(5,6,6)=1D0/36
       CNST(6,1,6)=1D0/72
       CNST(6,2,6)=1D0/72
       CNST(6,3,6)=1D0/36
       CNST(6,4,6)=1D0/36
       CNST(6,5,6)=1D0/36
       CNST(6,6,6)=1D0/18

       IF (IERROR .GT. 0) STOP 'ERROR IN GRID DEFINITION ' 

*....................RECOVER HCAL AND CCAL VECTORS USED TEMPORARI FOR SOTORING
*.............................................. FLOOR'S HEAD AND CONCENTRATION
       IF (IOFLLI.NE.0)THEN
         DO I=1,NUMNP
           HCAL(I)=0D0
         END DO
       ENDIF
       IF (IOTRLI.NE.0)THEN
         DO I=1,NUMNP
           CCAL(I)=0D0
         END DO
       ENDIF

* ..................... CORRIGES THICKNESS VALUES FOR THREE-DIMENSIONAL ELEMENTS
       IF (IOEQT.NE.1)THEN
         DO L=1,NUMEL
           IF(LTYPE(L).EQ.4.OR.LTYPE(L).EQ.6)ACTH(L)=1D0
         END DO
       ENDIF
c       do l=1,5
c       do is=1,6
c         write(5,*)' el bibi del elto',l
c         write(5,'(6f10.4)')(bibi((is-1)*6+kj,l),kj=1,6)
c       end do
c       end do
C______________________________Computes the projection of gravity vector
C______________________________over all the elements.

       IF (IOFLLI.EQ.1 .OR. IODENS.EQ.1) THEN

         CALL GRAV_PROJ
     ;(LMXNDL   ,NUMEL    ,NUMNP    ,COORD    ,GRAV     ,GRAVEL
     ;,KXX      ,LDIM)

       END IF
C______________________________ Prel. calc. related to matrix diffusion

          IF (IOEQT.NE.1. AND .NTDMT.NE.0) THEN

C______________________________ Computation of volume of nodes

            CALL COMP_VOL_NOD 
     ;(LMXNDL   ,NUMEL    ,NUMNP    ,ACTH   ,AREA     ,KXX      
     ;,LDIM     ,LNNDEL   ,VOLNOD)  

C______________________________ Store volumes associated to nodes
C______________________________ belonging to matrix zones.

            CALL COMP_PRE_DMT (VOLNOD,NUMNP)

          END IF       ! Is Matrix Diffusion considered?

       RETURN
       END      

********************************************************************************

       SUBROUTINE BIBITRANG3D (X,Y,Z,KXX,L,NUMEL,NUMNP,LMXNDL,BIBI,BTRA,
     ; IOEQT,HCAL,GRDFF,AREA,HBASE,ACTH,IOFLLI,IDIMBB,IODIM)

********************************************************************************
***     COMPUTES THE BASIS FUNCTION FOR THE TRIANGLE IN THREE DIMENSIONS     ***
********************************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION KXX(LMXNDL,NUMEL),X(NUMNP),Y(NUMNP),Z(NUMNP),
     ;       GRDFF(IODIM,LMXNDL,NUMEL),AREA(NUMEL),ACTH(NUMEL),XG(3),
     ;       YG(3),ZG(3),VI(3),VJ(3),VK(3),XLOC(3),YLOC(3),HBASE(NUMEL),
     ;     BTRA(NUMNP),BIBI(IDIMBB,NUMEL),HCAL(NUMNP)

             I1=KXX(1,L)
             I2=KXX(2,L)
             I3=KXX(3,L)
             XG(1)=X(I1)
             XG(2)=X(I2)
             XG(3)=X(I3)
             YG(1)=Y(I1)
             YG(2)=Y(I2)
             YG(3)=Y(I3)
             ZG(1)=Z(I1)
             ZG(2)=Z(I2)
             ZG(3)=Z(I3)
             CALL LOC_COORDINADES (XG(1),XG(2),XG(3),YG(1),YG(2),YG(3),
     ;                           ZG(1),ZG(2),ZG(3),VMINPI,
     ;                           VMINPJ,VMINPK,VMAXPI,VMAXPJ,VMAXPK,
     ;                           XLOCCERO,YLOCCERO,ZLOCCERO)

             DO I=1,3
               VI(I)=XG(I)-XLOCCERO
               VJ(I)=YG(I)-YLOCCERO
               VK(I)=ZG(I)-ZLOCCERO
               XLOC(I)=VI(I)*VMINPI+VJ(I)*VMINPJ+VK(I)*VMINPK
               YLOC(I)=VI(I)*VMAXPI+VJ(I)*VMAXPJ+VK(I)*VMAXPK
             END DO

             X1=XLOC(1)
             Y1=YLOC(1)
             X2=XLOC(2)
             Y2=YLOC(2)
             X3=XLOC(3)
             Y3=YLOC(3)

             B1=Y2-Y3
             B2=Y3-Y1
             C1=X3-X2
             C2=X1-X3
             B3=Y1-Y2
             C3=X2-X1
             IF (IOEQT.NE.1) 
     .       ACTH(L)=(BTRA(I1)+BTRA(I2)+BTRA(I3))/3.D+00
             IF (IOFLLI.NE.0)THEN
               HBASE(L)=(HCAL(I1)+HCAL(I2)+HCAL(I3))/3D0
             ENDIF

             AREA2=AREA(L)*4.D+00

*** bibi for 3-D triangles treated in its 2-D plane (treatment temporal)

             BIBI(1,L)=B1*B2/ AREA2               !xx
             BIBI(2,L)=C1*C2/ AREA2             !yy
             BIBI(3,L)=(B2*C1+B1*C2) / AREA2    !xy
             BIBI(4,L)=B3*B1/ AREA2               !xx
             BIBI(5,L)=C3*C1/ AREA2             !yy
             BIBI(6,L)=(B1*C3+B3*C1) / AREA2    !xy
             BIBI(7,L)=B2*B3/ AREA2              !xx
             BIBI(8,L)=C2*C3/ AREA2            !yy
             BIBI(9,L)=(B3*C2+B2*C3) / AREA2   !xy


             GRDFF(1,1,L)=B1/AREA(L)/2
             GRDFF(1,2,L)=B2/AREA(L)/2
             GRDFF(1,3,L)=-(B1+B2)/AREA(L)/2
             GRDFF(2,1,L)=C1/AREA(L)/2
             GRDFF(2,2,L)=C2/AREA(L)/2
             GRDFF(2,3,L)=-(C1+C2)/AREA(L)/2
             RETURN
             END        


       SUBROUTINE LOC_COORDINADES (X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,VMINPI,
     ;                           VMINPJ,VMINPK,VMAXPI,VMAXPJ,VMAXPK,
     ;                           XLOCCERO,YLOCCERO,ZLOCCERO)
********************************************************************************
**** COMPUTES THE DIRECTION OF MAXIMUN SLOPE FOR THE PLANE OF A FRACTURE     ***
********************************************************************************
             IMPLICIT REAL*8 (A-H,O-Z)
            
             
* assembles vector [1-3]
       UI=X3-X1
       UJ=Y3-Y1
       UK=Z3-Z1
       XMODU=DSQRT(UI**2+UJ**2+UK**2)

* assemles vector [1-2] (used as local x direction)
       VI=X2-X1
       VJ=Y2-Y1
       VK=Z2-Z1
       XMODV=DSQRT(VI**2+VJ**2+VK**2)

* assembles a vector perpendicular to the plane of the triangle, with
* positive k (asociated to z global coorddinate) component
       XNI=UJ*VK-UK*VJ
       XNJ=-UI*VK+UK*VI
       XNK=UI*VJ-UJ*VI
       XMODN=DSQRT(XNI**2+XNJ**2+XNK**2)
       CTE=1D0
       IF(XNK.LT.0D0)CTE=-1D0
       XNI=CTE*XNI/XMODN
       XNJ=CTE*XNJ/XMODN
       XNK=CTE*XNK/XMODN
       IF (DABS(XNK).LT.1E-10) THEN   !vertical plane
         IF (DABS(XNI).GT.1E-10)THEN
           IF (XNI.LT.0D0)THEN
             XNI=-1*XNI
             XNJ=-1*XNJ
             XNK=-1*XNK
           ENDIF
           VMAXPI=0D0
           VMAXPJ=0D0
           VMAXPK=1D0
           VMINPI=-XNJ
           VMINPJ=XNI
           VMINPK=0D0
         ELSE
           IF (XNJ.LT.0)THEN
             XNI=-1*XNI
             XNJ=-1*XNJ
             XNK=-1*XNK
           ENDIF
           VMAXPI=0D0
           VMAXPJ=0D0
           VMAXPK=1D0
           VMINPI=-XNJ
           VMINPJ=0D0
           VMINPK=0D0
         ENDIF
         XLOCCERO=X1
         YLOCCERO=Y1
         ZLOCCERO=Z1
       ELSE IF((DABS(XNI)+DABS(XNJ)).LT.1E-10)THEN  !horizontal plane
         VMAXPI=0D0
         VMAXPJ=1D0
         VMAXPK=0D0
         VMINPI=1D0
         VMINPJ=0D0
         VMINPK=0D0
         XLOCCERO=0D0
         YLOCCERO=0D0
         ZLOCCERO=0D0
       ELSE                                         !general case
* computes the coefficients of the equation of the plane 
         FACX=XNI
         FACY=XNJ
         FACZ=XNK
         TINDEP=-(XNI*X1+XNJ*Y1+XNK*Z1)

* computes a vector in direction of minimun slope
         IF (DABS(XNI).LT.1E-08)THEN  !plane paralell to axis x
           ZMIN=Z1
           YMIN=Y1
           XMIN=(DABS(X1)+DABS(X2)+DABS(X3))             !arbitrario
           VHORIZI=XMIN-X1
           IF(DABS(VHORIZI).LT.1E-08)VHORIZI=10.0
           VHORIZJ=0D0
           VHORIZK=0D0
c         ELSE
c           YMIN=Y2
c           XMIN=((Z1-ZMIN)*XNK+(Y1-YMIN)*XNJ)/XNI+X1
c           VHORIZI=XMIN-X1
c           VHORIZJ=YMIN-Y1
c         ENDIF
          ELSE
            XPLANO1=X1+X2+X3                                   !arbitrario
            YPLANO1=Y1+Y2+Y3                                   !arbitrario
            ZPLANO1=(-TINDEP-FACX*XPLANO1-FACY*YPLANO1)/FACZ
            YPLANO2=-YPLANO1                                   !arbitrario
            IF (YPLANO1.EQ.0D0) YPLANO2=10.00                  !arbitrario
            ZPLANO2=ZPLANO1                                !cond. horizontalidad
            XPLANO2=(-TINDEP-FACY*YPLANO2-FACZ*ZPLANO2)/FACX
            VHORIZI=XPLANO2-XPLANO1
            VHORIZJ=YPLANO2-YPLANO1
            VHORIZK=0D0
          ENDIF

* computes a vector in direction of maximun xlope
         VMAXPI=VHORIZJ*XNK-VHORIZK*XNJ
         VMAXPJ=-VHORIZI*XNK+VHORIZK*XNI
         VMAXPK=VHORIZI*XNJ-VHORIZJ*XNI
         CTE=1D0
         IF(VMAXPK.LT.0D0)CTE=-1D0
         XMODVMAXP=DSQRT(VMAXPI*VMAXPI+VMAXPJ*VMAXPJ+VMAXPK*VMAXPK)
         VMAXPI=VMAXPI*CTE/XMODVMAXP
         VMAXPJ=VMAXPJ*CTE/XMODVMAXP
         VMAXPK=VMAXPK*CTE/XMODVMAXP
       
* checks or corriges the direction of the minimun slope vector
         VMINPI=VMAXPJ*XNK-VMAXPK*XNJ
         VMINPJ=-VMAXPI*XNK+VMAXPK*XNI
         VMINPK=VMAXPI*XNJ-VMAXPJ*XNI
         XMODVMINP=DSQRT(VMINPI*VMINPI+VMINPJ*VMINPJ+VMINPK*VMINPK)
         VMINPI=VMINPI/XMODVMINP
         VMINPJ=VMINPJ/XMODVMINP
         VMINPK=VMINPK/XMODVMINP

         IF (DABS(XNI).GT.1.D-08 .AND. DABS(VHORIZI).GT.1.D-08 .AND.
     ;       DABS(VHORIZJ).GT.1.D-08) THEN

           IF (DABS(VMINPK).GT.1.D-05 .OR.
     ;         ((VMINPJ/VHORIZJ)-(VMINPI/VHORIZI)).GT.1.D-05) THEN

             WRITE(25,*)
             WRITE(25,*)'  ERROR DE PROGRAMACION EN LOC_COORDINADES 1'
             WRITE(25,*)'  XNI:',xni
             WRITE(25,*)'  VHORIZI:',vhorizi
             WRITE(25,*)'  VHORIZJ:',vhorizj
             WRITE(25,*)'  VMINPK:',vminpk
             WRITE(25,*)'  ((VMINPJ/VHORIZJ)-(VMINPI/VHORIZI)):',
     ;                     ((vminpj/vhorizj)-(vminpi/vhorizi))
             STOP 'LOC_COORDINADES 1'

           ENDIF

         ELSE

           IF (DABS(VMINPK).GT.1.D-05)THEN

             WRITE(25,*)
             WRITE(25,*)'  ERROR DE PROGRAMACION EN LOC_COORDINADES 2'
             WRITE(25,*)'  XNI:',xni
             WRITE(25,*)'  VHORIZI:',vhorizi
             WRITE(25,*)'  VHORIZJ:',vhorizj
             WRITE(25,*)'  VMINPK:',vminpk
             STOP ' LOC_COORDINADES 2'

           ENDIF

         ENDIF
       
         XLOCCERO=0D0
         YLOCCERO=0D0
         ZLOCCERO=-((XLOCCERO-X1)*XNI+(YLOCCERO-Y1)*XNJ)/XNK+Z1
       ENDIF
       RETURN
       END
       subroutine gradn_pt_3d (x1,y1,z1,x2,y2,z2,x3,y3,z3,
     .                         x4,y4,z4,x5,y5,z5,x6,y6,z6, 
     .                         gradn,peso,nval)

*................................................ triangular prism 
*                                                 (2, 3 or 6 integration points)
       implicit real*8 (a-h,o-z)
       dimension gradn(24,6),peso(6)


***************************************************************************
* Triangular prism (toblerone) with 2,3 or 6 integration points. Any form is
* permitted. Triangle 1,2,3 defines orientation towards triangle 4,5,6
*************************************************************************** 


          chi1=-1.d0
          chi2=+1.d0
 
          if (nval.eq.2) then 
           chi01=-.577350269189626             ! 2 integration points
           chi02=+.577350269189626
           al101= 0.333333333333333
           al201= 0.333333333333333
           al301= 0.333333333333333
           al102= 0.333333333333333
           al202= 0.333333333333333
           al302= 0.333333333333333
           weight01=1.0d0
           weight02=1.0d0
          else if (nval.eq.3) then
           chi01=-.774596669241483             ! 3 integration points
           chi02=+.000000000000000
           chi03=+.774596669241483
           al101= 0.333333333333333
           al201= 0.333333333333333
           al301= 0.333333333333333
           al102= 0.333333333333333
           al202= 0.333333333333333
           al302= 0.333333333333333
           al103= 0.333333333333333
           al203= 0.333333333333333
           al303= 0.333333333333333
           weight01=.5555555555555555
           weight02=.8888888888888888
           weight03=.5555555555555555
          else if (nval.eq.6) then 
           chi01=-.577350269189626             ! 6 integration points
           chi02=-.577350269189626             
           chi03=-.577350269189626             
           chi04=+.577350269189626
           chi05=+.577350269189626
           chi06=+.577350269189626
           al101= 0.5
           al201= 0.5
           al301= 0.0
           al102= 0.0
           al202= 0.5
           al302= 0.5
           al103= 0.5
           al203= 0.0
           al303= 0.5
           weight01=1.0d0/3.0
           weight02=1.0d0/3.0
           weight03=1.0d0/3.0
          end if

                               ! N1=1/2 (1+chi*chi1) L1
                               ! N2=1/2 (1+chi*chi1) L2
                               ! N3=1/2 (1+chi*chi1) L3
                               ! N4=1/2 (1+chi*chi2) L1
                               ! N5=1/2 (1+chi*chi2) L2
                               ! N6=1/2 (1+chi*chi2) L3
                               ! constraint: L3 = 1-L1-L2
       


       do ival=1,nval

        if (ival.eq.1) then 
           chi=chi01
           al1=al101
           al2=al201
           al3=al301
           weight=weight01
        else if (ival.eq.2) then
           chi=chi02
           al1=al102
           al2=al202
           al3=al302
           weight=weight02
        else if (ival.eq.3) then
           chi=chi03
           al1=al103
           al2=al203
           al3=al303
           weight=weight03
        else if (ival.eq.4) then 
           chi=chi04
           al1=al101
           al2=al201
           al3=al301
           weight=weight01
        else if (ival.eq.5) then
           chi=chi05
           al1=al102
           al2=al202
           al3=al302
           weight=weight02
        else if (ival.eq.6) then
           chi=chi06
           al1=al103
           al2=al203
           al3=al303
           weight=weight03
        end if                   

                                      ! local variable derivatives at gaus point
        dn1dchi= 0.5d0*al1*chi1
        dn1dl1 = 0.5d0*(1.0+chi*chi1)
        dn1dl2 = 0.d0

        dn2dchi= 0.5d0*al2*chi1
        dn2dl1 = 0.d0
        dn2dl2 = 0.5d0*(1.0+chi*chi1)

        dn3dchi=  0.5d0*al3*chi1
        dn3dl1 = -0.5d0*(1.0+chi*chi1)
        dn3dl2 = -0.5d0*(1.0+chi*chi1)

        dn4dchi= 0.5d0*al1*chi2
        dn4dl1 = 0.5d0*(1.0+chi*chi2)
        dn4dl2 = 0.d0

        dn5dchi= 0.5d0*al2*chi2
        dn5dl1 = 0.d0
        dn5dl2 = 0.5d0*(1.0+chi*chi2)

        dn6dchi=  0.5d0*al3*chi2
        dn6dl1 = -0.5d0*(1.0+chi*chi2)
        dn6dl2 = -0.5d0*(1.0+chi*chi2)
        


        a11= x1*dn1dchi+x2*dn2dchi+x3*dn3dchi+                 ! dx/dchi
     .       x4*dn4dchi+x5*dn5dchi+x6*dn6dchi
        a12= y1*dn1dchi+y2*dn2dchi+y3*dn3dchi+                 ! dy/dchi
     .       y4*dn4dchi+y5*dn5dchi+y6*dn6dchi
        a13= z1*dn1dchi+z2*dn2dchi+z3*dn3dchi+                 ! dz/dchi
     .       z4*dn4dchi+z5*dn5dchi+z6*dn6dchi

        a21= x1*dn1dl1 +x2*dn2dl1 +x3*dn3dl1+                  ! dx/dl1 
     .       x4*dn4dl1 +x5*dn5dl1 +x6*dn6dl1 
        a22= y1*dn1dl1 +y2*dn2dl1 +y3*dn3dl1+                  ! dy/dl1 
     .       y4*dn4dl1 +y5*dn5dl1 +y6*dn6dl1 
        a23= z1*dn1dl1 +z2*dn2dl1 +z3*dn3dl1+                  ! dz/dl1 
     .       z4*dn4dl1 +z5*dn5dl1 +z6*dn6dl1 

        a31= x1*dn1dl2 +x2*dn2dl2 +x3*dn3dl2+                  ! dx/dl2 
     .       x4*dn4dl2 +x5*dn5dl2 +x6*dn6dl2 
        a32= y1*dn1dl2 +y2*dn2dl2 +y3*dn3dl2+                  ! dy/dl2 
     .       y4*dn4dl2 +y5*dn5dl2 +y6*dn6dl2 
        a33= z1*dn1dl2 +z2*dn2dl2 +z3*dn3dl2+                  ! dz/dl2 
     .       z4*dn4dl2 +z5*dn5dl2 +z6*dn6dl2 

        
        det=a11*a22*a33+a12*a23*a31+a21*a32*a13
     .     -a13*a22*a31-a12*a21*a33-a23*a32*a11

        b11= (a22*a33-a23*a32)/det
        b12=-(a12*a33-a13*a32)/det 
        b13= (a12*a23-a13*a22)/det

        b21=-(a21*a33-a23*a31)/det
        b22= (a11*a33-a13*a31)/det
        b23=-(a11*a23-a13*a21)/det

        b31= (a21*a32-a22*a31)/det
        b32=-(a11*a32-a12*a31)/det
        b33= (a11*a22-a12*a21)/det

        peso(ival)=det/2.0*weight
             
        gradn(1,ival)=  dn1dchi*b11+dn1dl1 *b12+dn1dl2 *b13   ! dN1/dx
        gradn(2,ival)=  dn1dchi*b21+dn1dl1 *b22+dn1dl2 *b23   ! dN1/dy
        gradn(3,ival)=  dn1dchi*b31+dn1dl1 *b32+dn1dl2 *b33   ! dN1/dz

        gradn(4,ival)=  dn2dchi*b11+dn2dl1 *b12+dn2dl2 *b13   ! dN2/dx
        gradn(5,ival)=  dn2dchi*b21+dn2dl1 *b22+dn2dl2 *b23   ! dN2/dy
        gradn(6,ival)=  dn2dchi*b31+dn2dl1 *b32+dn2dl2 *b33   ! dN2/dz

        gradn(7,ival)=  dn3dchi*b11+dn3dl1 *b12+dn3dl2 *b13   ! dN3/dx
        gradn(8,ival)=  dn3dchi*b21+dn3dl1 *b22+dn3dl2 *b23   ! dN3/dy
        gradn(9,ival)=  dn3dchi*b31+dn3dl1 *b32+dn3dl2 *b33   ! dN3/dz

        gradn(10,ival)= dn4dchi*b11+dn4dl1 *b12+dn4dl2 *b13   ! dN4/dx
        gradn(11,ival)= dn4dchi*b21+dn4dl1 *b22+dn4dl2 *b23   ! dN4/dy
        gradn(12,ival)= dn4dchi*b31+dn4dl1 *b32+dn4dl2 *b33   ! dN4/dz

        gradn(13,ival)= dn5dchi*b11+dn5dl1 *b12+dn5dl2 *b13   ! dN5/dx
        gradn(14,ival)= dn5dchi*b21+dn5dl1 *b22+dn5dl2 *b23   ! dN5/dy
        gradn(15,ival)= dn5dchi*b31+dn5dl1 *b32+dn5dl2 *b33   ! dN5/dz

        gradn(16,ival)= dn6dchi*b11+dn6dl1 *b12+dn6dl2 *b13   ! dN6/dx
        gradn(17,ival)= dn6dchi*b21+dn6dl1 *b22+dn6dl2 *b23   ! dN6/dy
        gradn(18,ival)= dn6dchi*b31+dn6dl1 *b32+dn6dl2 *b33   ! dN6/dz

        gradn(19,ival)= 0.5*(1+chi*chi1)*al1
        gradn(20,ival)= 0.5*(1+chi*chi1)*al2
        gradn(21,ival)= 0.5*(1+chi*chi1)*al3
        gradn(22,ival)= 0.5*(1+chi*chi2)*al1
        gradn(23,ival)= 0.5*(1+chi*chi2)*al2
        gradn(24,ival)= 0.5*(1+chi*chi2)*al3


       end do

       return
       end



* desde aqui subrutinas de tetraedro OCD (provisional)   !provisional oct-30-96

       SUBROUTINE TETRAH_CONDUCTANCE_OCD_MATRIX 
     ; (BIBI,X1,X2,X3,X4,Y1,Y2,Y3,Y4,
     ; Z1,Z2,Z3,Z4,L,IDIMBB,NUMEL,LXTRA,IOCD)

********************************************************************************
* COMPUTES THE "ORTHOGONAL COLOCATION DOMAINN (OCD)" ij ELEMENTAL COEFFICIENTS,
* OF THE CONDUCTANCE FLOW MATRIX, IN TETRAHEDRON ELEMENTS
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION BIBI(IDIMBB,NUMEL),LXTRA(NUMEL)

       NZ=LXTRA(L)
       REWIND (IOCD)
       READ(IOCD,*)
 1000  READ (IOCD,*,END=1100)NZ1,DISTX,DISTY,DISTZ
       IF (NZ1.EQ.NZ) THEN
         CALL DISTORSION (DISTX,DISTY,DISTZ,X1,Y1,Z1)
         CALL DISTORSION (DISTX,DISTY,DISTZ,X2,Y2,Z2)
         CALL DISTORSION (DISTX,DISTY,DISTZ,X3,Y3,Z3)
         CALL DISTORSION (DISTX,DISTY,DISTZ,X4,Y4,Z4)
         GOTO 1100
       ENDIF
       GOTO 1000

 1100  CONTINUE

* CALCULATE COMPONENTS AND MODULE OF THE VECTORS FORMED BY THE NODES
* OF THE TETRAHEDRON.

       CALL VECTOR_I_J (R12,DELX12,DELY12,DELZ12,X1,Y1,Z1,X2,Y2,Z2)

       CALL VECTOR_I_J (R13,DELX13,DELY13,DELZ13,X1,Y1,Z1,X3,Y3,Z3)

       CALL VECTOR_I_J (R14,DELX14,DELY14,DELZ14,X1,Y1,Z1,X4,Y4,Z4)

       CALL VECTOR_I_J (R23,DELX23,DELY23,DELZ23,X2,Y2,Z2,X3,Y3,Z3)

       CALL VECTOR_I_J (R24,DELX24,DELY24,DELZ24,X2,Y2,Z2,X4,Y4,Z4)

       CALL VECTOR_I_J (R34,DELX34,DELY34,DELZ34,X3,Y3,Z3,X4,Y4,Z4)

c       write(5,*)' componentes y coordenadas del 1 y 4=',
c     ; r14,delx14,dely14,delz14,x1,y1,z1,x4,y4,z4

* COMPUTE THE "AREA VECTOR" COMPONENTS FOR EACH NODE OF THE TETRAHEDRON.
* "AREA VECTOR" OF A GENERIC NODE i IS DEFINED AS THE PERPENDICULAR
* VECTOR TO THE PLANE CONTAINING THE TRIANGLE IN FRONT OF THE NODE i.
* ALSO, THE MODULE OF THE "AREA VECTOR" IS EQUAL TO THE SURFACE OF 
* THE TRIANGLE REFERED.

       CALL VECTORS_PRODUCT (DELX24,DELY24,DELZ24,DELX23,DELY23,DELZ23,
     ;          XMOD243,AREA1X,AREA1Y,AREA1Z,SCALAR243,AREA1)

       CALL VECTORS_PRODUCT (DELX13,DELY13,DELZ13,DELX14,DELY14,DELZ14,
     ;          XMOD134,AREA2X,AREA2Y,AREA2Z,SCALAR134,AREA2)

       CALL VECTORS_PRODUCT (DELX14,DELY14,DELZ14,DELX12,DELY12,DELZ12,
     ;          XMOD142,AREA3X,AREA3Y,AREA3Z,SCALAR142,AREA3)

       CALL VECTORS_PRODUCT (DELX12,DELY12,DELZ12,DELX13,DELY13,DELZ13,
     ;          XMOD123,AREA4X,AREA4Y,AREA4Z,SCALAR123,AREA4)


* CALCULATE TERMS REQUIRED FOR COMPUTING THE VOLUME OF THE ELEMENT

       CALL VECTORS_PRODUCT 
     ; (AREA4X,AREA4Y,AREA4Z,DELX14,DELY14,DELZ14,XMOD1234,
     ; VCTLX12314,VCTLY12314,VCTLZ12314,SCALAR12314,AREA12314)

       VOLUME=SCALAR12314/6D0

* CALCULATE TERMS REQIRED FOR COMPUTING TERMS OF THE CONDUCTANCE MATRIX.
* THEY ARE NAMED SCALAR****. THEIR MATHEMATICAL FORM IS    -->    -->
*                                                           r_ij . r_kl
*  in this generic case, "****" is "IJKL"

       CALL VECTORS_PRODUCT 
     ; (DELX13,DELY13,DELZ13,DELX23,DELY23,DELZ23,XMOD1323,
     ; VCTLX1323,VCTLY1323,VCTLZ1323,SCALAR1323,AREA1323)

       CALL VECTORS_PRODUCT 
     ; (DELX14,DELY14,DELZ14,DELX24,DELY24,DELZ24,XMOD1424,
     ; VCTLX1424,VCTLY1424,VCTLZ1424,SCALAR1424,AREA1424)

       CALL VECTORS_PRODUCT 
     ; (DELX12,DELY12,DELZ12,-DELX23,-DELY23,-DELZ23,XMOD1232,
     ; VCTLX1232,VCTLY1232,VCTLZ1232,SCALAR1232,AREA1232)

       CALL VECTORS_PRODUCT 
     ; (DELX14,DELY14,DELZ14,DELX34,DELY34,DELZ34,XMOD1434,
     ; VCTLX1434,VCTLY1434,VCTLZ1434,SCALAR1434,AREA1434)

       CALL VECTORS_PRODUCT 
     ; (DELX12,DELY12,DELZ12,-DELX24,-DELY24,-DELZ24,XMOD1242,
     ; VCTLX1242,VCTLY1242,VCTLZ1242,SCALAR1242,AREA1242)

       CALL VECTORS_PRODUCT 
     ; (DELX13,DELY13,DELZ13,-DELX34,-DELY34,-DELZ34,XMOD1343,
     ; VCTLX1343,VCTLY1343,VCTLZ1343,SCALAR1343,AREA1343)

       CALL VECTORS_PRODUCT 
     ; (-DELX12,-DELY12,-DELZ12,-DELX13,-DELY13,-DELZ13,XMOD2131,
     ; VCTLX2131,VCTLY2131,VCTLZ2131,SCALAR2131,AREA2131)

       CALL VECTORS_PRODUCT 
     ; (DELX24,DELY24,DELZ24,DELX34,DELY34,DELZ34,XMOD2434,
     ; VCTLX2434,VCTLY2434,VCTLZ2434,SCALAR2434,AREA2434)

       CALL VECTORS_PRODUCT 
     ; (-DELX12,-DELY12,-DELZ12,-DELX14,-DELY14,-DELZ14,XMOD2141,
     ; VCTLX2141,VCTLY2141,VCTLZ2141,SCALAR2141,AREA2141)

       CALL VECTORS_PRODUCT 
     ; (DELX23,DELY23,DELZ23,-DELX34,-DELY34,-DELZ34,XMOD2343,
     ; VCTLX2343,VCTLY2343,VCTLZ2343,SCALAR2343,AREA2343)

       CALL VECTORS_PRODUCT 
     ; (-DELX13,-DELY13,-DELZ13,-DELX14,-DELY14,-DELZ14,XMOD3141,
     ; VCTLX3141,VCTLY3141,VCTLZ3141,SCALAR3141,AREA3141)

       CALL VECTORS_PRODUCT 
     ; (-DELX23,-DELY23,-DELZ23,-DELX24,-DELY24,-DELZ24,XMOD3242,
     ; VCTLX3242,VCTLY3242,VCTLZ3242,SCALAR3242,AREA3242)


* CALCULATE TERMS REQIRED FOR COMPUTING TERMS OF THE CONDUCTANCE MATRIX.
* THEY ARE NAMED "SCALAR****". THEIR MATHEMATICAL FORM IS:   -->   -->
*                                                           A_i . A_j
*  in this generic case, "****" is "AIAJ"


       CALL VECTORS_PRODUCT
     ; (AREA1X,AREA1Y,AREA1Z,AREA1X,AREA1Y,AREA1Z,XMODA1A1,
     ; VCTLXA1A1,VCTLYA1A1,VCTLZA1A1,SCALARA1A1,AREAA1A1)

       CALL VECTORS_PRODUCT
     ; (AREA2X,AREA2Y,AREA2Z,AREA2X,AREA2Y,AREA2Z,XMODA2A2,
     ; VCTLXA2A2,VCTLYA2A2,VCTLZA2A2,SCALARA2A2,AREAA2A2)

       CALL VECTORS_PRODUCT
     ; (AREA3X,AREA3Y,AREA3Z,AREA3X,AREA3Y,AREA3Z,XMODA3A3,
     ; VCTLAX3A3,VCTLYA3A3,VCTLZA3A3,SCALARA3A3,AREAA3A3)

       CALL VECTORS_PRODUCT
     ; (AREA4X,AREA4Y,AREA4Z,AREA4X,AREA4Y,AREA4Z,XMODA4A4,
     ; VCTLXA4A4,VCTLYA4A4,VCTLZA4A4,SCALARA4A4,AREAA4A4)

       CALL VECTORS_PRODUCT
     ; (AREA1X,AREA1Y,AREA1Z,AREA2X,AREA2Y,AREA2Z,XMODA1A2,
     ; VCTLXA1A2,VCTLYA1A2,VCTLZA1A2,SCALARA1A2,AREAA1A2)

       CALL VECTORS_PRODUCT
     ; (AREA1X,AREA1Y,AREA1Z,AREA3X,AREA3Y,AREA3Z,XMODA1A3,
     ; VCTLXA1A3,VCTLYA1A3,VCTLZA1A3,SCALARA1A3,AREAA1A3)

       CALL VECTORS_PRODUCT
     ; (AREA1X,AREA1Y,AREA1Z,AREA4X,AREA4Y,AREA4Z,XMODA1A4,
     ; VCTLXA1A4,VCTLYA1A4,VCTLZA1A4,SCALARA1A4,AREAA1A4)

       CALL VECTORS_PRODUCT
     ; (AREA2X,AREA2Y,AREA2Z,AREA3X,AREA3Y,AREA3Z,XMODA2A3,
     ; VCTLXA2A3,VCTLYA2A3,VCTLZA2A3,SCALARA2A3,AREAA2A3)

       CALL VECTORS_PRODUCT
     ; (AREA2X,AREA2Y,AREA2Z,AREA4X,AREA4Y,AREA4Z,XMODA2A4,
     ; VCTLXA2A4,VCTLYA2A4,VCTLZA2A4,SCALARA2A4,AREAA2A4)

       CALL VECTORS_PRODUCT
     ; (AREA3X,AREA3Y,AREA3Z,AREA4X,AREA4Y,AREA4Z,XMODA3A4,
     ; VCTLXA3A4,VCTLYA3A4,VCTLZA3A4,SCALARA3A4,AREAA3A4)

* .................... COMPUTE EVERY ij TERM OF THE (OCD) CONDUCTANCE MATRIX

       CALL ISOTR_STIFF_TERM_OCD
     ; (SCALAR1323,SCALAR1424,SCALARA3A4,SCALARA4A4,SCALARA3A3,TERM12,
     ; VOLUME)

       CALL ISOTR_STIFF_TERM_OCD
     ; (SCALAR1232,SCALAR1434,SCALARA2A4,SCALARA4A4,SCALARA2A2,TERM13,
     ; VOLUME)

       CALL ISOTR_STIFF_TERM_OCD
     ; (SCALAR1242,SCALAR1343,SCALARA2A3,SCALARA3A3,SCALARA2A2,TERM14,
     ; VOLUME)

       CALL ISOTR_STIFF_TERM_OCD
     ; (SCALAR2131,SCALAR2434,SCALARA1A4,SCALARA4A4,SCALARA1A1,TERM23,
     ; VOLUME)
      
       CALL ISOTR_STIFF_TERM_OCD
     ; (SCALAR2141,SCALAR2343,SCALARA1A3,SCALARA3A3,SCALARA1A1,TERM24,
     ; VOLUME)

       CALL ISOTR_STIFF_TERM_OCD
     ; (SCALAR3141,SCALAR3242,SCALARA1A2,SCALARA2A2,SCALARA1A1,TERM34,
     ; VOLUME)

* ... STORAGE COMPONENTS OF THE CONDUCTANCE MATRIX IN THE "BIBI" MATRIX

       CALL ASSIGN_BIBI (BIBI,IDIMBB,NUMEL,L,3,TERM12)

       CALL ASSIGN_BIBI (BIBI,IDIMBB,NUMEL,L,9,TERM13)

       CALL ASSIGN_BIBI (BIBI,IDIMBB,NUMEL,L,15,TERM14)

       CALL ASSIGN_BIBI (BIBI,IDIMBB,NUMEL,L,21,TERM23)

       CALL ASSIGN_BIBI (BIBI,IDIMBB,NUMEL,L,27,TERM24)

       CALL ASSIGN_BIBI (BIBI,IDIMBB,NUMEL,L,33,TERM34)

       RETURN
       END

       SUBROUTINE VECTOR_I_J 
     ; (RIJ,DELXIJ,DELYIJ,DELZIJ,XI,YI,ZI,XJ,YJ,ZJ)

********************************************************************************
* COMPUTE THE DIRECTIONAL COMPONENTS (DEL*IJ) OF THE VECTOR FORMED BY 
*            NODES I AND J, AND ITS MODULE (RIJ)
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       DELXIJ=XJ-XI
       DELYIJ=YJ-YI
       DELZIJ=ZJ-ZI
       RIJ=(DELXIJ**2+DELYIJ**2+DELZIJ**2)**(0.5)
 
c       write(5,*)' ****************************** todos ',
c     ; xj,yj,zj,xi,yi,zi,delxij,delyij,delzij,rij

       RETURN
       END


       SUBROUTINE VECTORS_PRODUCT 
     ; (DELXV1,DELYV1,DELZV1,DELXV2,DELYV2,DELZV2,XMODV1V2,
     ; VCTLV1V2X,VCTLV1V2Y,VCTLV1V2Z,SCALARV1V2,AREAV1V2)

********************************************************************************
* COMPUTE THE COMPONENTS OF THE VECTORIAL
* PRODUCT BETWEEN TWO VECTORS. ALSO, THE SCALAR PRODUCT AND THE AREA
* BETWEEN THEY ARE CALCULATED
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       VCTLV1V2X=DELYV1*DELZV2-DELZV1*DELYV2
       VCTLV1V2Y=DELZV1*DELXV2-DELXV1*DELZV2
       VCTLV1V2Z=DELXV1*DELYV2-DELYV1*DELXV2

       XMODV1V2=DSQRT(VCTLV1V2X**2+VCTLV1V2Y**2+VCTLV1V2Z**2)
       SCALARV1V2=DELXV1*DELXV2+DELYV1*DELYV2+DELZV1*DELZV2
       AREAV1V2=XMODV1V2/2D0
 
       RETURN
       END
       

       SUBROUTINE ISOTR_STIFF_TERM_OCD
     ; (PROD1,PROD2,AREAS1,AREAS2,AREAS3,TERM,VOLUME)

********************************************************************************
* CALCULATE THE GENERIC TERM OF THE OCD ELEMENT CONDUCTANCE MATRIX 
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       TERM= -1/(48D0*VOLUME)/3D0*
     ; ( 2*PROD1*PROD2+ AREAS1* ( PROD1**2/AREAS2 + PROD2**2/AREAS3) )

       RETURN 
       END       


       SUBROUTINE ASSIGN_BIBI (BIBI,IDIMBB,NUMEL,L,NSUP,TERM)

********************************************************************************
** STORAGE THE CONDUCTANCE MATRIX ELEMENT TERMS IN THE BIBI MATRIX 
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION BIBI(IDIMBB,NUMEL)

       BIBI(NSUP-2,L)=TERM 
       BIBI(NSUP-1,L)=TERM 
       BIBI(NSUP,L)=TERM 

       RETURN
       END

       SUBROUTINE DISTORSION (DISTX,DISTY,DISTZ,X,Y,Z)
       IMPLICIT REAL*8 (A-H,O-Z)

********************************************************************************
***   AFFECT THE COORDINATES OF A POINT BY A DIAGONAL DISTORSION TENSOR      ***
********************************************************************************

       X=X*DISTX
       Y=Y*DISTY
       Z=Z*DISTZ

       RETURN
       END
