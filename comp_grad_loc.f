      SUBROUTINE COMP_GRAD_LOC
     &          (GRADLOC  ,IODIM    ,LMXNDL   ,LTYP     ,NP       ,XLOC
     &          ,YLOC     ,ZLOC)
     
****************************************************************************
*     PURPOSE
*         To compute the local gradient of shape functions im a certain
*         number of points.
*
*     DESCRIPTION
*         Computes the local gradient in the given points acoording to the
*         element type. The points coordinates are XLOC, YLOC and ZLOC.
*         Local gradient is stored in the variable GRADLOC.
*
*         GRADLOC: Contains the gradient of shape functions
*                  evaluated in the given points. It is a matrix of
*                  (Num. of shape funcs. x Num. of given points x
*                  Num. of dims.) elements.
*                  The storage is similar to the one used for full matrices.
*                  GRADLOC(i,j,k) contains the gradient of shape function j
*                  in direction i evaluted in node k.
*
*                  The notation for the example above in the inline comments
*                  is: @Nj/@Xi (k)
*                  
*                  Where Xi = X, Y, Z (I=1,2,3, respectively)
*                  i.e. Derivative of shape function "j", with respect
*                  direction "Xi", evaluates in node "k".
*
*                  One more example. 
*                  If the gradient of the shape functions of a triangle wants 
*                  to be compute in two points, the output will be a matrix 
*                  of shape (dims., 3,NP)
*                  Then GRADLOC will be:
*
*                                   | @N1/@X(p1)   @N2/@X(p1)   @N3/@X(p1) |
*                                   |                                      |
*                  GRADLOC(:,:,1) = | @N1/@Y(p1)   @N2/@Y(p1)   @N3/@Y(p1) | = 
*                                   |                                      |
*                                   | @N1/@Z(p1)   @N2/@Z(p1)   @N3/@Z(p1) |
*
*
***************************************************************************

      IMPLICIT NONE
      
      INTEGER*4::IODIM,LMXNDL,LTYP,NP
      REAL*8::XLOC(NP),YLOC(NP),ZLOC(NP)
      REAL*8::GRADLOC(IODIM,LMXNDL,NP)


C-------------------- Initialization
      GRADLOC = 0D0
      
C-------------------- Computes local gradient of shape functions
C-------------------- in the given points.

      SELECT CASE(LTYP)

          CASE(1)  ! 1D elements

              GRADLOC(1,1,1:NP) = -1D0                           ! @N1/@X (1:NP)
              GRADLOC(1,2,1:NP) =  1D0                           ! @N2/@X (1:NP)
      

          CASE(2,5) ! Linear 2D triangles, 3 nodes.    
                
              GRADLOC(1,1,1:NP) = -1D0                           ! @N1/@X (1:NP)
              GRADLOC(2,1,1:NP) = -1D0                           ! @N1/@Y (1:NP)
              GRADLOC(1,2,1:NP) =  1D0                           ! @N2/@X (1:NP)
              GRADLOC(2,2,1:NP) =  0D0                           ! @N2/@Y (1:NP)
              GRADLOC(1,3,1:NP) =  0D0                           ! @N3/@X (1:NP)
              GRADLOC(2,3,1:NP) =  1D0                           ! @N3/@Y (1:NP)


          CASE(3) ! Quadrilaterales

          
              GRADLOC(1,1,1:NP) = -1D0 + YLOC(1:NP)              ! @N1/@X (1:NP)
              GRADLOC(2,1,1:NP) = -1D0 + XLOC(1:NP)              ! @N1/@Y (1:NP)
              GRADLOC(1,2,1:NP) =  1D0 - YLOC(1:NP)              ! @N2/@X (1:NP)
              GRADLOC(2,2,1:NP) = -XLOC(1:NP)                    ! @N2/@Y (1:NP)
              GRADLOC(1,3,1:NP) =  YLOC(1:NP)                    ! @N3/@X (1:NP)
              GRADLOC(2,3,1:NP) =  XLOC(1:NP)                    ! @N3/@Y (1:NP)
              GRADLOC(1,4,1:NP) = -YLOC(1:NP)                    ! @N4/@X (1:NP)
              GRADLOC(2,4,1:NP) =  1D0 - XLOC(1:NP)              ! @N4/@Y (1:NP)


          CASE(4)  ! Tetrahedricons

              GRADLOC(1,1,1:NP) = -1D0                           ! @N1/@X (1:NP)
              GRADLOC(2,1,1:NP) = -1D0                           ! @N1/@Y (1:NP)
              GRADLOC(3,1,1:NP) = -1D0                           ! @N1/@Z (1:NP)
              GRADLOC(1,2,1:NP) =  1D0                           ! @N2/@X (1:NP)
              GRADLOC(2,2,1:NP) =  0D0                           ! @N2/@Y (1:NP)
              GRADLOC(3,2,1:NP) =  0D0                           ! @N2/@Z (1:NP)
              GRADLOC(1,3,1:NP) =  0D0                           ! @N3/@X (1:NP)
              GRADLOC(2,3,1:NP) =  1D0                           ! @N3/@Y (1:NP)
              GRADLOC(3,3,1:NP) =  0D0                           ! @N3/@Z (1:NP)
              GRADLOC(1,4,1:NP) =  0D0                           ! @N4/@X (1:NP)
              GRADLOC(2,4,1:NP) =  0D0                           ! @N4/@Y (1:NP)
              GRADLOC(3,4,1:NP) =  1D0                           ! @N4/@Z (1:NP)


          CASE(6) ! Prismatics


              GRADLOC(1,1,1:NP) = -1D0 + YLOC(1:NP) + ZLOC(1:NP) ! @N1/@X (1:NP)
              GRADLOC(2,1,1:NP) = -1D0 + XLOC(1:NP)              ! @N1/@Y (1:NP)
              GRADLOC(3,1,1:NP) = -1D0 + XLOC(1:NP)              ! @N1/@Z (1:NP)
              GRADLOC(1,2,1:NP) = -1D0*YLOC(1:NP)                ! @N2/@X (1:NP)
              GRADLOC(2,2,1:NP) =  1D0 - XLOC(1:NP)              ! @N2/@Y (1:NP)
              GRADLOC(3,2,1:NP) =  0D0                           ! @N2/@Z (1:NP)
              GRADLOC(1,3,1:NP) = -1D0*ZLOC(1:NP)                ! @N3/@X (1:NP)
              GRADLOC(2,3,1:NP) =  0D0                           ! @N3/@Y (1:NP)
              GRADLOC(3,3,1:NP) =  1D0 - XLOC(1:NP)              ! @N3/@Z (1:NP)
              GRADLOC(1,4,1:NP) =  1D0 - YLOC(1:NP) - ZLOC(1:NP) ! @N4/@X (1:NP)
              GRADLOC(2,4,1:NP) = -1D0*XLOC(1:NP)                ! @N4/@Y (1:NP)
              GRADLOC(3,4,1:NP) = -1D0*XLOC(1:NP)                ! @N4/@Z (1:NP)
              GRADLOC(1,5,1:NP) =  YLOC(1:NP)                    ! @N5/@X (1:NP)
              GRADLOC(2,5,1:NP) =  XLOC(1:NP)                    ! @N5/@Y (1:NP)
              GRADLOC(3,5,1:NP) =  0D0                           ! @N5/@Z (1:NP)
              GRADLOC(1,6,1:NP) =  ZLOC(1:NP)                    ! @N6/@X (1:NP)
              GRADLOC(2,6,1:NP) =  0D0                           ! @N6/@Y (1:NP)
              GRADLOC(3,6,1:NP) =  XLOC(1:NP)                    ! @N6/@Z (1:NP)

      END SELECT   


      END SUBROUTINE COMP_GRAD_LOC