      SUBROUTINE AREA_ELEM
     ;(IERROR   ,IOWAR    ,LMXNDL   ,MAINF    ,NUMEL    ,NUMNP
     ;,AREA     ,KXX      ,LTYPE    ,X        ,Y        ,Z
     ;,FILENAME)

********************************************************************************
*
* PURPOSE Computes the size (lenght, area or volume for 1,2 or 3-D problems,
*         respectively)
*
* DESCRIPTION Assigns the element type and performs different computations,
*             depending on element shape. The elements considered are:
*                     LTYPE(L)=1 : 1-D Element (line in space)
*                     LTYPE(L)=2 : 2-D Element (triangular)
*                     LTYPE(L)=3 : 2-D Element (cuadrangular)
*                     LTYPE(L)=4 : 3-D Element (tetrahedron)
*                     LTYPE(L)=5 : 2-D Element (triangle in space)
*                     LTYPE(L)=6 : 3-D Element (prismatic)
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LTYPE                  Vector containing the type of each element            
*  X                      X-coord for a given node                              
*  Y                      Y-coord for a given node                              
*  Z                      Z-coord for a given node                              
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  GRADN_PT_3D                                                                  
*  LOC_COORDINADES                                                              
*
* HISTORY
*
*     AAR        01-2001     First coding
*
********************************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
     
      DIMENSION XG(3),YG(3),ZG(3),VI(3),VJ(3),VK(3),XLOC(3),YLOC(3)
     ;      ,AREA(NUMEL),PESO(6),GRADN(24,6),LTYPE(NUMEL)
     ;      ,X(NUMNP),Y(NUMNP),Z(NUMNP),KXX(LMXNDL,NUMEL)

      CHARACTER FILENAME(20)*20

C______________________________ Starts loop over elements

      DO L=1,NUMEL

C______________________________ Modifies LTYPE array (COHERENCY WITH PRODAT)

        LTYPE2=0
        IF (LTYPE(L).EQ.3) THEN                       ! Triangular element (2-D)
          LTYPE2=2
        ELSE IF (LTYPE(L).EQ.5) THEN                ! Cuadrangular element (2-D)
          LTYPE2=3
        ELSE IF (LTYPE(L).EQ.9) THEN                 ! Tetrahedron element (3-D)
          LTYPE2=4
        ELSE IF (LTYPE(L).EQ.10) THEN                 ! Triangular element (3-D)
          LTYPE2=5
        ELSE IF (LTYPE(L).EQ.11)THEN                           ! Toblerone (3-D)
          LTYPE2=6
        END IF

C______________________________ Identifies some useful variables

        IF (LTYPE2.NE.0) THEN
          LTYPE1=LTYPE2
        ELSE
          LTYPE1=LTYPE(L)
        END IF

        I1=KXX(1,L)                                             ! Connectivities
        I2=KXX(2,L)                
        X1=X(I1)                                                 ! X-coordinates
        X2=X(I2)
        Y1=Y(I1)                                                 ! Y-coordinates
        Y2=Y(I2)
        Z1=Z(I1)                                                 ! Z-coordinates
        Z2=Z(I2)

C______________________________ Computations depending on element shape

        IF (LTYPE1.EQ.1) THEN                                     ! 1-D elements

           AREA1=DSQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)
     ;                +(Z2-Z1)*(Z2-Z1))

        ELSE IF (LTYPE1.EQ.2) THEN                   ! Triangular elements (2-D)

           I3=KXX(3,L)
           X3=X(I3)
           Y3=Y(I3)
           B2=Y3-Y1
           C2=X1-X3
           B3=Y1-Y2
           C3=X2-X1
           AREA1=(C3*B2-C2*B3)/2.D+00

        ELSE IF (LTYPE1.EQ.3) THEN                 ! Cuadrangular elements (2-D)

          I4=KXX(4,L)
          X4=X(I4)
          Y4=Y(I4)
          AREA1=DABS(X4-X2)*DABS(Y4-Y2)

        ELSE IF (LTYPE1.EQ.4) THEN                  ! Tetrahedron elements (3-D)

          I3=KXX(3,L)
          I4=KXX(4,L)
          X3=X(I3)
          Y3=Y(I3)
          Z3=Z(I3)
          X4=X(I4)
          Y4=Y(I4)
          Z4=Z(I4)

          A1=X2*Y3*Z4+X3*Y4*Z2+X4*Y2*Z3
     .      -X4*Y3*Z2-X3*Y2*Z4-X2*Y4*Z3
          A2=-X1*Y3*Z4-X3*Y4*Z1-X4*Y1*Z3
     .      +X4*Y3*Z1+X3*Y1*Z4+X1*Y4*Z3
          A3=+X1*Y2*Z4+X2*Y4*Z1+X4*Y1*Z2
     .      -X4*Y2*Z1-X2*Y1*Z4-X1*Y4*Z2
          A4=-X1*Y2*Z3-X2*Y3*Z1-X3*Y1*Z2
     .      +X3*Y2*Z1+X2*Y1*Z3+X1*Y3*Z2

          AREA1=(A1+A2+A3+A4)/6.0D0

        ELSE IF (LTYPE1.EQ.5) THEN                   ! Triangular elements (3-D)

          I3=KXX(3,L)
          X3=X(I3)
          Y3=Y(I3)
          Z3=Z(I3)

          XG(1)=X1
          XG(2)=X2
          XG(3)=X3
          YG(1)=Y1
          YG(2)=Y2
          YG(3)=Y3
          ZG(1)=Z1
          ZG(2)=Z2
          ZG(3)=Z3

C______________________________ Computes the maximum slope vector. It will be 
C______________________________ used to compute element area

          CALL LOC_COORDINADES 
     ; (XG(1),     XG(2),     XG(3),     YG(1),     YG(2),     YG(3),
     ;  ZG(1),     ZG(2),     ZG(3),    VMINPI,    VMINPJ,    VMINPK,
     ; VMAXPI,    VMAXPJ,    VMAXPK,  XLOCCERO,  YLOCCERO,   ZLOCCERO)

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

          B2=Y3-Y1
          C2=X1-X3
          B3=Y1-Y2
          C3=X2-X1

          AREA1=(C3*B2-C2*B3)/2.D+00

        ELSE IF (LTYPE1.EQ.6) THEN                    ! Prismatic elements (3-D)
                                                        ! Six integration points

          I3=KXX(3,L)
          I4=KXX(4,L)
          I5=KXX(5,L)
          I6=KXX(6,L)
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

          CALL GRADN_PT_3D 
     ;(X1    ,Y1    ,Z1    ,X2    ,Y2    ,Z2
     ;,X3    ,Y3    ,Z3    ,X4    ,Y4    ,Z4
     ;,X5    ,Y5    ,Z5    ,X6    ,Y6    ,Z6
     ;,GRADN ,PESO  ,6)

          AREA1=0.D0
          DO I=1,6
            AREA1=AREA1+PESO(I)
          END DO
          
        END IF

C______________________________ Checks if element size is smaller than tolerance

        IF (DABS(AREA1) .LT. 1.D-20)
     ;     CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;       ' WARNING: SIZE LOWER THAN 1.D-20. ELEMENT ',L,0,2,3,0.00)

C______________________________ Final check. Negative size

        IF (AREA1 .LT. 0.D+00) THEN
           WRITE (MAINF,100) L
 100       FORMAT(/,' CRITICAL STOP: NEGATIVE AREA OF ELEMENT ',I5,'.'
     ;           ,/,' PROBABLY, THAT IS DUE TO A VERY SHARP SHAPE.'
     ;          ,//,' NODAL COORDINATES:'
     ;           ,/,' ===== ============'
     ;          ,//,'  NODE          X          Y          Z'
     ;           ,/,' ===== ========== ========== ==========')
           DO J=1,LMXNDL
             INODE=KXX(J,L)
             IF (INODE.GT.0) 
     ;          WRITE(MAINF,200) INODE,X(INODE),Y(INODE),Z(INODE)
 200         FORMAT(1X,I5,3(1X,E10.3))
           END DO
           
           STOP ' FATAL ERROR. CHECK FILE RES.OUT'

        END IF

C______________________________ If element size is not negative, it is stored
C______________________________ in array AREA

        AREA(L)=AREA1

      END DO                                                      ! Next element
      RETURN
      END
