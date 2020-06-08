      SUBROUTINE BASISFUNC_OBS
     ;(LMXNDL   ,MAINF    ,NEL      ,NUMEL    ,NUMNP    ,XOB
     ;,YOB      ,ZOB      ,AREA     ,BF       ,KXX      ,LTYPE
     ;,X        ,Y        ,Z)       

***********************************************************************
* PURPOSE
*
* Determines basis function values.
*
* DESCRIPTION
*
* This subroutine is called once for every basic unit. Finds the
* element to which a point belongs and determines the value of the
* basis function for the nodes belonging to this element. It returns
* the basis function values BF.
* The subroutine carries out more or less what INTER_CAL and TRBASIS
* carried out in TRANSIN III.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  BF                     Basis function values
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LTYPE                  Vector containing the type of each element            
*  X                      X-coord for a given node
*                         (e.g.: X(KXX(2,5)) is the X-coord.
*                         of the 2nd node of the 5th element)
*  Y                      Y-coord for a given node                              
*  Z                      Z-coord for a given node                              
*
* EXTERNAL VARIABLES: SCALARS
*
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NEL                    Element number at which basic unit point belongs to
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  XOB                    X-coordinate of basic unit
*  YOB                    Y-coordinate of basic unit                            
*  ZOB                    Z-coordinate of basic unit                            
*
* INTERNAL VARIABLES: SCALARS
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  TRANGLE3DBASIS
*  TOBLEOBS
*
* HISTORY
*
*     AMS      4-1998     Revision and common elimination (INTER_CAL)
*     CK      11-1999     Modification of INTER_CAL and TRBASIS
*     AAR     03-2001     Revision 
*
***********************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       DIMENSION KXX(LMXNDL,NUMEL),X(NUMNP),Y(NUMNP),AREA(NUMEL),          
     ; BF(6),LTYPE(NUMEL),Z(NUMNP)

       ZERO=-0.0000001

* ............................................... LOOK FOR THE ELEMENT

       DO 20 NEL=1,NUMEL

C______________________________ Modifies LTYPE array (COHERENCY WITH PRODAT)

         LTYPE2=0
         IF (LTYPE(NEL).EQ.3) THEN                       ! Triangular element (2-D)
           LTYPE2=2
         ELSE IF (LTYPE(NEL).EQ.5) THEN                ! Cuadrangular element (2-D)
           LTYPE2=3
         ELSE IF (LTYPE(NEL).EQ.9) THEN                 ! Tetrahedron element (3-D)
           LTYPE2=4
         ELSE IF (LTYPE(NEL).EQ.10) THEN                 ! Triangular element (3-D)
           LTYPE2=5
         ELSE IF (LTYPE(NEL).EQ.11)THEN                           ! Toblerone (3-D)
           LTYPE2=6
         END IF

C______________________________ Identifies some useful variables

         IF (LTYPE2.NE.0) THEN
           LTYPE1=LTYPE2
         ELSE
           LTYPE1=LTYPE(NEL)
         END IF

* .......................................RECTANGLES

         IF(LTYPE1.EQ.3) THEN
           AR=-AREA(NEL)
           KY=4
           KX=2
           DO 30 I=1,4
             AR=-AR
             BF(I)=(Y(KXX(KY,NEL))-YOB)*(X(KXX(KX,NEL))-XOB)/AR
             KY=KY-1
             KX=KX-1
             IF(KX.EQ.0)KX=4
   30      CONTINUE             
           IF (BF(1).GE.ZERO.AND.BF(2).GE.ZERO.AND.BF(3).GE.ZERO
     ;                      .AND.BF(4).GE.ZERO) RETURN    ! Point belongs to NEL

* ................................ TRIANGLES DEFINED ON 2-D SPACE

         ELSE IF(LTYPE1.EQ.2)THEN
           AR=2*AREA(NEL)
           X1=X(KXX(1,NEL))
           Y1=Y(KXX(1,NEL))           
           X2=X(KXX(2,NEL))
           Y2=Y(KXX(2,NEL))
           X3=X(KXX(3,NEL))
           Y3=Y(KXX(3,NEL))

           BF(1)=(X2*Y3-X3*Y2+(Y2-Y3)*XOB+(X3-X2)*YOB)/AR               
           BF(2)=(X3*Y1-X1*Y3+(Y3-Y1)*XOB+(X1-X3)*YOB)/AR               
           BF(3)=(X1*Y2-X2*Y1+(Y1-Y2)*XOB+(X2-X1)*YOB)/AR 
           IF (BF(1).GE.ZERO.AND.BF(2).GE.ZERO.AND.BF(3).GE.ZERO) RETURN !OK

*....................................... LINEAR SEGMENTS (1D ELEMENTS)

         ELSE IF(LTYPE1.EQ.1)THEN
           I1=KXX(1,NEL)
           I2=KXX(2,NEL)
           X1=X(I1)
           Y1=Y(I1)
           X2=X(I2)
           Y2=Y(I2)
           D12=DSQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
           D1P=DSQRT((X1-XOB)*(X1-XOB)+(Y1-YOB)*(Y1-YOB))
           D2P=DSQRT((X2-XOB)*(X2-XOB)+(Y2-YOB)*(Y2-YOB))
           BF(1)=1-D1P/D12
           BF(2)=D1P/D12
           IF(DABS(D12-D1P-D2P).LE.-ZERO) RETURN         ! Point belongs to NEL

* .......................................................... TETRAHEDRA

         ELSE IF(LTYPE1.EQ.4)THEN
           I1=KXX(1,NEL)
           I2=KXX(2,NEL)
           I3=KXX(3,NEL)
           I4=KXX(4,NEL)

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

           B1=-(Y3*Z4-Y4*Z3)+(Y2*Z4-Y4*Z2)-(Y2*Z3-Y3*Z2)
           B2=+(Y3*Z4-Y4*Z3)-(Y1*Z4-Y4*Z1)+(Y1*Z3-Y3*Z1)
           B3=-(Y2*Z4-Y4*Z2)+(Y1*Z4-Y4*Z1)-(Y1*Z2-Y2*Z1)
           B4=+(Y2*Z3-Y3*Z2)-(Y1*Z3-Y3*Z1)+(Y1*Z2-Y2*Z1)

           C1=+(X3*Z4-X4*Z3)-(X2*Z4-X4*Z2)+(X2*Z3-X3*Z2)
           C2=-(X3*Z4-X4*Z3)+(X1*Z4-X4*Z1)-(X1*Z3-X3*Z1)
           C3=+(X2*Z4-X4*Z2)-(X1*Z4-X4*Z1)+(X1*Z2-X2*Z1)
           C4=-(X2*Z3-X3*Z2)+(X1*Z3-X3*Z1)-(X1*Z2-X2*Z1)

           D1=-(X3*Y4-X4*Y3)+(X2*Y4-X4*Y2)-(X2*Y3-X3*Y2)
           D2=+(X3*Y4-X4*Y3)-(X1*Y4-X4*Y1)+(X1*Y3-X3*Y1)
           D3=-(X2*Y4-X4*Y2)+(X1*Y4-X4*Y1)-(X1*Y2-X2*Y1)
           D4=+(X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1)
  
           A1=X2*Y3*Z4+X3*Y4*Z2+X4*Y2*Z3
     .     -X4*Y3*Z2-X3*Y2*Z4-X2*Y4*Z3
           A2=-X1*Y3*Z4-X3*Y4*Z1-X4*Y1*Z3
     .     +X4*Y3*Z1+X3*Y1*Z4+X1*Y4*Z3
           A3=+X1*Y2*Z4+X2*Y4*Z1+X4*Y1*Z2
     .     -X4*Y2*Z1-X2*Y1*Z4-X1*Y4*Z2
           A4=-X1*Y2*Z3-X2*Y3*Z1-X3*Y1*Z2
     .     +X3*Y2*Z1+X2*Y1*Z3+X1*Y3*Z2

           SEISVOL=(A1+A2+A3+A4)
           VOLELM=SEISVOL/6.0
           IF(DABS(VOLELM-AREA(NEL)).GT.0.00001)THEN
             WRITE(MAINF,*)
     .       ' DISCREPANCIAS EN EL CALCULO DE VOLUMENES'
             STOP
           ENDIF

           BF(1)=(A1+B1*XOB+C1*YOB+D1*ZOB)/SEISVOL
           BF(2)=(A2+B2*XOB+C2*YOB+D2*ZOB)/SEISVOL
           BF(3)=(A3+B3*XOB+C3*YOB+D3*ZOB)/SEISVOL
           BF(4)=(A4+B4*XOB+C4*YOB+D4*ZOB)/SEISVOL
           IF (BF(1).GE.ZERO.AND.BF(2).GE.ZERO.AND.BF(3).GE.ZERO
     ;                      .AND.BF(4).GE.ZERO) RETURN    ! Point belongs to NEL

* ................................ TRIANGLES DEFINED ON 3-D SPACE

         ELSE IF(LTYPE1.EQ.5)THEN
           I1=KXX(1,NEL)
           I2=KXX(2,NEL)
           I3=KXX(3,NEL)

           X1=X(I1)
           Y1=Y(I1)
           Z1=Z(I1)
              
           X2=X(I2)
           Y2=Y(I2)
           Z2=Z(I2)

           X3=X(I3)
           Y3=Y(I3)
           Z3=Z(I3)
           AR=2*AREA(NEL)

           CALL TRANGLE3DBASIS
     ;     (X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,XOB,YOB,BF,
     ;      ZOB,MAINF,NEL,AR,INDIC)

           IF(INDIC.EQ.0) THEN
             GOTO 20                        ! Go on to next element
           ELSE
             RETURN                         ! Point belongs to element NEL
           END IF
* .......................................................... TOBLERONES

         ELSE
           I1=KXX(1,NEL)
           I2=KXX(2,NEL)
           I3=KXX(3,NEL)
           I4=KXX(4,NEL)
           I5=KXX(5,NEL)
           I6=KXX(6,NEL)

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

           CALL TOBLEOBS (X1,X2,X3,X4,X5,X6,Y1,Y2,Y3,Y4,Y5,Y6,
     ;                    Z1,Z2,Z3,Z4,Z5,Z6,XOB,YOB,ZOB,MAINF,
     ;                    NEL,BF,INDIC)

           IF(INDIC.EQ.0) THEN
             GOTO 20                        ! Go on to next element
           ELSE
             RETURN                         ! Point belongs to element NEL
           END IF

         ENDIF

*.............. Current point does not belong to this element, look for next one
 
 20    CONTINUE

       RETURN

       END

********************************************************************************
********************************************************************************
********************************************************************************

       SUBROUTINE TRANGLE3DBASIS (X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,XOB,YOB,BF,
     ;                           ZOB,MAINF,L,AR,INDIC)
********************************************************************************
***  THIS ROUTINE COMPUTES THE INTERPOLATION FUNCTIONS AT AN OBSERVATION     ***
***            POINT FOR TRIANGLES DEFINED ON 3-D SPACE                      ***
********************************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION BF(4)
       ZERO=-0.0000001

* computes the distance from the observation point to the plane of the triangle


* assemles vector [1-3]
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

* assemles vector [1-observation point] and computes the distence from the
* observation point, to the plane.
       PI=XOB-X1
       PJ=YOB-Y1
       PK=ZOB-Z1
       XMODP=DSQRT(PI**2+PJ**2+PK**2)
       DIST=(PI*XNI+PJ*XNJ+PK*XNK)

* defines the proyection of the observation poin on the plane
       EPSILON=0.01
       IF (DABS(DIST).LE.EPSILON)THEN
         IF (DABS(XNK).LT.1E-08)THEN   !vertical plane (particular case)
           CALL VERTICAL_PLANE (XNI,XNJ,XNK,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,
     ;     Z3,XOB,YOB,ZOB,XL1,YL1,XL2,YL2,XL3,YL3,XLO,YLO)
           GOTO 27
         ENDIF
         XOB=XOB-DIST*XNI             ! SI ESTA POR ENCIMA DEL PLANO 
         YOB=YOB-DIST*XNJ             ! (EN RELACION A Z),
         ZOB=ZOB-DIST*XNK             ! LA DISTANCIA ES POSITIVA Y SE RESTA

* assemles vector [1-point of proyection] 
         PPLI=XOB-X1
         PPLJ=YOB-Y1
         PPLK=ZOB-Z1 
         COMP=PPLI*XNI+PPLJ*XNJ+PPLK*XNK         !checks ortogonality
         IF (DABS(COMP).GT.0.0001)THEN
           WRITE(MAINF,*)' PROBLABLY THE ELEMENT',L,
     ;     ' IS BAD CONDITIONED GEOMETRICALY, I HAVE PROBLEMS TO',
     ;     ' DEFINE ITS BASIS FUNCTIONS, PLEACE CHECK IT'
           STOP
         ENDIF

* defines the local coordinates of the point 1
         XL1=0D0
         YL1=0D0

* defines the local coordinates of the point 2
         XL2=XMODV    
         YL2=0D0      
         
* assembles a vector in local y coordinate
         VYI=VJ*XNK-VK*XNJ
         VYJ=-VI*XNK+VK*XNI
         VYK=VI*XNJ-VJ*XNI
         XMODVY=DSQRT(VYI**2+VYJ**2+VYK**2)
         VYI=-VYI/XMODVY
         VYJ=-VYJ/XMODVY
         VYK=-VYK/XMODVY

* checks ortogonality
         PPUNTO=VYI*VI+VYJ*VJ+VYK*VK
         IF (PPUNTO.GT.0.0001)THEN
           WRITE(MAINF,*)' PROBLABLY THE ELEMENT',L,
     ;     ' IS BAD CONDITIONED GEOMETRICALY, I HAVE PROBLEMS TO',
     ;     ' DEFINE ITS BASIS FUNCTIONS, PLEACE CHECK IT---'
           STOP
         ENDIF

* defines the local coordinates of the point 3
         XL3=(UI*VI+UJ*VJ+UK*VK)/XMODV
         YL3=UI*VYI+UJ*VYJ+UK*VYK

* defines the local coordinates of the observation point
         XLO=(PPLI*VI+PPLJ*VJ+PPLK*VK)/XMODV
         YLO=PPLI*VYI+PPLJ*VYJ+PPLK*VYK
 27      BF(1)=(XL2*YL3-XL3*YL2+(YL2-YL3)*XLO+(XL3-XL2)*YLO)/AR
         BF(2)=(XL3*YL1-XL1*YL3+(YL3-YL1)*XLO+(XL1-XL3)*YLO)/AR
         BF(3)=(XL1*YL2-XL2*YL1+(YL1-YL2)*XLO+(XL2-XL1)*YLO)/AR 
         IF (BF(1).LT.ZERO. OR .BF(2).LT.ZERO. OR .BF(3).LT.ZERO)THEN
           INDIC=0
         ELSE
           INDIC=1
         ENDIF
       ELSE
         INDIC=0
         RETURN
       ENDIF
       END

********************************************************************************
********************************************************************************
********************************************************************************

       SUBROUTINE VERTICAL_PLANE (XNI,XNJ,XNK,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,
     ;                         Z3,XOB,YOB,ZOB,X1LOC,Y1LOC,X2LOC,Y2LOC,
     ;                         X3LOC,Y3LOC,XOBLOC,YOBLOC)
   
********************************************************************************
*** CHECKS WETHER THE OBSERVATION POINT BELONGS TO THE TRIANGLE WHICH IS     ***
***                       INSIDE A VERTICAL PLANE                            ***
********************************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)

       VYLOCI=0D0                    !defines a vector in Y local direction
       VYLOCJ=0D0
       VYLOCK=1D0
       IF (DABS(XNI).GT.1E-08) THEN  !component x(global) no equal to zero
         IF (XNI.LT.0D0) THEN            !prescribes component x positive
           XNI=-XNI
           XNJ=-XNJ
           XNK=-XNK
         ENDIF
       ELSE                          !component x(global) equal to zero
         IF (XNJ.LT.0D0)THEN             !prescribes component y positive
           XNI=-XNI
           XNJ=-XNJ
           XNK=-XNK
         ENDIF 
       ENDIF
       VXLOCI=-XNJ           !defines a vector in X local direction
       VXLOCJ=XNI
       VXLOCK=0D0

*                            !defines the 0,0 of the local coordenades
       DELX=XOB-X1           
       DELY=YOB-Y1
       DELZ=ZOB-Z1
       IF (DSQRT(DELX*DELX+DELY*DELY+DELZ*DELZ).LT.1E-08)THEN !definition of
         XCERO=X2       !vector o-po requires that the relative position of po 
         YCERO=Y2       !be different as o
         ZCERO=Z2 
       ELSE
         XCERO=X1
         YCERO=Y1
         ZCERO=Z1 
       ENDIF

       VOBI=XOB-XCERO  !defines a vector in direction o-po (o=0,0 local and
       VOBJ=YOB-YCERO  !po= observation point)
       VOBK=ZOB-ZCERO
       XOBLOC=VOBI*VXLOCI+VOBJ*VXLOCJ+VOBK*VXLOCK
       YOBLOC=VOBI*VYLOCI+VOBJ*VYLOCJ+VOBK*VYLOCK

       V1I=X1-XCERO    !defines vector o-1 and the local coord. of 1
       V1J=Y1-YCERO
       V1K=Z1-ZCERO
       X1LOC=V1I*VXLOCI+V1J*VXLOCJ+V1K*VXLOCK
       Y1LOC=V1I*VYLOCI+V1J*VYLOCJ+V1K*VYLOCK

       V2I=X2-XCERO    !defines vector o-2 and the local coord. of 2
       V2J=Y2-YCERO
       V2K=Z2-ZCERO
       X2LOC=V2I*VXLOCI+V2J*VXLOCJ+V2K*VXLOCK
       Y2LOC=V2I*VYLOCI+V2J*VYLOCJ+V2K*VYLOCK

       V3I=X3-XCERO    !defines vector o-3 and the local coord. of 3
       V3J=Y3-YCERO
       V3K=Z3-ZCERO
       X3LOC=V3I*VXLOCI+V3J*VXLOCJ+V3K*VXLOCK
       Y3LOC=V3I*VYLOCI+V3J*VYLOCJ+V3K*VYLOCK
       RETURN
       END

********************************************************************************
********************************************************************************
********************************************************************************

       SUBROUTINE TOBLEOBS (X1,X2,X3,X4,X5,X6,Y1,Y2,Y3,Y4,Y5,Y6,
     ;                      Z1,Z2,Z3,Z4,Z5,Z6,XOB,YOB,ZOB,MAINF,
     ;                      L,BF,INDIC)
********************************************************************************
*** COMPUTE THE BASIS FUNCTION AT OBSERVATION POINTS FOR TOBLERONE ELEMENTS  ***
********************************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION U(3,6),Y0(3),Y(3),A(3,3),XJACOB(3,3),B(3),
     ;           F(3),DELTA(3),ABANDA(3,5),WORK(9),BF(6)

       Y(1)=0.33
       Y(2)=0.33
       Y(3)=0.33

C       write(5,*)' aqui',X1,X2,X3,X4,X5,X6,Y1,Y2,Y3,Y4,Y5,Y6,
C     ;                      Z1,Z2,Z3,Z4,Z5,Z6,XOB,YOB,ZOB

       U(1,1)=X1
       U(1,2)=X2
       U(1,3)=X3
       U(1,4)=X4
       U(1,5)=X5
       U(1,6)=X6

       U(2,1)=Y1
       U(2,2)=Y2
       U(2,3)=Y3
       U(2,4)=Y4
       U(2,5)=Y5
       U(2,6)=Y6

       U(3,1)=Z1
       U(3,2)=Z2
       U(3,3)=Z3
       U(3,4)=Z4
       U(3,5)=Z5
       U(3,6)=Z6

       Y0(1)=XOB
       Y0(2)=YOB
       Y0(3)=ZOB

***  PRELIMINARY CHECK OF POINT WITHIN ELEMENT   !revisar
       INDIC=0          
       DO I=1,3
          UMX=U(I,1)
          UMN=U(I,1)
          DO J=2,6
             IF (U(I,J).GT.UMX) UMX=U(I,J)
             IF (U(I,J).LT.UMN) UMN=U(I,J)
          ENDDO
          IF (Y0(I).LT.UMN .OR. Y0(I).GT.UMX) RETURN ! point clearly outside elt
       ENDDO                                   

*** POINT MAY BELONG TO ELEMNT, START ITERATIVE PROCESS !REVISAR HASTA AQUI

       ITER=0             !iban antes de definir la U
 10    ITER=ITER+1        !revisar que este bien

       DO I=1,3
         A(I,1)=(1D0-Y(3))*(U(I,1)-U(I,3)) + 
     ;          (1D0+Y(3))*(U(I,4)-U(I,6))

         A(I,2)=(1D0-Y(3))*(U(I,2)-U(I,3)) + 
     ;          (1D0+Y(3))*(U(I,5)-U(I,6))
                 
         A(I,3)=U(I,6)-U(I,3)
         B(I)=-2*Y0(I)+U(I,3)+U(I,6)
       END DO

       DO I=1,3
         XJACOB(I,1)=A(I,1)
         XJACOB(I,2)=A(I,2)
         XJACOB(I,3)=A(I,3)+(-U(I,1)+U(I,4)+U(I,3)-U(I,6))*Y(1)
     ;                    +(-U(I,2)+U(I,5)+U(I,3)-U(I,6))*Y(2)
       END DO

       FMAXIM=0D0
       DO I=1,3
         F(I)=B(I)
         DO J=1,3
           F(I)=F(I)+A(I,J)*Y(J)
         END DO
         IF (DABS(F(I)).GT.FMAXIM) THEN
           IRES=I
           FMAXIM=DABS(F(I))
         ENDIF
       END DO

       DO I=1,3
         DO J=1,5
           ABANDA(I,J)=0D0
         END DO
       END DO
       ABANDA(1,1)=0D0
       ABANDA(1,2)=0D0
       ABANDA(1,3)=XJACOB(1,1)
       ABANDA(1,4)=XJACOB(1,2)
       ABANDA(1,5)=XJACOB(1,3)
       ABANDA(2,1)=0D0
       ABANDA(2,2)=XJACOB(2,1)
       ABANDA(2,3)=XJACOB(2,2)
       ABANDA(2,4)=XJACOB(2,3)
       ABANDA(2,5)=0D0
       ABANDA(3,1)=XJACOB(3,1)
       ABANDA(3,2)=XJACOB(3,2)
       ABANDA(3,3)=XJACOB(3,3)
       ABANDA(3,4)=0D0
       ABANDA(3,5)=0D0

       DO I=1,3
         DELTA(I)=-F(I)
       END DO

       CALL LEQT1B (ABANDA,3,2,2,
     .               3,DELTA,1,3,0,WORK,IER)

       DELMAXIM=0D0
       DO I=1,3
         IF (DABS(DELTA(I)).GT.DELMAXIM) THEN
           IDELTA=I
           DELMAXIM=DABS(DELTA(I))
         ENDIF
       END DO

C       WRITE(5,*) 
C       WRITE(5,*) 
C       WRITE(5,*) 
C     ;   ' ITERACION --NOD RES   --RESIDUO ----NOD F ---FUNCION'
C       WRITE(5,'(I6,4X,I6,4X,E12.6,2X,I6,2X,E12.6)')ITER,IRES,FMAXIM,
C     ;                                       IDELTA,DELMAXIM

       DO I=1,3
         Y(I)=Y(I)+DELTA(I)
       END DO
       F1=.5*(1D0-Y(3))*Y(1)
       F2=.5*(1D0-Y(3))*Y(2)
       F3=.5*(1D0-Y(3))*(1D0-Y(1)-Y(2))
       F4=.5*(1D0+Y(3))*Y(1)
       F5=.5*(1D0+Y(3))*Y(2)
       F6=.5*(1D0+Y(3))*(1D0-Y(1)-Y(2))
       SUMA=F1+F2+F3+F4+F5+F6

c       write(5,*)'****************************************************'   
c       write(5,*)
c       write(5,*)
c     ; '     f1     f2     f3     f4     f5     f6    suma'
c       write(5,'(7(g13.5,1x))')f1,f2,f3,f4,f5,f6,suma

       EPSILON=1E-10
       EPS=-1.e-4
       IF (DELMAXIM.LT.EPSILON. AND .FMAXIM.LT.EPSILON)THEN
         IF (F1.GE.EPS. AND . F2.GE.EPS. AND .F3.GE.EPS. AND .
     ;     F4.GE.EPS. AND . F5.GE.EPS. AND .F6.GE.EPS)THEN
           BF(1)=F1         
           BF(2)=F2         
           BF(3)=F3         
           BF(4)=F4         
           BF(5)=F5         
           BF(6)=F6         
           INDIC=1
         ELSE
           INDIC=0
         ENDIF
         RETURN
       ELSE
         IF (ITER.GT.30) THEN
           WRITE(MAINF,*)' I HAVE PROBLEMS SOLVING BASIS FUNCTIONS'
           WRITE(MAINF,*)' AT ELEMENT ',L,' SUBROUTINE TOBLOBS'
           STOP
         ENDIF
         GOTO 10
       ENDIF
       END

***%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*            SUBROUTINAS QUE NO SE USAN!!!!!!!!! (DE MOMENTO)

       SUBROUTINE  TRBASIS_MULTICAPA (KXX,X,Y,Z,XOBS,YOBS,ZOBS,LOBS,
     ;                      BFOBS,AREA,LNNDEL)

***********************************************************************
***    COMPUTES ELEMENT AND INTERPOLATION FUNCTION VALUES AT 
***    OBSERVATION POINTS
***********************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       INCLUDE 'COMMON.FOR'

       DIMENSION KXX(LMXNDL,NUMEL),X(NUMNP),Y(NUMNP),AREA(NUMEL),          
     ; XOBS(NUOBS),YOBS(NUOBS),LOBS(NUOBS),BFOBS(LMXNDL,NUOBS),
     ; LNNDEL(NUMEL),BF(4),Z(NUMNP),ZOBS(NUMNP)

       ZERO=-0.0000001

***  CROSS ALL OBSERVATION POINTS

       DO 10 NOB=1,NUOBS
          XOB=XOBS(NOB)
          YOB=YOBS(NOB)
          ZOB=ZOBS(NOB)

***  FIRST SEARCH AROUND 1-D ELEMENTS

          DO 20 NEL=1,NUMEL
             AR=AREA(NEL)
             NNUD=LNNDEL(NEL)
             IF (NNUD.EQ.2) THEN
                I1=KXX(1,NEL)
                I2=KXX(2,NEL)
                X1=X(I1)
                Y1=Y(I1)
                Z1=Z(I1)
                X2=X(I2)
                Y2=Y(I2)
                Z2=Z(I2)
                XNORMA_VDIR2= (X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)
     ;                               +(Z1-Z2)*(Z1-Z2)
                D2=(X1-XOB)*(X1-XOB)+(Y1-YOB)*(Y1-YOB)
     ;                               +(Z1-ZOB)*(Z1-ZOB)
                PESC= (X2-X1)*(XOB-X1) + (Y2-Y1)*(YOB-Y1)
     ;                                 + (Z2-Z1)*(ZOB-Z1)
                XX=PESC/XNORMA_VDIR2
                DISTANCE2= D2- XX*PESC
                TOLER=XNORMA_VDIR2/1.D6                 ! 1 PER MIL
                IF (XX.GT.1-ZERO .OR. XX.LT.ZERO .OR.   ! OBSERVATION POINT
     ;                    DISTANCE2.GT.TOLER) GOTO 20   ! TOO FAR
                BF(1)=1-XX
                BF(2)=XX

***  STORES INTERPOLATION FUNCTIONS VALUES AND ELEMENT NUMBER

                LOBS(NOB)=NEL
                DO 40 I=1,NNUD
 40                BFOBS(I,NOB)=BF(I)
                GOTO 10
             END IF

***  THE POINT DOESN'T BELONG TO THIS 1-D ELEMENT, GO TO THE NEXT

 20       CONTINUE

          ZAUX=1.D20
          NELAUX=0
          DO 50 NEL=1,NUMEL
             AR=AREA(NEL)
             NNUD=LNNDEL(NEL)
             IF (NNUD.EQ.2) GOTO 50
             I1=KXX(1,NEL)
             I2=KXX(2,NEL)
             I3=KXX(3,NEL)
             CALL FUNFOR (XOB,YOB,NEL,KXX,
     ;                    X,Y,NNUD,AR,LMXNDL,NUMEL,BF,NUMNP)

*** RECTANGULAR ELEMENT

             IF (NNUD.EQ.4) THEN
                IF (BF(1).LT.ZERO .OR. BF(2).LT.ZERO .OR. 
     ;                      BF(3).LT.ZERO.OR. BF(4).LT.ZERO) GOTO 50

                TOLERZ=DSQRT(AR)/100D0
                I4=KXX(4,NEL)

                ZL=Z(I1)*BF(1)+Z(I2)*BF(2)+Z(I3)*BF(3)+Z(I4)*BF(4)

*** TRIANGULAR ELEMENT 

             ELSE IF (NNUD.EQ.3) THEN
                IF (BF(1).LT.ZERO .OR. BF(2).LT.ZERO 
     ;                    .OR. BF(3).LT.ZERO) GOTO 50

                TOLERZ=(DSQRT(AR))/25D0

                ZL=Z(I1)*BF(1)+Z(I2)*BF(2)+Z(I3)*BF(3)
             END IF
             ZNEW=DABS(ZOB-ZL)
             IF (ZNEW.LE.TOLERZ) GOTO 60   ! THE OBSERVATION POINT BELONGS TO 
                                           ! THIS ELEMENT
             IF (ZNEW.GT.TOLERZ .AND. ZAUX.GT.ZNEW) THEN
                ZAUX=ZNEW
                NELAUX=NEL
             END IF

 50       CONTINUE                         ! NEXT ELEMENT

          IF (NELAUX.NE.0) THEN       ! THE NODE IS BETWEEN TWO GRID LAYERS

***  STORES INTERPOLATION FUNCTIONS VALUES AND ELEMENT NUMBER

             WRITE(MAINF,1000) NOB,NELAUX
 1000        FORMAT (1X,'WARNING: THE OBSERVATION POINT',I5,
     ;               ' DOES NOT BELONG TO A 2-D ELEMENT.',
     ;               ' IT IS ASSIGNED TO THE NEAREST ELEMENT:',I5)
             NEL=NELAUX
             NNUD=LNNDEL(NEL)
             CALL FUNFOR (XOB,YOB,NEL,KXX,
     ;                 X,Y,NNUD,AREA(NEL),LMXNDL,NUMEL,BF,NUMNP)
 60          LOBS(NOB)=NEL
             DO 70 I=1,NNUD
 70             BFOBS(I,NOB)=BF(I)
          ELSE
             WRITE(MAINF,2000) NOB
 2000        FORMAT (1X,'ERROR: THE OBSERVATION POINT',I5,
     ;               ' DOES NOT BELONG TO THE GRID.',
     ;               ' CHECK ITS COORDINATES',I5)
             STOP
          END IF
 10    CONTINUE                            ! NEXT OBSERVATION POINT


***   WRITE THE INTERPOLATION FUNCTION VALUES AND ELEMENT FOR EVERY OBS. P.


*       IF (INPWR.NE.0) THEN         !PROV: DE MOMENTO SE ESCRIBE SIEMPRE
          WRITE(MAINF,100)
 100      FORMAT(///,
     ;          ' VALUES OF BASIS FUNCTIONS AT THE OBSERVATION POINTS'
     ;          ,//,' NODE   ELEMENT        F1        F2',
     ;          '        F3        F4')
          DO 80 NOB=1,NUOBS
             NEL=LOBS(NOB)
             NNUD=LNNDEL(NEL)
             WRITE(MAINF,200) NOB,NEL,(BFOBS(I,NOB),I=1,NNUD)
 200         FORMAT(I5,I10,4F10.3)
 80       CONTINUE
*       END IF

       RETURN 
       END

************************************************************************
************************************************************************

       SUBROUTINE FUNFOR
     ;      (XOB,YOB,NEL,KXX,X,Y,NNUD,AR1,LMXNDL,NUMEL,BF,NUMNP)

******************************************************************
***         COMPUTES INTERPOLATION FUNCTIONS
******************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION KXX(LMXNDL,NUMEL),BF(4),X(NUMNP),Y(NUMNP)

       IF (NNUD.EQ.3)THEN
         AR=2*AR1
         X1=X(KXX(1,NEL))
         Y1=Y(KXX(1,NEL))           
         X2=X(KXX(2,NEL))
         Y2=Y(KXX(2,NEL))
         X3=X(KXX(3,NEL))
         Y3=Y(KXX(3,NEL))
         BF(1)=(X2*Y3-X3*Y2+(Y2-Y3)*XOB+(X3-X2)*YOB)/AR               
         BF(2)=(X3*Y1-X1*Y3+(Y3-Y1)*XOB+(X1-X3)*YOB)/AR               
         BF(3)=(X1*Y2-X2*Y1+(Y1-Y2)*XOB+(X2-X1)*YOB)/AR 
       ELSE
         AR=-AR1
         KY=4
         KX=2
         DO 30 I=1,4
           AR=-AR
           BF(I)=(Y(KXX(KY,NEL))-YOB)*(X(KXX(KX,NEL))-XOB)/AR
           KY=KY-1
           KX=KX-1
           IF(KX.EQ.0)KX=4
 30      CONTINUE             
       ENDIF
       END

************************************************************************

