      SUBROUTINE COMP_BUOYANCY_CF
     &          (BETAC      ,BUOYANCY   ,CREF       ,CAUX1
     &          ,COORD      ,DBUOYANCY  ,LTYP       ,GRADLOC
     &          ,GRAVEL     ,IDIM       ,IPARTNER   ,IOCALCDEVT
     &          ,KXX        ,IODIM      ,L          ,LMXNDL
     &          ,MAXPG      ,NNUD       ,NUMEL      ,NUMNP)

********************************************************************************
*
*     PURPOSE
*         Calculate the coefficients needed to calculate
*         bouyancy term consistently
*
*      IPARTNER: contains for every node the local node nr of the node that
*               should be used with this node in the integral term.
*
*      GRADLOC: GRADLOC(j,i,k) contains the gradient in local coordinates 
*                 of shape function j in direction i evaluted in node k.
*                 The notation in inline comments is:
*                  @Ni/@Xj (k)
*                 Where Xj = X, Y, Z (J=1,2,3, respectively)
*                 i.e. Derivative of shape function "i", with respect
*                 direction "Xj", evaluates in node "k".
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

C--------------------  EXTERNAL VARIABLES: SCALARS

      REAL*8::BETAC,CREF

      INTEGER*4::LTYP,IDIM,IODIM,L,LMXNDL,NNUD,NUMEL,NUMNP
     &          ,IOCALCDEVT,MAXPG

C--------------------  EXTERNAL VARIABLES: ARRAYS

      REAL*8::CAUX1(NUMNP),BUOYANCY(IODIM,LMXNDL,NUMEL),GRAVEL(NUMEL,3)
     &       ,DBUOYANCY(IODIM,LMXNDL*LMXNDL,NUMEL),COORD(NUMNP,3)
     &       ,GRADLOC(IODIM,LMXNDL,MAXPG)
      INTEGER*4::IPARTNER(6,3,6),KXX(LMXNDL,NUMEL)

C--------------------  INTERNAL VARIABLES: SCALARS

      REAL*8::BETAC2,BYNC,DBYNCII,DBYNCIJ,EXP_I,EXP_J,EXPNT
     &       ,FRAC,DENSINT,CI,CI2,CJ,CJ2,FRAC1,FRAC2,DLDCI,DLDCJ

      INTEGER*4::I_DIM,II_POS,IJ_POS,INODE,J_DIM,JNODE,KNODE

C--------------------  INTERNAL VARIABLES: ARRAYS

      REAL*8::XLOC(NNUD),YLOC(NNUD),ZLOC(NNUD)


      BUOYANCY(:,:,L) = 0D0
      DBUOYANCY(:,:,L) = 0D0

C-------------------- Initializes nodes local coordinates
C-------------------- and local gradient in the nodes

      XLOC(1:NNUD) = COORD(KXX(1:NNUD,L),1)
      YLOC(1:NNUD) = COORD(KXX(1:NNUD,L),2)
      ZLOC(1:NNUD) = COORD(KXX(1:NNUD,L),3)

C-------------------- Computes local gradient of shape functions
C-------------------- in the element nodes

      CALL COMP_GRAD_LOC
     &    (GRADLOC  ,IODIM    ,LMXNDL   ,LTYP     ,NNUD     ,XLOC
     &    ,YLOC     ,ZLOC)

C-------------------- and the mask for null buoyancy terms.
C-------------------- (if -1, this coefficient is always null)

      SELECT CASE(LTYP)

          CASE(1)  ! 1D elements

              BUOYANCY(1,1,L) = -1D0

          CASE(2,5) ! Linear 2D triangles, 3 nodes.    
                
              BUOYANCY(1,1,L) = -1D0   
              BUOYANCY(2,1,L) = -1D0
              BUOYANCY(2,2,L) = -1D0
              BUOYANCY(1,3,L) = -1D0

          CASE(3) ! Quadrilaterales
 
              BUOYANCY(1,1,L) = -1D0
              BUOYANCY(2,1,L) = -1D0
              BUOYANCY(2,2,L) = -1D0
              BUOYANCY(1,4,L) = -1D0

          CASE(4)  ! Tetrahedricons

              BUOYANCY(1,1,L) = -1D0
              BUOYANCY(2,1,L) = -1D0
              BUOYANCY(3,1,L) = -1D0
              BUOYANCY(2,2,L) = -1D0
              BUOYANCY(3,2,L) = -1D0
              BUOYANCY(1,3,L) = -1D0
              BUOYANCY(3,3,L) = -1D0
              BUOYANCY(1,4,L) = -1D0
              BUOYANCY(2,4,L) = -1D0

          CASE(6) ! Prismatics

              BUOYANCY(1,1,L) = -1D0
              BUOYANCY(2,1,L) = -1D0
              BUOYANCY(3,1,L) = -1D0
              BUOYANCY(1,2,L) = -1D0
              BUOYANCY(3,2,L) = -1D0
              BUOYANCY(1,3,L) = -1D0
              BUOYANCY(2,3,L) = -1D0
              BUOYANCY(2,4,L) = -1D0
              BUOYANCY(3,4,L) = -1D0
              BUOYANCY(3,5,L) = -1D0
              BUOYANCY(2,6,L) = -1D0

      END SELECT   
      

C-------------------- Computations

      DO INODE=1,NNUD
         
          CI = CAUX1(KXX(INODE,L))
             
          DO I_DIM=1,IDIM

C-------------------- If this component is not always equal to zero

              IF (BUOYANCY(I_DIM,INODE,L).GE.0) THEN

	            BYNC = 0D0
             
C-------------------- if we are not dealing with a rectangular prism !?

                  JNODE = IPARTNER(INODE,I_DIM,LTYP)
                  CJ = CAUX1(KXX(JNODE,L))          

C-------------------- Calculate directional derivative of gravity in the I_DIM
C-------------------- direction

                  DO KNODE=1,NNUD  !loop over shape functions

                      DO J_DIM=1,IDIM    !loop over dimensions
                          
                          BYNC = BYNC
     &                         - GRAVEL(L,J_DIM)
     &                           *COORD(KXX(KNODE,L),J_DIM)
     &                           *GRADLOC(I_DIM,KNODE,INODE)
                      END DO !DIMJ 
                  END DO !KNODE

 
                  IF (IOCALCDEVT .EQ.1) THEN
                      DBYNCII = BYNC 
                      DBYNCIJ = BYNC 
                  END IF
              
C-------------------- Computation of density integral.
C-------------------- The integral is made in the interval (0,1). It terms of
C-------------------- concentration it means betwen the interval (C0,C1),i.e.(CJ,CI),
C-------------------- since the node whose coefficient is non-zero is always
C-------------------- located at one in coordinate 1.
C-------------------- The analitical expresion for the integral is
C-------------------- [1/BETA*(C0-C1)]*[EXP(BETA*(C0-CREF) - EXP(BETA*(C1-CREF)]
C-------------------- In terms of C_I and C_J it becames
C-------------------- [1/BETA*(CJ-CI)]*[EXP(BETA*(CJ-CREF) - EXP(BETA*(CI-CREF)]
C-------------------- However this is only important when BETA*(C0-C1)<<1, since
C-------------------- the expresion above is symmetrical, i.e., f(CI,CJ) = f(CJ,CI).

                  EXPNT = BETAC*(CJ-CI)

C-------------------- If there is no density gradient in the current direction:

                  IF (DABS(EXPNT).LT.1D-8) THEN 

                      EXP_I = EXP(BETAC*(CI - CREF))

                      DENSINT = EXP_I*(1+EXPNT/2.+EXPNT*EXPNT/6.) - 1D0
                          
                      BYNC=BYNC*DENSINT

C-------------------- Calculate derivatives of integral term 

                      IF(IOCALCDEVT.EQ.1) THEN

C-------------------- Derivative w. r. t. CI

                          BETAC2 = BETAC*BETAC
                          CJ2 = CJ*CJ
                          CI2 = CI*CI

                          DENSINT = BETAC*EXP_I
     &                              *(BETAC2*(CJ2+CI2) + BETAC*(CJ-CI)
     &                                - 2D0*BETAC2*CI*CJ + 3D0)/6D0

                          DBYNCII = DBYNCII*DENSINT

C-------------------- Derivative w. r. t. CJ

	                    DENSINT = BETAC*EXP_I*(2D0*BETAC*(CJ-CI)+3D0)
     &                             /6D0
                          DBYNCIJ = DBYNCIJ*DENSINT

                      END IF 

                  ELSE

C-------------------- If there is a density gradient in the current direction
                     
                      EXP_I = EXP(BETAC*(CI - CREF))
                      EXP_J = EXP(BETAC*(CJ - CREF))
                      FRAC = 1D0/(BETAC*(CJ - CI))
                      DENSINT = FRAC*(EXP_J - EXP_I)

                      BYNC = BYNC*(DENSINT - 1D0)
                    
C-------------------- Calculate derivatives of integral term 

                      IF (IOCALCDEVT.EQ.1) THEN
                    
                          FRAC1 = 1D0/(BETAC*((CI-CJ)*(CI-CJ)))
                          FRAC2 = 1D0/(CI-CJ)

                          DLDCI = -FRAC1*(EXP_I-EXP_J)+FRAC2*EXP_I
                          DLDCJ =  FRAC1*(EXP_I-EXP_J)-FRAC2*EXP_J 

                          DBYNCII = DBYNCII*DLDCI

                          DBYNCIJ = DBYNCIJ*DLDCJ
                    
                      END IF  

                  END IF !(IF (ABS(CI-CJ) .LT. 1E-6) THEN )
              
                  II_POS = (INODE - 1)*NNUD + INODE
                  IJ_POS = (INODE - 1)*NNUD + JNODE

                  DBUOYANCY(I_DIM,II_POS,L) = DBYNCII
                  DBUOYANCY(I_DIM,IJ_POS,L) = DBYNCIJ
              ELSE

	            BYNC = 0D0

              END IF !BUOYANCY.GE.0 (BUOYANCY(I_DIM,INODE,L))

              BUOYANCY(I_DIM,INODE,L) = BYNC              

          END DO ! Loop over dimensions
      END DO     ! loop over nodes    


      END SUBROUTINE COMP_BUOYANCY_CF
