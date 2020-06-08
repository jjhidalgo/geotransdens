      SUBROUTINE GRAV_PROJ
     ;(LMXNDL   ,NUMEL    ,NUMNP    ,COORD    ,GRAV     ,GRAVEL
     ;,KXX      ,LDIM)

********************************************************************************
*
* PURPOSE  To compute the projection of the gravity vector over all the elements
*          of the mesh.It is required in order to compute the gravity flow 
*          between the nodes of the element.
*
* DESCRIPTION  A projection of the 3D global absolute gravity vector over the 
*              plane or line containing the element is done for 1D and 2D
*              elements. For 3D elements, the global gravity one is the 
*              searched vector.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COORD                  Nodal coordinates                                     
*  GRAV                   Gravity array direction                               
*  GRAVEL                                                                       
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LTYPE                  Vector containing the type of each element            
*
* EXTERNAL VARIABLES: SCALARS
*
*  LMXNDL                 Maximum number of nodes per element                   
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  LOC_COORDINADES                                                              
*
* HISTORY: First coding: German Galarza (July-97)
*          Revision: Andres Alcolea (Oct-98)
*
********************************************************************************

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       DIMENSION
     ;      KXX(LMXNDL,NUMEL)   ,GRAV(3)       ,GRAVEL(NUMEL,3)  
     ;      ,COORD(NUMNP,3)     ,LDIM(NUMEL)


       DO L=1,NUMEL
  
C______________________________Accesing to global coordinates of the
C______________________________element nodes and other constants.

         LD=LDIM(L)

         I1=KXX(1,L)
         X1=COORD(I1,1)
         Y1=COORD(I1,2)
         Z1=COORD(I1,3)
  
         I2=KXX(2,L)
         X2=COORD(I2,1)
         Y2=COORD(I2,2)
         Z2=COORD(I2,3)
         
C______________________________1D elements

         IF (LD.EQ.1) THEN

           XMOD12=(X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)+(Z2-Z1)*(Z2-Z1)
           DELX=X2-X1
           DELY=Y2-Y1
           DELZ=Z2-Z1
           
C______________________________Computes the projection of the gravity vector
C______________________________over the current element

           ALPHA=GRAV(1)*DELX+GRAV(2)*DELY+GRAV(3)*DELZ 

           GRAVEL(L,1)=ALPHA*DELX/XMOD12
           GRAVEL(L,2)=ALPHA*DELY/XMOD12
           GRAVEL(L,3)=ALPHA*DELZ/XMOD12
  
C______________________________2D elements

         ELSE IF(LD.EQ.2) THEN

           I3=KXX(3,L)
           X3=COORD(I3,1)
           Y3=COORD(I3,2)
           Z3=COORD(I3,3)

           U1_X=X2-X1
           U1_Y=Y2-Y1
           U1_Z=Z2-Z1

C-------------------------Computes a orthogonal vector
C-------------------------to U1 in the plane of the element
  
           U2_X=X3-X1
           U2_Y=Y3-Y1
           U2_Z=Z3-Z1

C------------------------- First computes an orthogonal vector
C------------------------- to the plane of the element
  
           PRODVECTX=U1_Y*U2_Z-U1_Z*U2_Y
           PRODVECTY=U1_Z*U2_X-U1_X*U2_Z
           PRODVECTZ=U1_X*U2_Y-U1_Y*U2_X

C------------------------- Now makes the vector product of U1 and 
C------------------------- PRODVECT to compute an orthogonal vector
C------------------------- to U1 in the plane of the element

           U3_X=U1_Y*PRODVECTZ-U1_Z*PRODVECTY
           U3_Y=U1_Z*PRODVECTX-U1_X*PRODVECTZ
           U3_Z=U1_X*PRODVECTY-U1_Y*PRODVECTX

C------------------------- Now makes the orthogonal projection 
C------------------------- (GRAV**U1 *U1/norm square of U1 +
C-------------------------  GRAV**U3 *U3/norm square of U3) 
C------------------------- ** is scalar product, and * is a number times
C------------------------- a vector

           XMODU3=U3_X*U3_X+U3_Y*U3_Y+U3_Z*U3_Z
           XMODU1=U1_X*U1_X+U1_Y*U1_Y+U1_Z*U1_Z

           ALPHA=GRAV(1)*U1_X+GRAV(2)*U1_Y+GRAV(3)*U1_Z
           BETA =GRAV(1)*U3_X+GRAV(2)*U3_Y+GRAV(3)*U3_Z
 

           GRAVEL(L,1)=ALPHA*U1_X/XMODU1+BETA*U3_X/XMODU3
           GRAVEL(L,2)=ALPHA*U1_Y/XMODU1+BETA*U3_Y/XMODU3
           GRAVEL(L,3)=ALPHA*U1_Y/XMODU1+BETA*U3_Y/XMODU3

C______________________________For 3D elements, the global gravity
C______________________________vector is that read in the input data.

         ELSE

           GRAVEL(L,1)=GRAV(1)
           GRAVEL(L,2)=GRAV(2)
           GRAVEL(L,3)=GRAV(3)

         ENDIF
  
       END DO

       RETURN
       END
