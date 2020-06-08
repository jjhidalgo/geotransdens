      SUBROUTINE RHS_NOLI_IN 
     &          (DERA     ,DERB     ,FACTOR   ,IODERB   ,KXX
     &          ,LMXNDL   ,LNNDEL   ,NPA      ,NUMEL    ,NUMNP
     &          ,RHSNEW   ,RHSOLD)

*****************************************************************************
*
* PURPOSE
*
*     Computes the contribution of the nonlinear part to the right hand side
*     of the inverse problem
*
* DESCRIPTION
*     Computes the contribution of the nonlinear part to the right hand side
*     of the inverse problem
*
* EXTERNAL VARIABLES: ARRAYS
*
*  DERADFT                Contains the derivatives of AFLU (or ATRA) and 
*                         DFLU (or DTRA) matrices with respect to head (or 
*                         concentrations). This term appears on the left
*                         hand side of both flow and transport direct and 
*                         inverse problems, as well as on the right hand side
*                         of flow/transport inverse problem.
*  DERHC                  Nodal head (or concentrations) derivatives with 
*                         respect to estimated parameters.
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  FACTOR                 Equal to EPSFLU (EPSTRA) or 1-EPSFLU (1-EPSTRA)
*  IDIM                   Second dimension of matrix DERADFT
*  INEW                   Index that controls the location of the right hand 
*                         side for the computation of nodal derivatives at 
*                         the current computation time
*  IOLD                   Index that controls the location of the nodal 
*                         derivatives at the previous computation time inside
*                         array DERH
*  LMXNDL                 Maximum number of nodes per element                   
*  NPA                    Number of paramaters to be estimated                                                      
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  NNUD                   Number of nodes of the current element                
*  NP                     Counter
*
* HISTORY
*
*     GGL        1988     First coding
*     AMS     10-1998     Revision
*     JHG      5-2004     IOCPAR eliminated.
*
*****************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IODERB,LMXNDL,NPA,NUMEL,NUMNP

      REAL*8::FACTOR

      INTEGER*4::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

      REAL*8::DERA(NUMEL,LMXNDL*LMXNDL),DERB(NUMNP)
     &       ,RHSNEW(NUMNP,NPA),RHSOLD(NUMNP,NPA)

C------------------------- Internal

      INTEGER*4::I,IJPOS,INODE,J,JNODE,L,NNUD

      REAL*8::DERBFACTOR,DERI
      
C------------------------- First executable statement

      DO L=1,NUMEL

          NNUD = LNNDEL(L)

          DO I=1,NNUD

              INODE = KXX(I,L)

              DO J=1,NNUD

                  JNODE = KXX(J,L)
                  IJPOS = (I-1)*NNUD + J
                  
	            DERI = FACTOR*DERA(L,IJPOS)


                  RHSNEW(INODE,1:NPA) = RHSNEW(INODE,1:NPA)
     &                                 + DERI*RHSOLD(JNODE,1:NPA)

              END DO ! J=1,NNUD

          END DO ! I=1,NNUD

      END DO ! L=1,NUMEL


	IF (IODERB.NE.0) THEN

          IF (IODERB.LT.0) THEN

	        DERBFACTOR = -FACTOR

	    ELSE

              DERBFACTOR = FACTOR

	    END IF !IODERB.LT.0


	    DO I=1,NUMNP

              DERI = DERBFACTOR*DERB(I)

              RHSNEW(I,1:NPA) = RHSNEW(I,1:NPA) - DERI*RHSOLD(I,1:NPA)

	    END DO !I=1,NUMNP

	END IF !IODERB.EQ.1

      END SUBROUTINE RHS_NOLI_IN 
