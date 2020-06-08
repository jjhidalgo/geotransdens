      SUBROUTINE ASSEMBLE_SYM_NO_DIAG_INTO_SPARSE
     ;(A           ,A_DSC       ,FACTOR       ,KXX          ,IA_COLS   
     ;,IAD_S       ,IADD_S      ,IADN_S       ,LMXNDL       ,LNNDEL   
     ;,NB          ,NN          ,NUMEL)
     
************************************************************************
*   PURPOSE
*   To add the content of a symetric matrix to an assembled matrix in
*   sparse storage. The diagonal elements of the matrix are equal to minus
*   the sum of the elements on tis row and therefore are not stored in the
*   incomming array.
*
*   
*
*   DESCRIPTION
*   Step 1: loop over elements
*   Step 2: loop over element entries
*   Step 3: placing element entries
*
*
***********************************************************************



      IMPLICIT NONE
* EXTERNAL VARIABLES: SCALARS
      INTEGER*4 NUMEL,IA_COLS, LMXNDL,NB,NN
	REAL*8 FACTOR

* EXTERNAL VARIABLES: ARRAYS
      INTEGER  KXX(LMXNDL,NUMEL),LNNDEL(NUMEL),IAD_S(NB,NN),IADN_S(NN)
     ;        ,IADD_S(NN)
      REAL*8 A(NUMEL,IA_COLS)
     ;               ,A_DSC(NB,NN)

* INTERNAL VARIABLES: SCALARS
      REAL*8  VALUE                     
	INTEGER*4 I,J,L,INODE,JNODE,ICOL ,NNUD,IPOS

      
C____________________________________Step 1: loop over elements
	DO L=1,NUMEL

	   NNUD=LNNDEL(L)

C____________________________________Step 2: loop over element entries
         DO I=1,NNUD-1
	      INODE = KXX(I,L)
	      DO J= I+1,NNUD
	         JNODE = KXX(J,L)
	          
               IPOS =(I - 1)*NNUD + J - I *(I +1)/2

	      
			  !identifying nodes
			  VALUE=A(L,IPOS)*FACTOR

C____________________________________Step 3: placing element entries
			  !placing element (gnode1,gnode2)
			  CALL FIND(INODE,JNODE,ICOL,IAD_S,IADN_S,NB,NN)
			  A_DSC(ICOL,INODE)=A_DSC(ICOL,INODE)+VALUE
             
			  !updating diagonal
			  ICOL=IADD_S(INODE)
			  A_DSC(ICOL,INODE)=A_DSC(ICOL,INODE)-VALUE

			  !placing element (gnode2,gnode1)
			  CALL FIND(JNODE,INODE,ICOL,IAD_S,IADN_S,NB,NN)
			  A_DSC(ICOL,JNODE)=A_DSC(ICOL,JNODE)+VALUE
             
			  !updating diagonal
			  ICOL=IADD_S(JNODE)
			  A_DSC(ICOL,JNODE)=A_DSC(ICOL,JNODE)-VALUE

             ENDDO

        ENDDO          ! loop over element entries
      ENDDO            ! loop over elements
          
	RETURN
	END