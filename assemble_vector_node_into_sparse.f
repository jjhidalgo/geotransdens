      SUBROUTINE ASSEMBLE_VECTOR_NODE_INTO_SPARSE
     &(A_DSC    ,FACTOR   ,IADD_S     ,NB    ,NN     ,VECTOR)
           


      IMPLICIT NONE
* EXTERNAL VARIABLES: SCALARS
      INTEGER*4  NB,NN           
      REAL*8 FACTOR

* EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 IADD_S(NN)
      REAL*8 VECTOR(NN),A_DSC(NB,NN)

* INTERNAL VARIABLES: SCALARS
      INTEGER*4 INODE
      REAL*8 VALUE

* INTERNAL VARIABLES: ARRAYS
      
      DO INODE = 1,NN
         VALUE=VECTOR(INODE)
         VALUE = VALUE * FACTOR
         A_DSC(IADD_S(INODE),INODE)= A_DSC(IADD_S(INODE),INODE)+VALUE
      ENDDO
      
      RETURN
      END    

*********************************************************************
*********************************************************************

      SUBROUTINE ASSEMBLE_VECTOR_ELEM_INTO_SPARSE
     ;(FACTOR,LMXNDL,NB,NN,NUMEL,A,A_DSC,IADD_S,KXX,LNNDEL)

	IMPLICIT NONE

* EXTERNAL VARIABLES: SCALARS
      INTEGER*4 NB,NN,NUMEL,LMXNDL
      REAL*8 FACTOR

* EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 IADD_S(NN),KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)
      REAL*8 A(NUMEL,LMXNDL),A_DSC(NB,NN)

* INTERNAL VARIABLES: SCALARS
      INTEGER*4 L,I,INODE,NNUD
      REAL*8 VALUE

* INTERNAL VARIABLES: ARRAYS


	DO L=1,NUMEL
	   NNUD = LNNDEL(L)
	   DO I=1, NNUD
	      INODE = KXX(I,L)
	      VALUE = A(L,I)
            VALUE = VALUE*FACTOR
	      A_DSC(IADD_S(INODE),INODE)=A_DSC(IADD_S(INODE),INODE)+VALUE
	   ENDDO
	ENDDO
      
	RETURN
	END	   