	SUBROUTINE ASSEMBLE_CROSSTERM_INTO_SPARSE
     ;(FACTOR     ,IA_COLS    ,NB         ,NUMEL      ,NN      
     ;,LMXNDL     ,A          ,A_DSC      ,IAD_S      ,IADD_S
     ;,IADN_S     ,KXX        ,LNNDEL)

  
	IMPLICIT NONE
! EXTERNAL VARIABLES, SCALARS
      INTEGER*4 IA_COLS,LMXNDL,NB,NN,NUMEL
      REAL*8 FACTOR 

! EXTERNAL VARIABLES, ARRAYS
      INTEGER*4 IAD_S(NB,NN),IADN_S(NN),KXX(LMXNDL,NUMEL)
     ;         ,LNNDEL(NUMEL),IADD_S(NN)
      REAL*8 A(NUMEL,IA_COLS),A_DSC(NB,NN) 

! INTERNAL VARIABLES, SCALARS
      INTEGER*4 INODE, JNODE,L,I,J,NNUD,KIJ,KJI,KII
	REAL*8 ELEM_IJ,ELEM_JI,ELEM_II
      
! INTERNAL VARIABLES, ARRAYS
	REAL*8, DIMENSION(:), ALLOCATABLE:: ROWVECTOR 
	

	!allocate rowvector
	ALLOCATE(ROWVECTOR(LMXNDL))

      !loop over elements
	DO L=1,NUMEL
		NNUD=LNNDEL(L)
 
		!loop over nodes of element
 	    !fill local vector with row and calculate missing element
	    ROWVECTOR=0D0
		DO I=2,NNUD
			ROWVECTOR(1) = ROWVECTOR(1)-A(L,I-1)
              ROWVECTOR (I) = A(L,I-1)
	    ENDDO
		
		!assemble matrix to nondiagonal terms	
          DO I=1,NNUD-1
			INODE = KXX(I,L)
			DO J=I,NNUD
				JNODE = KXX(J,L)
				
				!search the values to be assembled at ij and ji
				ELEM_IJ = ROWVECTOR(J)*FACTOR
				ELEM_JI = ROWVECTOR(I)*FACTOR
				
	            !getting column index of both
                  CALL FIND(INODE,JNODE,KIJ,IAD_S,IADN_S,NB,NN)
			    CALL FIND(JNODE,INODE,KJI,IAD_S,IADN_S,NB,NN)

                  !and assemble 
                  A_DSC(KIJ,INODE)= A_DSC(KIJ,INODE) + ELEM_IJ 
                  A_DSC(KJI,JNODE)= A_DSC(KJI,JNODE) + ELEM_JI

			enddo
	    enddo
		
		!assemble matrix to diagonal terms
		DO I=1,NNUD
             ELEM_II= ROWVECTOR(I)*FACTOR
             INODE = KXX(I,L)
		   KII = IADD_S(INODE)
		   A_DSC(KII,INODE)= A_DSC(KII,INODE)+ ELEM_II 
		ENDDO  

	ENDDO !loop over elements
      
	!Finally destroy local array
      DEALLOCATE(ROWVECTOR)

      RETURN
	END