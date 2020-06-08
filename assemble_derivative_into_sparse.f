	SUBROUTINE ASSEMBLE_DERIVATIVE_INTO_SPARSE
     ;(FACTOR     ,IA_COLS    ,NB         ,NUMEL      ,NN      
     ;,LMXNDL     ,A          ,A_DSC      ,IAD_S      ,IADN_S     
     ;,KXX        ,LNNDEL)

  
	IMPLICIT NONE
! EXTERNAL VARIABLES, SCALARS
      INTEGER*4 IA_COLS,LMXNDL,NB
     ;         ,NN,NUMEL
      
      REAL*8 FACTOR 

! EXTERNAL VARIABLES, ARRAYS
      INTEGER*4 IAD_S(NB,NN),IADN_S(NN),KXX(LMXNDL,NUMEL)
     ;         ,LNNDEL(NUMEL)
      REAL*8 A(NUMEL,IA_COLS),A_DSC(NB,NN) 

! INTERNAL VARIABLES, SCALARS
      INTEGER*4 INODE, JNODE,L,I,J,NNUD,KIJ,KJI,KII
      REAL*8   DI_DJ
     


      DO L=1,NUMEL

         NNUD =LNNDEL(L)
c________________________________Nondiagonal elements

	   DO I= 1,NNUD-1
	      INODE = KXX(I,L)
	      DI_DJ= FACTOR * A(L,I)
               
		  DO J=I+1,NNUD
	         JNODE = KXX(J,L)
	          
	         !getting column index
               CALL FIND(INODE,JNODE,KIJ,IAD_S,IADN_S,NB,NN)
			 CALL FIND(JNODE,INODE,KJI,IAD_S,IADN_S,NB,NN)

	         !summing contribution
               A_DSC(KIJ,INODE) = A_DSC(KIJ,INODE) + DI_DJ
	         A_DSC(KJI,JNODE) = A_DSC(KJI,JNODE) + DI_DJ

	     ENDDO !J
	  ENDDO    !I

c________________________________Diagonal elements
        DO I = 1,NNUD
		 INODE = KXX(I,L)
	     DI_DJ= FACTOR * A(L,I)
	     CALL FIND(INODE,INODE,KII,IAD_S,IADN_S,NB,NN)
	     A_DSC(KII,INODE)=A_DSC(KII,INODE)+DI_DJ
	  ENDDO


	ENDDO      !L


	RETURN
	END