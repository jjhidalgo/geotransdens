	SUBROUTINE ASSEMBLE_DERIVATIVE_INTO_SHUFFLED_SPARSE
     ;(FACTOR     ,IA_COLS    ,IA_ROWS    ,LMXNDL     ,LBLOCK_NR  
     ;,NB         ,NN         ,NUMEL      ,A          ,A_DSC      
     ;,IAD_D      ,IADN_D     ,KXX        ,LNNDEL)

  
	IMPLICIT NONE

! EXTERNAL VARIABLES, SCALARS
      INTEGER*4 IA_COLS,IA_ROWS,NUMEL,LMXNDL,LBLOCK_NR,NB,NN    
      REAL*8 FACTOR 

! EXTERNAL VARIABLES, ARRAYS
      INTEGER*4 IAD_D(NB, NN), IADN_D(NN),KXX(LMXNDL,NUMEL)
     ;         ,LNNDEL(NUMEL)
      REAL*8 A(IA_ROWS,IA_COLS),A_DSC(NB,NN) 

! INTERNAL VARIABLES, SCALARS
      INTEGER*4 IIROW,IICOL,IJROW,IJCOL,JIROW,JICOL,INODE, JNODE
     ;         ,KII, KIJ, KJI,L,I,J,NNUD
      REAL*8   DI_DJ,DJ_DI



       
      DO L=1,NUMEL                           !Loop over elements

         NNUD =LNNDEL(L)

c_________________________________________1) nondiagonal elements

	   DO I= 1,NNUD-1                        !loop over nodes of element
	      INODE = KXX(I,L)

	      DI_DJ= FACTOR * A(L,I)               

		  DO J=I+1,NNUD                      !loop over nodes of element
	         JNODE = KXX(J,L)
	         DJ_DI =FACTOR * A(L,J) 

               SELECT CASE(LBLOCK_NR)

			      CASE (1)  !block dfludflu
	                  IJROW = (INODE-1)*2+1
	                  IJCOL = (JNODE-1)*2+1
	                  JIROW = (JNODE-1)*2+1
	                  JICOL = (INODE-1)*2+1

	              CASE (2)  !block dfludtra
	                  IJROW = (INODE-1)*2+1
	                  IJCOL = JNODE*2  
	                  JIROW = (JNODE-1)*2+1
	                  JICOL = INODE*2


	              CASE (3)  !block dtradflu
	                  IJROW = INODE*2   
	                  IJCOL = (JNODE-1)*2+1
	                  JIROW = JNODE*2
	                  JICOL = (INODE-1)*2+1


	              CASE(4)   !block dtradtra
	                  IJROW = INODE*2
	                  IJCOL = JNODE*2
	                  JIROW = JNODE*2
	                  JICOL = INODE*2


           	 END SELECT
 
               !sum the element of the input matrix to Aij
               CALL FIND( IJROW,IJCOL,KIJ,IAD_D,IADN_D,NB,NN)
               A_DSC(KIJ,IJROW) = A_DSC(KIJ,IJROW) + DI_DJ
               

			 CALL FIND( JIROW,JICOL,KJI,IAD_D,IADN_D,NB,NN)
			 A_DSC(KJI,JIROW) = A_DSC(KJI,JIROW) + DJ_DI

	     ENDDO !J
	  ENDDO    !I


c_________________________________________2) diagonal elements

	   DO I= 1,NNUD              !loop over nodes of element
	      INODE = KXX(I,L)

	      DI_DJ= FACTOR * A(L,I)               

            SELECT CASE(LBLOCK_NR)

			 CASE (1)  !block dfludflu
	              IIROW = (INODE-1)*2+1
	              IICOL = (INODE-1)*2+1

	          CASE (2)  !block dfludtra
	              IIROW = (INODE-1)*2+1
	              IICOL = INODE*2  


	          CASE (3)  !block dtradflu
	              IIROW = INODE*2   
	              IICOL = (INODE-1)*2+1

	          CASE(4)   !block dtradtra
	              IIROW = INODE*2
	              IICOL = INODE*2

           END SELECT

           !sum the element of the input matrix to Aij
           CALL FIND( IIROW,IICOL,KII,IAD_D,IADN_D,NB,NN)
           A_DSC(KII,IIROW) = A_DSC(KII,IIROW) + DI_DJ
           

	 ENDDO !I

	ENDDO !L


	RETURN
	END