      SUBROUTINE ASSEMBLE_SYM_ND_INTO_SHUFFLED_SPARSE
     ;(FACTOR    ,IA_COLS   ,IA_ROWS    ,LBLOCK_NR ,LNMXDL    
     ;,NB        ,NN        ,NUMEL      ,A         ,A_DSC    
     ;,KXX       ,LNNDEL    ,IAD_D      ,IADN_D    ,IADD_D) 
   
  

*******************************************************************
*  PURPOSE
*  To assemble an element-wise stored matrix without diagonal into a sparse matrix,
*  the rows and columns of which have been shuffled to decrease bandwidth.
*
*   EXTERNAL VARIABLES, SCALARS
* NUMEL
* LNMXDL
* LBLOCK_NR: the block of the non-shufled matrix where the element-wise matrix
*            should go. 
*            if 1, the block of derivatives of flow eqn to head
*            if 2, the block of derivatives of flow eqn to mass fraction
*            if 3, the block of derivatives of tpt eqn to head
*            if 4, the block of derivatives of tpt eqn to mass fraction
* 
* IDIM
******************************************************************
 
      IMPLICIT NONE

! EXTERNAL VARIABLES, SCALARS
      INTEGER*4 IA_COLS,IA_ROWS,LBLOCK_NR ,LNMXDL ,NB,NN ,NUMEL
      REAL*8 FACTOR

! EXTERNAL VARIABLES, ARRAYS
      INTEGER*4 KXX(LNMXDL,NUMEL),LNNDEL(NUMEL),IAD_D(NB,NN),IADD_D(NN)
     ;        ,IADN_D(NN)
      REAL*8 A(IA_ROWS,IA_COLS),A_DSC(NB,NN) 

! INTERNAL VARIABLES, SCALARS
      INTEGER*4 IJROW,JIROW,IJCOL,JICOL, KIJ,INODE,JNODE
     ;          ,J,L,NNUD,KJI,I,IPOS
      REAL*8  VALUE  

      
      DO L=1,NUMEL                !loop over elements
	   NNUD=LNNDEL(L)             

          DO I=1,NNUD-1
              
			INODE= KXX(I,L)  
                    
              DO J=I+1,NNUD
	           
                  JNODE= KXX(J,L)	      

                  IPOS = (I - 1)*NNUD + J - I*(I+1)/2
		        VALUE=A(L,IPOS)
                  VALUE=VALUE*FACTOR

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


			 CALL FIND( IJROW,IJCOL,KIJ,IAD_D,IADN_D,NB,NN)
			 CALL FIND( JIROW,JICOL,KJI,IAD_D,IADN_D,NB,NN)          
 

			 !ading element ij
			 A_DSC(KIJ,IJROW)= A_DSC(KIJ,IJROW) + VALUE

			 !adding element ji
			 A_DSC(KJI,JIROW)= A_DSC(KJI,JIROW) + VALUE

			 !updating diagonal
			 A_DSC(IADD_D(IJROW),IJROW) = 
     ;                A_DSC(IADD_D(IJROW),IJROW)- VALUE
			 A_DSC(IADD_D(JIROW),JIROW) = 
     ;                A_DSC(IADD_D(JIROW),JIROW)- VALUE   
           ENDDO 
        ENDDO  !J
      ENDDO    !L
	
	RETURN
	END  
            




