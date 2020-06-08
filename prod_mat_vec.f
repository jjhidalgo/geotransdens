       SUBROUTINE PROD_MAT_VEC
     ;(FACTOR        ,IAD            ,IADN           
     ;,IDIMMAT_COLS  ,IDIMMAT_ROWS   ,IDIMVEC_COLS   
     ;,IDIMVEC_ROWS  ,I_TYPE_MAT     ,I_TYP_VEC      ,LMXNDL         
     ;,NUMEL         ,NUMNP          ,KXX            ,LNNDEL         
     ;,MATRIX        ,TARGETVEC      ,VECTOR)

******************************************************************
*
* PURPOSE: 
* Used to multiply a matrix with a vector 
* 
*
* ACRONYM: 
* Adds the PRODuct of a MATrix and a VECtor to a vector
* 
* DESCRIPTION:
* The subroutine works for matrices that are stored elementwise and
* vectors that are stored nodewise.
* Step 1: Multiply Condenced (diagonal) matrix with vector
* Step 2: Multiply symetric matrix with vector
* Step 3: Multiply full matrix with vector
*
* VARIABLES
*   I_TYPE_MAT   type of matrix
*                             1 -->   nodewise Vector.
*                             2 -->   elementwise vector
*                             3 -->   derivative type matrix
*                             4 -->   Full matrix (elementwise).
*                             5 -->   Symmetric matrix (elementwise).
*                             6 -->   Symmetric matrix without diagonal (elementwise).
*                                     (It is supposed tha the sum of the terms
*                                     out of the diagonal is equal to the diagonal with
*                                     the opposite sign).
*                             7 -->   Symmetric banded matrix.
*                             8 -->   Non symmetric banded matrix.
*                             9 -->   Sparse matrix (as the one used by WatSolve).
*                            10 -->   Elementwise matrix whose rows are identical and the first
*                                     element equal to minus sum of the rest
* I_TYP_VEC     type of vector: types 1 and 2 of i_typ_mat
* 
* HISTORY: Programmed by G.Galarza at november of 1997.
******************************************************************


      IMPLICIT NONE

C_____External variables: scalars
      INTEGER*4 IDIMMAT_COLS,IDIMMAT_ROWS,IDIMVEC_COLS,IDIMVEC_ROWS
     ;         ,I_TYPE_MAT,I_TYP_VEC,LMXNDL,NUMEL,NUMNP  
      REAL*8   FACTOR        

C_____External variables: arrays
      INTEGER*4 KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)
     ;         ,IAD(IDIMMAT_ROWS,IDIMMAT_COLS)
     ;         ,IADN(IDIMMAT_COLS)
	REAL*8 MATRIX(IDIMMAT_ROWS,IDIMMAT_COLS)
     ;      ,TARGETVEC(NUMNP),VECTOR(IDIMVEC_ROWS,IDIMVEC_COLS)
     ;      

C_____Internal variables: scalars
      INTEGER*4 L,I,J,POSIJ,NNUD,NODEI,NODEJ,NODE1,
     &          NBAND,NBAND1,NBAND2,INI,IFIN,KI

 

**********************************NODEWISE VECTOR******************************     
      IF (I_TYP_VEC .EQ.1) THEN	


		SELECT CASE (I_TYPE_MAT)

			CASE(2)

C___________________________________Part 1:Condenced scheme for storage matrix 
      			DO L=1,NUMEL                 !Loop over elements
 				   NNUD=LNNDEL(L)            
				   DO I=1,NNUD              !loop over nodes of element
					  NODEI=KXX(I,L)
					  TARGETVEC(NODEI)=TARGETVEC(NODEI)+MATRIX(L,I)
     ;                    *VECTOR(1,NODEI)*FACTOR
				   ENDDO
				ENDDO


		  CASE(4) 
C___________________________________Part 2 :full  matrix 
				DO L=1,NUMEL                 !Loop over elements
					NNUD=LNNDEL(L)
					DO I=1,NNUD              !loop over nodes of elements
					   NODEI=KXX(I,L)
					   DO J=1,NNUD           !loop over nodes of elements
						   NODEJ= KXX(J,L)
						   POSIJ= (I-1)*NNUD + J
						   TARGETVEC(NODEI)=TARGETVEC(NODEI)
     ;                         +FACTOR*MATRIX(L,POSIJ)*VECTOR(1,NODEj)
					   ENDDO
					ENDDO
				ENDDO

		   CASE(5)

C___________________________________Part 3 :symetric  matrix 
				DO L=1,NUMEL                 !Loop over elements
					NNUD=LNNDEL(L)
					DO I=1,NNUD              !loop over nodes of elements
					   NODEI=KXX(I,L)
					   DO J=1,I              !loop over nodes of elements
						   NODEJ= KXX(J,L)
						   POSIJ= ((I-1)*I)/2+J
						   TARGETVEC(NODEI)=TARGETVEC(NODEI)
     ;                          +FACTOR*MATRIX(L,POSIJ)*VECTOR(1,NODEI)
						   TARGETVEC(NODEJ)=TARGETVEC(NODEJ)
     ;                          +FACTOR*MATRIX(L,POSIJ)*VECTOR(1,NODEJ)
					   ENDDO
					ENDDO
				ENDDO

             CASE (6) 

				DO L=1,NUMEL                 !Loop over elements
					NNUD=LNNDEL(L)
					DO I=1,NNUD-1              !loop over nodes of elements
					   NODEI=KXX(I,L)
					   DO J=I+1,NNUD              !loop over nodes of elements
						   NODEJ= KXX(J,L)
						   POSIJ= (I - 1)*NNUD + J - I*(I+1)/2
						   TARGETVEC(NODEI)=TARGETVEC(NODEI)
     ;         +FACTOR*MATRIX(L,POSIJ)*(VECTOR(1,NODEJ)-VECTOR(1,NODEI))
						   TARGETVEC(NODEJ)=TARGETVEC(NODEJ)
     ;         +FACTOR*MATRIX(L,POSIJ)*(VECTOR(1,NODEI)-VECTOR(1,NODEJ))
					   ENDDO
					ENDDO
				ENDDO

             CASE (8)
                  NBAND2=IDIMMAT_COLS
	            NBAND=(NBAND2-1)/2
	            NBAND1=NBAND+1
				DO I=1,IDIMMAT_ROWS
					INI=MAX(1,NBAND1+1-I)
					KI=MAX(1,I-NBAND)
					IFIN=MIN(NBAND2,NBAND1+NUMNP-I)
				  DO J=INI,IFIN
					TARGETVEC(I)=TARGETVEC(I)
     &                            +MATRIX(I,J)*VECTOR(1,KI)*FACTOR
					KI=KI+1
	                !MATRIX(IDIMMAT_ROWS,IDIMMAT_COLS)
			      END DO
                 END DO

		   CASE(9)  !Sparse matrix
                DO I=1,IDIMMAT_COLS
				DO J = 1,IADN(I)
					NODEJ = IAD(J,I)
	                TARGETVEC(I)=TARGETVEC(I)+
     ;                      MATRIX(J,I)*VECTOR(1,NODEJ)
	            ENDDO
	          ENDDO

		  CASE (10) !Elementwise matrix whose rows are identical and the first
	                !element equal to minus sum of the rest


			DO L=1,NUMEL

				NNUD=LNNDEL(L)
	            NODE1=KXX(1,L)

				DO I=1,NNUD

					NODEI=KXX(I,L)

					DO J=2,NNUD

						NODEJ=KXX(J,L)

						TARGETVEC(NODEI) = TARGETVEC(NODEI)
     &	                + MATRIX(L,J-1)*FACTOR*
     &                     (VECTOR(1,NODEJ)-VECTOR(1,NODE1))

					END DO !I
	            END DO !J
	        END DO !L


C__________________________________Part 4. symetric matrix without diagonal
                   
		END SELECT





*******************************ELEMENTWISE VECTOR*******************************      
	ELSEIF (I_TYP_VEC .EQ.2 ) THEN
		
		SELECT CASE (I_TYPE_MAT)

			CASE(2)

C___________________________________Part 5:Condenced scheme for storage matrix 
      			DO L=1,NUMEL                 !Loop over elements
 				   NNUD=LNNDEL(L)            
				   DO I=1,NNUD              !loop over nodes of element
					  NODEI=KXX(I,L)
					  TARGETVEC(NODEI)=TARGETVEC(NODEI)+MATRIX(L,1)
     ;                  *FACTOR*VECTOR(L,I)
				   ENDDO
				ENDDO


		  CASE(4) 
C___________________________________Part 6 :full  matrix 
				DO L=1,NUMEL                 !Loop over elements
					NNUD=LNNDEL(L)
					DO I=1,NNUD              !loop over nodes of elements
					   NODEI=KXX(I,L)
					   DO J=1,NNUD           !loop over nodes of elements
						   NODEJ= KXX(J,L)
						   POSIJ= (I-1)*NNUD + J
						   TARGETVEC(NODEI)=TARGETVEC(NODEI)
     ;                         +FACTOR*MATRIX(L,POSIJ)*VECTOR(L,J)
					   ENDDO
					ENDDO
				ENDDO

		   CASE(5)

C___________________________________Part 7 :symetric  matrix 
				DO L=1,NUMEL                 !Loop over elements
					NNUD=LNNDEL(L)
					DO I=1,NNUD              !loop over nodes of elements
					   NODEI=KXX(I,L)
					   DO J=1,I              !loop over nodes of elements
						   NODEJ= KXX(J,L)
						   POSIJ= ((I-1)*I)/2+J
						   TARGETVEC(NODEI)=TARGETVEC(NODEI)
     ;                         +FACTOR*MATRIX(L,POSIJ)*VECTOR(L,I)
						   TARGETVEC(NODEJ)=TARGETVEC(NODEJ)
     ;                         +FACTOR*MATRIX(L,POSIJ)*VECTOR(L,J)
					   ENDDO
					ENDDO
				ENDDO


             CASE(6)  

				DO L=1,NUMEL                 !Loop over elements
					NNUD=LNNDEL(L)
					DO I=1,NNUD-1              !loop over nodes of elements
					   NODEI=KXX(I,L)
					   DO J=I+1,NNUD              !loop over nodes of elements
						   NODEJ= KXX(J,L)
						   POSIJ= ((I-1)*I)/2+J
						   TARGETVEC(NODEI)=TARGETVEC(NODEI)
     ;         +FACTOR*MATRIX(L,POSIJ)*(VECTOR(L,I)-VECTOR(L,J))
						   TARGETVEC(NODEJ)=TARGETVEC(NODEJ)
     ;         +FACTOR*MATRIX(L,POSIJ)*(VECTOR(L,J)-VECTOR(L,I))
					   ENDDO
					ENDDO
				ENDDO

C__________________________________Part 8. symetric matrix without diagonal
		END SELECT	
	 
	ENDIF
      RETURN                           
      END
