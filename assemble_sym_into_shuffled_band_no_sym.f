	SUBROUTINE ASSEMBLE_SYM_INTO_SHUFFLED_BAND_NO_SYM
     ;(FACTOR    ,IA_COLS   ,IA_ROWS    ,IDSC_COLS  ,IDSC_ROWS
     ;,LMXNDL    ,NBAND     ,LBLOCK_NR  ,A          ,A_DSC
     ;,LNNDEL    ,KXX)


********************************************************************************
*
* PURPOSE  To assemble a symetric matrix into a non-symmetric banded matrix.
*
*
* DESCRIPTION The subroutine assembles the elements (multiplied by a factor)
*			of a symmetric matrix into symmetric banded one.
*
*
*             The symmetric matrix must be stored in a vector in the following
*             fashion:
*
*                      | 11  12  13 |
*                      | --  22  23 |   =>  |11, 12, 13, 22, 23, 33|
*                      | --  --  33 |
*
*             Then, The element (I,J) is located the position 
*             (I - 1)*NNUD + J - I*(I-1)/2 in the vector.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A						Origin matrix (full)
*  A_DSC					Destiny matrix (non-symmetric banded).
*  LBLOCK_NR				The block of the non-shufled matrix where the element-wise matrix
*						should go. 
*							if 1, the block of derivatives of flow eqn to head
*							if 2, the block of derivatives of flow eqn to mass fraction
*							if 3, the block of derivatives of tpt eqn to head
*							if 4, the block of derivatives of tpt eqn to mass fraction
*
* EXTERNAL VARIABLES: SCALARS
*
*  FACTOR					Factor that multiplies the matrix to be assembled.
*						For time-weighted schemes pourposes (theta-scheme).
*
*  IA_COLS				Number of columns (first dimension) of A matrix.
*  IA_ROWS				Number of rows (second dimension) of A matrix.
*  IDSC_COLS				Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS				Number of rowss (second dimension) of A_DSC matrix.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY: First coding: JHG (7-2003)
*
********************************************************************************



		IMPLICIT NONE
	
		INTEGER:: IA_COLS, IA_ROWS, IDSC_COLS,IDSC_ROWS

		INTEGER::NNUD,I,J,KEXT,AIJ,L,LMXNDL,NODE1
     &			 ,NBAND,IDIAG,LBLOCK_NR,LMDIAG,IJ_BAND,JI_BAND
     &             ,IJ_COL,JI_COL,IJ_ROW,JI_ROW


		REAL*8::  A(IA_ROWS,IA_COLS),
     &			  A_DSC(IDSC_ROWS,IDSC_COLS), 
     &			  FACTOR

		INTEGER*4::LNNDEL(IA_ROWS),KXX(LMXNDL,IA_ROWS)
     &               ,NODEI,NODEJ,KINT


	 LMDIAG = 2*NBAND +2
	 DO L=1,IA_ROWS

		NNUD=LNNDEL(L)

		DO I=1,NNUD-1 
              NODEI= KXX(I,L)  
		
			DO J=I+1,NNUD
				NODEJ=KXX(J,L)

				SELECT CASE (LBLOCK_NR)
					CASE(1)
						IJ_ROW = (NODEI - 1) * 2 + 1
     					    IJ_COL = (NODEJ - 1) * 2 + 1
						JI_ROW=  (NODEJ - 1) * 2 + 1
	                    JI_COL = (NODEI - 1) * 2 + 1

					CASE(2)
						IJ_ROW = (NODEI - 1) * 2 + 1
						IJ_COL = NODEJ * 2 
	                    JI_ROW = (NODEJ - 1) * 2 + 1
                          JI_COL = NODEI * 2 

					CASE(3)
						IJ_ROW =  NODEI * 2 
						IJ_COL=  (NODEJ - 1) * 2 + 1
	                    JI_ROW = (NODEJ * 2)
						JI_COL = (NODEI - 1) * 2 + 1

					CASE(4)
						IJ_ROW = NODEI * 2
						IJ_COL = NODEJ * 2
						JI_ROW = NODEJ * 2
						JI_COL = NODEI * 2
				END SELECT


                  IF (NODEI .LT. NODEJ) THEN									   
				   KINT = IJ_ROW
	               KEXT = JI_ROW
                     IJ_BAND = LMDIAG + IABS(IJ_ROW - IJ_COL)
				   JI_BAND = LMDIAG - IABS(JI_ROW - JI_COL) 
				ELSE
				   KINT = JI_ROW
	               KEXT = IJ_ROW
                     IJ_BAND = LMDIAG + IABS(JI_ROW - JI_COL) 
				   JI_BAND = LMDIAG - IABS(IJ_ROW - IJ_COL) 
                  ENDIF
				 				    
				AIJ = (I - 1)*NNUD + J - I*(I-1)/2

				
			   A_DSC(KINT,IJ_BAND)=A_DSC(KEXT,IJ_BAND)+A(L,AIJ)*FACTOR
                 A_DSC(KEXT,JI_BAND)=A_DSC(KEXT,JI_BAND)+A(L,AIJ)*FACTOR

			END DO ! J=I+1,NNUD

		END DO ! I=1,NNUD-1


C------------------------- Diagonal.		
		DO I=1,NNUD

			SELECT CASE (LBLOCK_NR)
				CASE(1)
					NODE1 =	(KXX(I,L) - 1) * 2 + 1
	                IDIAG = LMDIAG
	            CASE(2) 
					NODE1 = (KXX(I,L) - 1) * 2 + 1
	                IDIAG = LMDIAG + 1
				CASE(3)
					NODE1 =	KXX(I,L) * 2 
	                IDIAG = LMDIAG - 1
				CASE(4)
					NODE1 =	KXX(I,L) * 2 
	                IDIAG = LMDIAG
			END SELECT

			
			AIJ = (I - 1)*NNUD + I
			A_DSC(NODE1,IDIAG) = A_DSC(NODE1,IDIAG)+A(L,AIJ)*FACTOR

		END DO !I=1,NNUD

	END DO ! L=1,IA_ROWS

	END SUBROUTINE ASSEMBLE_SYM_INTO_SHUFFLED_BAND_NO_SYM