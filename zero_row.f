	SUBROUTINE ZERO_ROW
     ;(IDIMMAT1,IDIMMAT2,ISPARSE,MAXNN,NROW,IADN,MATRIX)
      
*************************************************************
*
* PURPOSE
* To put to zero all the elements of a given row. 
*
*************************************************************

	IMPLICIT NONE

C EXTERNAL VARIABLES: SCALARS
      INTEGER*4 MAXNN, IDIMMAT1,IDIMMAT2,ISPARSE,NROW

C EXTERNAL VARIABLES: ARRAYS
      INTEGER*4  IADN(MAXNN)
	REAL*8 MATRIX(IDIMMAT1, IDIMMAT2)

C INTERNAL VARIABLES: SCALARS
      INTEGER*4 I

      IF (ISPARSE .EQ.0) THEN 
	    DO I=1,IDIMMAT2                                    
              MATRIX(NROW,I)=0D0
          END DO
      ELSEIF (ISPARSE .EQ. 1) THEN
		DO I= 1,IADN(NROW)
	       MATRIX(I,NROW)=0D0
	    ENDDO
	ENDIF 

	RETURN
	END