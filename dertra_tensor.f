      SUBROUTINE DERTRA_TENSOR(DERTRA,ISOZ,IANISTRP,LDIM)

***********************************************************************
* PURPOSE
*
*   Builds the structure of the derivative of transmissivity tensor
*   w. r. t. parameter in a given anisotropy direction 
*
*
* 
* EXTERNAL VARIABLES: ARRAYS
*
*
*  DERTRA structure of the derivative of transmissivity tensor
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  ISOZ       Degree of anistropy of the element
*  IANISTRP   Anistropy direction with respect transmissivity tensor
*             is derivated
*  LDIM       Dimension of the element
*
*
***********************************************************************

      IMPLICIT NONE

      REAL*8::DERTRA(6)
      INTEGER*4::ISOZ,IANISTRP,LDIM


C------------------------- First executable statement

      DERTRA(:) = 0D0

      SELECT CASE (LDIM)

C------------------------- 1d element

	    CASE(1)

              DERTRA(1) = 1D0


C------------------------- 2d Element

	    CASE(2)

              SELECT CASE (ISOZ)

	            CASE(1)

                      DERTRA(1:2) = 1D0

	            CASE(2:3)

                        DERTRA(IANISTRP) = 1D0
          
	        END SELECT !ISOZ

C------------------------- 3d element
	    CASE(3)
      
	        SELECT CASE (ISOZ)

	            CASE(1)

	                DERTRA(1:3) = 1D0

	            CASE(2)

	                IF(IANISTRP.EQ.1) THEN
	                    DERTRA(1) = 1D0
	                ELSE IF (IANISTRP.EQ.2) THEN
                          DERTRA(2) = 1D0
	                    DERTRA(3) = 1D0
	                END IF !IANISTRP.EQ.1

	            CASE(3:6)

                        DERTRA(IANISTRP) = 1D0
	
	        END SELECT !ISOZ

	END SELECT !LDIM

      END SUBROUTINE DERTRA_TENSOR