      SUBROUTINE INCREMENT_ITER
     &          (IOITERVEND ,ITERTOTV ,ITERV    ,ITERVMX)

***********************************************************************
* PURPOSE
*
*	Generic subroutine to increment iterations counters and control
*     if maximum nuber of iterations has been surpassed.
*
* DESCRIPTION
*
*	Generic subroutine to increment iterations counters.
*
* EXTERNAL VARIABLES: ESCALARS
*
*  IOITERVEND   1 if maximum number of iterations has been surpassed.
*  ITERVTOT     Cumulative iterations
*  ITERV        Iterations in current time-step
*  ITERVMX      Maximum number of iterations.
*
* HYSTORY
*
*  JHG	12/2003	First coding.
***********************************************************************

      IMPLICIT NONE

      INTEGER*4::IOITERVEND ,ITERTOTV,ITERV,ITERVMX


      ITERV = ITERV +1

      ITERTOTV = ITERTOTV +1

	IF (ITERV.GT.ITERVMX) THEN

		IOITERVEND = 1

	ELSE

		IOITERVEND = 0

	END IF


      END SUBROUTINE INCREMENT_ITER
