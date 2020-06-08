       SUBROUTINE CHECK_DIVERGENCE(DELVMAX    ,DELVMAXOLD ,IONEWT
     &                            ,IREDTIMV   ,MAXNUMDIV  ,NUMDIV
     &                            ,RESVMAX    ,RESVMAXOLD)

*********************************************************
*
* PURPOSE   Checks the divergence criterion
*
* DESCRIPTION This subroutine checks if a reduction time step must
*             be carried out and if code must stop for sveral convergence 
*             problems
*
* EXTERNAL VARIABLES: SCALARS
*
*  DELVMAX                Maximum change in the state variable                                                              
*  DELVMAXOLD             Maximum change in the state variable in last
*                         iteration                                                      
*  EQUATION               String containing the name of the equation 
*                         to be solved (FLOW or TRANSPORT)                                                      
*  IREDTIMV               Indicator variable (=1 --> reduction time must be
*                         carried out
*  MAXNUMDIV              Maximum number of consecutive divergences.
*                         This criterion allows avoiding divergence when
*						there are oscilations in the state variable
*                         increment or in the residue.
*  NUMDIV                 Number of consecutrive divergences.
*  RESVMAX                Maximum residuo otained during the system solution                                                      
*  RESVMAXOLD             Maximum residuo otained during the system solution
*                         in previous iteration                                                                                                            
*
* HISTORY: First coding : German Galarza
*          Revision:      Andres Alcolea
*          Modified:      JHG 
*                        (Time increment reduction moved to UPDATE_TIME_INCR
*                         Writing of information moved to WRITE_NO_LIN).
********************************************************************************

      IMPLICIT NONE
      
      INTEGER*4::IREDTIMV,IONEWT,MAXNUMDIV,NUMDIV

      REAL*8::DELVMAX, DELVMAXOLD, RESVMAX, RESVMAXOLD


C_____________________________initialize
      IREDTIMV = 0

c________________________________Newton Rhapson: check residuals and increments
	IF(IONEWT.EQ.1) THEN 
	   IF (DELVMAX.GT.DELVMAXOLD .OR. RESVMAX.GT.RESVMAXOLD) THEN

		 NUMDIV = NUMDIV + 1

	     IF (NUMDIV.GT.MAXNUMDIV) THEN
			IREDTIMV = 1
           END IF

	   END IF

	ELSE
c________________________________Picard: check increments only
	   IF (DELVMAX.GT.DELVMAXOLD) THEN

	      NUMDIV = NUMDIV + 1

	      IF (NUMDIV.GT.MAXNUMDIV) THEN

	        IREDTIMV = 1
            END IF
	  ENDIF

	ENDIF !(IONEWT.EQ.1)

      RETURN

      END SUBROUTINE CHECK_DIVERGENCE
