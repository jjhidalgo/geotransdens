      REAL*8 FUNCTION WEIGHTRESID
     ;(COLUMN, IDIMCOV,NBANDCOV,NUMTOBS,COVINV,RESID)

*******************************************************************************
*
* PURPOSE
*
*   This function calc. the product of the row vector of residuals with the 
*   i-th column of the covariance matrix, to obtain the i-th weighted residual.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix
*  RESID                  Vector of unweighted residuals    
*
* EXTERNAL VARIABLES: SCALARS
*
*  COLUMNS                Column number of covariance matrix
*  IDIMCOV                Used to dimension array COVINV
*  NUMTOBS                Total number of observations      

* INTERNAL VARIABLES: SCALARS
*
*  I                      Dummy counter
*  POS                    Index of element in COVINV
*  POSITION               external function. Index of element in COVINV
*  WEIGHT                 Element of COVINV
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  POSITION               Gets the index of element in COVINV
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan, 2003): Revision and formatting
*
*******************************************************************************

      IMPLICIT NONE
                                                  ! External variables: scalars
      INTEGER*4 IDIMCOV,NUMTOBS,COLUMN,nbandcov
                                                   ! External variables: arrays
      REAL*8 COVINV(IDIMCOV),RESID(NUMTOBS) 
                                                  ! Internal variables: scalars
      INTEGER*4 I
      REAL*8 WEIGHT, getvalue2
                                                   ! Internal variables: arrays
      
      WEIGHTRESID=0.D0
      DO I=1,NUMTOBS                                   ! Loop over observations
         WEIGHT = DSQRT(GETVALUE2
     ;(COLUMN  ,IDIMCOV  ,NBANDCOV  ,I  ,NUMTOBS  ,COVINV))
         WEIGHTRESID=WEIGHTRESID+RESID(I)* WEIGHT      
      ENDDO

      RETURN
      END 
