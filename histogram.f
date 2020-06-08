      INTEGER*4 FUNCTION HISTOGRAM(DEVICESTDDEV,MEANRESID )
      IMPLICIT NONE
************************************************
* PURPPOSE
* This function calculates the histogramm class.
* The histogram class means: How many times the 
* standard deviation is the mean residual different 
* from the expected mean residual (which is zero)
* In this histogram there are 6 categories:
*
* mean residual smaller than -2* stddev
* mean residual smaller than -1* stddev and larger than -2* stddev
* mean residual smaller than 0 and larger than -1* stddev
* mean residual larger than 0 and smaller than 1* stddev
* mean residual larger than 1 * stddev and smaller than 2* stddev
* mean residual larger than 2* stddev
*
*  EXTERNAL VARIABLES: SCALARS
* meanresid: the mean of the dataset of which we want 
*     to calculate the histogram class 
* devicestddev: the standard deviation of the device that 
*              made the measurements. This standard deviation
*              reflects the a priori estimate of reliablility
*              and should not be confused with the standard deviation
*              of the dataset that the device measured.
* INTERNAL VARIABLES: SCALARS
* histogram: an integer which can take the following values
*histogram = -3:: mean residual < -2* stddev
*histogram = -2:: mean residual < -1* stddev and >= -2* stddev
*histogram = -1:: mean residual < 0 and >= -1* stddev
*histogram = +1:: mean residual >= 0 and < 1* stddev
*histogram = +2:: mean residual >= 1 * stddev and < 2* stddev
*histogram = +3:: mean residual >= 2* stddev
           
      REAL*8 MEANRESID, DEVICESTDDEV
******************************************************
      IF (MEANRESID .LT. (-2.*DEVICESTDDEV)) THEN 
          HISTOGRAM = -3
      ENDIF
******************************************************
      IF ((MEANRESID .LT. (-1.*DEVICESTDDEV)) .AND.(MEANRESID .GE. (-2.*
     &     DEVICESTDDEV))) THEN 
           HISTOGRAM = -2
      ENDIF
******************************************************
      IF ((MEANRESID .LT. 0) .AND. (MEANRESID .GE. (-1.*DEVICESTDDEV))) 
     &    THEN 
           HISTOGRAM = -1
      ENDIF
******************************************************
      IF ((MEANRESID .GE. 0) .AND. (MEANRESID .LT. DEVICESTDDEV)) 
     &THEN 
          HISTOGRAM = 1
      ENDIF
******************************************************
      IF ((MEANRESID .GE. (1.*DEVICESTDDEV)) .AND. (MEANRESID .LT. (2.*
     &    DEVICESTDDEV))) THEN 
           HISTOGRAM = 2
      ENDIF
******************************************************
      IF (MEANRESID .GE. (2.*DEVICESTDDEV)) THEN 
          HISTOGRAM = 3
      ENDIF

      RETURN 
      END
