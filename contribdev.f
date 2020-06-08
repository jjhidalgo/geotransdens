      REAL*8 FUNCTION CONTRIBDEV 
     ;(IDIMCOV  ,NBANDCOV ,NOF      ,NOFOBS   ,NUMTOBS
     ;,COVINV   ,RESID)

********************************************************************************
*
* PURPOSE
*
* This function calculates the contribution to the penalty criterion of one 
* device. 
*
* THEORY 
*
* The contribution to the penalty criterion of one observation is given
* by the product      
*                     t   -1
*                   J  * V   * J
*
* where the only nonzero elements in J are the residuals corresponding to.
* the device whose contribution is being calculated. 
* So, J is NOT the jacobian.
* V^(-1) is the inverse of the covariance matrix corresponding to the 
* measurement type that the device has been measuring.
*
* DESCRIPTION
*
* Use has been made of the fact that any element ov [V^-1](j,k) is nonzero only
* if the absolute value of (j-k) is smaller than (nbandcov-1) where nbandcov is 
* the number of nonequal bands of the inverse covariance matrix.
*
*  Step 1: Loop over observations in device
*     Step 2: Loop over bands of covariance matrix
*        Step 3:Adds the contribution of current device
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix                      
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCOV                Used to dimension COVINV
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NOF                    First observation defining a particular device.
*  NOFOBS                 Number of observations defining a given dev.
*  NOL                    Last observation defining a particular device.
*  NUMTOBS                Total number of observations                          
*
* INTERNAL VARIABLES: SCALARS
*
*  BAND                   Band counter
*  COL                    Column counter
*  J                      Dummy counter                      
*  K                      Dummy counter                                         
*  POS                    Pointer to COVINV
*  ROW                    Row counter                                        
*  VJK                    Value of the covariance matrix at row j and column k
*  RESID_J                J-esimo residuo del device 
*  RESID_K                K-esimo residuo del device 
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  POSITION               Gets the right pointer to COVINV  
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan. 2003): Revision and formatting
*
********************************************************************************

C_______________________________________________Step 0: Declaration of variables

      IMPLICIT NONE
      INTEGER*4 J, K , ROW, COL, NOFOBS, NBANDCOV,POSITION,NUMTOBS
     ;          ,IDIMCOV,NOF,BAND,POS
      REAL*8 VJK,COVINV(IDIMCOV),RESID(NUMTOBS)
      REAL*8 RESID_J, RESID_K



      CONTRIBDEV = 0.D0                   ! Init. contribution of current device

C_________________________________ Step 1: Loop over observations in this device

      DO J= 1, NOFOBS     

C__________________________________Step 2: Loop over bands of covariance matrix

        DO K = J-(NBANDCOV-1),J+(NBANDCOV-1)           ! Loop over non-void pos.
       
                     ! To prevent that non existing array elements are asked for 

          IF (K .GT. 0 .AND. K .LE. NOFOBS) THEN 
             ROW = NOF+J-1
             COL = NOF+K-1  

             IF (COL .GT. ROW) THEN 
               ROW = K
               COL = J
             ENDIF

                                                       ! Find out the right band

             BAND=ABS(ROW - COL) + 1
             IF (BAND .LE. NBANDCOV) THEN
               POS = POSITION(COL,NUMTOBS,ROW)
               VJK = COVINV(POS)
             ELSE
               VJK=0.D0
             ENDIF

C_________________________________Step 3:Adds the contribution of current device

             RESID_J= RESID(NOF+J-1) 
             RESID_K= RESID(NOF+K-1)
             CONTRIBDEV = CONTRIBDEV +  RESID_J * VJK * RESID_K
              
          ENDIF                                   ! K .GT. 0 .AND. K .LE. NOFOBS
        ENDDO                                ! K = J-(NBANDCOV-1),J+(NBANDCOV-1)
      ENDDO                                                       ! J= 1, NOFOBS

      RETURN
      END
