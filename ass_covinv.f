      SUBROUTINE ASS_COVINV
     ;(CORRCOEF,IDIMCOV,NDEVS,STDEV,COVINV,IODEVICE)

***********************************************************************
* PURPOSE
*
* Assembly of the inverse correlation matrix (stored in vector COVINV)
*
* DESCRIPTION
*
* The inverse correlation matrix, is assembled depending on the value
* of IODEVICE(NODEV,6) which for each device contains the information
* on the structure of the inverse correlation matrix for the device.
* The matrix is stored as a vector in COVINV.
* 
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix                      
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Status for calc. of obs.                    
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  CORRCOEF               Correlation coefficient within device
*  IDIMCOV                Used to dimension array COVINV                 
*  NDEVS                  Number of devices                        
*  STDEV                  Default standard deviation for device
*
* INTERNAL VARIABLES: SCALARS
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR     03-2001     First revision
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION COVINV(IDIMCOV),IODEVICE(NDEVS+1,10)

* This subroutine has to be checked. Only the diagonal V is correct!

C______________________________ Loop over devices

      DO ND=1,NDEVS

        NOF=IODEVICE(ND,8)                        ! First obs. of current device
        NOL=IODEVICE(ND+1,8)-1                     ! Last obs. of current device

*_____________Independent observations (diagonal covariance matrix)

        IF (IODEVICE(ND,6).EQ.1) THEN

          DO I=NOF,NOL
            COVINV(I)=1D0/COVINV(I)/COVINV(I)
            A=12.0
          ENDDO

*_____________Autoregressive errors (only if all obs. have same st. dev.)

C I am not 100% certain that the inverse covariance matrix is correct. The matrix
C described in (Carrera and Walthers, 1996) is erroneous or at least ambiguous......

        ELSE IF (IODEVICE(ND,6).EQ.2) THEN

          AUX1=1D0/STDEV/STDEV/(1-CORRCOEF*CORRCOEF)      !(?)
          AUX2=-CORRCOEF*AUX1
          DO I=NOF,NOL
            COVINV(I)=AUX1              !Diagonal
            COVINV(NOL+I)=AUX2          !Band next to diagonal
          ENDDO
          COVINV(NOL)=AUX1              !Lower right corner

*_____________Exponential

        ELSE IF (IODEVICE(ND,6).EQ.3) THEN

*_____________Gaussian

        ELSE IF (IODEVICE(ND,6).EQ.4) THEN

*_____________Linear

        ELSE IF (IODEVICE(ND,6).EQ.5) THEN

C       And so on.........

        ENDIF

      ENDDO

C It should be verified that the previously supplied value of NBANDCOV (the
C bandwidth of the inverse covariance matrix) is correct. The value cannot be
C determined in ASS_COVINV because it is previously needed to dimension COVINV.

C If the correlation between two obs. measured over an time interval depends on
C the temporal difference between the two measurements, this should be
C determined as: dt = 0.5 * ((TOBS(1,2)-TOBS(1,1))-(TOBS(2,2)-TOBS(2,1)).

C If the user wants to enter the covariance matrix (and not its inverse), maybe
C a call to an inversion routine should be added....?

      RETURN
      END
