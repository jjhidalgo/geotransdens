      SUBROUTINE ESTIMCOV
     ;(IDIMHESS   ,IOCORMAT   ,IODETHESSOK  ,IOEIGENVEC  ,ITYPEPAR
     ;,MAINF      ,NPAR       ,NSTAT        ,NTYPAR      ,COVAR
     ;,COREL      ,EIGENVAL   ,EIGENVEC     ,RESIDPAR    ,STAT
     ;,TYPENAME)

********************************************************************************
*
* PURPOSE
*   
*   This subroutine calculate the lower and upper confidence bounds of the 
* parameters and optionally the eigenvalues and eigenvectors of the parameter, 
* the covariance correlation matrix
*
* THEORY
*
*   Correlation of paramterers i and j is given by 
*
*                      covariance of parameters i and j
*   corr(i,j) = -----------------------------------------------------------
*                standard dev of par i* standard dev of par j
*
* while the upper and lower boundaries of a 95% confidence interval of a 
* parameter i are given by:
*
* lower boundary: estimate - 2 * standard deviation
* upper boundary: estimate + 2 * standard deviation
*
* The eigenvalues and eigenvectors are calculated in the subroutine 
* eigrs of file AUXXX.f
*
* PROCEDURE 
*
* - Step 1: Write a header for the outputfile
* - Step 2: Calculate the 95% confidence interval boundaries
* - Step 3: Write confidence intervals to file   
* - Step 4: Write covariance matrix to file
* - Step 5: Calculate the correlation matrix if asked for
* - Step 6: Write the correlation matrix to file
* - Step 8: Call the subroutine to calculate eigenvalues and -vectors 
*         if necessary. In this case, calculate and write them
*
* EXTERNAL VARIABLES  SCALARS
*
*  IDIMHESS               Used to dimension array HESS. It is equal to          
*                         NPAR*(NPAR+1)/2     
*  IOCORMAT               Option for writing correl. matrix of estim. param.
*  IODETHESSOK            Indicates if problems were encountered with the determinant of
*                         the hessian.
*                         If 0: the determinant of the hessian is acceptable.
*                         If 1: the hessian is singular
*  IOEIGENVEC             Option for writing eigenvectors and eigenvalues
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NPAR                   Total number of parameters to be estimated            
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NZPAR                  Total number of zones for all nodal and element       

*
* EXTERNAL VARIABLES  ARRAYS
*
*  COVAR                  called hessinv by statout.
*  COREL                  the same as hessinv. this is because EIGRS destroys
*                         its input. Called "work" by statout
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  PARC                   Vector containing calculated values for all           
*                         parameters      
*  STPAR                  Vector containing standard deviation errors of        
*                         all parameters prioo information   
*  TYPENAME               Array containing  the names of the state var. and 
*                         parameter types in the same order as OBSCLASS and STAT
*
* INTERNAL VARIABLES  SCALARS
*
*  I, J                   dummy counters
*  IER                    error parameter of egrs
*  IZ                     Input row dimension of matrix z                                       
*  JOB                    variable telling EIGRS what to do (see eigrs)
*  NPAR1                  equal to npar(to prevent npar from being 
*                         destroyed by EIGRS)
*  PLACE                  the array index in COVAR of the diagonal element 
*                         of the inverse hessian matrix which contains the 
*                         variance of a given parameter estimate
*  PLACEII                index of element (i,i) in the covariance matrix
*
*  PLACEIJ                index of element (i,j) in the covariance matrix
*  PLACEJJ                index of element (j,j) in the covariance matrix
*  STDDEV                 standard deviation of parameter estimate
*  TYPECOUNT              used to see to what parameter type a parameter 
*                         belongs
*  WK                     see eigrs
*  ZONE                   zone number

*
* INTERNAL VARIABLES  ARRAYS
*
*  EIGENVAL              contains the eigenvalues of the covariance matrix
*  EIGENVEC              contains the eigenvectors of the covariance matrix
*
* SUBROUTINES REFERENCED: 
*
*  EIGRS                  calculates eigenvectosr of a symmetric matrix
*  OUTPUT2                this subroutine writes a matrix in symmetric band
*                         storage mode to file. It is written for the 
*                         correlation and covariance matrices. 
*  OUTPUT3                this subroutine writes a matrix in full storage mode
*                         to file. It is written for the eigenvector components
*                         matrix.
* 
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan. 2003): Revision and formatting
*
*************************************************

C______________________________________________ Step 0: Declaration of variables
      IMPLICIT NONE 
                                                   ! External variables: scalars
      INTEGER*4 IOEIGENVEC, IOCORMAT, IDIMHESS, NPAR
     ;         ,NTYPAR,MAINF,NSTAT,IODETHESSOK

                                                   ! External variables: arrays
      REAL*8 COVAR(IDIMHESS),COREL(IDIMHESS)
     ;      ,EIGENVEC(NPAR, NPAR), EIGENVAL(NPAR),STAT(40,11)
     ;      ,RESIDPAR(NPAR,2),ITYPEPAR(NPAR,2)
      CHARACTER*10 TYPENAME(NSTAT+NTYPAR)
                                                  ! Internal variables: scalars
      INTEGER*4 JOB, IZ,  IER, PLACEII, PLACEJJ
     ;         ,I, J,NPAR1, PLACEIJ
      REAL*8 TAUW1

    
C___________________________ Step 4: If so desired, writes the covariance matrix 

      IF (IODETHESSOK.EQ.0) THEN
         IF (IOCORMAT .GE. 1) THEN 
            WRITE(MAINF,3000)
 3000       FORMAT(//,36X,'COVARIANCE MATRIX',/,36X,'---------- ------'
     ;            ,//)

                                         !scaling to get the inverse fisher info 
            I=1
            DO WHILE  (STAT(I,1).EQ. 0D0)
              I=I+1
            ENDDO
            TAUW1=STAT(I,1)/STAT(I,11)

            DO I=1,IDIMHESS
                COVAR(I)=COVAR(I)*TAUW1
            ENDDO
            CALL OUTPUT2
     ;(NPAR,IDIMHESS,MAINF,NSTAT,ITYPEPAR,COVAR,TYPENAME)
C________________________ Step 5: Calculate the correlation matrix if so desired

            DO I=1,NPAR
                DO J=1,I
                   PLACEIJ=(I*(I+1))/2-(I-J)
                   PLACEII=(I*(I+1))/2
                   PLACEJJ=(J*(J+1))/2
                   COREL(PLACEIJ)=COVAR(PLACEIJ)/(DSQRT(COVAR(PLACEII))
     ;              *DSQRT(COVAR(PLACEJJ)))
                 ENDDO 
             ENDDO                       
     
      
C___________________________ Step 6: Write the correlation matrix, if so desired

           WRITE(MAINF,4000)
 4000      FORMAT(//,36X,'CORRELATION MATRIX',/,36X,'----------- ------'
     ;            ,//)

           CALL OUTPUT2
     ;(NPAR,IDIMHESS,MAINF,NSTAT,ITYPEPAR,COREL,TYPENAME)

         ENDIF        !IF (IOCORMAT .GT. 1) THEN


C____________________ Step 7: Check if eigenvalues or eigenvectors are asked for
C____________________         In this case, calculate and write them
 
         IF (IOEIGENVEC .GT. 0) THEN
            JOB=1                       ! Calculate eigenvalues and eigenvectors

            NPAR1=NPAR
            IZ=NPAR
            CALL EIGRS  (COVAR,NPAR1,JOB,EIGENVAL,EIGENVEC,IZ
     ;                  ,RESIDPAR,IER) 

            WRITE(MAINF,5000)
 5000       FORMAT(//,34X,'EIGENVECTORS',/,34X,'------------',/)

            CALL OUTPUT3(EIGENVEC, NSTAT,TYPENAME, NPAR
     ;,EIGENVAL,MAINF,itypepar)

 
         ENDIF                                      !IF (IOEIGENVEC .GT. 0) THEN

      ENDIF                                        !IF (IODETHESSOK .EQ. 0) THEN
     
      RETURN
      END
