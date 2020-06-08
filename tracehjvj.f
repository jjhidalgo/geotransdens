      REAL*8 FUNCTION TRACEHJVJ
     ;(IDIMCOV  ,IDIMHESS ,NBANDCOV ,NDEVS    ,NPAR     ,NSTAT
     ;,NUMTOBS  ,OBSTYPE  ,COVINV   ,HESSINV  ,IODEVICE ,OBSCLASS
     ;,VJAC     ,OUT)

*********************************************************
*
* PURPOSE
*
* This function calculates the trace of the matrix which is defined by the
* product of

*         -1    t     -1
*        H  *  J  *  V   *  J
*
* where H is the Hessian, J is the jacobian matrix, v is the covariance matrix.
* Since J(t) * V^-1 * J is the partial objective function for this type of
* observation, the product comes down upon the inverse of the Hessian 
* multiplied by a scalar. 
*The trace of a matrix is the sum of the diagonal elements
*
* DESCRIPTION
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix                      
*  HESSINV                Inverse of hessian matrix
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Status for calc. of obs.                    
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*  OBSCLASS               See description in STAT_OUTP
*  OUT                    Matrix that is the result of multiplying J*V*J
*  VJAC                   Jacobian matrix                                       
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCOV                Used to dimension array COVINV
*  IDIMHESS               Used to dimension array HESS. It is equal to          
*                         NPAR*(NPAR+1)/2                                       
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NDEVS                  Number of devices
*  NPAR                   Total number of parameters to be estimated            
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)               
*  NUMTOBS                Total number of observations                          
*  OBSTYPE                See descrition at STAT_OUTP
*
* INTERNAL VARIABLES: SCALARS
*
*  COVINVLK               Element (l,k) at COVINV
*  HESSINVIK              Element (l,k) at HESSINV
*  I                      Dummy counter
*  J                      Dummy counter                                        
*  JVJKI                  Element (k,i) at Jt*V*J matrix
*  K                      Dummy counter
*  L                      Dummy counter                                         
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  GETVALUE               Gets matrices value at right pointers
*  JVJ                    Calcs. the product JtVJ
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan. 2003): Revision and formatting
*
********************************************************************************

      IMPLICIT NONE
 
                                                   ! External variables: scalars

      INTEGER*4 NUMTOBS,NPAR,OBSTYPE,NDEVS,NSTAT,IDIMHESS,NBANDCOV
     ;          ,IDIMCOV

                                                   ! External variables: arrays

       INTEGER*4 OBSCLASS(NSTAT,NDEVS+1),IODEVICE(NDEVS+1,9)
      REAL*8 VJAC(NUMTOBS,NPAR), HESSINV(IDIMHESS)
     ;      , COVINV(IDIMCOV)

                                                   ! Internal variables: scalars

      INTEGER*4 I,K, POS1, POS2
      REAL*8 HESSINVIK,JVJKI

                                                   ! Internal variables: arrays

      REAL*8 OUT(IDIMHESS)
      CALL ZERO_ARRAY(OUT,IDIMHESS)
      TRACEHJVJ=0D0
C_______________________________________________ Step 1: multiply Jtransp * V*J

      CALL JVJ
     ;(IDIMCOV  ,IDIMHESS ,NBANDCOV ,NDEVS    ,NPAR     ,NSTAT
     ;,NUMTOBS  ,OBSTYPE  ,COVINV   ,IODEVICE ,OBSCLASS ,OUT
     ;,VJAC)

C_______________ Step 2:multiply jvj with hessinv. Calculate diagonal elements

      DO I=1,NPAR
        DO K=1,NPAR

             !calculate position i,k
          IF (I .GT. K) THEN
              POS1=I*(I+1)/2 -(I-K)
              POS2=K*(K+1)/2 -(K-I)
          ELSE
              POS1 = K*(K+1)/2 -(K-I)
              POS2=I*(I+1)/2 -(I-K)
          ENDIF

          HESSINVIK=HESSINV(POS1)
          JVJKI=OUT(POS2)
          TRACEHJVJ=TRACEHJVJ+HESSINVIK*JVJKI
        ENDDO
      ENDDO

      RETURN
      END
