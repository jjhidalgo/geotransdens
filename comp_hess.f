      SUBROUTINE COMP_HESS
     ;(IDIMCOV  ,NBANDCOV ,NPAR     ,NSTAT    ,NUMTOBS
     ;,COVINV   ,FOBJ_WGT ,HESS     ,MEASTYP  ,VJAC)

********************************************************************************
*
* PURPOSE
*
* Computes an approximation of the Hessian (the second order derivative
* of the observations with respect to the parameters to be estimated)
*
* DESCRIPTION
*
*                        t    -1
* Computes:   lambda   *J   *V  *J
*                   obs  obs      obs
*
*
*                 where: lambda        : weighting factor for each
*                              obs       observation type
*
*                        J             : Jacobian (gradient of obs.
*                         obs            with respect to the parameters
*                                        to be estimated)
*
* The subroutine actually computes the Hessian divided by 2. As also
* the gradient is reduced by a factor two, the parameter change,
* determined as the fraction gradient/Hessian remains unchanged.
*
* For each value in the Hessian matrix, the calculation is organised
* according to the inverse covariance matrix. As this matrix is stored
* in bands, the calculation sequence is, for each element in the
* Hessian matrix, initiated with calculation for the diagonal
* followed by calculations for the non-diagonal bands.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix                      
*  FOBJ_WGT               Array containing all objective function weights for   
*                         state variables (heads, concentrations, etc)          
*  HESS                   Hessian matrix of objective function.                 
*  VJAC                   Jacobian matrix                                       
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCOV                Used to dimension array COVINV
*  IDIMHESS               Used to dimension array HESS. It is equal to          
*                         NPAR*(NPAR+1)/2                                       
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NPAR                   Total number of parameters to be estimated            
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)               
*  NUMTOBS                Total number of observations                          
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR      7-2001     Revision and header inclusion
*
********************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION COVINV(IDIMCOV),VJAC(NUMTOBS,NPAR),HESS(NPAR*(NPAR+1)/2)
     ;         ,FOBJ_WGT(NSTAT),MEASTYP(NUMTOBS)

      NPOS=0
      DO NP1=1,NPAR                                 ! Row in Hessian matrix
         DO NP2=1,NP1                               ! Column in Hessian matrix
            SUM=0.D0

            NPOS=NPOS+1                             ! Position in Hessian vector

            DO I=1,NUMTOBS                          ! Loop over diagonal of COVINV

               SUM=SUM+VJAC(I,NP1)*COVINV(I)*VJAC(I,NP2)*
     ;                     FOBJ_WGT(MEASTYP(I))

               IF (NBANDCOV.GE.2) THEN              ! If COVINV is not diagonal
                  K=NUMTOBS
                  DO NB=2,NBANDCOV                  ! Loop for bands of COVINV
                     DO J=1,NUMTOBS-NB+1            ! Loop for elements in non-diagonal bands
                        K=K+1
                        SUM=SUM+COVINV(K)*(
     ;                          VJAC(J,NP1)*VJAC(NB+J-1,NP2)+
     ;                                VJAC(NB+J-1,NP1)*VJAC(J,NP2)  )*
     ;                                            FOBJ_WGT(MEASTYP(J))
                     ENDDO
                  ENDDO
               ENDIF

            ENDDO      ! I=1,NUMTOBS 

            HESS(NPOS)=SUM

         ENDDO  !  NP2=1,NP1
      ENDDO     !  NP1=1,NPAR

      RETURN
      END
