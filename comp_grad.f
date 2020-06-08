       SUBROUTINE COMP_GRAD
     ; (IDIMCOV  ,NBANDCOV ,NDEVS    ,NPAR     ,NSTAT    ,NUMTOBS
     ; ,COVINV   ,FOBJ_WGT ,GRAD     ,MEASTYP  ,VJAC     ,VOBS
     ; ,VOBSC)


***********************************************************************
* PURPOSE
*
* Computes the gradient of the objective function with respect to the
* parameters to be estimated
*
* DESCRIPTION
*
*                        t    -1         *
* Computes:   lambda   *J   *V   *(V   -V )
*                   obs  obs        obs  obs
*
*
*                 where: lambda        : weighting factor for each
*                              obs       observation type
*
*                        J             : Jacobian (gradient of obs.
*                         obs            with respect to the parameters
*                                        to be estimated)
*                         -1
*                        V             : inverse covariance matrix
*
*                        V             : observation
*                         obs
*
* The subroutine actually computes (a part of) the gradient divided by 2. As
* also the Hessian is reduced by a factor two, the parameter change, determined
* as the fraction gradient/Hessian remains unchanged.
*
* The calculation is organised according to the inverse covariance matrix. As
* this matrix is stored in bands, the calculation sequence is initiated with
* calculation for the diagonal followed by calculations for the non-diagonal
* bands.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix                      
*  GRAD                   Vector containing objective function gradient         
*  VJAC                   Jacobian matrix                                       
*  VOBS                   Observation value                                     
*  VOBSC                  Value of simulated value corresponding to observation 
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCOV                                                                      
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NPAR                   Total number of parameters to be estimated            
*  NUMTOBS                Total number of observations                          
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR      7-2001     Revision and header inclusion
*     AMS      1-2004     Inclusion of objective function weights
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION COVINV(IDIMCOV),VOBS(NUMTOBS),VOBSC(NUMTOBS+NDEVS),
     ;          VJAC(NUMTOBS,NPAR),GRAD(NPAR),MEASTYP(NUMTOBS),
     ;          FOBJ_WGT(NSTAT)

C---------------------- Computes gradient

      DO NPOS=1,NPAR                    ! Loop over positions in gradient vector
         SUM=0.D0

         DO I=1,NUMTOBS                            ! Loop over diagonal of COVINV

            SUM=SUM + 2*COVINV(I)*(VOBSC(I)-VOBS(I))*VJAC(I,NPOS)*
     ;                          FOBJ_WGT(MEASTYP(I))

C---------------------- Non diagonal covariance matrix

            IF (NBANDCOV.GE.2) THEN           ! Temp. corr. between obs. with device
               K=NUMTOBS
               DO NB=2,NBANDCOV                                  ! Non-diagonal bands
                  DO J=1,NUMTOBS-NB+1                 ! Elements in non-diagonal bands
                     K=K+1
                     SUM=SUM+COVINV(K)*(  
     ;                  (VOBSC(J)-VOBS(J))*VJAC(NPOS,NB+J-1)+
     ;                       (VOBSC(NB+J-1)-VOBS(NB+J-1))*VJAC(NPOS,J)
     ;                                 )*FOBJ_WGT(MEASTYP(J))
                  ENDDO
               ENDDO
            ENDIF

         ENDDO   ! I=1,NUMTOBS

         GRAD(NPOS)=SUM

      ENDDO                                                     ! Next component

      RETURN

      END
