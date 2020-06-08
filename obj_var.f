      SUBROUTINE OBJ_VAR
     ;(FNEW     ,FOLD     ,IDIMCOV  ,ISUMFO   ,NDEVS
     ;,NSTAT    ,NUMTOBS  ,OBJCON   ,OBJFLO   ,OBJGEO,  OBJHED   ,OBJHUM
     ;,COVINV   ,FOBJ_WGT ,VOBS     ,VOBSC    ,MEASTYP)   

***********************************************************************
* PURPOSE
*
* Computes the observation part of the objective function
*
* DESCRIPTION
*                               *t   -1        *
* Computes:   lambda   *((V   -V)) *V  *(V   -V)
*                   obs    obs  obs       obs  obs
*
*
*                 where: lambda        : weighting factor for each
*                              obs       observation type
*                         -1
*                        V             : inverse covariance matrix
*
*                        V             : observation 
*                         obs
*
* The calculation is organised according to the inverse covariance
* matrix. As this matrix is stored in bands, the calculation sequence
* is initiated with calculation for the diagonal followed by
* calculations for the non-diagonal bands. As the observations are
* assumed to be independent between devices, the calculations are
* carried out within a loop for devices.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix                      
*  FOBJ_WGT               Array containing all objective function weights for   
*                         state variables (heads, concentrations, etc)          
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Status for calc. of obs.                    
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*  VOBS                   Observation value                                     
*  VOBSC                  Value of simulated value corresponding to observation 
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  FNEW                   Objective function value in the current iteration     
*  FOLD                   Objective function computed value in last iteration   
*  IDIMCOV                                                                      
*  ISUMFO                                                                       
*  NDEVS                                                                        
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)               
*  NUMTOBS                Total number of observations                          
*  OBJCON                 Concentration contribution to objective function      
*  OBJFLO                                                                       
*  OBJHED                 Head contribution to objective function               
*  OBJHUM
*  OBJGEO                 Contribution of geophysics      
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR      7-2001     Revision. Small changes
*     JHG      5-2020     Added geohysical observations
*     
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION COVINV(IDIMCOV),VOBS(NUMTOBS),VOBSC(NUMTOBS+NDEVS),
     ;          FOBJ_WGT(NSTAT),MEASTYP(NUMTOBS)

C------------------------- Initialization of obj. functions

       OBJHED=0.D0
       OBJCON=0.D0
       OBJHUM=0.D0
       OBJFLO=0.D0
       OBJGEO=0.D0

       DO I=1,NUMTOBS                        ! Loop over diagonal of COVINV
          FOBJ_I=COVINV(I)*(VOBSC(I)-VOBS(I))*(VOBSC(I)-VOBS(I))*
     ;                         FOBJ_WGT(MEASTYP(I))

          IF (MEASTYP(I).EQ.1) THEN
            OBJHED=OBJHED+FOBJ_I                   ! Contr. of heads
          ELSE IF (MEASTYP(I).EQ.2) THEN
            OBJCON=OBJCON+FOBJ_I                   ! Contr. of conc.
          ELSE IF (MEASTYP(I).EQ.3) THEN
            OBJHUM=OBJHUM+FOBJ_I                   ! Contr. of humid.
          ELSE IF (MEASTYP(I).EQ.4) THEN
            OBJFLO=OBJFLO+FOBJ_I                   ! Contr. of flow.
          ELSE IF (MEASTYP(I).EQ.6) THEN
            OBJGEO=OBJGEO+FOBJ_I                   ! Contr. of geophysics.
          ENDIF

       ENDDO

C------------------------- Computes total objective function

       FNEW=FNEW+OBJHED+OBJCON+OBJHUM+OBJFLO+OBJGEO

      
       IF (FNEW.GT.FOLD) THEN
          ISUMFO=1
       ELSE
          ISUMFO=0
       ENDIF

       RETURN
       END
