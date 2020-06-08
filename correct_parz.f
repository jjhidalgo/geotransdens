      SUBROUTINE CORRECT_PARZ
     ;(MXGRPZN  ,NTYPAR    ,NUMITER   ,NZPAR   ,NZTRA  ,INORPAR
     ;,IOPT_GS  ,ISOZ      ,IVPAR     ,LDIM_GS ,PARZ)

********************************************************************************
*
* PURPOSE
*
* This routine corrects the value of PARZ, due to a log-transformation of the
* GEOSTATISTICALLY estimated parameters that have varied (1st iteration=ALL 
* COMPONENTS or variable pilot points=ONLY THOSE GROUPS WHERE PILOT POINTS MAY 
* VARY). 
*
* DESCRIPTION
*
* EXTERNAL VARIABLES: ARRAYS
*
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IOPT_GS                General options for inverse problem. 
*  ISOZ                   Anisotropy of every transmissivity zone               
*  IVPAR                  Array containing estimation indexes
*  LDIM_GS                Array containing dimension of transmisivity zones
*  PARZ                   Vector containing calculated values for all parameters                                            
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  MXGRPZN                Maximum number of groups of zones
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMITER                Inverse problem iteration number
*  NZPAR                  Total number of zones for all nodal and element       
*  NZTRA                  Number of transmissivity zones
*
* INTERNAL VARIABLES: SCALARS
*
*  IDIM                   Dummy counter of dimensions
*  IGROUP                 Group to which zone belongs to
*  ITYPE                  Type of parameter of zone IZPAR
*  IZPAR                  Dummy counter of zones
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY: AAR. First coding (Aug-2003)
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NZPAR,MXGRPZN,NUMITER,NZTRA,NTYPAR
     ;         ,IOPT_GS(MXGRPZN,20),IVPAR(NZPAR,4),ISOZ(NZTRA)
     ;         ,LDIM_GS(NZTRA),INORPAR(NTYPAR)
                                                                 ! Real external
      REAL*8 PARZ(NZPAR)
                                                              ! Integer internal
      INTEGER*4 IZPAR,IGROUP,ITYPE,IDIM

C_______________________ Step 1: Corrects zonal parameters that are estimated 
C_______________________         geostatistically and logarithmically and may 
C_______________________         have a variation due to variable pilot points

      DO IZPAR=1,NZPAR

        IF (IVPAR(IZPAR,1).NE.0) THEN

           IGROUP=IVPAR(IZPAR,3)
           ITYPE=IOPT_GS(IGROUP,1)    ! Type of parameter
           IF (IOPT_GS(IGROUP,2).EQ.1 .AND. IVPAR(IZPAR,4).NE.0) THEN

              IF (NUMITER.EQ.1 .OR. IOPT_GS(IGROUP,6).NE.0) THEN
                  PARZ(IZPAR)=10.0D0**PARZ(IZPAR)
                  IF (ITYPE.EQ.1 .AND. ISOZ(IZPAR).EQ.1) THEN
                     DO IDIM=1,LDIM_GS(IZPAR)
                       PARZ(INORPAR(IDIM)+IZPAR)=PARZ(IZPAR)
                     END DO ! IDIM=1,LDIM(NZ)
                  END IF
              END IF

           END IF ! IOPT_GS(IGROUP,2).EQ.1 .AND. IVPAR(IZPAR,4).NE.0

        END IF ! IVPAR(IZPAR,1).NE.0

      END DO
        
      RETURN
      END

