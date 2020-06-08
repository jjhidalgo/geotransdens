      SUBROUTINE UPDATE_PARZ
     ;(IDIMWGT    ,NPAR    ,NZPAR    ,DLT_PAR    ,IPNT_PAR    
     ;,IVPAR      ,PARZ    ,WGT_PAR)

********************************************************************************
*
* PURPOSE Updates vector of zonal parameters from iteration K to iteration K+1
*
* DESCRIPTION 
*
* - Step 0: Declaration of variables
* - Step 1: Loop over zones of parameters
*   - Step 1.1: Checks if zonal parameter varies along the optimization process
*     - Step 1.1.A: Loop over parameterization components
*     - Step 1.1.B: Adds updated value, depending on the log. estimation index
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  DLT_PAR                Vector containing increments of unknown parameters
*                         (solution of Marquardt's process)
*  IPNT_PAR               Array containing pointers to IPNT_PAR
*  IVPAR                  Array containing estimation indexes
*                           - Column 1: First useful position at IPNT_PAR
*                           - Column 2: Last useful position at IPNT_PAR
*                           - Column 3: Group of zones
*                           - Column 4: Type of parameter
*  PARZ                   Array containing zonal parameters
*  WGT_PAR                Array containing the weights of the parameterization
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMWGT                Used to dimension IPNT_PAR, WGT_PAR
*  NPAR                   Total number of parameters to be estimated            
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*
* INTERNAL VARIABLES: ARRAYS
*
*
* INTERNAL VARIABLES: SCALARS
*
*  IPOS                   Position at IPNT_PAR
*  IZPAR                  Dummy counter of zonal parameters
*  UPDATE                 Quantity to be added to PARZ
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY
*
*     AAR        9-2003   First coding
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NZPAR,IDIMWGT,NPAR
     ;         ,IVPAR(NZPAR,4),IPNT_PAR(IDIMWGT*NZPAR)
                                                                 ! Real external
      REAL*8 WGT_PAR(IDIMWGT*NZPAR),DLT_PAR(NPAR),PARZ(NZPAR)
                                                              ! Integer internal
      INTEGER*4 IZPAR,IPOS
                                                                 ! Real internal
      REAL*8 UPDATE

C_______________________ Step 1: Loop over zones of parameters

      DO IZPAR=1,NZPAR

C_______________________ Step 1.1: Checks if zonal parameter varies along the
C_______________________           optimization process

         IF (IVPAR(IZPAR,1).NE.0) THEN

            UPDATE=0.0D0       ! Initializes updating value

C_______________________ Step 1.1.A: Loop over parameterization components

            DO IPOS=IVPAR(IZPAR,1),IVPAR(IZPAR,2)
               UPDATE=UPDATE+WGT_PAR(IPOS)*DLT_PAR(IPNT_PAR(IPOS))
            END DO

C_______________________ Step 1.1.B: Adds updated value, depending on the log. 
C_______________________             estimation index

            IF (IVPAR(IZPAR,4).EQ.0) THEN   ! Aritm. estimation
               PARZ(IZPAR)=PARZ(IZPAR)+UPDATE
            ELSE                            ! Log. estimation
               PARZ(IZPAR)=10.0D0**( DLOG10( PARZ(IZPAR) ) + UPDATE )
            END IF

         END IF ! IVPAR(IZPAR,1).NE.0

      END DO ! IZPAR=1,NZPAR

      RETURN
      END
