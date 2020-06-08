      SUBROUTINE LIKELIHOOD
     ;(DETHESS  ,EXPS1    ,EXPS2    ,EXPS3    ,IDIMCOV  ,IDIMHESS
     ;,LIKS1    ,LIKS2    ,LIKS3    ,NDEVS    
     ;,NPAR     ,NSTAT    ,NTYPAR   ,NUMTOBS  ,OBJF
     ;,COVINV   ,COVPAR   ,EIGENV   ,INORPAR2
     ;,OBSCLASS ,STAT     ,ALPHANEG)


********************************************************************************
*
* PURPOSE
*
* This subroutine calculates the likelihood function S and  its expected value.
* Both are calculated using three sets of statistical parameters. These 
* statistical parameters are the two set of alphas calculated in the subroutine 
* ALPHA, plus a special case where all alphas are set to 1.
*
* REFERENCE
*
* The equations implemented in this subroutine are taken from 
*     A.Medina, J.Carrera, 2001:'Geostatistical inversion of coupled
*     problems: dealing with computational problems and different types
*     of data', Barcelona, Spain 
*
* THEORY
*
* As described in the alpha subroutine there are two ways to calculate 
* sets of the statistical parameters alpha. Which set of alphas is used 
* controls the shape of the expected likelihood equation. 
* Likelihood function (and its expected value) calculated with the set of 
* alphas according to form. 1 is called here S1, while the one calc. according 
* to form. 2 is S2. The additional one (alphas set to 1) is called S3.
* 
* DESCRIPTION
*
* - Step 1: calculate part of expected likelihood function with set of alpha 1
* - Step 2: calculate part of the expected likelihood function using alpha 2
* - Step 3: add terms that need number of observations or of zones
* - Step 4: calculate the expected likelihood function using alpha = 1
* - Step 5: calculate the likelihood function using alpha1, alpha2 and
*          alpha3; neglecting all constant terms
*
*  EXTERNAL VARIABLES: SCALARS
*
*  DETHESS                Determinant of the hessian matrix
*  IDIMCOV                Used to dimension array COVINV
*  IDIMHESS               Used to dimension array HESS. It is equal to          
*                         NPAR*(NPAR+1)/2         
*  NDEVS                  Number of devices      
*  NPAR                   Total number of parameters to be estimated            
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning) 
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMTOBS                Total number of observations 
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                 
*   
*  EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix of measurements   
*  COVPAR                 Inverse of the covariance matrix of parameters  
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR  
*  IVPAR                  Vector containing estimation index for all            
*                         parameters 
*  OBSCLASS               Array containing number of devices where a particular
*                         type of state var. was measured
*                         COLUMN 1 contains number of devices, while COLUMNS 2-?
*                         contain identifiers to those devices
*                         ROW is related to state var. type
*  STAT                   Contains statistical parameters and properties
*                         for all parameter types and for all observation types. 
*
*  INTERNAL VARIABLES: SCALARS
*
*  EXPS1                  Expected likelihood acc. to form 1
*  EXPS2                  Expected likelihood acc. to form 2
*  EXPS3                  Expected likelihood acc. to form 3
*  I,J                    dummys ofor do loops                
*  LIKS1                  Log-likelihood function acc. to form. 1
*                         Neglecting constant terms             
*  LIKS2                  Log-likelihood function acc. to form. 2             
*                         Neglecting constant terms
*  LIKS3                  Log-likelihood function acc. to form. 3
*                         Neglecting constant terms
*  NOF                    index of starting location of data of device in vobs 
*  NOFOBSTYPE             number of observations of measurement type
*  NOFZON                 number of zones of parameter type
*  NOL                    index of last entry of data of device in vobs
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ZERO_ARRAY             Fills an array with zeros
*
*  HISTORY: LJS (Dec. 2002): First coding
*           AAR (Jan. 2003): Revision and formatting
*
********************************************************************************

C______________________________________________ Step 0: Declaration of variables

      IMPLICIT NONE

                                                   ! EXTERNAL VARIABLES: SCALARS
      INTEGER*4 NDEVS, NUMTOBS, NPAR, NTYPAR,NSTAT,IDIMHESS
     ;         ,IDIMCOV
      REAL*8 DETHESS, EXPS1, EXPS2, EXPS3, LIKS1, LIKS2, LIKS3, OBJF
                                                    ! EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 OBSCLASS(NSTAT, NDEVS+1), INORPAR2(NTYPAR+1)
       
      REAL*8 STAT(40,11),COVINV(IDIMCOV)
     ;      ,COVPAR(IDIMHESS), EIGENV(NPAR*NPAR)
                                                   ! INTERNAL VARIABLES: SCALARS
      INTEGER*4 I, J,NOFZON, LENGTH, POS1
     ;          ,POS2,ROW, KOL, STARTPOS, NOFOBS
      REAL*8 DET,PI

      LOGICAL ALPHANEG


C_______________________________________Initialize the different functions
      EXPS1=0.D0
      EXPS2=0.D0
      EXPS3=0.D0
      LIKS1=0.D0
      LIKS2=0.D0
      LIKS3=0.D0
      PI = 2D0*DASIN(1D0)

C___________________________calculate the determinants of the covariance matrices
C____________________________STATE VARIABLES 
                                            
      DO I = 1, NSTAT                                 !loop over state variables

         STARTPOS=0
         DET = 0
         NOFOBS=INT(STAT(I,11)       )!nr of observations of current st.var.type

         IF (NOFOBS .GT. 0 ) THEN                     !if there are observations
         
             DO J=1,I-1                     !getting array index of observations
                STARTPOS=STARTPOS+INT(STAT(J,11))
             ENDDO
      
             DO J=STARTPOS+1,STARTPOS+INT(STAT(I,11))
                DET=DET-DLOG(COVINV(J))   !-1 to get ln of the inverse of covinv
             ENDDO    
             STAT(I,10)=DET
         ENDIF                                           !IF (NZON .GT. 0 ) THEN
      ENDDO      
    


C___________________________calculate the determinants of the covariance matrices
C___________________________PARAMETERS

      DO I = 1, NTYPAR                                     !loop over parameters

         STARTPOS=0
         NOFZON=INT(STAT(I+NSTAT,11))
         IF (NOFZON .GT. 0 ) THEN               !if there are parameter zones of 
            CALL ZERO_ARRAY(EIGENV,NPAR*NPAR)                ! current parameter

            DO J=1,I-1   !getting array index of zones of current parameter type 
                STARTPOS=STARTPOS+INT(STAT(J+NSTAT,11))
            ENDDO

            DO ROW = STARTPOS+1,STARTPOS+NOFZON
                DO KOL= STARTPOS+1,ROW
                 POS1=ROW*(ROW+1)/2  - (ROW-KOL) !POSITION IN COVPAR
                 POS2= (ROW-STARTPOS)*(ROW-STARTPOS+1)/2 - (ROW-KOL)
                 EIGENV(POS2)=COVPAR(POS1)
              ENDDO
            ENDDO

            LENGTH=NOFZON*(NOFZON+1)/2
            CALL DETERMINANT(DET,LENGTH,NOFZON,NOFZON,EIGENV)
            STAT(I+NSTAT,10)=-1D0*DLOG(DET) 
         ENDIF  !IF (NZON .GT. 0 ) THEN
      ENDDO 
      

      
C_______________________________________________calculate the log likelihood function
      CALL COMP_LIKELIHOOD
     ;(LIKS1    ,LIKS2    ,LIKS3    ,NDEVS    
     ;,NSTAT    ,NTYPAR   ,NUMTOBS  ,OBJF
     ;,INORPAR2 ,OBSCLASS ,STAT)


C_____________________________________________Calculate the expected value of log-likelihood
C_____________________________________________using optimal alfas. If some of the alfa are 
C_____________________________________________negative, then use the alternative calculation
C_____________________________________________method, based on the ln of the exponential of 
C_____________________________________________expected likelihood
      IF (ALPHANEG .EQV. .FALSE.) THEN  
         CALL  COMP_EXP_LIKELIHOOD_1
     ;(DETHESS  ,EXPS1    ,NDEVS    ,NSTAT    ,NTYPAR 
     ;,INORPAR2 ,OBSCLASS ,STAT)

      ENDIF


C______________________________________calculate the expected value of log likelihood
C______________________________________using alphas consistent with lambda ( EXPS2 ) and
c______________________________________using all alpha equal to 1 (EXPS3)
      CALL COMP_EXP_LIKELIHOOD_2_3
     ;(DETHESS  ,EXPS2    ,EXPS3   ,NDEVS    ,NPAR     
     ;,NSTAT    ,NTYPAR   ,NUMTOBS ,INORPAR2 ,OBSCLASS
     ;,STAT)


      END 
********************************************************************
********************************************************************
********************************************************************
      SUBROUTINE COMP_LIKELIHOOD
     ;(LIKS1    ,LIKS2    ,LIKS3    ,NDEVS    
     ;,NSTAT    ,NTYPAR   ,NUMTOBS  ,OBJF
     ;,INORPAR2 ,OBSCLASS ,STAT)
     ;
C______________________________________________ Step 0: Declaration of variables

      IMPLICIT NONE

                                                   ! EXTERNAL VARIABLES: SCALARS
      INTEGER*4 NDEVS, NUMTOBS, NTYPAR,NSTAT
      REAL*8 LIKS1, LIKS2, LIKS3, OBJF
                                                    ! EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 OBSCLASS(NSTAT, NDEVS+1), INORPAR2(NTYPAR+1)
       
      REAL*8 STAT(40,11)
                                                   ! INTERNAL VARIABLES: SCALARS
      INTEGER*4 I, NOFZON, NOFOBSTYPE
     ;         
      REAL*8 PI,NTO,NPA

      PI =2 * ASIN(1D0)

C________________ Step 3: add terms that need number of observations or of zones

      IF (STAT(1,2).NE.0) THEN
          
                                                ! Loop for counting observations
         LIKS2=LIKS2+ STAT(1,11) * DLOG(STAT(1,2)) !later this term  is substracted

      ELSE
         WRITE(25,*) 'PROBLEMAS EN EL CALCULO DE LOG(0)'
         RETURN
      END IF
      
      DO I=1,NSTAT                                   ! Loop over state variables
        IF (OBSCLASS(I,1) .NE. 0) THEN                       ! Class is not void
           NOFOBSTYPE=INT(STAT(I,11))     
           LIKS2 = LIKS2 - REAL(NOFOBSTYPE) * DLOG(STAT(I,2))
        ENDIF    ! Void class?
      ENDDO     ! Next state variable

      DO I=1,NTYPAR                                ! Loop over parameters types
         IF (INORPAR2(I) .NE. INORPAR2(I+1) .AND.STAT(I+NSTAT,2) .NE.0)
     ;   THEN
           NOFZON=INT(STAT(I+NSTAT,11))                     ! -n*ln(alpha1/tauw)
           LIKS2 = LIKS2 - REAL(NOFZON) * DLOG(STAT(I+NSTAT,2))
          ENDIF
      ENDDO





C____________ Step 5: calculate the likelihood function using TAU1, TAU2 and
C_____________       TAU3
      
      ! to prevent taking into account terms with zero weight
      NTO=0
      DO I= 1,NSTAT
         IF (STAT(I,2) .NE. 0D0 .AND. STAT(I,11) .GT. 0D0) THEN
            NTO = NTO +STAT(I,11)
         ENDIF
      ENDDO
      
      
      NPA=0
      DO I= 1,NTYPAR
          IF (STAT(I+NSTAT,2) .NE. 0D0 .AND. STAT(I+NSTAT,11) .GT. 0D0)
     ;    THEN
             NPA = NPA +STAT(I+NSTAT,11)
          ENDIF
      ENDDO


      LIKS1=LIKS1+REAL(NUMTOBS+NPA)*(1D0+DLOG(2D0*PI))
      LIKS2=LIKS2+REAL(NUMTOBS+NPA)*(1D0+DLOG(2D0*PI))
      LIKS2=LIKS2+REAL(NUMTOBS+NPA)*DLOG(OBJF/REAL(NUMTOBS+NPA))
      LIKS3=LIKS3+REAL(NUMTOBS+NPA)*DLOG(2D0*PI)

      DO I = 1,NSTAT                             ! Loop over state variables
        IF (OBSCLASS(I,1) .GT.0) THEN                        ! Class is not void
            LIKS1=LIKS1+STAT(I,10)
            LIKS1=LIKS1+STAT(I,11)*DLOG(STAT(I,1)/STAT(I,11))
            LIKS2=LIKS2+STAT(I,10)
            LIKS3=LIKS3+STAT(I,10)
            LIKS3=LIKS3+STAT(I,1)
        ENDIF
      ENDDO



      DO I = 1,NTYPAR                               ! Loop over parameter types
        IF (INORPAR2(I).NE.INORPAR2(I+1) .AND. STAT(I+NSTAT,2).NE.0D0
     ;        .AND. STAT(I+NSTAT,1).NE.0)
     ;  THEN
            LIKS1=LIKS1+STAT(I+NSTAT,10)
            LIKS1=LIKS1+STAT(I+NSTAT,11)*DLOG(STAT(I+NSTAT,1)
     ;                                       /STAT(I+NSTAT,11))
            LIKS2=LIKS2+STAT(I+NSTAT,10)
            LIKS3=LIKS3+STAT(I+NSTAT,10)
            LIKS3=LIKS3+STAT(I+NSTAT,1)
        ENDIF
      ENDDO


      END 
********************************************************************
********************************************************************
********************************************************************

       SUBROUTINE COMP_EXP_LIKELIHOOD_1
     ;(DETHESS  ,EXPS1    ,NDEVS    ,NSTAT    ,NTYPAR 
     ;,INORPAR2 ,OBSCLASS ,STAT)
C______________________________________________ Step 0: Declaration of variables

      IMPLICIT NONE

                                                   ! EXTERNAL VARIABLES: SCALARS
      INTEGER*4 NDEVS,NTYPAR,NSTAT
      REAL*8 EXPS1, DETHESS
                                                    ! EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 OBSCLASS(NSTAT, NDEVS+1), INORPAR2(NTYPAR+1)
       
      REAL*8 STAT(40,11)
                                                   ! INTERNAL VARIABLES: SCALARS
      INTEGER*4 I, NOFZON, NOFOBSTYPE
   


C___________ Step 1: calculate part of expected likelihood function with alpha 1

      DO I = 1,NSTAT                                  !Loop over state variables
        IF (OBSCLASS(I,1) .NE. 0) THEN
           EXPS1 = EXPS1 + (1./STAT(I,3))*STAT(I,1)  ! (1/alpha) * penalty crit.
        ENDIF
      ENDDO

      DO I = 1, NTYPAR                                    ! Loop over parameters
                                                 ! (1/alpha) * penalty criterion
        IF (INORPAR2(I) .NE.INORPAR2(I+1) .AND. STAT(I+NSTAT,2) .NE. 0) 
     ;      EXPS1 = EXPS1 + (1./STAT(I+NSTAT,3))*STAT(I+NSTAT,1) 
      ENDDO

      EXPS1 = EXPS1 + DETHESS


C________________ Step 3: add terms that need number of observations or of zones
 
                                                ! Loop for counting observations
      
      DO I=1,NSTAT                                   ! Loop over state variables
        IF (OBSCLASS(I,1) .NE. 0) THEN                       ! Class is not void
           NOFOBSTYPE=INT(STAT(I,11))     
           EXPS1 = EXPS1 - REAL(NOFOBSTYPE)*DLOG(1./STAT(I,3)) !-n*ln(1/alpha)
        ENDIF    ! Void class?
      ENDDO     ! Next state variable

      DO I=1,NTYPAR                                ! Loop over parameters types
         IF (INORPAR2(I) .NE. INORPAR2(I+1) .AND.STAT(I+NSTAT,2) .NE.0)
     ;   THEN
           NOFZON=INT(STAT(I+NSTAT,11))                     ! -n*ln(alpha1/tauw)
           EXPS1 = EXPS1 - REAL(NOFZON)* DLOG(STAT(1,3)/STAT(I+NSTAT,3))
          ENDIF
      ENDDO

      RETURN
      END

********************************************************************
********************************************************************
********************************************************************

      SUBROUTINE COMP_EXP_LIKELIHOOD_2_3
     ;(DETHESS  ,EXPS2    ,EXPS3   ,NDEVS    ,NPAR     
     ;,NSTAT    ,NTYPAR   ,NUMTOBS ,INORPAR2 ,OBSCLASS 
     ;,STAT)


C______________________________________________ Step 0: Declaration of variables

      IMPLICIT NONE

                                                   ! EXTERNAL VARIABLES: SCALARS
      INTEGER*4 NDEVS, NUMTOBS, NPAR, NTYPAR,NSTAT
      REAL*8 DETHESS, EXPS2, EXPS3
                                                    ! EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 OBSCLASS(NSTAT, NDEVS+1), INORPAR2(NTYPAR+1)
       
      REAL*8 STAT(40,11)
                                                   ! INTERNAL VARIABLES: SCALARS
      INTEGER*4 I,NOFZON, NOFOBSTYPE


C___________ Step 1: calculate part of expected likelihood function with alpha 1

      DO I = 1,NSTAT                                  !Loop over state variables
        IF (OBSCLASS(I,1) .NE. 0) THEN
           EXPS3 = EXPS3 +  STAT(I,1)
        ENDIF
      ENDDO

      DO I = 1, NTYPAR                                    ! Loop over parameters
                                                 ! (1/alpha) * penalty criterion
        EXPS3=EXPS3 +  STAT(I+NSTAT,1)
 
      ENDDO

      EXPS3 = EXPS3 + DETHESS

C______ Step 2: calculate part of the expected likelihood function using alpha 2

      EXPS2 = DETHESS +  NUMTOBS + NPAR
      IF (STAT(1,4).NE.0D0) THEN
         EXPS2 = EXPS2 + (NUMTOBS + NPAR) * DLOG(STAT(1,4))
      ELSE
         WRITE(25,*)'PROBLEMAS EN EL CALCULO DE LOG(0)'
         RETURN
      END IF

C________________ Step 3: add terms that need number of observations or of zones
 
                                                ! Loop for counting observations
 
      
      DO I=1,NSTAT                                   ! Loop over state variables
        IF (OBSCLASS(I,1) .NE. 0) THEN                       ! Class is not void
           NOFOBSTYPE=INT(STAT(I,11))     
           EXPS2 = EXPS2 - REAL(NOFOBSTYPE) * DLOG(STAT(I,2))
        ENDIF    ! Void class?
      ENDDO     ! Next state variable

      DO I=1,NTYPAR                                ! Loop over parameters types
         IF (INORPAR2(I) .NE. INORPAR2(I+1) .AND.STAT(I+NSTAT,2) .NE.0)
     ;   THEN
           NOFZON=INT(STAT(I+NSTAT,11))                     ! -n*ln(alpha1/tauw)
           EXPS2 = EXPS2 - REAL(NOFZON) * DLOG(STAT(I+NSTAT,2))
          ENDIF
      ENDDO


      RETURN
      END





