      SUBROUTINE RESIDUALS
     ;(IDIMCOV  ,IOINV      ,NBANDCOV   ,NDEVS        
     ;,NSTAT    ,NTYPAR     ,NUMTOBS    ,NPAR      ,TOTALMAX
     ;,TOTALMEAN,TOTALMIN   ,TOTALSTDEV ,COVINV    ,COVPAR   ,DEVICESTAT 
     ;,IODEVICE ,ITYPEPAR   ,OBSCLASS ,RESID    ,RESIDPAR
     ;,STAT     ,VOBSC      ,WORK)

********************************************************************************
*
* PURPOSE
*
* This subroutine calculates the residuals and their statistics (mean, std, dev.
* contribution to the objective function of each device and the autocorrelation)
* These statistics are first calculated per device, then per 
* measurement type and then for parameters and then for all residuals.
*
* DESCRIPTION
*
*  - Step 1: Define RESID to contain resids and VOBSC with weighted residuals
*  - Step 2: Loop over meas. types 
*      - Step 2.1: Loop over the devices defining current type
*           - Step 2.1.a: Fill array with residuals & get min and max
*           - Step 2.1.b: Calculate statistics of this device
*           - Step 2.1.c: Update mean and std dev. of this type of measurement
*  - Step 3: Loop over parameter types
*      - Step 3.1: Loop over the zones defining current type
*           - Step 3.1.a: Fill array with residuals & get min and max
*           - Step 3.1.b: Calculate statistics of this zone
*           - Step 3.1.c: Update mean and std dev. of this type of parameter
*  - Step 4: Calculate final stats for all residuals
* 
* EXTERNAL VARIABLES: SCALARS
*
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NUMTOBS                Total number of observations
*  NDEVS                  Number of devices    
*  NPAR                   Number of estimated parameters
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)  
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)     
*
* EXTERNAL VARIABLES: ARRAYS
*    
*  COVINV                 Inverse of the covariance matrix of observations
*  COVPAR                 Inverse of the covariance matrix  of parameters
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR 
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Status for calc. of obs.                    
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*                         Column 6: Covariance matrix type
*                         Column 9: Related problem
*  OBSCLASS               See stat_outp                  
*  STAT                   Contains statistical parameters and properties
*                         for all parameter types and for all observation types. 
*  TOBS                   Time of observation            
*  VOBS                   On incoming: observation value.
*                         Is redefined to carry residuals
*  VOBSC                  On incoming:simulated value corresponding to obs. 
*                         Is redefined to carry weighted residuals
*
* INTERNAL VARIABLES: SCALARS
*
*  AUTOCO1, AUTOCO2, AUTOCO3; dummy variables to calculate the autocorrelation
*  AUTOCORRELATION       the autocorrelation of a device
*  HISTOGRAM             How many times the standard deviation is the mean residual 
*                        different from the expected mean residual
*  ITYPMEAS,ND ,NO,NZO, ITYPPAR  do loop variable
*  MAX                   maximum of weighted residuals of observations of a 
*                        device
*  MEAN                  mean of weighted residuals of observations of a device 
*  MEANYI                mean of all residuals but last residual
*  MEANYI1               mean of all residuals but first residual
*  MIN                   minimum of weighted resids of observations of a device  
*  NOF                   index of starting location of data of device in vobs 
*  NOFOBS                number of obsevations of device
*  NOL                   index of last entry of data of device in vobs
*  STDDEV                standard deviation of weighted residuals of observations 
*                        of a device 
*  TOTALMAXIMUM          maximum of all residuals
*  TOTALMEAN             mean of all residuals
*  TOTALMINIMUM          minimum of all rsiduals
*  TOTALNOBS             total number of observations
*  TOTALSTDDEV           standard deviation of all residuals
*  TYPEMAX, TYPEMIN, TYPESTDDEV  same as max, min, stddev, but for a type
*                                of parameter or state variable
*  TYPEMEAN              mean of weighted residuals of a type of state
*                        variables or parameters
*  TYPENOBS              counts the number of observations that are available 
*                        for a state variable type
*  WEIGHT                1/stddev of prior information 
*  WEIGHTEDRESID         a weighted residual
*
* INTERNAL VARIABLES: ARRAYS
*
*  DEVICESTAT             row 1: calculated mean residual of this device              
*                         row 2: calculated standard deviation
*                         row 3: minimum residual of device 
*                         row 4: maximum residual of device 
*                         row 5: autocorrelation  of device 
*                         row 6: histogram class of device
*                         row 7: contribution to objective function
*                         column i: device numnmber i 
*
* FUNCTIONS AND SUBROUTINES REFERENCED 
*
*  WEIGHTEDRESID          Calcs weighted residuals         
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan. 2003): Revision and formatting
*
********************************************************************************
          
C______________________________________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                   ! External variables: scalars
      INTEGER*4 NDEVS,NUMTOBS, NBANDCOV, NTYPAR,NSTAT,IDIMCOV
     ;         ,NPAR,IOINV
                                                    ! External variables: arrays
      INTEGER*4 OBSCLASS(NSTAT,NDEVS+1),IODEVICE(NDEVS+1,9)
     ;         ,ITYPEPAR(NPAR,2)


      REAL*8 VOBSC(NUMTOBS), STAT(40,11)
     ;      ,COVINV(IDIMCOV),RESID(NUMTOBS)
     ;      ,DEVICESTAT(7,NDEVS),RESIDPAR(NPAR,2)
     ;      ,COVPAR(NPAR*(NPAR+1)/2), WORK(NPAR*(NPAR+1)/2)

                                                   ! Internal variables: scalars
      REAL*8 TOTALMEAN, TOTALSTDEV,TOTALMIN, TOTALMAX, MEANYI
     ;      ,MEAN, STDDEV, MIN, MAX, AUTOCORRELATION, WEIGHTEDRESID
     ;      ,TYPEMEAN,TYPEMAX, TYPEMIN, TYPESTDEV, MEANYI1
     ;      ,AUTOCO1, AUTOCO2, AUTOCO3,WEIGHTRESID

      INTEGER*4 ITYPMEAS,ND ,NO, ITYPPAR,NOF, NOL, NOFOBS
     ;         ,TOTALNOBS, TYPENOBS, i,IPAR

      TOTALMEAN = 0.D0                     ! Initialization of resid. statistics
      TOTALSTDEV = 0.D0
      TOTALMIN = 1.D50
      TOTALMAX = -1.D50
      TOTALNOBS = 0 

C_____ Step 1: Redefine VOBS to contain resids and VOBSC with weighted residuals
         


      

      DO I=1,NUMTOBS                      !filling vobsc with weighted residuals
         VOBSC(I)= WEIGHTRESID(I,IDIMCOV,NBANDCOV,NUMTOBS,COVINV,RESID)
      ENDDO

C_________________________________________________ Step 2: Loop over meas. types 

      DO ITYPMEAS = 1,NSTAT
   
        TYPEMEAN = 0.D0             ! Initializ. stats of residuals of this type
        TYPESTDEV = 0.D0
        TYPEMIN = 1.D50
        TYPEMAX = -1.D50
        TYPENOBS = 0

       IF (OBSCLASS(ITYPMEAS,1) .NE.0) THEN    ! If there are meas. of this type

C_________________________ Step 2.1: Loop over the devices defining current type

        DO ND = 2, OBSCLASS(ITYPMEAS,1)+1                     !loop over devices 
    
          MEAN = 0D0                           !initializing stat of this device
          STDDEV = 0D0
          MIN = 1D50
          MAX = -1D50
          AUTOCO1 = 0D0
          AUTOCO2 = 0D0
          AUTOCO3 = 0D0
          AUTOCORRELATION = 0D0
   
C_______________________ Step 2.1.a: Fill array with residuals & get min and max

          NOF = IODEVICE(OBSCLASS(ITYPMEAS,ND),8)    ! First meas. within device
          NOL = IODEVICE(OBSCLASS(ITYPMEAS,ND)+1,8) - 1   ! Last
          NOFOBS = NOL-NOF + 1

          DO NO=NOF,NOL               ! loop over the observations of the device

            MEAN = MEAN + VOBSC(NO)                   ! Updates Device mean
            TOTALMEAN = TOTALMEAN + VOBSC(NO)         ! Updates Total Mean
            TOTALSTDEV=TOTALSTDEV+VOBSC(NO)*VOBSC(NO) ! Updates total stdev
            TYPEMEAN = TYPEMEAN + VOBSC(NO)           ! Updates Meas. Type Mean
            TYPENOBS = TYPENOBS +1                    ! Upd. type number of obs.
            TOTALNOBS = TOTALNOBS +1                  ! Upd. Total number of obs

                                      ! Upd. device maximum and minimum residual
            IF (VOBSC(NO)  .GT. MAX) MAX = VOBSC(NO)
            IF (VOBSC(NO) .LT. MIN) MIN = VOBSC(NO)

                                  ! Upd. meas. type maximum and minimum residual
            IF (VOBSC(NO) .GT. TYPEMAX) TYPEMAX = VOBSC(NO)
            IF (VOBSC(NO) .LT. TYPEMIN) TYPEMIN = VOBSC(NO)

                                       ! Upd. total maximum and minimum residual
            IF (VOBSC(NO) .GT. TOTALMAX) TOTALMAX = VOBSC(NO)
            IF (VOBSC(NO) .LT. TOTALMIN) TOTALMIN = VOBSC(NO)

         ENDDO   ! DO NO=NOF,NOL  
    
C_______________________________ Step 2.1.b: Calculate statistics of this device

         MEANYI= MEAN-VOBSC(NOF + NOFOBS - 1)     ! Mean of all but last observ.
         MEANYI1= MEAN-VOBSC(NOF)                ! Mean of all but first observ.

         MEAN = MEAN/REAL(NOFOBS)                  
         MEANYI=MEANYI/REAL(NOFOBS-1)         
         MEANYI1=MEANYI1/REAL(NOFOBS-1)

         DEVICESTAT(1,OBSCLASS(ITYPMEAS,ND))= MEAN      ! Store stats of devices
         DEVICESTAT(3,OBSCLASS(ITYPMEAS,ND))= MIN
         DEVICESTAT(4,OBSCLASS(ITYPMEAS,ND))= MAX

C__________________ Calculate the standard deviation & autocorrelation of device

          DO NO=1, NOFOBS              
             STDDEV = STDDEV + 
     ;              (VOBSC(NOF+NO-1) - MEAN) * (VOBSC(NOF+NO-1) - MEAN)
             IF (NO .NE. NOFOBS) THEN
                 AUTOCO1= AUTOCO1 + 
     ;              (VOBSC(NOF+NO)-MEANYI1) * (VOBSC(NOF+NO-1)- MEANYI)

                 AUTOCO2 = AUTOCO2 +  
     ;          (VOBSC(NOF+NO-1) - MEANYI) * (VOBSC(NOF+NO-1) - MEANYI)

                 AUTOCO3 = AUTOCO3 + 
     ;              (VOBSC(NOF+NO)-MEANYI1) * (VOBSC(NOF+NO)-MEANYI1)
             ENDIF
          ENDDO                        

          STDDEV = DSQRT(STDDEV / DFLOAT(NOFOBS))

          IF (NOFOBS .GE. 3) THEN                   ! Autocorr can be calculated
              AUTOCORRELATION = AUTOCO1/(DSQRT(AUTOCO2)*DSQRT(AUTOCO3))
          ELSE
              AUTOCORRELATION=-2.D0         ! ERROR indicating that the autocorr.
          ENDIF                             ! could not be calculated

          DEVICESTAT(2,OBSCLASS(ITYPMEAS,ND))= STDDEV
          DEVICESTAT(5,OBSCLASS(ITYPMEAS,ND))= AUTOCORRELATION 

        ENDDO   ! DO NO=NOF,NOL                      
 
C______________ Step 2.1.c: Update mean and std dev. of this type of measurement

        TYPEMEAN = TYPEMEAN /TYPENOBS              
        DO ND = 2, OBSCLASS(ITYPMEAS,1)+1                   ! Loop over devices 
      
          NOF = IODEVICE(OBSCLASS(ITYPMEAS,ND),8)    
          NOL = IODEVICE(OBSCLASS(ITYPMEAS,ND)+1,8) - 1
          NOFOBS = NOL-NOF + 1                    ! Number of obs. within device
      
          DO NO=1, NOFOBS              
             TYPESTDEV = TYPESTDEV + 
     ;      (VOBSC(NO-1+NOF) - TYPEMEAN) * (VOBSC(NO-1+NOF) - TYPEMEAN)
          ENDDO
      
        ENDDO

        TYPESTDEV= DSQRT(TYPESTDEV/FLOAT(TYPENOBS))

        STAT(ITYPMEAS,5) = TYPEMEAN           ! Store stats of measurement types
        STAT(ITYPMEAS,7) = TYPEMIN
        STAT(ITYPMEAS,8) = TYPEMAX
        STAT(ITYPMEAS,6) = TYPESTDEV

        ENDIF ! IF (OBSCLASS(ITYPMEAS,1) .NE.0)  (Non void class???)
      ENDDO  ! DO ITYPMEAS = 1,NSTAT



C___ Step 2.1: For each parameter to be estimated calculates residuals , related 
C___           weight, component at PARC and identifies parameter type


      IF (IOINV.EQ.0) RETURN !in this case there are no paramete residuals

        
	DO I= 1, NPAR*(NPAR+1)/2
         IF (COVPAR(I).GE.0D0) WORK(I)=DSQRT(COVPAR(I))
         IF (COVPAR(I).LT.0D0) WORK(I)=(-1D0)*DSQRT(ABS(COVPAR(I)))
      ENDDO     



C_____________________ Step 2.3: Calculates final product and updates obj. func.

      DO ITYPPAR=1,NTYPAR
        IF (STAT(ITYPPAR+NSTAT,2).NE.0) THEN
          TYPEMEAN = 0              ! Initializing stats of residuals of this type
          TYPESTDEV = 0.D0
          TYPEMAX = -1D50
          TYPEMIN = 1D50
          TYPENOBS = 0

          DO IPAR=1,NPAR
            IF (ITYPPAR.EQ.ITYPEPAR(IPAR,1)) THEN

                WEIGHTEDRESID=RESIDPAR(IPAR,2)
                TYPENOBS=TYPENOBS+1
                TOTALNOBS =TOTALNOBS +1
                                                    ! Updates param. type stats.
                IF (WEIGHTEDRESID  .GT. TYPEMAX) TYPEMAX = WEIGHTEDRESID
                IF (WEIGHTEDRESID .LT. TYPEMIN) TYPEMIN = WEIGHTEDRESID
                IF (WEIGHTEDRESID .GT. TOTALMAX)TOTALMAX = WEIGHTEDRESID
                IF (WEIGHTEDRESID .LT. TOTALMIN)TOTALMIN = WEIGHTEDRESID
                TYPEMEAN = TYPEMEAN + WEIGHTEDRESID
                TYPESTDEV=TYPESTDEV+ WEIGHTEDRESID * WEIGHTEDRESID

                                                         ! Updates global stats
                TOTALMEAN = TOTALMEAN + WEIGHTEDRESID
                TOTALSTDEV=TOTALSTDEV+ WEIGHTEDRESID * WEIGHTEDRESID

            END IF   ! ITYPPAR.EQ.ITYPEPAR(IPAR)
        END DO  ! IPAR=1,NPAR        

        TYPEMEAN=TYPEMEAN/FLOAT(TYPENOBS)
        TYPESTDEV=DSQRT(TYPESTDEV/FLOAT(TYPENOBS)-TYPEMEAN*TYPEMEAN)
        STAT(ITYPPAR+NSTAT,5)= TYPEMEAN        ! Store stats for parameter types
        STAT(ITYPPAR+NSTAT,6)= TYPESTDEV
        STAT(ITYPPAR+NSTAT,7)= TYPEMIN
        STAT(ITYPPAR+NSTAT,8)= TYPEMAX    
       ENDIF
      END DO ! ITYPPAR=1,NTYPAR

      TOTALMEAN=TOTALMEAN/FLOAT(TOTALNOBS)
      TOTALSTDEV=DSQRT(TOTALSTDEV/FLOAT(TOTALNOBS)-TOTALMEAN*TOTALMEAN)
       

      RETURN
      END
