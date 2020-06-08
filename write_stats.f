      SUBROUTINE OUTPUT
     ;(ALPHANEG ,CPERO       ,DOF        ,EXPS1     ,EXPS2     ,EXPS3
     ;,FOVERDOF ,IDIMCOV     ,IODETHESSOK,IOINV     ,IORESLIST ,LIKS1      
     ;,LIKS2    ,LIKS3       ,MAINF      ,NDEVS     ,NPAR     
     ;,NSTAT    ,NTYPAR      ,NUMTOBS    ,NZPAR     ,OBJF      ,TOTALMAX 
     ;,TOTALMEAN,TOTALMIN    ,TOTALSTDDEV,COVINV    ,DEVICESTAT,DEVNAME  
     ;,INORPAR2 ,IODEVICE    ,IOLG_PAR   ,ITYPEPAR  ,IVPAR     
     ;,OBSCLASS ,PARC        ,RESID      ,RESIDPAR  ,STAT     
     ;,TOBS     ,TYPENAME    ,VOBSC
     ;!NUEVOS
     ;,IDIMWGT  ,IPNT_PAR    ,COVAR      ,IDIMHESS  ,STPAR)

********************************************************************************
*
* PURPOSE
*   
*   Echoes some statistics to the main outputfile
*
* DESCRIPTION
*
*   - Step 1: headers & general statistics
*   - Step 2: likelihood functions
*   - Step 3: statistics by data type
*   - Step 4: weighted residuals per device
*   - Step 5: weighted residuals per parameter type
*   - Step 6: statistics for all weighted residuals
*   - Step 7: list of all residuals
*   - Step 8: list of residuals of the parameters 
*
* EXTERNAL VARIABLES: SCALARS
*
*  ALPHANEG               Logical=true if some of the calculated alphas<0
*  CPERO                  Contribution per observation. Calculated as 
*                             objF/(#observation + #prior information)
*  DOF                    Degrees Of Freedom                     
*  EXPS1                  Expected likelihood acc. to form 1
*  EXPS2                  Expected likelihood acc. to form 2
*  EXPS3                  Expected likelihood acc. to form 3
*  FOVERDOF               objF/DOF
*  IDIMCOV                Used to dimension array COVINV
*  IODETHESSOK            Indicates if problems were encountered with the determinant of
*                         the hessian.
*                         If 0: the determinant of the hessian is acceptable.
*                         If 1: the hessian is singular
*  LIKS1                  Log-likelihood function acc. to form. 1             
*  LIKS2                  Log-likelihood function acc. to form. 2             
*  LIKS3                  Log-likelihood function acc. to form. 3             
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NDEVS                  Number of devices         
*  NSTAT                  Maximum number of state variables whose data is used  
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMTOBS                Total number of observations                          
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  OBJF                   Objective function          
*  TOTALMAX               Maximum residual
*  TOTALMEAN              Mean residual                      
*  TOTALMIN               Minimum residual
*  TOTALSTDDEV            Residual standard deviation
*                         for calibration (used for dimensioning)               
*   
*  EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix                      
*  DEVICESTAT             row 1: calculated mean residual of this device              
*                         row 2: calculated standard deviation
*                         row 3: minimum residual of device 
*                         row 4: maximum residual of device 
*                         row 5: autocorrelation  of device 
*                         row 6: histogram class of device
*                         row 7: contribution to objective function
*                         column i: device numnmber i 
*  DEVNAME                Device name                                           
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
*  IOLG_PAR               Array containing all logarithmic options of           
*                         estimated parameters                                  
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  OBSCLASS               Array containing number of devices where a particular
*                         type of state var. was measured
*                         COLUMN 1 contains number of devices, while COLUMNS 2-?
*                         contain identifiers to those devices
*                         ROW is related to state var. type
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*  PARM                   Vector containing measured values for all             
*  STAT                   Contains statistical parameters and properties
*                         for all parameter types and for all observation types. 
*  STPAR                  Vector containing standard deviation errors of        
*                         all parameters prioo information                      
*  TOBS                   Time of observation                                   
*  TYPENAME               Array containing  the names of the state var. and 
*                         parameter types in the same order as OBSCLASS and STAT
*  VOBS                   Observation value                                     
*  VOBSC                  Value of simulated value corresponding to observation 
*                         parameters                                            
*
* INTERNAL VARIABLES: SCALARS
*
*  COBJF                  Fraction of tatal objective function contributed
*                         by partial objective function
*  I                      Dummy variable for do loops
*  LASTVAL                The last nonzero value of the array inorpar
*  NUM                    The number of a device
*  WOBJF                  Lambda*sum(weighted residuals)
*  ESTIMPAR               Logical, = true when a parameter is estimated
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  CHECK_CONSIST          To check wether all estimated parameters have
*                         estimation weights.
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan. 2003): Revision and formatting
*
********************************************************************************

      IMPLICIT NONE

                                                   ! External variables: scalars
      INTEGER*4 NTYPAR, NDEVS,MAINF,NSTAT,NUMTOBS,NZPAR
     ;         ,IDIMCOV ,IODETHESSOK,IORESLIST,NPAR,IOINV,IDIMWGT
     ;         ,IDIMHESS
      REAL*8 OBJF, DOF, CPERO, FOVERDOF,LIKS1, LIKS2,LIKS3
     ;      ,EXPS1, EXPS2, EXPS3, TOTALMEAN, TOTALMAX, TOTALMIN
     ;      ,TOTALSTDDEV
      LOGICAL ALPHANEG
                                                    ! External variables: arrays
      CHARACTER*10 TYPENAME(NSTAT+NTYPAR) ,DEVNAME(NDEVS)*10
      REAL*8 STAT(40,11),DEVICESTAT(7,NDEVS)
     ;       ,VOBSC(NUMTOBS) ,TOBS(NUMTOBS),PARC(NPAR)
     ;       ,RESIDPAR(NPAR,2),COVINV(IDIMCOV),RESID(NUMTOBS)
     ;       ,COVAR(IDIMHESS),STPAR(NZPAR)
      INTEGER*4 OBSCLASS(NSTAT, NDEVS+1), INORPAR2(NTYPAR+1)
     ;         ,IOLG_PAR(NTYPAR),IODEVICE(NDEVS+1,9),IVPAR(NZPAR,4)
     ;         ,ITYPEPAR(NPAR,2),IPNT_PAR(NZPAR*IDIMWGT)
                                                   ! Internal variables: scalars
      REAL*8 COBJF,WOBJF,STDEV,STDDEV,PARAMVAL
      INTEGER*4 I,J,NUM,NOL,NOF,NOFOBS,K,HISTOGRAM, ZONENR,IPAR,IPLACE
     ;         ,IPOS
     ;         
      LOGICAL LPARESTIM

C__________________________________________ Step 1: WARNINGS AND  & general stats.

      WRITE(MAINF,1000)
 1000 FORMAT(//,36X,'OUTPUT STATISTICS',/,36X,'====== ==========',//)

       CALL CHECK_CONSIST
     ;(MAINF,NTYPAR,NZPAR,INORPAR2,IVPAR,STAT,TYPENAME)
      

      IF (IODETHESSOK.EQ. 1) 
     ;   WRITE(MAINF,2500)
 2500 FORMAT(//,' WARNING: HESSIAN MATRIX IS NOT POSITIVE DEFINITE',/
     ;         ,'          FOLLOWING STATS. COULD NOT BE CALCULATED:',/
     ;         ,'          - ALPHA-PARAMETERS',/
     ;         ,'          - A POSTERIORI COVARIANCE MATRIX',/
     ;         ,'          - CORRELATION MATRIX',/
     ;         ,'          - LIKELIHOOD FUNCTIONS',/
     ;         ,'          - A POSTERIORI STANDARD DEVIATION',/
     ;         ,'          - CONFIDENCE INTERVALS OF ESTIMATED PARAM.',/
     ;         ,'          - MODEL SELECTION CRITERIA',/)



      IF (ALPHANEG .EQV. .TRUE.)
     ;   WRITE (MAINF,3000) 
 3000 FORMAT(//,' WARNING: SOME OF THE ALPHA-PARAMETERS ARE NEGATIVE.'/
     ;         ,'          THE EXPECTD LOG LIKELIHOOD FUNCTION BASED',/
     ;         ,'          ON ALPHAS EQN 18-19-20 COULD NOT BE ',/
     ;         ,'          CALCULATED  ',/)

        

      WRITE(MAINF,4000) OBJF,INT(DOF),CPERO,FOVERDOF
 4000 FORMAT(36X,'GENERAL STATISTICS',/,36X,'------- ----------',//
     ;      ,10X,'OBJECTIVE FUNCTION                     = ',E13.5,/
     ;      ,10X,'DEGREES OF FREEDOM                     = ',I13,/
     ;      ,10X,'MEAN CONTRIBUTION PER OBSERVATION      = ',E13.5,/
     ;      ,10X,'OBJECTIVE FUNCTION /DEGREES OF FREEDOM = ',E13.5,/)


C__________________________________________________ Step 2: Likelihood functions

      IF (IODETHESSOK .EQ.0 .AND. IOINV.GT.0)

     ; WRITE(MAINF,5000) LIKS1,LIKS2,LIKS3
 5000 FORMAT(//,39X,'LIKELIHOOD',/,39X,'----------',//
     ;         ,10X,'TAU(I) EQUAL TO  N(I)/F(I)           = ',E12.5,/
     ;         ,10X,'TAU(I) EQUAL TO 1/LAMBDA(I)          = ',E12.5,/
     ;         ,10X,'TAU SET TO ONE                       = ',E12.5,//)

      IF ((IODETHESSOK .EQ.0 .AND. ALPHANEG .EQV. .FALSE.)
     ;    .AND. IOINV.GT.0)
     ; WRITE(MAINF,5001) EXPS1,EXPS2,EXPS3
 5001 FORMAT(//,35X,'EXPECTED LIKELIHOOD',/,
     ;          35X,'-------- ----------',//
     ;         ,10X,'ALPHAS ACCORDING TO EQNS 18-19-20    = ',E12.5,/
     ;         ,10X,'ALPHAS ACCORDING TO EQNS 21-22-23    = ',E12.5,/
     ;         ,10X,'ALPHAS SET TO ONE                    = ',E12.5,//)

      IF ((IODETHESSOK .EQ.0 .AND. ALPHANEG .EQV. .TRUE.)
     ;    .AND. IOINV.GT.0)
     ; WRITE(MAINF,5002) EXPS2,EXPS3
 5002 FORMAT(//,35X,'EXPECTED LIKELIHOOD',/,
     ;          35X,'-------- ----------',//
     ;         ,10X,'ALPHAS OF  EQNS 18-19-20: WAS  NOT CALCULATED ',/
     ;         ,10X,'ALPHAS ACCORDING TO EQNS 21-22-23    = ',E12.5,/
     ;         ,10X,'ALPHAS SET TO ONE                    = ',E12.5,//)


C_______________________________________________ Step 3: statistics by data type

      WRITE(MAINF,6000)
 6000 FORMAT(//,32X,'STATISTICS PER DATA TYPE',/
     ;         ,32X,'---------- --- ---- ----',//
     ;         ,3X,'NAME      ',3X,'OBJ. FUNC.',3X,'  * LAMBDA'
     ;         ,3X,'* LAMBDA %',3X,'    ALFA 1',3X,'    ALFA 2'
     ;         ,3X,' LOG. OPT.'
     ;       ,/,3X,'----      ',3X,'---- -----',3X,'  - ------'
     ;         ,3X,'- ------ -',3X,'    ---- -',3X,'    ---- -'
     ;         ,3X,' ---- ----')

      DO I=1,NSTAT                            ! loop over state variable types
        IF (OBSCLASS(I,1) .NE. 0) THEN                       ! Class is not void
            WOBJF = STAT(I,2)*STAT(I,1)                           ! lambda*obj f
            COBJF = (STAT(I,2)*STAT(I,1))/OBJF           ! contribution to obj f
            WRITE(MAINF,7000) TYPENAME(I), STAT(I,1), WOBJF, 100*COBJF,
     ;      STAT(I,3),  STAT(I,4),0
         ENDIF
      ENDDO
 
      DO I=1,NTYPAR * MIN(IOINV,1)                         ! loop over parameter types
       IF (STAT(I+NSTAT,2) .GT. 0) THEN                    ! Type of  parameters estimated?
           WOBJF = STAT(I+NSTAT,2)*STAT(I+NSTAT,1)         !lambda*obj f
           COBJF = (STAT(I+NSTAT,2)*STAT(I+NSTAT,1))/OBJF  !contribution to obj f
           WRITE(MAINF,7000) TYPENAME(I+NSTAT),STAT(I+NSTAT,1)
     ;    ,WOBJF,100*COBJF  ,STAT(I+NSTAT,3),STAT(I+NSTAT,4)
     ;    ,IOLG_PAR(I)             
       ENDIF
      ENDDO

 7000 FORMAT(3X,A10,5(3X,E10.3),4X,I5)

 
C___________________________________ Step 4: Weighted residuals device to device

      WRITE(MAINF,8000)
 8000 FORMAT(//,23X,'ANALYSIS OF WEIGHTED RESIDUALS PER DEVICE',/
     ;         ,23X,'-------- -- -------- --------- --- ------',//)

      DO I=1,NSTAT                        ! Loop over state variable types 
         IF (OBSCLASS(I,1) .NE. 0) THEN                         ! Non-void class

            WRITE(MAINF,9000) TYPENAME(I)
 9000       FORMAT(37X,'----------',/,37X,A10,/,37X,'----------',/)

            WRITE(MAINF,9010) 
 9010       FORMAT(2X,'      NAME',2X,'      MEAN',
     ;             2X,'   MINIMUM',2X,'   MAXIMUM',2X,'     STDEV',
     ;             2X,'     FOBJ.',2X,' AUTOCORR.',/
     ;             2X,'      ----',2X,'      ----',
     ;             2X,'   -------',2X,'   -------',2X,'     -----',
     ;             2X,'     -----',2X,' ---------')

            DO J=2, OBSCLASS(I,1)+1            ! Loop over devices of this class

              NUM=OBSCLASS(I,J)                                  ! Device number

              IF (DEVICESTAT(5,NUM) .NE. -2) THEN         ! Autocorr. calculated

                 WRITE(MAINF,9020) DEVNAME(NUM)
     ;          ,DEVICESTAT(1,NUM),DEVICESTAT(3,NUM), DEVICESTAT(4,NUM) 
     ;          ,DEVICESTAT(2,NUM),DEVICESTAT(7,NUM), DEVICESTAT(5,NUM)

              ELSE                   

                 WRITE(MAINF,9021) DEVNAME(NUM)
     ;          ,DEVICESTAT(1,NUM),DEVICESTAT(3,NUM), DEVICESTAT(4,NUM) 
     ;          ,DEVICESTAT(2,NUM),DEVICESTAT(7,NUM),' NOT CALC.'

              ENDIF

 9020         FORMAT(2X,A10,6(2X,E10.3))
 9021         FORMAT(2X,A10,5(2X,E10.3),2X,A10)

            ENDDO    ! Next device of this class

            WRITE(MAINF,9030) TYPENAME(I),STAT(I,5),STAT(I,6)*STAT(I,6)
 9030       FORMAT(//,10X,'GLOBAL STATISTICS FOR MEASUREMENT TYPE: ',
     ;          A10,/,10X,'------ ---------- --- ----------- -----',
     ;             //,5X,' MEAN WEIGHTED RESIDUAL        = ',E10.3,/
     ;               ,5X,' VARIANCE OF WEIGHTED RESIDUAL = ',E10.3,/)

         ENDIF                                         ! Measurements available?
      ENDDO                                              ! Next measurement type

C_________________________________ Step 5: Weighted residuals per parameter type 

      WRITE(MAINF,9040)
 9040 FORMAT(//,20X,'ANALYSIS OF WEIGHTED RESIDUALS PER PARAMETER TYPE',
     ;     /,20X,'-------- -- -------- --------- --- --------- ----',//)


            WRITE(MAINF,9050) TYPENAME(I)
 9050       FORMAT(35X,'----------',/,35X,A10,/,35X,'----------',/)

            WRITE(MAINF,9060) 
 9060       FORMAT(2X,'NAME      ',2X,'      MEAN',2X,'   MINIMUM',
     ;             2X,'   MAXIMUM',2X,'     STDEV',2X,'     FOBJ.',/
     ;             2X,'----      ',2X,'      ----',2X,'   -------',
     ;             2X,'   -------',2X,'     -----',2X,'     -----',2X)

      DO I=1,NTYPAR * MIN(IOINV,1)                          ! Loop over types of param.

         NOF= INORPAR2(I)+1
         NOL= INORPAR2(I+1)
         LPARESTIM = .FALSE.

         DO J=NOF,MIN(NOL,NZPAR)
            IF(IVPAR(J,1) .NE. 0) THEN
		     DO K= IVPAR(J,1) ,IVPAR(J,2)
			    IF(IPNT_PAR(K).EQ.I) THEN
				    LPARESTIM = .TRUE. ! Some of this class are estim.
	            ENDIF
	         ENDDO
	      ENDIF
         ENDDO 

         IF (LPARESTIM .EQV. .TRUE.)
     ;       WRITE(MAINF,9070) TYPENAME(I+NSTAT),STAT(I+NSTAT,5)
     ;       ,STAT(I+NSTAT,7),STAT(I+NSTAT,8),STAT(I+NSTAT,6)
     ;       ,STAT(I+NSTAT,2)*STAT(I+NSTAT,1)

      ENDDO  ! Next class

 9070 format(2X,a10,5(2X,E10.3))

C_________________________________ Step 6: statistics for all weighted residuals

      WRITE(MAINF,9080) TOTALMEAN,TOTALSTDDEV*TOTALSTDDEV,TOTALMIN
     ;                 ,TOTALMAX

 9080 FORMAT(//,20X,'STATISTICS FOR ALL WEIGHTED RESIDUALS',/
     ;         ,20X,'---------- --- --- -------- ---------',//
     ;         ,10X,'MEAN WEIGHTED RESIDUAL        = ',E10.3,/
     ;         ,10X,'VARIANCE OF WEIGHTED RESIDUAL = ',E10.3,/
     ;         ,10X,'MINIMUM WEIGHTED RESIDUAL     = ',E10.3,/
     ;         ,10X,'MAXIMUM WEIGHTED RESIDUAL     = ',E10.3,/)

C_____________________________________________________ Step 7: List of residuals

      
      IF (IORESLIST .GE. 1) THEN                              ! If so desired...

         WRITE(MAINF,9090)  
 9090    FORMAT(/,36X,'LIST OF RESIDUALS',/,36X,'---- -- ---------',//,
     ;          37X,'---------------',/,
     ;          37X,'STATE VARIABLES',/,
     ;          37X,'---------------',//,
     ;           2X,' DEV. NAME',2X,'   TYPE   ',2X,'MEAS. TIME',
     ;           2X,'    RESID.',2X,'WGT. RESID',2X,' HISTOGRAM',/,
     ;           2X,' ---- ----',2X,'   ----   ',2X,'----- ----',
     ;           2X,'    ------',2X,'---- -----',2X,' ---------')

         DO I=1,NSTAT                    ! Loop over state variables

            IF (OBSCLASS(I,1) .GT. 0) THEN    ! Class is not void

               DO J=2,OBSCLASS(I,1)+1
                 NOF=IODEVICE(OBSCLASS(I,J),8)
                 NOL=IODEVICE(OBSCLASS(I,J)+1,8)-1
                 NOFOBS=NOL-NOF+1
                 DO K=NOF,NOL
                      STDEV=(1.D0/dSQRT(COVINV(K)))
                      WRITE(MAINF,9100) DEVNAME(OBSCLASS(I,J)),
     ;                  TYPENAME(IODEVICE(OBSCLASS(I,J),1)),TOBS(K)
     ;                 ,RESID(K),VOBSC(K),HISTOGRAM(STDEV,RESID(K)) 
                 ENDDO
              ENDDO                                      !DO J=2,OBSCLASS(I,1)+1
            ENDIF                                !IF (OBSCLASS(I,1) .GT. 0) THEN
         ENDDO                                                     !DO I=1,NSTAT
 9100    FORMAT(2X,A10,2X,A10,3(2X,E10.3),I5)

C___________________________________ Step 8: List of residuals of the parameters 

         WRITE(MAINF,9110)  
 9110    FORMAT(//,39X,'----------',/,39X,'PARAMETERS',/,
     ;             39X,'----------',//,
     ;           2X,' NUM. ZONE',2X,'   TYPE   ',2X,
     ;           2X,'    RESID.',2X,'WGT. RESID',2X,/,
     ;           2X,' ---- ----',2X,'   ----   ',2X,
     ;           2X,'    ------',2X,'---- -----')

          DO IPAR=1, NPAR                               ! Loop over parameter zones 
	       
                ZONENR=ITYPEPAR(IPAR,2)
                WRITE(MAINF,9101) ZONENR
     ;          ,TYPENAME(ITYPEPAR(IPAR,1)+10),RESIDPAR(IPAR,1)
     ;          ,RESIDPAR(IPAR,2)            				
			
          ENDDO

         

           
      ENDIF                                         !IF (RESLIST .GE. 1) THEN  
 9101 FORMAT(3X,I5,8X,A10,3(2X,E10.3))

C____________________________________ Step 9: Write a header for the outputfile

      WRITE(MAINF,10000)
10000 FORMAT(//,25X,'A POSTERIORI STATISTICS OF ESTIMATED PARAMETERS',/
     ;         ,25X,'= ========== ========== == ========= ==========',//
     ;          ,2X,'      TYPE',2X,' NUM. ZONE',2X,'  ESTIMATE'
     ;          ,2X,'PRIOR STD.',2X,'POST. STD.',2X,'LOW. CONF.'
     ;          ,2X,'UPP. CONF.',/
     ;          ,2X,'      ----',2X,' ---- ----',2X,'  --------'
     ;          ,2X,'----- ----',2X,'----- ----',2X,'---- -----'
     ;          ,2X,'---- -----',/)


C______________________ Step 2: calculate the 95% confidence interval boundaries
      IPAR = 0
      DO IPAR=1,NPAR                               ! Loop over parameter zones
	 IPLACE = (IPAR*(IPAR+1))/2
         STDDEV = DSQRT(COVAR(IPLACE))
         ZONENR= ITYPEPAR(IPAR,2)
         IPOS= INORPAR2(ITYPEPAR(IPAR,1))+ ITYPEPAR(IPAR,2)
         PARAMVAL=PARC(IPAR)
         
                
C___________________________________ Step 3: Echoes confidence intervals to file 

         WRITE(MAINF,10010) TYPENAME(ITYPEPAR(IPAR,1)+NSTAT),ZONENR
     ;             ,PARAMVAL,STPAR(IPOS),STDDEV,PARAMVAL-2*STDDEV
     ;             ,PARAMVAL+2*STDDEV

                                      ! Estimated

      ENDDO                                                         ! I=1, NZPAR    

10010 FORMAT(7X,A10,2X,I5,5(2X,E10.3)) 

!=================================LIST OF INTERPOLATED PARAMETERS
      
*10030 FORMAT( //,'VALUES OF CALCULATED PARAMETERS',
*     ;       /,'      TYPE        ZONENR       VALUE  GROUP',/
*     ;        ,'      ----        ------       -----  -----')
*
*	IF (NPAR.NE.NPARDET) THEN
*	   WRITE(MAINF,10030)

	! Count groups
*         NGROUPS = 0
*	   DO I=1,NZPAR
*	      NGROUPS= MAX(NGROUPS,IVPAR(I,3))
*	   ENDDO

*         DO IGROUP =1,NGROUPS
*	      DO I=1,NZPAR
*	         IF (IVPAR(I,3).EQ.IGROUP) THEN
*	            LPARESTIM = .TRUE.
*    	            IF (IVPAR(I,1).NE.0) THEN
*			       DO K= IVPAR(I,1),IVPAR(I,2)
*			       	 IF (IPNT_PAR(K).NE.I .AND. IPNT_PAR(K).NE.0)
*     ;                 LPARESTIM = .FALSE.
*			       ENDDO
*		         ENDIF

*		         IF (LPARESTIM.EQ. .FALSE.) THEN
*			        J=1
*			        DO WHILE(INORPAR2(J+1).LT.I) 
*				       J=J+1
*			        ENDDO

*			        ZONENR=I-INORPAR2(J)

*			        WRITE(MAINF,10020) TYPENAME(J+NSTAT),ZONENR
*    ;                    , PARZ(I),IGROUP
*		         ENDIF
*	         ENDIF
*	      ENDDO
*        ENDDO
*	ENDIF 
*10020 FORMAT(7X,A10,2X,I5,2X,E10.3,I5)

      RETURN
      END

********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************

      SUBROUTINE OUTPUT2
     ;(NPAR,IDIMMATRIX,MAINF,NSTAT,ITYPEPAR,MATRIX,TYPENAME)

********************************************************************************
*
* PURPOSE
*
* This subroutine writes a matrix in symmetric band storage mode to file. It is 
* used to write the correlation and covariance matrices.
* 
*
* It makes the following format if there are M parameters:
*             typename1  typename1  typename1  typename2  typename2
*               zonenr     zonenr      zonenr    zonenr    zonenr
* type1,zone1  NUMBER        /          /          /        / 
*       zone2  NUMBER      NUMBER       /          /        /
* typeM,zonenr NUMBER      NUMBER      NUMBER    NUMBER    NUMBER      
*
*
*             typenameN  typename?  typename?  typename?  typenameM
*               zonenr     zonenr      zonenr    zonenr    zonenr
* typeN,zonenr NUMBER         /          /        /         /    
*     ...      NUMBER      NUMBER        /        /         /
* typeM,zonenr NUMBER      NUMBER      NUMBER    NUMBER    NUMBER   
*
* DESCRIPTION
*
*   - Step 1: Loop over all parameters
*        - Step 2: If parameter is estimated then
*              - Step 2.1: Save name and index of parameter. 
*              - Step 2.2: If 5 new columns have been gathered or if all 
*                          parameters have been checked
*                  - Step 2.2.1: Write column headers 
*                  - Step 2.2.2: Loop over rows, starting with the first row 
*                                that can have an entry
*                       - Step 2.2.1.a: Write row to file
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  COVINV                 Inverse of the covariance matrix    
*  MAINF                  Unit number of the main output file (RES.OUT)     
*  NPAR                   Total number of parameters to be estimated 
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning) 
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component. 
*  SIZE                   Number of elements of the matrix considered. 
*                         Here: idimhess
*
* EXTERNAL VARIABLES: ARRAYS
*
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IVPAR                  Vector containing estimation index for all            
*                         parameters  
*  MATRIX                 The matrix that is to be printed
*  TYPENAME               Array containing  the names of the state var. and 
*                         parameter types in the same order as OBSCLASS and STAT
*
* INTERNAL VARIABLES: SCALARS
*
*  APP                    Variable counting the number of columns that have
*                         been detected that need to be printed. When it reaches 
*                         5, the 5 columns are printed and the variable is
*                         set to 0.
*  CPTYPE                 acronym for column-parameter type. It is a scalar that
*                         has as value the INORPAR index number of the parameter
*                         type to which currently printed column is associated 
*  GV                     Dummys counter for do loops.
*  INDEX                  The index of a parameter in ivpar.
*  KLO                    Dummy counter for do loops
*  MATCOL                 The number of the matrix column that is being printed
*  MATROW                 The number of the matrix row that is being printed
*  PRINTCOL               Dummys counter for do loops.
*  ROW                    Dummys counter for do loops.
*  RPARNAME               Acronym for row-parameter name. The name of the
*                         parameter type to which the row currently being
*                         printed belongs.
*  RPTYPE                 Acronym for row-parameter type. It is a scalar that 
*                         has asvalue the INORPAR index number of the parameter
*                         type to which the currently printed row is associated 
*  ZNO                    Scalar to get the zone number of the parameter to 
*                         which the row currently being printed belongs.
*
* INTERNAL VARIABLES: ARRAYS
*
* MATCOLEL                Contains the column numbers of the columns currently
*                         being printed of matrix to be printed.
* PARNAME                 The parameter names associated to the columns being
*                         printed
* PARNO                   Array of IVPAR indexes of the parameters to which the
*                         columns currently being printed belongs
* STARTPOS                Array of indexes of position in IVPAR where a
*                         parameter type begins. Used to get the zonenr                 
* VALUE                   The value of the elements on the row being printed
*
* ZONENR                  Array of zone numbers of the parameters to which the
*                         columns currently being printed belongs
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* GETVALUE                To retrieve an element positioned at a given row and
*                         column of matrix in symmetric triangle storage mode. 
*
* HISTORY: LJS (Dec. 2002): First coding
*
***************************************************************
*                                    !step 0 : variable declaration
      IMPLICIT NONE

* external variables: scalars       
      INTEGER*4 NPAR,IDIMMATRIX,MAINF,NSTAT

* external variables: arrays
      INTEGER*4 ITYPEPAR(NPAR,2)
	REAL*8 MATRIX(IDIMMATRIX)
      CHARACTER*10 TYPENAME(40)

* internal variables: scalars
      integer*4 ipar,jpar,kpar,npar_printed,npartoprint
     ;          ,n,ipos

* internal variables: arrays
      INTEGER*4 ZONENR(5)
	REAL*8 VALUES(5)
	CHARACTER*10 TITLES(5), TITLE

100   FORMAT(18X,5(5X,A10))
200   FORMAT(13X,5(10X,I5))
300   FORMAT(A10,I5,5(5X,E10.3))

      !INITIALIZE
      NPAR_PRINTED = 0  
      IPAR = 1
    
	DO WHILE (IPAR .LE. NPAR)  !LOOP OVER PARAMETERS
           

C                          !INITIALIZE
	    NPARTOPRINT = MIN(5,NPAR-IPAR+1)
	    DO N=1,5
	       TITLES(N)= '          '
	       ZONENR(N)=0
	    ENDDO
	    
c           !gather column titles
	    DO N=1,NPARTOPRINT
	       TITLES(N)= TYPENAME(NSTAT+ITYPEPAR(IPAR+N-1,1))
	       ZONENR(N)= ITYPEPAR(IPAR+N-1,2)
	    ENDDO
	    

		WRITE(MAINF,100) TITLES(1),TITLES(2),TITLES(3),TITLES(4)
     ;                    ,TITLES(5)
		WRITE(MAINF,200) ZONENR(1),ZONENR(2),ZONENR(3),ZONENR(4)
     ;                    ,ZONENR(5)

          
          DO JPAR = IPAR, NPAR
	      
		   DO N=1 ,5 
	          VALUES(N)=0D0
	       END DO 
		  
		   DO N= 1,NPARTOPRINT
		      KPAR = IPAR + N - 1
		      IF(JPAR.GT.KPAR) THEN
                   IPOS = JPAR*(JPAR+1)/2 - (JPAR-KPAR)
	          ELSE
                   IPOS = KPAR*(KPAR+1)/2 - (KPAR-JPAR)
	          ENDIF
               
	          VALUES(N)= MATRIX(IPOS)
	          TITLE=TYPENAME(NSTAT+ITYPEPAR(JPAR,1))
		                       
	       ENDDO   !N
		   
		   WRITE (MAINF,300) TITLE,ITYPEPAR(JPAR,2)
     ;,VALUES(1),VALUES(2),VALUES(3),VALUES(4),VALUES(5)     

	    ENDDO  !JPAR
	    
		IPAR = IPAR + NPARTOPRINT

	ENDDO      !IPAR   
	         
             

      return 
	end


********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************

      SUBROUTINE OUTPUT3(EIGENVEC, NSTAT,TYPENAME, NPAR
     ;,EIGENVAL,MAINF,itypepar)

      

****************************************************************
* PURPOSE
* 
* this subroutine writes a matrix in full storage mode to file. 
* It is used to write the eigenvector components matrix.
*
* It makes the following format if there are M parameters:
*  eigenvalue1
*            typename1  typename1  typename1  typename2  typename2
*               zonenr     zonenr      zonenr    zonenr    zonenr
*              NUMBER      NUMBER      NUMBER    NUMBER    NUMBER    
*     ....     
*            typenameM  typenameM  typenameM  typenameM  typenameM
*               zonenr     zonenr      zonenr    zonenr    zonenr
*              NUMBER      NUMBER      NUMBER    NUMBER    NUMBER  
*  
* eigenvalue2
*            typename1  typename1  typename1  typename2  typename2
*               zonenr     zonenr      zonenr    zonenr    zonenr
*              NUMBER      NUMBER      NUMBER    NUMBER    NUMBER    
*     ....     
*            typenameM  typenameM  typenameM  typenameM  typenameM
*               zonenr     zonenr      zonenr    zonenr    zonenr
*              NUMBER      NUMBER      NUMBER    NUMBER    NUMBER  
*
* EXTERNAL VARIABLES: SCALARS
*
*  NPAR                   Total number of parameters to be estimated            
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)  
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                 
*
* EXTERNAL VARIABLES: ARRAYS
*
*  EIGENVAL               The eigenvalues of the covariance matrix
*  EIGENVEC               The eigenvectors of the covariance matrix
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IVPAR                  Vector containing estimation index for all            
*                         parameters  
*  TYPENAME               Array containing  the names of the state var and 
*                         parameter types in the same order as STAT
*
* INTERNAL VARIABLES  ARRAYS
*
*  VALUE                   Contains the next 5 components of an eigenvectors that will 
*                         be printed
*
* INTERNAL VARIABLSE: SCALARS:
* 
*  PARAMS                  Number of estimated zones
*  TYPES                   Identifier of parameter type
*  DUMMY,I,J,L             Dummy counters
*  PARSTOGO                Number of unprinted eigenvectors
*
* EXTERNALS
*
*  GETVALUE                Function to retrieve an element of a matrix in 
*                         symetric band storage mode.
*  GETPARZNNR              Function to get the zone number of an estimated par. 
*  GETPARNAME              Function to get the name of an estimated parameter
*
***************************************************************

      IMPLICIT NONE

C                                                   External variables: scalars
      INTEGER*4  NSTAT, NPAR,MAINF
C                                                    External variables: arrays
      REAL*8 EIGENVEC(NPAR, NPAR), EIGENVAL(NPAR)
      INTEGER*4 ITYPEPAR(NPAR,2)
      CHARACTER*10 TYPENAME(40)
C                                                    Internal variables: arrays
      REAL*8  VALUE(5)
      INTEGER*4 ZONE(5)
      CHARACTER*10 NAME(5)
C                                                   Internal variablse: scalars
      REAL*8 DUMMY
      INTEGER*4 I,J,K,KK,PARSTOGO,L,M

C_____________________________________________________________format statements
50    FORMAT(/,A11,I3,A3,G10.4)
100   FORMAT(A10,' ',A10,' ',A10,' ',A10,' ',A10)
200   FORMAT(I3,' ',I10,' ',I10,' ',I10,' ',I10)
300   FORMAT(G10.4,' ',G10.4,' ',G10.4,' ',G10.4,' ',G10.4,/,/)
C___________________________________________________Step 1 sort the eigenvalues
    
      DO I=1,NPAR
        DO J=I,NPAR
           IF (EIGENVAL(I) .LT. EIGENVAL(J)) THEN
               DUMMY=EIGENVAL(J)
               EIGENVAL(J)=EIGENVAL(I)
               EIGENVAL(I)=DUMMY
               DO K=1, NPAR                          !switching rows of eigenvec
                  DUMMY = EIGENVEC(K,I)
                  EIGENVEC(K,I)=EIGENVEC(K,J)
                  EIGENVEC(K,J)=DUMMY
               ENDDO
           ENDIF                         !IF (EIGENVAL(I) .LT. EIGENVAL(J)) THEN
       ENDDO
      ENDDO

C______________________Step 3: print eigenvalue and contributions to eigenvector

      DO I=1,NPAR
         WRITE(MAINF,50) 'EIGENVALUE ',I,': ', EIGENVAL(I)
         WRITE(MAINF,*) 'CONTRIBUTIONS:'
         J=0
         DO WHILE (J .LT. NPAR)
             PARSTOGO=NPAR-J
             
		   DO M=1,MIN(PARSTOGO,5)
                NAME(M) =TYPENAME(NSTAT+ITYPEPAR(J+M,1))
                ZONE(M) = ITYPEPAR(J+M,2)
             ENDDO
   


**************************In case there are still 5 or more values to be printed
           IF (PARSTOGO .GE. 5) THEN 
               DO KK=1,5
                 VALUE(KK)=EIGENVEC(J+KK,I)
               ENDDO
               WRITE(MAINF,100)  NAME(1),NAME(2),NAME(3),NAME(4),NAME(5)
               WRITE(MAINF,200) ZONE(1),ZONE(2),ZONE(3),ZONE(4),ZONE(5)
               WRITE(MAINF,300) VALUE(1),VALUE(2),VALUE(3),VALUE(4)
     ;,VALUE(5)
               J=J+5
            ENDIF
*************************************In case there are only 4 values to be printed
            IF (PARSTOGO .EQ. 4) THEN 
               DO L=1,4
                  VALUE(L)=EIGENVEC(J+L,I)
               ENDDO
               WRITE(MAINF,100) NAME(1),NAME(2),NAME(3),NAME(4)
               WRITE(MAINF,200)  ZONE(1),ZONE(2),ZONE(3),ZONE(4)
               WRITE(MAINF,300) VALUE(1),VALUE(2),VALUE(3),VALUE(4)
               J=J+4
            ENDIF
***********************************In case there are only 3 values to be printed
            IF (PARSTOGO .EQ. 3) THEN 
               DO L=1,3
                 VALUE(L)=EIGENVEC(J+L,I)
               ENDDO
               WRITE(MAINF,100) NAME(1),NAME(2),NAME(3)
               WRITE(MAINF,200) ZONE(1),ZONE(2),ZONE(3)
               WRITE(MAINF,300) VALUE(1),VALUE(2),VALUE(3)
               J=J+3
            ENDIF
***********************************In case there are only 2 values to be printed
            IF (PARSTOGO .EQ. 2) THEN 
               DO L=1,2
                  VALUE(l)=EIGENVEC(J+L,I)
               ENDDO
               WRITE(MAINF,100) NAME(1),NAME(2)
               WRITE(MAINF,200) ZONE(1),ZONE(2)
               WRITE(MAINF,300) VALUE(1),VALUE(2)
               J=J+2
            ENDIF
***********************************In case there are only 1 values to be printed
            IF (PARSTOGO .EQ. 1) THEN 
               VALUE(1)=EIGENVEC(J+1,I)
               WRITE(MAINF,100)  NAME(1)
               WRITE(MAINF,200)  ZONE(1)
               WRITE(MAINF,300) VALUE(1)
               J=J+1
            ENDIF

         ENDDO
      ENDDO
      RETURN
      END

