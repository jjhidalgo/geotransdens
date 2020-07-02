       SUBROUTINE STAT_OUTPUT
     ;(IDIMCOV    ,IDIMHESS   ,IDIMWORK  ,IODIM     ,IOVAR     ,IOPRHED
     ;,ISOT       ,MAINF      ,NBANDCOV  ,NDEVS     ,NFLAGS    ,NPAR
     ;,NSTAT      ,NTYPAR     ,NUMTOBS   ,NZPAR     ,COVINV    ,COVPAR
     ;,DEVICESTAT ,DEVNAME    ,EIGENVAL  ,EIGENVEC  ,FOBJ_WGT  ,HESSINV
     ;,IFLAGS     ,INORPAR    ,IODEVICE  ,IOLG_PAR  ,ITYPEPAR  ,IVPAR
     ;,NZONE_PAR  ,OBSCLASS   ,PARC      ,PARM      ,PAR_WGT
     ;,RESID      ,STPAR      ,TOBS      ,RESIDPAR  ,VJAC      ,VOBS
     ;,VOBSC      ,WORK       ,IOINV     ,MEASTYP
     ;!NUEVOS
     ;,WGT_UNK    ,IPNT_PAR  ,IDIMWGT   ,IOPT_GS    ,MXGRPZN
     ;,NPARDET)

********************************************************************************
*
* PURPOSE
*
* This subroutine arranges all statistical output related to the inverse problem
* It calculates the necessary statistics and generates the output file
*
* This subroutine does not calculate anything itself but calls the appropriate 
* subroutines and passes the appropriate variables from one to the other.
*
* If the hessian is not positive definite a number of operations cannot be
* carried out. These are then skipped.
*
* DESCRIPTION
*
* Step 0: Declaration of variables 
* Step 1: Initialises local arrays
* Step 1.1: Redefine array inorpar so that zero-elements become equal to 
*           the last nonzero entry
* Step 2: Fill array obsclass. Groups devices by data type.
* Step 3: The partial objective functions are calculated,the contribution 
*         of each device and the total objective function
* Step 4: Calculate general problem statistics
* Step 5: Calculation and inversion of the hessian matrix
* Step 6: The alphas are calculated
* Step 7: The expected likelihood criterions are calculated
* Step 8: The statistics of the residuals are calculated
* Step 9: The estimate covariance matrix is calculated
* Step 10: The model selection criterions are calculated
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix of measurements   
*  COVPAR                 Inverse of the covariance matrix of parameters                     
*  DEVICESTAT             Array containing different stats. for each devices.
*                         COLUMN is related to a particular device.
*                         ROW description:
*                              1: Calculated mean residual of this device
*                              2: Calculated standard deviation
*                              3: Minimum residual of device
*                              4: Maximum residual of device
*                              5: Autocorrelation  of device
*                              6: Histogram class of device
*                              7: Contribution to objective function
*  DEVNAME                Device name                                           
*  EIGENVEC               Array containing the eigenvalues and eigenvectors of 
*                         the covariance matrix
*  FOBJ_WGT               Array containing all objective function weights for   
*                         state variables (heads, concentrations, etc)          
*  HESSINV                                                                      
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
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters      
*  OBSCLASS               Array containing number of devices where a particular
*                         type of state var. was measured
*                         COLUMN 1 contains number of devices, while COLUMNS 2-?
*                         contain identifiers to those devices
*                         ROW is related to state var. type
*                          1:  heads
*                          2:  concentrations
*                          3:  humidity
*                          4:  measured flux
*                          5-10:   measurement types that can be added later
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*  PARM                   Vector containing measured values for all             
*                         parameters                                            
*  PAR_WGT                Array containing objective function weights for       
*                         all estimated parameters                              
*  STPAR                  Vector containing standard deviation errors of        
*                         all parameters prioo information                      
*  TOBS                   Time of observation                                   
*  VJAC                   Jacobian matrix                                       
*  VOBS                   Observation value                                     
*  VOBSC                  Value of simulated value corresponding to observation 
*  WORK                   Workspace array used in some subroutines.             
*
* INTERNAL VARIABLES: ARRAYS
*
*  STAT                   A given row is related to a parameter or state 
*                         variable type. A given column is a particular 
*                         statistic value
*                         ROWS:
*                              1: Head level/pressure
*                              2: Concentrations
*                              3: Humidity
*                              4: flow
*                              5-10: Other types of state var. (not used by now)
*                              11: Txx
*                              12: Tyy
*                              13: Txy
*                              14: Tzz
*                              15: Txz
*                              16: Tyz
*                              17: Storage coefficient
*                              18: Recharge
*                              19: Presc. head
*                              20: stats for Presc. flow
*                              21: Leakage
*                              22: Longit.  disp.
*                              23: Transv. disp.
*                              24: Molecular diff.
*                              25: Porosity
*                              26: First order decay
*                              27: Retardation
*                              28: External conc.
*                              29: Generic param.
*                              30: Age coefficient
*                              31: Formation Factor
*                              32-40: Make your own parameter! 
*                             
*                         COLUMNS:
*                              1: partial objective function
*                              2: lambda
*                              3: alpha_1
*                              4: alpha_2
*                              5: mean weighted residual
*                              6: standard deviation of weighted residual
*                              7: minimum weighted residual
*                              8: maximum weighted residual
*                              9: correlation with normal distrib. of weighted
*                                 squared residuals
*                             10: determinant of covariance matrix
*                             11: number of observations or number of estimated
*                                 zones
*  TYPENAME               Array containing  the names of the state var. and 
*                         parameter types in the same order as OBSCLASS and STAT
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCOV                Used to dimension array COVINV
*  IDIMHESS               Used to dimension array HESS. It is equal to          
*                         NPAR*(NPAR+1)/2                                       
*  IOPRHED                Indicates whether the flow state variable state is    
*                         preasure (set to 1) or head (set to 0)                
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NBANDCOV               Band width of the inverse covariance matrix           
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
* INTERNAL VARIABLES: SCALARS
*
*  AIC1                   AIC model selection criteria according to form. 1
*  AIC2                   AIC model selection criteria according to form. 2
*  AIC3                   AIC model selection criteria according to form. 3
*  ALPHANEG               Logical=true if some of the calculated alphas<0
*  ARMA1                  ARMA selec. criteria according to form. 1
*  ARMA2                  ARMA selec. criteria according to form. 2
*  ARMA3                  ARMA selec. criteria according to form. 3
*  BIC1                   BIC  model selection criteria according to form. 1
*  BIC2                   BIC  model selection criteria according to form. 2
*  BIC3                   BIC  model selection criteria according to form. 3
*  CPERO                  Contribution per observation. Calculated as 
*                             objF/(#observation + #prior information)
*  D1                     Used to calculate hessian determinant
*  D2                     Used to calculate hessian determinant                 
*  DC1                    Kashyap model sel. crit. acc. to form. 1
*  DC2                    Kashyap model sel. crit. acc. to form. 2
*  DC3                    Kashyap model sel. crit. acc. to form. 3
*  DETHESS                Determinant of the hessian matrix
*  DOF                    Degrees Of Freedom                     
*  EXPS1                  Expected likelihood acc. to form 1
*  EXPS2                  Expected likelihood acc. to form 2
*  EXPS3                  Expected likelihood acc. to form 3
*  FOVERDOF               objF/DOF
*  I                      Dummy counter
*  IER                    Error index at LIN1VP                          
*  IOCORRMAT              Option for writing correl. matrix of estim. param.
*  IODETHESSOK            Indicates if problems were encountered with the determinant of
*                         the hessian.
*                         If 0: the determinant of the hessian is acceptable.
*                         If 1: the hessian is singular
*  IOEIGENVEC             Option for writing eigenvectors and eigenvalues
*  IORESLIST              Option for writing list of all residuals  
*  J                      Dummy counter
*  LIKS1                  Log-likelihood function acc. to form. 1             
*  LIKS2                  Log-likelihood function acc. to form. 2             
*  LIKS3                  Log-likelihood function acc. to form. 3             
*  MAXNOFOBS              Maximum number of observations in any device
*  ND                     Dummy counter of devices
*  NOF                    First observation defining a particular device.
*  NOFOBS                 Number of observations defining a given dev.        
*  NOL                    Last observation defining a particular device.
*  OBJF                   Objective function          
*  OBSTYPE                Observation type                                    
*  TOTALMAX               Maximum residual
*  TOTALMEAN              Mean residual                      
*  TOTALMIN               Minimum residual
*  TOTALSTDDEV            Residual standard deviation
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ALPHA                  Calculates alpha values according to different form.
*  DEFINENAMES            Fills array TYPENAME  
*  ESTIMCOV               Calculates covariance matrix of param. estimates
*  GENSTAT                Calculates general statistics      
*  LIKELIHOOD             Calculates likelihood according to different form.
*  LINV1P                 IMSL routine used to invert hessian matrix    
*  MODELSELEC             Calculates model selection criteria            
*  OUTPUT                 Echoes general output                   
*  RESIDUALS              Calculates residuals (weighted/non weighted)
*  SO_OBJ                 Calculates obj. func. and the contrib. of each device
*  ZERO_ARRAY                            
*  ZERO_ARRAY_I                                                                 
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan. 2003): Revision and formatting
*
********************************************************************************

C______________________________________________ Step 0: Declaration of variables

       IMPLICIT NONE

                                                   ! External variables: scalars

       INTEGER*4 NDEVS,NTYPAR, NUMTOBS, NZPAR, NBANDCOV, NPAR,MAINF
     ;          ,IDIMHESS,NSTAT, IOVAR ,IDIMCOV,IOPRHED,IODIM,ISOT
     ;          ,IDIMWORK,NFLAGS, IOINV ,IDIMWGT,MXGRPZN,NPARDET

                                                    ! External variables: arrays

       INTEGER*4 IODEVICE(NDEVS+1,9), INORPAR(NTYPAR), IVPAR(NZPAR,4)
     ;          ,OBSCLASS(NSTAT,NDEVS+1), IOLG_PAR(NTYPAR)
     ;          ,NZONE_PAR(NTYPAR),ITYPEPAR(NPAR,2) ,IOPT_GS(MXGRPZN,20)
     ;          ,IFLAGS(NFLAGS),MEASTYP(NUMTOBS),IPNT_PAR(NZPAR*IDIMWGT)

       REAL*8 VOBS(NUMTOBS), VOBSC(NUMTOBS),TOBS(NUMTOBS,2)
     ;        ,PARM(NPAR), PARC(NPAR), COVINV(IDIMCOV)
     ;        ,STPAR(NZPAR),PAR_WGT(NTYPAR),WORK(IDIMWORK)
     ;        ,FOBJ_WGT(NSTAT),VJAC(NUMTOBS,NPAR),COVPAR(IDIMHESS)
     ;        ,RESIDPAR(NPAR,2),WGT_UNK(NPAR)

       REAL*8  DEVICESTAT(7,NDEVS),EIGENVEC(NPAR,NPAR)
     ;        ,HESSINV(IDIMHESS), EIGENVAL(NPAR),RESID(NUMTOBS)

       CHARACTER DEVNAME(NDEVS)*10

                                                   ! Internal variables: scalars

       REAL*8 TOTALMEAN, TOTALSTDEV,DETHESS, OBJF, D1, D2, DOF,
     ;        CPERO, FOVERDOF,LIKS1, LIKS2, LIKS3, EXPS1, EXPS2, EXPS3,
     ;        TOTALMIN, TOTALMAX

       INTEGER*4 IER, I,ND,OBSTYPE,MAXNOFOBS,NOF,NOL
     ;          ,NOFOBS, J, IOEIGENVEC, IORESLIST
     ;          ,IOCORMAT, NOFPAR,IODETHESSOK
     ;          ,K

       LOGICAL ALPHANEG

     
       REAL*8 STAT(40,11)
       CHARACTER*10 TYPENAME(40)
       INTEGER*4 INORPAR2(NTYPAR+1)

C_________________________Step 1: Initialises local arrays and get users desires
      
      STAT = 0.0D0
      INORPAR2 = 0D0
      IORESLIST = 0
      IOEIGENVEC = 0
      IOCORMAT = 0
      IF (IOVAR .GE. 1000) THEN          !if the list of all residuals is wanted
         IORESLIST = 1
         IOVAR = IOVAR - 1000
      ENDIF
      IF (IOVAR .GE. 100) THEN                    !if the eigenvectors are wanted 
         IOEIGENVEC = 1
         IOVAR = IOVAR - 100
      ENDIF
      IF (IOVAR .GE. 10) THEN               !if the correlation matrix is wanted
          IOCORMAT = 1
      ENDIF

*________ Options checking
  
      IF (IOINV.GT.0) THEN
        IF (IOEIGENVEC.EQ.1 .AND. IOCORMAT.EQ.0)THEN
           IOCORMAT=1
           WRITE(MAINF,666)
 666       FORMAT(//,' WARNING: EIGENVECTORS ANALYSIS REQUESTED AND '
     ;            'A POSTERIORI COVARIANCE MATRIX IS NOT CALCULATED.'
     ;         ,/,' TRANSIN HAS ACTIVATED THIS OPTION AUTOMATICALLY',/)
        ENDIF
      ELSE   
        IF (IOEIGENVEC.EQ.1) WRITE(MAINF,667)
 667    FORMAT(//,' WARNING: EIGENVEC. ANALYSIS REQUESTED AND INVERSE '
     ;            'PROBLEM OPTION IS NOT ACTIVATED. ',/,
     ;            'EIGENVECTORS WILL NOT BE CALCULATED')

        IF (IOCORMAT.EQ.1) WRITE(MAINF,668)
 668    FORMAT(//,' WARNING: A POSTERIORI CORRELATION MATRIX REQUESTED '
     ;            'AND INVERSE PROBLEM OPTION IS NOT ACTIVATED. ',/,
     ;            'CORRELATION MATRIX WILL NOT BE CALCULATED')

        IOEIGENVEC=0
        IOCORMAT=0
      END IF
                                        !cloning inorpar to work easier with it 
      INORPAR2(1)=0     

      DO I=2,MAX (ISOT,IODIM)+1
          INORPAR2(I)= INORPAR2(I-1)+NZONE_PAR(1)
      ENDDO
      DO I= MAX (ISOT,IODIM)+2,6
         INORPAR2(I)= INORPAR2(I-1)
      ENDDO

      INORPAR2(7)=INORPAR2(6)
      IF (ISOT .EQ. 6) INORPAR2(7)=INORPAR2(7)+NZONE_PAR(1)

      DO I= 8,NTYPAR+1
         INORPAR2(I)=INORPAR2(I-1)+NZONE_PAR(I-6)
      ENDDO
      
     
C______________________ Step 2: Fills array OBSCLASS. Groups dev. by data type.

      MAXNOFOBS = 0

      OBSCLASS = 0D0

      DO ND= 1 , NDEVS                                       ! Loop over devices
              
         NOF = IODEVICE(ND,8)              ! First measurement at current device
         NOL = IODEVICE(ND+1,8) - 1         ! Last measurement at current device
         NOFOBS= NOL - NOF +1            ! Number of measurements at curr. dev.
         STAT(IODEVICE(ND,1),11) = STAT(IODEVICE(ND,1),11)+NOFOBS
         IF (NOFOBS .GT. 0) THEN

                    ! Identifies meas. data type (current row at OBSCLASS array)
             OBSTYPE = IODEVICE(ND,1) 
                                 ! Updates number of devs.of the same meas. type
             OBSCLASS(OBSTYPE,1) = OBSCLASS(OBSTYPE,1) + 1 
                            ! Stores number of current dev. at the app. position
             OBSCLASS(OBSTYPE,OBSCLASS(OBSTYPE,1)+1) = ND  

                                                             ! Updates MAXNOFOBS
             IF (NOFOBS .GT. MAXNOFOBS)  MAXNOFOBS = NOFOBS

          ENDIF                                   !if measurements are available
      ENDDO                                                        ! Next device


                                                          !  define array stats
      
      DO I=1,NSTAT                   !fill stat with lambdas of state variables
         STAT(I,2)=FOBJ_WGT(I)      
      ENDDO

      DO I=1, 6 !fill stat with lambdas of T and anisotropy component
         NOFPAR = 0
         DO J=INORPAR2(I)+1,INORPAR2(I+1)
             IF (IVPAR(J,1) .NE.0)THEN
	          DO K=IVPAR(J,1),IVPAR(J,2)
	             IF (IPNT_PAR(K).EQ.J) THEN
                       STAT(I+NSTAT,2)=PAR_WGT(1)
                       NOFPAR= NOFPAR +1
	             ENDIF
	          ENDDO
             ENDIF
         ENDDO
         STAT(I+NSTAT,11)=NOFPAR
      ENDDO

                                          ! Fill STAT with lambdas of parameters
      DO I=NSTAT+7,NTYPAR+NSTAT       ! -5 Anisot. components are not in PAR_WGT 
         NOFPAR = 0
         NOF = INORPAR(I-NSTAT)+1
         NOL = INORPAR(I-NSTAT)+NZONE_PAR(I-NSTAT-5)

         DO J=NOF,NOL
             IF (IVPAR(J,1) .NE. 0)THEN
	          DO K= IVPAR(J,1),IVPAR(J,2)
	             IF (IPNT_PAR(K).EQ.J) THEN
                          STAT(I+NSTAT,2)=PAR_WGT(I-NSTAT-5)
                          NOFPAR= NOFPAR +1
	             ENDIF
	          ENDDO
             ENDIF
         ENDDO
         STAT(I,11)=NOFPAR 
      ENDDO

                        ! Fill array with names of state variable and parameters

      CALL DEFINENAMES (IOPRHED,TYPENAME)

      DO I=1,NUMTOBS                                
         RESID(I)=VOBSC(I)-VOBS(I)
      ENDDO
 

      !fill array itypepar with parameter types
	
	IF (NPARDET.NE. 0) THEN
           DO K=1,NTYPAR-1
             
             IF (INORPAR2(K)+1 .LE. INORPAR2(K+1) ) THEN

                DO I=INORPAR2(K)+1, INORPAR2(K+1)
                   IF (IVPAR(I,1).NE.0) THEN
                      DO J=IVPAR(I,1),IVPAR(I,2)
                         NOFPAR=IPNT_PAR(J)
                         ITYPEPAR(NOFPAR,1)=K
                         ITYPEPAR(NOFPAR,2)=I-INORPAR2(K)
                      ENDDO
                   ENDIF
                ENDDO

             ENDIF
          ENDDO 
	ENDIF
	IF (NPARDET.NE. NPAR) THEN
	    DO I = 1,MXGRPZN
			IF (IOPT_GS(I,2) .GT.0) THEN
	            DO J= NOFPAR +1,NOFPAR+IOPT_GS(I,5)
				   ITYPEPAR(J,1) = IOPT_GS(I,1)
	               ITYPEPAR(J,2)= J-NPARDET
                  ENDDO
	            NOFPAR =NOFPAR + IOPT_GS(I,5)
	         ENDIF
         ENDDO
	ENDIF
C_______________________ Step 3: Calculate total and partial objective functions

      CALL SO_OBJ 
     ;(IDIMCOV  ,NBANDCOV   ,NDEVS      ,NPAR     ,NSTAT      ,NTYPAR
     ;,NUMTOBS  ,OBJF       ,COVINV     ,COVPAR   ,DEVICESTAT ,IODEVICE
     ;,ITYPEPAR ,OBSCLASS   ,PARC       ,PARM     ,RESID      ,RESIDPAR
     ;,STAT)

 
C__________________________________ Step 4: Calculate general problem statistics

      CALL GENSTAT
     ;(CPERO    ,DOF      ,FOVERDOF ,NPAR     ,NUMTOBS  ,OBJF)

C________________________Step 5: Calculation and inversion of  the hessin matrix

      CALL HESSIANO
     ;(IDIMCOV  ,0          ,MAINF    ,NBANDCOV  ,NFLAGS   ,NPAR     
     ;,NSTAT    ,NUMTOBS    ,COVINV   ,COVPAR    ,FOBJ_WGT ,HESSINV     
     ;,IFLAGS   ,VJAC       ,WGT_UNK  ,MEASTYP)

      CALL EQUAL_ARRAY(WORK,HESSINV,IDIMHESS)   
      
!                        On OUTPUT, who knows what is WORK, but HESSINV is H\-1

      CALL LINV1P(WORK, NPAR,HESSINV ,D1,D2,IER)

      DETHESS = LOG(D1)+D2*LOG(2D0)   !Logarithm of the determinant of the hessian

* ojo2
*      write(6667,*) 'lnH=',dethess
* fin del ojo2

      IODETHESSOK = 0
      IF (IER .EQ. 129) IODETHESSOK = 1       ! check if hessian is not singular

                    ! Hessian matrix is positive definite. Calculates everything

      IF (IODETHESSOK.EQ.0 .AND. IOINV.GT.0) THEN  

C______________________________________________________ Step 6: Calculate alphas

        CALL ALPHA
     ;(IDIMCOV  ,IDIMHESS ,NBANDCOV ,NDEVS    ,NPAR     ,NSTAT
     ;,NTYPAR   ,NUMTOBS  ,OBJF     ,COVINV   ,COVPAR
     ;,HESSINV  ,INORPAR2 ,IODEVICE ,ITYPEPAR ,OBSCLASS
     ;,STAT     ,VJAC     ,WORK)    

        ALPHANEG = .FALSE.

        DO I=1,40
           IF (STAT(I,3).LT. 0.D0 .OR. STAT(I,4) .LT.0.D0) 
     ;         ALPHANEG = .TRUE.
        ENDDO

C______________________________ Step 7: Calculates xpected likelihood criterions

        CALL LIKELIHOOD
     ;(DETHESS  ,EXPS1    ,EXPS2    ,EXPS3    ,IDIMCOV  ,IDIMHESS
     ;,LIKS1    ,LIKS2    ,LIKS3    ,NDEVS    ,NPAR     
     ;,NSTAT    ,NTYPAR   ,NUMTOBS  ,OBJF     ,COVINV   ,COVPAR   
     ;,EIGENVEC ,INORPAR2 ,OBSCLASS ,STAT     ,ALPHANEG)



      ENDIF   ! Hessian matrix is positive definite
      
C_________________________________________ Step 8 Calculates residual statistics

      CALL RESIDUALS
     ;(IDIMCOV  ,IOINV      ,NBANDCOV   ,NDEVS        
     ;,NSTAT    ,NTYPAR     ,NUMTOBS    ,NPAR      ,TOTALMAX
     ;,TOTALMEAN,TOTALMIN   ,TOTALSTDEV ,COVINV    ,COVPAR   ,DEVICESTAT 
     ;,IODEVICE   ,ITYPEPAR ,OBSCLASS ,RESID    ,RESIDPAR
     ;,STAT     ,VOBSC      ,WORK)
C___________________________________________________________ Echoes some results

      CALL OUTPUT
     ;(ALPHANEG ,CPERO       ,DOF        ,EXPS1     ,EXPS2     ,EXPS3
     ;,FOVERDOF ,IDIMCOV     ,IODETHESSOK,IOINV     ,IORESLIST ,LIKS1      
     ;,LIKS2    ,LIKS3       ,MAINF      ,NDEVS     ,NPAR      ,NSTAT    
     ;,NTYPAR   ,NUMTOBS     ,NZPAR      ,OBJF      ,TOTALMAX ,TOTALMEAN
     ;,TOTALMIN ,TOTALSTDEV  ,COVINV     ,DEVICESTAT,DEVNAME   ,INORPAR2 
     ;,IODEVICE ,IOLG_PAR    ,ITYPEPAR   ,IVPAR     ,OBSCLASS  ,PARC         
     ;,RESID    ,RESIDPAR    ,STAT       ,TOBS      ,TYPENAME  ,VOBSC
     ;!NUEVOS
     ;,IDIMWGT  ,IPNT_PAR    ,HESSINV    ,IDIMHESS  ,STPAR)


C_____________________________ Step 9: Calculates paremeter estimates statistics
      
      IF (IOINV .GT.0) THEN   
         CALL ESTIMCOV 
     ;(IDIMHESS   ,IOCORMAT   ,IODETHESSOK  ,IOEIGENVEC  ,ITYPEPAR
     ;,MAINF      ,NPAR       ,NSTAT        ,NTYPAR      ,HESSINV
     ;,WORK       ,EIGENVAL   ,EIGENVEC     ,RESIDPAR    ,STAT
     ;,TYPENAME)

C__________________________________ Step 10: Calculates model selection criteria

         IF (IODETHESSOK.EQ.0) CALL MODELSELEC
     ;      (DETHESS,LIKS1,LIKS2,LIKS3,MAINF,NPAR,NUMTOBS
     ;       ,STAT)
      ENDIF

      RETURN
      END

