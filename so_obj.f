      SUBROUTINE SO_OBJ 
     ;(IDIMCOV  ,NBANDCOV   ,NDEVS      ,NPAR     ,NSTAT      ,NTYPAR
     ;,NUMTOBS  ,OBJF       ,COVINV     ,COVPAR   ,DEVICESTAT ,IODEVICE
     ;,ITYPEPAR ,OBSCLASS   ,PARC       ,PARM     ,RESID      ,RESIDPAR
     ;,STAT)

********************************************************************************
*
* PURPOSE
*
* This subroutine calculates the parameters objective function (contrib. for
* each param. type), state variable/s objective function for each observation 
* type,the contribution to it device to device and the total objective function.
*
* REFERENCES:
*
* the equations implemented in this subroutine are taken from
*     A.Medina, J.Carrera, 2001:'Geostatistical inversion of coupled
*     problems: dealing with computational problems and different types
*     of data', Barcelona, Spain
*
*     J.Carrera, S.P.Neumann, 1986:'Estimation of aquifer parameters
*     under transient and steady state conditions: 1', water resources
*     research , vol 22, no 2, pp 199-210
*
*
* THEORY:
*
*  The penalty criterion for a type of state variables is defined as
*           t      -1
*         J   *  V    *   J
*           m      m       m
*
* where Jm is the vector of residuals of state variable type m, and Vm is a
* matrix proportional to the covariance matrix.
*
* The penalty criterion for a tye of parameters p is defined as
*           t      -1
*         J   *  V    *   J
*           p      p       p
*
*
* The contribution of 1 device consist of the product of the vector of
* residuals of this device times V^-1 times the vector of residuals of this
* device.
*
* The objective function is defined as the sum of all penalty criterions
* multplied by their respective lambdas (where lambda is the approximated
* ratio between the variance of a variable type and the variance of heads).
*
* DESCRIPTION:
*
*  The general structure of this subroutine is as follows:
*
*  - Step 1: LOOP over all the types of state variables
*     - Step 1.1: LOOP over the devices that measure this state variable
*       - Step 1.1.1: Calculate the contribution of device to penalty criterion
*       - Step 1.1.2: Add this contribution to the current value of the
*                   penalty criterion
*
*  - Step 2: LOOP over all the estimated parameters
*    - Step 2.1: Calculates vector of residuals, identifies WEIGHT,IZPAR
*                & ITYPEPAR
*    - Step 2.2: Calculate the product COVPAR*RESID
*    - Step 2.3: Add this contribution to the current value of the partial obj.
*                function
*
*  - Step 3: calculate total objective function
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
*                              8: A priori standard deviation associated to
*                                 device. This entry is filled in entdatobs.f
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Status for calc. of obs.                    
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*  IOLG_PAR               Array containing all logarithmic options of           
*                         estimated parameters                                  
*  ITYPEPAR               Array containing type of a given parameter, as set 
*                         array PAR
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  IZPAR                  Similar to ITYPEPAR, but containing component at PARC
*                         PARM or IVPAR
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
*  RESID                  Array containing residuals of observations
*  RESIDPAR               Array containing residuals of parameters and the 
*                         product of them times COVPAR
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
*                              31-40: Make your own parameter!
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
*  WEIGHT                Similar to ITYPEPAR, but containing related weigths for
*                        objective function
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCOV                Used to dimension COVINV
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NDEVS                  Number of observation devices
*  NPAR                   Total number of parameters to be estimated            
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)               
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMTOBS                Total number of observations                          
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  OBJF                   Total objective function                         
*
* INTERNAL VARIABLES: SCALARS
*
*  IPAR                   Dummy counter of estimated parameters 
*  ITYPMEAS               Dummy counter for state var. types
*  ITYPPAR                Dummy counter for parameter types                    
*  ND                     Dummy counter for devices
*  NOF                    First observation defining a particular device.
*  NOFOBS                 Number of observations defining a given dev.
*  NOL                    Last observation defining a particular device.
*  ZONECONTRIB            Contribution to obj. func. of a given zone
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  CONTRIBDEV             Calculates contribution to obj. fun. of a given dev.
*  MUL_SYMMAT_VEC         Multiplies a symm. matrix times a vector
*  POS_WGT_LOG            Identifies ITYPEPAR,IZPAR and WEIGHT. Also calculates 
*                         the residual
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan. 2003): Revision and formatting
*
********************************************************************************

C______________________________________________ Step 0: Declaration of variables

      IMPLICIT NONE

                                                   ! External variables: scalars
      INTEGER*4 NBANDCOV,NUMTOBS,NDEVS,NTYPAR,NSTAT,IDIMCOV,NPAR
     ;        

      REAL*8 OBJF

                                                    ! External variables: arrays
      INTEGER*4 IODEVICE(NDEVS+1,9),OBSCLASS(NSTAT, NDEVS+1)
     ;         ,ITYPEPAR(NPAR,2)
     ;         

      REAL*8  PARM(NPAR), PARC(NPAR), COVINV(IDIMCOV)
     ;       ,DEVICESTAT(7,NDEVS) ,RESID(NUMTOBS) ,STAT(40,11)
     ;       ,COVPAR(NPAR*(NPAR+1)/2),RESIDPAR(NPAR,2)

                                                   ! Internal variables: scalars
      INTEGER*4 NOF,NOL,NOFOBS,ND,ITYPMEAS,ITYPPAR,IPAR

      REAL*8 CONTRIBDEV,ZONECONTRIB,XLAMBDA

C____________________________ Step 1: LOOP over all the types of state variables

       DO ITYPMEAS= 1,NSTAT
         STAT(ITYPMEAS,1) = 0.D0      ! Init. contrib. to obj. func of this type

         IF (OBSCLASS(ITYPMEAS,1) .GT. 0) THEN          ! If meas. are available

C_____________ Step 1.1: LOOP over the devices that measure this state variable

            DO ND = 2, OBSCLASS(ITYPMEAS,1) + 1

               NOF = IODEVICE(OBSCLASS(ITYPMEAS,ND),8) !first measurement of device
               NOL = IODEVICE(OBSCLASS(ITYPMEAS,ND)+1,8) - 1 !last MEASUREMENT OF DEVICE
               NOFOBS= NOL - NOF + 1

C_____________________________ Calculates non-weighted residual

              

C_____ Step 1.1.1: Calculate the contrib. of device to penalty criterion 

               DEVICESTAT(7,OBSCLASS(ITYPMEAS,ND)) = CONTRIBDEV
     ;(IDIMCOV  ,NBANDCOV ,NOF      ,NOFOBS   ,NUMTOBS
     ;,COVINV   ,RESID)
                      
                
C_________________________________ Step 1.1.2: Add contribution to obj. function

               STAT(ITYPMEAS,1) = STAT(ITYPMEAS,1)+
     ;                            DEVICESTAT(7,OBSCLASS(ITYPMEAS,ND))

            ENDDO ! ND = 2, OBSCLASS(ITYPMEAS,1) + 1
         ENDIF  ! OBSCLASS(ITYPMEAS,1) .GT. 0
        ENDDO   ! ITYPMEAS= 1,10

********************************************************************************
*                                                                              *
***************** PARAMETER RESIDUALS (i.e. when i .ge. 10)*********************
*                                                                              *
********************************************************************************

C___________________________ Step 2: Calculates objective function of parameters

      CALL ZERO_ARRAY(RESIDPAR,2*NPAR)



C___ Step 2.1: For each parameter to be estimated calculates residuals , related 
C___           weight, component at PARC and identifies parameter type

      DO IPAR=1,NPAR                      ! Loop over parameters being estimated
            RESIDPAR(IPAR,1)=PARC(IPAR)-PARM(IPAR)
	ENDDO
            
C__________________________________ Step 2.2: Calculates product of COVPAR*RESID

      CALL MUL_SYMMAT_VEC (NPAR,COVPAR,RESIDPAR(1,2),RESIDPAR(1,1))



C_____________________ Step 2.3: Calculates final product and updates obj. func.

      DO IPAR=1,NPAR
	  ZONECONTRIB=RESIDPAR(IPAR,1)*RESIDPAR(IPAR,2)
	 
        STAT(ITYPEPAR(IPAR,1)+NSTAT,1)=
     ;       STAT(ITYPEPAR(IPAR,1)+NSTAT,1)+ ZONECONTRIB       
        XLAMBDA = STAT(ITYPEPAR(IPAR,1)+NSTAT,2)
	  RESIDPAR(IPAR,2)= XLAMBDA*RESIDPAR(IPAR,2)
  
      END DO

C________________________________________ Step  3: calculate total obj. function

      OBJF = 0.D0

      DO ITYPMEAS= 1,NSTAT                      !loop over all state variables
          OBJF = OBJF + STAT(ITYPMEAS,2) * STAT(ITYPMEAS,1)
      ENDDO             ! objF = objF + lambda(i) * [partial objective function]

      DO ITYPPAR= 1,NTYPAR                          !loop over all parameters
        IF (STAT(ITYPPAR+NSTAT,2) .NE. 0D0 .AND. 
     ;  STAT(ITYPPAR+NSTAT,1).NE. 0D0)
     ; OBJF = OBJF + STAT(ITYPPAR+NSTAT,2) * STAT(ITYPPAR+NSTAT,1)
      ENDDO              !objF = objF + lambda(i) * [partial objective function]

      RETURN
      END
