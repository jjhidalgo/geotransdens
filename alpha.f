      SUBROUTINE ALPHA
     ;(IDIMCOV  ,IDIMHESS ,NBANDCOV ,NDEVS    ,NPAR     ,NSTAT
     ;,NTYPAR   ,NUMTOBS  ,OBJF     ,COVINV   ,COVPAR
     ;,HESSINV  ,INORPAR2 ,IODEVICE ,ITYPEPAR ,OBSCLASS
     ;,STAT     ,VJAC     ,WORK)    

********************************************************************************
*
* PURPOSE
*
* This subroutine calculates the values of alpha, where alpha is a scalar 
* multiplying the inverse (proportional) covariance matrix. Two different sets 
* of alphas are calculated, according to the 2 formulations as described in the 
* reference
* It also calculates the tau1 that is needed to claculate the log-likelihood
* function
*
* REFERENCE
* The equations implemented in this subroutine are taken from
*     A.Medina, J.Carrera, 2001:'Geostatistical inversion of coupled
*     problems: dealing with computational problems and different types
*     of data', Barcelona, Spain
*
* THEORY:
*
* In the ln(expected likelihood) equation appear a set of scalars, each one 
* related to a type of state variable or a prior information type.  These 
* scalars are called alpha�s and are a "measurement" of the absolute degree of 
* uncertainty of measurements of the type they belong to, contrarily to lambda's
* which are a "measurement" of the relative certainty of different meas. types
* Intepretations of the alpha values in the S equation:
*  1) alphas are chosen such that S is minimized
*  2) alphas are chosen consistent with lambdas.
* The alpha-parameters related to prior information which are called 'touw'
* in the reference are here also called alpha
*
* DESCRIPTION
* - Step 1:Calculation of alphas accordingly to form. 1
*   - Step 1.1: calculate alpha 1 (head levels)
*   - Step 1.2: calculate alpha 2 to 10  (related to the different measurement 
*               types; 2 means concentrations)
*   - Step 1.3: calculate alpha 11 to 30 (related to the different parameter 
*               types)
* - Step 2:Calculation of alphas accordingly to form. 2
*   - Step 2.1: calculate alpha 1
*   - Step 2.2: calculate alpha 2 to 30  (related to the different measurement 
*               and parameterstypes)
* FOOTNOTE: Alpha related to head level is calculated separately, as the rest of
*           alpha parameters are calculated as a function of it. 
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix                      
*  HESSINV                Array containing inverse of hessian matrix
*  INORPAR2               Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Flag for inclusion of data
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*                         Column 9: Related flow/tpt. problem
*  ITYPEPAR               For all estimated parameters, ITYPEPAR contains the 
*                         param. type they belong to
*  IVPAR                  Vector containing estimation index for all            
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
*                         A complete description can be found at STAT_OUTP
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
*                         all parameters prioo information                      
*  VJAC                   Jacobian matrix                                       
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCOV                Used to dimension array COVINV
*  IDIMHESS               Used to dimension array HESS. It is equal to          
*                         NPAR*(NPAR+1)/2                                       
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NDEVS                  Number of devices            
*  NPAR                   Total number of parameters to be estimated            
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)               
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMTOBS                Total number of observations                          
*  OBJF                   Total objective function
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Dummy counter                                    
*  J                      Dummy counter                                         
*  NOF                    First observation defining a particular device.
*  NOFOBSTYPE             Total number of meas. of state var. of a given type
*  NOL                    Last observation defining a particular device.
*  NOFZON                 Total number of zones of a given param. type   
*  TRACE                  Trace of H*J\t*V\-1*J or H*V          
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  TRACEHJVJ              Calcs. trace of H*J\t*V\-1*J                          
*  ZERO_ARRAY                                                                   
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan. 2003): Revision and formatting
*
********************************************************************************

C____________________________________________  Step 0: Declaration of variables 

      IMPLICIT NONE
                                                   ! External variables: scalars

      INTEGER*4 NUMTOBS,NPAR, NDEVS,NTYPAR,IDIMHESS,NSTAT
     ;         ,IDIMCOV,NBANDCOV
      REAL*8 OBJF
                                                   ! External variables: arrays

      REAL*8 STAT(40,11),HESSINV(IDIMHESS),VJAC(NUMTOBS,NPAR)
     ;      ,COVINV(IDIMCOV),WORK(IDIMHESS),COVPAR(IDIMHESS)

      INTEGER*4 INORPAR2(NTYPAR+1),ITYPEPAR(NPAR,2) 
     ;         ,OBSCLASS(NSTAT, NDEVS+1),IODEVICE(NDEVS+1,9)

                                                   ! Internal variables: scalars

      INTEGER*4 I, NOFOBSTYPE, NOF, NOL,OBSTYPE,IPAR,K
     ;        ,POS1,POS2
      REAL*8  TRACE, TRACEHJVJ, NOFZON
     ;        



C__________________________________________________  Step 1.1: Calculates alpha1

                     ! Finds the first state variable type that has observations
      OBSTYPE = 1
      DO WHILE (OBSCLASS(OBSTYPE,1) .EQ. 0)
         OBSTYPE = OBSTYPE + 1
      ENDDO

      NOFOBSTYPE = INT(STAT(OBSTYPE,11))!the number of observations of this type

  
      TRACE=TRACEHJVJ     
     ;(IDIMCOV  ,IDIMHESS ,NBANDCOV ,NDEVS    ,NPAR     ,NSTAT
     ;,NUMTOBS  ,OBSTYPE ,COVINV   ,HESSINV  ,IODEVICE ,OBSCLASS
     ;,VJAC     ,WORK)

      STAT(OBSTYPE,3) = STAT(OBSTYPE,1)/(NOFOBSTYPE-TRACE)   ! AlPHA1 by form. 1




C_________________ Step 1.2: Calculate alpha concerning rest ov state var. types

      DO I = OBSTYPE+1,NSTAT                              ! Loop over var. types
         
         IF (OBSCLASS(I,1) .GT. 0) THEN       ! Measur. available for this class

           NOFOBSTYPE=INT(STAT(I,11))

           TRACE=TRACEHJVJ
     ;(IDIMCOV  ,IDIMHESS ,NBANDCOV ,NDEVS    ,NPAR     ,NSTAT
     ;,NUMTOBS  , I       ,COVINV   ,HESSINV  ,IODEVICE ,OBSCLASS
     ;,VJAC     ,WORK)

           STAT(I,3) = (STAT(I,1)+STAT(OBSTYPE,3)*TRACE)
     ;                                 /REAL(NOFOBSTYPE)            
         ENDIF                                            ! OBSCLASS(I,1) .GT. 0
      ENDDO                                                  ! I=OBSTYPE+1,NSTAT

C______________________________________ Step 1.3 Calculate alphas for parameters

      DO I = 1,NTYPAR                                   ! Loop over param. types
         TRACE=0.D0
         
         IF (STAT(I+NSTAT,2) .NE. 0) THEN              ! Param. type estimated ?

            DO IPAR=1,NPAR
               IF (I.EQ.ITYPEPAR(IPAR,1)) THEN
                  DO K=1,NPAR
                     
                           !getting array indexes of elements ipar,k and k, ipar
                      IF (I .GT. K) THEN   
                         POS1=IPAR*(IPAR+1)/2 -(IPAR-K)  
                         POS2=K*(K+1)/2 -(K-IPAR)
                      ELSE
                        POS1 = K*(K+1)/2 -(K-IPAR)
                        POS2=IPAR*(IPAR+1)/2 -(IPAR-K)
                      ENDIF
                    
                      
                      TRACE=TRACE+HESSINV(POS1)*COVPAR(POS2)

                  END DO
               END IF
            END DO


           NOFZON = INT(STAT(I+NSTAT,11))
           
           STAT(I+NSTAT,3) = (STAT(I+NSTAT,1) 
     ;            + STAT(OBSTYPE,3)*TRACE)/NOFZON
        ENDIF                                       !if parameter type estimated
      ENDDO   ! Next param. type


C_______________________________________________ALPHA ACCORDING TO FORMULATION 2
C____________________________________ Step 2.1 Calculate alpha for heads form. 2

      STAT(OBSTYPE,4) = OBJF/(NUMTOBS+NPAR)

C_____________________________Step 2.2 Calculate alphas for other param.and vars.

      DO I=OBSTYPE+1,NSTAT                        ! Rest of state variable types
         IF (OBSCLASS(I,1) .GT. 0) THEN
            STAT(I,4)=STAT(OBSTYPE,4)/STAT(I,2)
         ENDIF
      ENDDO
 
 
      DO I=1, NTYPAR                                ! Loop over parameter types
          NOF = INORPAR2(I)+1
          NOL= INORPAR2(I+1)
          IF ( NOF .LE. NOL .AND.  STAT(I+NSTAT,2) .NE. 0.D0)
     ;      STAT(I+NSTAT,4)=STAT(OBSTYPE,4)/STAT(I+NSTAT,2)
      ENDDO

      RETURN
      END
