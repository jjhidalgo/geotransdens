      SUBROUTINE INIT_GEOEST
     ;(DRELMX_ORIG  ,DABSMX_ORIG  ,ICALL      ,IDIMQ    ,IERROR
     ;,IODIM        ,IOEQT        ,IORTS      ,IOTRS    ,MAINF
     ;,MAXITER_ORIG ,NPARALG      ,NUMEL      ,NUMNP    ,NWRITE
     ;,XMARQ_ORIG   ,ACTH         ,CAUDAL     ,CCAL     ,HCAL         
     ;,IOWRITE      ,IPAR_INV     ,LDIM       ,LTYPE    ,PAR_DIR      
     ;,PAR_INV      ,QXYZ         ,VD         ,XNORVD   ,NFL_SIM      
     ;,NTP_SIM      ,PAREL        ,PARNP      ,NPPEL    ,NPPNP
     ;,IOVAR        ,IOVAR_ORIG)

********************************************************************************
*
* PURPOSE Due to to changes in some scalars related to linear/non linear
*         inverse problem trhough the algorithm, some of them have to be
*         saved at the very beginning of the calibration process (only if more 
*         than one conditional simulation is being updated). Thus, first time,
*         some values are saved. After first simulation, before starting the
*         calib. process, these original values are recovered and changing
*         values are reassigned. Also, reads again initial conditions if that
*         is necessary
*
* DESCRIPTION
*
* EXTERNAL VARIABLES: ARRAYS
*
*  IPAR_INV               Array containing all integer inverse problem          
*                         parameters                                            
*  PAR_DIR                Array containing all real direct problem              
*                         parameters                                            
*  PAR_INV                Array containing all real inverse problem             
*                         parameters                                            
*
* EXTERNAL VARIABLES: SCALARS
*
*  DABSMX_ORIG            Non linear problem convergence criteria
*  DRELMX_ORIG            Non linear problem convergence criteria
*  ICALL                  Index of calling
*  MAXITER_ORIG           Maximum number of iterations, as read in *DIM.DAT
*  NPARALG                Maximum number of algorithm parameters, including     
*                         minimization and simulation (used for dimensioning)   
*  XMARQ_ORIG             Initial value of Marquardt's param. as read origin.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY    AAR (Jan-2003): First coding
*
********************************************************************************
 
      IMPLICIT NONE

C______________________________________________________ Declaration of variables

      INTEGER*4 MAXITER_ORIG,NPARALG,ICALL,IOTRS,IORTS,IERROR,NWRITE
     ;         ,MAINF,NUMNP,NUMEL,IODIM,IDIMQ,IOEQT,NFL_SIM
     ;         ,NTP_SIM,NPPEL,NPPNP,IOVAR,IOVAR_ORIG  
     ;         ,IPAR_INV(NPARALG),LTYPE(NUMEL),LDIM(NUMEL)
     ;         ,IOWRITE(NWRITE)
      REAL*8 XMARQ_ORIG,DRELMX_ORIG,DABSMX_ORIG
     ;      ,PAR_INV(NPARALG),PAR_DIR(NPARALG),CCAL(NUMNP),HCAL(NUMNP)
     ;      ,ACTH(NUMEL),XNORVD(NUMEL),VD(IODIM,NUMEL),CAUDAL(NUMNP)
     ;      ,QXYZ(IDIMQ,NUMEL),PAREL(NUMEL,NPPEL),PARNP(NUMNP,NPPNP)

      INTEGER*4 NROW
      CHARACTER FILENAME(20)*20

C______________________________ Initializes or returns to initial value inverse
C______________________________ problem related variables. Only those that might
C______________________________ change during inversion process.

      IF (ICALL.EQ.1) THEN  ! First entry. Nothing is initialized. Clones values

         MAXITER_ORIG=IPAR_INV(4)
         XMARQ_ORIG=PAR_INV (1)
         DRELMX_ORIG=PAR_DIR (4)
         DABSMX_ORIG=PAR_DIR (5)
         IOVAR_ORIG=IOVAR
          
      ELSE            ! 2nd. entry. Recovers values and reads initial conditions

        NROW=0 
        REWIND(15)         ! Rewinds file with initial conditions
        IF (IOTRS.EQ.1.OR.IORTS.EQ.1) CALL ENDATINICOND 
     ; (IERROR   ,IOWRITE(1),IORTS    ,IOTRS    ,IOWRITE(2) ,15
     ; ,MAINF    ,NFL_SIM   ,NTP_SIM  ,NUMNP    ,CCAL       ,FILENAME 
     ; ,HCAL)

        IF (IOEQT.EQ.2) THEN
          CALL ENTVEL 
     ; (IDIMQ    ,IERROR   ,IOWRITE(1),IODIM   ,IOWRITE(2),15
     ; ,MAINF    ,NROW     ,NUMEL    ,ACTH     ,FILENAME ,LDIM
     ; ,LTYPE    ,QXYZ     ,VD       ,XNORVD)

          CALL ENTFLOW
     ; (IERROR   ,IOWRITE(1),IOWRITE(2),15  ,MAINF    ,NROW
     ; ,NUMNP    ,CAUDAL   ,FILENAME)

        ENDIF

! Recovers original values

        IPAR_INV(4) = MAXITER_ORIG
        PAR_INV (1) = XMARQ_ORIG
        PAR_DIR (4) = DRELMX_ORIG 
        PAR_DIR (5) = DABSMX_ORIG 
        IOVAR = IOVAR_ORIG

C------------------------- Initializes PAREL and PARNP

        CALL ZERO_ARRAY(PAREL,NUMEL*NPPEL)
        CALL ZERO_ARRAY(PARNP,NUMNP*NPPNP)

      END IF

      RETURN
      END

