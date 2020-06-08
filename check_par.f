       SUBROUTINE CHECK_PAR
     ;(ERNUM    ,IERROR  ,IGRP1     ,IOINV   ,IOPBLI   ,IOTIM    
     ;,IOWAR    ,IUPAR   ,IVVAR1    ,MAINF   ,MXGRPZN  ,NFNL    
     ;,NFNLVAR1 ,NFTVAR1 ,NGROUP_ZN ,NPFNL   ,STVAR1   ,VAR      
     ;,FILENAME  ,NROW)

*****************************************************************************
*
* PURPOSE
*
*      Checks zonal data. Estimation indexes, prior info. time func, etc.
*
* DESCRIPTION
*
* - Step 0: Declaration of variables
* - Step 1: Checks time function number. Echoes a warning if problem is 
*           in steady state and time function was defined
* - Step 2: Checks non linear function number. Problem has to be under 
*           transient conditions
* - Step 3: Echoes a warning if direct problem is solved and estim.
*           option or standard deviation are not zero
* - Step 4: Checks group of zones exceeding maximum allowed. If it is zero, it
*           is set to a number of group of zones not defined in GEO file
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  ERNUM                  Starting error number for the current parameter
*  IERROR                 Current number of errors on input data                
*  IGRP1                  Group of zones to which actual zone belongs to (READ)
*  IOINV                  Inverse problem option                                
*  IOPBLI                 If zero, linear  problem, otherwise nonlinear.
*  IOTIM                  Problem regime (0, SS, 1, transient with given
*                         initial cond., 2 transient with SS initial cond.)
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUPAR                  Unit number of file PAR                               
*  IVVAR1                 Read estimation index                                                      
*  MAINF                  Unit number of the main output file (RES.OUT)       
*  MXGRPZN                Maximum number of groups of zones allowed (PARAMETER)  
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NFNLVAR1               Read nonlinear function number
*  NFTVAR1                Read time function number
*  NGROUP_ZN              Number of groups of zones reads in GEO file
*  NPFNL                  Counts the number of nonlinear functions read
*  STVAR1                 Read standard deviation of a parameter                                                      
*  VAR                    String. Name of the current parameter
*
* INTERNAL VARIABLES: SCALARS
*
*  NROW                   Current record number                                 
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*
* HISTORY
*
*     SCR      5-1997     First version
*     AMS      1-1998     Revision. Addition of header
*     AAR      7-2003     Revision
*
*****************************************************************************

C_______________________ Step 0: Declaration of variables

       IMPLICIT NONE
       
                                                             ! Integer external
       INTEGER*4 NFNL,NFTVAR1,NFNLVAR1,IOTIM,IERROR
     ;          ,IOWAR,MAINF,IUPAR,IOPBLI,NPFNL,IOINV,IVVAR1
     ;          ,IGRP1,MXGRPZN,NGROUP_ZN
                                                                ! Real external
       REAL*8 STVAR1
       REAL*4 ERNUM
                                                             ! Integer internal
       INTEGER*4 NROW
                                                                   ! Characters
       CHARACTER*20 FILENAME(20),VAR

C_______________________ Step 1: Checks time function number. Echoes a warning
C_______________________         if problem is in steady state and time funct.
C_______________________         was defined

       IF (IOTIM.EQ.0 .AND. NFTVAR1.NE.0) 
     ;     CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;     ' STEADY STATE PROBLEM AND TIME FUNCTION IS NOT ZERO AT ,//
     ;     ('//VAR//')',NROW,5,IUPAR,0,0.00)

C_______________________ Step 2: Checks non linear function number. Problem has
C_______________________         to be under transient conditions

       IF (IOPBLI.NE.0. AND .NFNLVAR1.NE.0)THEN

          IF (IOTIM.NE.0)THEN
             IF (NFNLVAR1.LE.NFNL. AND. NFNLVAR1.GT.0)THEN
                NPFNL=NPFNL+1
             ELSE
                CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     :           ' NON LINEAR FUNCTION OUT OF RANGE AT '//
     ;           VAR,NROW,6,IUPAR,1,ERNUM+0.01)
             ENDIF                   
          ELSE
             CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;       'STEADY STATE PROBLEM DO NOT ALLOW NON LINEAR CONDITIONS',
     ;       NROW,6,IUPAR,1,ERNUM+0.02)
          ENDIF

       ELSE IF (IOPBLI.EQ.0 .AND. NFNLVAR1.NE.0) THEN
           CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;     ' LINEAR FLOW OR TRANSPORT DEFINED AND NON LINEAR FUNCTION'//
     ;     ' IS NOT ZERO AT ('//VAR//')',NROW,6,IUPAR,0,ERNUM+0.03)
       ENDIF

C_______________________ Step 3: Echoes a warning if direct problem is solved
C_______________________         and estim.option or standard deviation are not
C_______________________         zero. 

       IF (IOINV.LE.0) THEN

           IF (IVVAR1.NE.0) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;     ' DIRECT PROBLEM SOLVED AND ESTIM. OPTION IS NOT ZERO AT ,//
     ;     ('//VAR//')',NROW,3,IUPAR,0,ERNUM+0.04)

           IF (STVAR1.NE.0D0) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;     ' DIRECT PROBLEM SOLVED AND ST. DEVIATION IS NOT ZERO AT ,//
     ;     ('//VAR//')',NROW,3,IUPAR,0,ERNUM+0.05)

       END IF ! IOINV.LE.0


C_______________________ Step 4: Checks group of zones exceeding max. allowed. 

       IF (IGRP1.GT.MXGRPZN) THEN 
          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;     ' NUMBER OF GROUP OF ZONES EXCEEDS MAXIMUM ALLOWED AT '//
     ;     ' AT ('//VAR//')',NROW,6,IUPAR,1,ERNUM+0.06)
          STOP ' CHECK FILE RES.OUT'
       END IF

       IF (IGRP1.EQ.0) IGRP1=NGROUP_ZN+1

       RETURN 
       END
