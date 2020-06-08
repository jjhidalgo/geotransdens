**************************************************************************
****  ONLY MACHINE DEPENDENT ROUTINES ARE INCLUDED IN THIS FILE
**************************************************************************

       SUBROUTINE FECHA (MAINF,TITEL)

******************************************************************
***    IT WRITES MODEL'S DATE AND TITLE
******************************************************************

       CHARACTER*24 DATIM,FDATE,TITEL*80

       DATIM=FDATE()
       WRITE(MAINF,10) TITEL,DATIM(1:3)//'  '//DATIM(9:10)//
     ;         '-'//DATIM(5:7)//'-'//DATIM(21:24)//'  '//DATIM(12:19)
 10    FORMAT(///////,
     .        15X,'M O D E L  N A M E',/,
     .        '+',14X,18('_'),//,1X,A80,///,
     .        15X,'D A T E ........: ',A25)

       RETURN
       END

************************************************************************
************************************************************************

       SUBROUTINE CPUTIME (MAINF)

******************************************************************
*** IT WRITES THE TOTAL CPUTIME SINCE THE BEGINNING OF THE PROGRAM
******************************************************************

       REAL*4 DTIME,TARRAY(2),WW,XX

       WW=DTIME(TARRAY)
       XX=TARRAY(1)
       IHORAS=INT(XX/3600)
       XX=XX-IHORAS*3600.
       IMINUT=INT(XX/60)
       XX=XX-IMINUT*60.
       ISEGU=INT(XX+0.5)
       IF (IHORAS.NE.0) THEN
          WRITE(MAINF,1) IHORAS,IMINUT,ISEGU
       ELSE
          WRITE(MAINF,2) IMINUT,ISEGU
       END IF

 1     FORMAT(//,2X,'TOTAL CPUTIME :',I3,' HOURS',I3,
     ;        ' MINUTES',I3,' SECONDS')
 2     FORMAT(//,2X,'TOTAL CPUTIME :',I3,' MINUTES',I3,' SECONDS')

       RETURN
       END

************************************************************************
************************************************************************

       SUBROUTINE INIT_CPUTIME(I)

       RETURN
       END

************************************************************************
       SUBROUTINE OPENFILES
     ; (IOBALC   ,IOBALGC  ,IOBALGH  ,IOBALH   ,IOCMC    ,IOCMH
     ; ,IOEQT    ,IOFLLI   ,IOINV    ,IOINV_GS ,IOMHC    ,IOMHH    
     ; ,IOPINITC ,IOPINITH ,IOPLC    ,IOPLH    ,IORTS    ,IOSEC    
     ; ,IOSEH    ,IOTRLI   ,IOTRS    ,IOWIT    ,IOWPI    ,IPROCESS 
     ; ,MAINF    ,NBLCVP   ,NUMNP    ,FILENAME)


*****************************************************************************
* PURPOSE 
*     Opens all input, output and auxiliary files 
*
* DESCRIPTION                                    
*     
*     This subroutine opens the following files in the same order as they are
*     refered below.
*     
*     Input data files:
*
*         FILENAME(2): Grid definition. 
*         FILENAME(3): Parameters definition
*         FILENAME(4): Time steps and time functions
*         FILENAME(5): Measured data
*         FILENAME(6): Initial conditions
*         FILENAME(7): Geostat. inv. prob. information
*         FILENAME(8): A priori covariance matrix of parameters to be estim.
*
*      Output files:
*                   
*         FILENAME(10):  PLT file. Temporal variations of computed and measured 
*                                 heads and concentrations at all obs. times
*         FILENAME(11): MSH file. All information relative to the finish 
*                                 element mesh       
*         FILENAME(13): MHH file. Nodal heads at every IOMSH observation time,
*                                 including first and last.
*         FILENAME(14): MSC file. Similar to MHH for concentrations
*         FILENAME(16): MCC       Idem to MHH for concentrations
*                                                               
*         FILENAME(17): SEC file. Same format as PLT FILE, contains heads and 
*                                 concentrations along sectiond of the model
*         FILENAME(18): CVM file. Same format as PLT FILE, 
*                                     X axis:head or concentrations
*                                     Y axis:measured values 
*         FILENAME(19): BLH file. Flow mass balance 
*         FILENAME(20): BLC file. Transport mass balance
*
*      Auxiliary files (SCRATCH):
*
*      UNITS 
*
*         50              Direct access. Heads at computed times in the 
*                         current/previous optimization iteration
*         51              Direct access. Heads at computed times in the 
*                         current/previous optimization iteration. Units 50 
*                         and 51 interchange their roles.
*         52              Seq. access. Computed times in the current/previous
*                         optimization iteration
*         53              Seq. access. Computed times in the current/previous 
*                         optimization iteration. Units 52 and 53
*                         interchange their roles
*         54              Same as 50 for concentrations. 
*         55              Same as 51 for concentrations.
*         56              Same as 52 for concentrations.
*         57              Same as 53 for concentrations.
*         70              Estimated parameters through the optimization 
*                         process
*         71              History of the optimization iterations.
*         81              State variables computed at computation times and 
*                         observation devices
*         93              Head at times defined by variables TIME and IOMHH
*         94              Concentration at times defined by variables TIME 
*                         and IOMHC
*       
* EXTERNAL VARIABLES: ARRAYS 
*
*       FILENAME          Array containing all the input and output filenames
*
* EXTERNAL VARIABLES: SCALARS
*
*  IOBALC                 Zonal transport mass balance at observation times     
*  IOBALGC                Global transport  mass balance (integrated in time)   
*  IOBALGH                Global flow mass balance (integrated in time)         
*  IOBALH                 Zonal flow mass balance at observation times          
*  IOCMC                  If 1, computed vs. measured values of concent. at all 
*                         observation points (steady state) are written in      
*                         CVM file                                              
*  IOCMH                  If 1, computed vs. measured values of heads at all    
*                         observation points (steady state) are written in      
*                         CVM file                                              
*  IOEQT                  Type of problem to be solved                          
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1
*  IOINV                  Inverse problem option                                
*  IOMHC                  If non zero, computed values of concent. at all nodes 
*                         are writen in MCC file every IOMCC observation times  
*  IOMHH                  If non zero, computed values of heads at all nodes    
*                         are writen in MHH file every IOMHH observation times  
*  IOPINITC               Option for the extrapolation of concentrations
*                         at the next time step in the Newton process.
*  IOPINITH               Option for the extrapolation of heads or pressures
*                         at the next time step in the Newton process.
*  IOPLC                  Controls when computed and/or measured concent. are   
*                         written in PLT file                                   
*  IOPLH                  Controls when computed and/or measured heads are      
*                         written in PLT file                                   
*  IORTS                  Transport regime                                      
*  IOSEC                  If non zero, time evolution of computed concent.      
*                         at IOSEC 1-D sections is written in SEC file          
*  IOSEH                  If non zero, time evolution of computed heads         
*                         at IOSEH 1-D sections is written in SEC file          
*  IOTRLI                 If zero, linear transport problem, otherwise set to 1
*  IOTRS                  Flow regime                                           
*  IOWIT                  If 1, estimated parameters history through the
*                         optimization problem are written
*  IPROCESS               Controls if it is the first time that subroutine      
*                         LECDIM is called.                                     
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NUMNP                  Number of nodes

* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  OPEN_ONE_FILE          Opens one specified file, and stops program if
*                         it does not exist. 
*
* HISTORY
*
*     AMS        1988     First coding
*     SCR      4-1997     Revision and verification
*     AMS      1-1998     Common elimination. One GOTO sentence is also elim. 
*     AMS      5-1998     Revision of unit numbers. 
*
*****************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER FILENAME(20)*20

C------------------------- FIRST EXECUTABLE STATEMENT.
C------------------------- The first time input data file have to be open

       IF (IPROCESS.EQ.0) THEN

C------------------------- OPENS INPUT FILES

C------------------------- Opens GRI file

          CALL OPEN_ONE_FILE(FILENAME(2),MAINF,11)

C------------------------- Opens PAR file

          CALL OPEN_ONE_FILE(FILENAME(3),MAINF,12)

C------------------------- Opens TIM file when transient conditions apply

          IF (IOTRS+IORTS.NE.0) 
     ;       CALL OPEN_ONE_FILE(FILENAME(4),MAINF,13)
          
C------------------------- Opens OBS file (measurements and observation points 
C------------------------- coordinates)

          IF (IOINV.NE.0.OR.IOPLH.NE.0.OR.IOPLC.NE.0
     ;                  .OR.IOCMH.NE.0.OR.IOCMC.NE.0) 
     ;       CALL OPEN_ONE_FILE(FILENAME(5),MAINF,14)
          
C------------------------- Opens INI file (initial conditions and flow and 
C------------------------- velocities if necessary)

          IF ((IOTRS.EQ.1 .OR. IORTS.EQ.1) .OR. (IOEQT .EQ. 2) ) 
     ;       CALL OPEN_ONE_FILE(FILENAME(6),MAINF,15)

C------------------------- Opens GEO file (geost. inv. prob. information)

          IF (IOINV_GS.NE.0) CALL OPEN_ONE_FILE (FILENAME(7),MAINF,16)

C------------------------- Opens COV file (A PRIORI COV. MATRIX OF PARAMETERS)

          NBLCVP=0
          IF (NBLCVP.GT.0 .AND. IOINV.GT.0) 
     ;       CALL OPEN_ONE_FILE (FILENAME(8),MAINF,17)

C------------------------- OUTPUT FILES


C------------------------- PLT file (computed vs. measured heads and/or 
C------------------------- concentrations at observation points)

          IF (IOPLH+IOPLC .NE. 0)
     ;       OPEN(26,FILE=FILENAME(10),STATUS='UNKNOWN')
 
C------------------------- Opens MHH and MSH files (mesh and nodal values of 
C------------------------- heads)

          IF (IOMHH .NE. 0) THEN
             OPEN(UNIT=27,FILE=FILENAME(11),STATUS='UNKNOWN')
             OPEN(UNIT=29,FILE=FILENAME(13),STATUS='UNKNOWN')
          ENDIF

C------------------------- Opens MSC and MCC files (mesh and nodal values of 
C------------------------- concentrations))

          IF (IOMHC .NE. 0 ) THEN
             OPEN(UNIT=30,FILE=FILENAME(14),STATUS='UNKNOWN')
             OPEN(UNIT=32,FILE=FILENAME(16),STATUS='UNKNOWN')
          ENDIF

C------------------------- Opens SEC file (1-D sections of computed heads 
C------------------------- and/or concentrations)
      
          IF (IOSEH+IOSEC .GT. 0) 
     ;       OPEN(UNIT=33,FILE=FILENAME(17),STATUS='UNKNOWN')  

C------------------------- Opens CVM file (stationary measured vs. computed 
C------------------------- heads and/or concentrations at all observ. points)

          IF (IOCMH+IOCMC .GT. 0) 
     ;       OPEN(UNIT=34,FILE=FILENAME(18),STATUS='UNKNOWN')  

C------------------------- Opens BLH and BLC files (flow/transport mass balance)

          IF (IOBALH.GT.0 .OR. IOBALGH.GT.0)
     &       OPEN(UNIT=35,FILE=FILENAME(19),STATUS='UNKNOWN')

          IF (IOBALC.GT.0 .OR. IOBALGC.GT.0)
     &       OPEN(UNIT=36,FILE=FILENAME(20),STATUS='UNKNOWN') 

       ENDIF ! IPROCESS.EQ.0

C------------------------- OPENS AUXILIARY STORAGE FILES (SCRATCH)
 
       IF (IOEQT.NE.2) THEN

C------------------------- Opens file to store heads at all nodal points and 
C------------------------- at all computed times (direct access)

          OPEN(UNIT=50,STATUS='SCRATCH',FORM='UNFORMATTED',
     ;                 ACCESS='DIRECT',RECL=8*NUMNP)

C------------------------- File to store computation times

          OPEN(UNIT=52,STATUS='SCRATCH',FORM='UNFORMATTED')

C------------------------- If flow is nonlinear and heads at each time step 
C------------------------- are initialized with the heads at previous 
C------------------------- optimization iteration, two additional files are 
C------------------------- needed (heads and computation times)

          IF (IOFLLI.NE.0 .AND. IOPINITH.NE.0) THEN
             OPEN(UNIT=51,STATUS='SCRATCH',FORM='UNFORMATTED',
     ;                 ACCESS='DIRECT',RECL=8*NUMNP)

             OPEN(UNIT=53,STATUS='SCRATCH',FORM='UNFORMATTED')
          ENDIF

          IF (IOMHH.NE.0) OPEN (94,STATUS='SCRATCH',FORM='UNFORMATTED')
       END IF         ! IOEQT.NE.2

       IF (IOEQT.NE.1) THEN

C------------------------- Opens file to store concentrations at all nodal 
C------------------------- points and at all computed times

          OPEN(UNIT=54,STATUS='SCRATCH',FORM='UNFORMATTED',
     ;                 ACCESS='DIRECT',RECL=8*NUMNP)

C------------------------- File to store computation times

          OPEN(UNIT=56,STATUS='SCRATCH',FORM='UNFORMATTED')

C------------------------- If transport is nonlinear and conc. at each time 
C------------------------- step are initialized with the conc. at previous 
C------------------------- optimization iteration, two additional files are 
C------------------------- needed (conc. and computation times)

          IF (IOTRLI.NE.0 .AND. IOPINITC.NE.0) THEN
             OPEN(UNIT=55,STATUS='SCRATCH',FORM='UNFORMATTED',
     ;                 ACCESS='DIRECT',RECL=8*NUMNP)

             OPEN(UNIT=57,STATUS='SCRATCH',FORM='UNFORMATTED')
          ENDIF

          IF (IOMHC.NE.0) OPEN (93,STATUS='SCRATCH',FORM='UNFORMATTED')
       END IF    ! IOEQT.NE.1

C------------------------- Opens file to store heads and/or concentrations 
C------------------------- at all observation devices and at all computed 
C------------------------- times (used for file PLT)

       IF (IOPLH.GT.1 .OR. IOPLC.GT.1) 
     ;                OPEN(81,STATUS='SCRATCH',FORM='UNFORMATTED')


C------------------------- Opens files to write estimated parameters through
C------------------------- the optimization process and inverse problem values

       IF (IOWIT.NE.0) OPEN(70,STATUS='SCRATCH',FORM='UNFORMATTED')
       IF (IOWPI.NE.0) OPEN(71,STATUS='SCRATCH',FORM='UNFORMATTED')
       RETURN
       END

