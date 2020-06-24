        SUBROUTINE LECDIM 
     ;(IDIMFNT  ,IERROR   ,IOBALC   ,IOBALDC  ,IOBALDH  ,IOBALGC
     ;,IOBALGH  ,IOBALH   ,IOCNSF   ,IOCNST   ,IOCRITRAP
     ;,IODIM    ,IOEQT    ,IOFLLI   ,IOFLSAT  ,IOFMLF
     ;,IOFMLT   ,IOINV    ,IOPART   ,IOPINITC 
     ;,IOPINITH           ,IORDCH   ,IORTS    ,IOSMFL   ,IOSMTP   
     ;,IOSUCHUM ,IOTRLI   ,IOTRS    ,IOVAR    ,IPROCESS 
     ;,ISOT     ,LMXNDL   ,MAINF
     ;,NBAND    ,NBANDCOV ,NBLCVP   
     ;,NDEVS
     ;,NFLAGS   ,NFNL     ,NFNT     ,NINT     ,NMAXF
     &,NMAXT   
     ;,NOPTS    ,NPAR     ,NPARALG  ,NPARF    ,NPARFPRG 
     ;,NPARPRG  ,NPBFL    ,NPBTP    ,NSTAT    ,NTDMT    
     ;,NTYPAR   ,NUMEL    ,NUMNP    ,NUMTIT   ,NUMTNOD  ,NUMTOBS  
     ;,NWRITE   ,FOBJ_WGT ,IFLAGS   ,IOLG_PAR ,IOPTS    
     ;,IOWRITE  ,IPAR_DIR ,IPAR_INV ,NZONE_PAR,PAR_DIR  
     ;,PAR_INV  ,PAR_WGT  ,FILENAME
     ;!NUEVOS 
     ;,IODENS_INI      ,IODIRECT    ,IOSPARSE  ,ITPTVAR  ,LINMET
     &,IOCONSRC
     ;,IOINV_GS ,MXGRPZN   ,MXLINCMB ,NGROUP_ZN ,IO_KG_GS ,IOPT_GS)

********************************************************************************
*
* PURPOSE Reads dimensions, options and scalar param. of the current problem
*
* DESCRIPTION 
*
* This subroutine reads the dimension and options information from DIM file,
* and writes it as read to check it.
*
* It proceeds as follows:
*
*     - Opens DIM FILE and reads TITLE and input filenames, if these are not
*       specified in interactive mode. The program will be stopped if any of
*       these input files are not found.
*
*     - Reads the groups cards A3 to A9.
*
*         * First, it looks for the current card label on DIM FILE through
*           SRC_NCARD subroutine.
*         * Reads the current card trough the character function LEEL,which
*           returns the LEAUX string
*         * Load the current card variables by reading LAUX string
*         * Checks the read data and writes error messages on MAIN FILE
*           calling subroutine ERROR
*         * Write in MAIN FILE interpreted data, after each card is
*           read.
*
*     - Calls subroutine OPENFILES to open the rest of input, output and auxili`
*       files
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  FOBJ_WGT               Array containing all objective function weights for   
*                         state variables (heads, concentrations, etc)          
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  IOLG_PAR               Array containing all logarithmic options of           
*                         estimated parameters                                  
*  IOPTS                  Array containing different problem options
*  IOPT_GS                General options for inverse problem. Each row contains 
*                         information of a given group of zones
*  IO_KG_GS               Kriging dimensions and options for geost. inv. prob. 
*                         Each row contains information of a given group of zones 
*                         Only sense if group is estimated geost. On each row:
*  IOWRITE                Array containing all output options                   
*  IPAR_DIR               Array containing all integer direct problem           
*                         parameters                                            
*  IPAR_INV               Array containing all integer inverse problem          
*                         parameters                                            
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PAR_DIR                Array containing all real direct problem              
*                         parameters                                            
*  PAR_INV                Array containing all real inverse problem             
*                         parameters                                            
*  PAR_WGT                Array containing objective function weights for       
*                         all estimated parameters                              
*
* INTERNAL VARIABLES: ARRAYS
*
*  FILETYPE               Contains the filenames of all input and output files  
*  NCARD                  Vector. Contains the label cards                      
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMFNT                First dimension of array FNT, it coincides with       
*  IERROR                 Current number of errors on input data                
*  IOBALC                 Zonal transport mass balance at observation times     
*  IOBALDC                Detailed (nodal) flow mass balance   
*  IOBALDH                Detailed (nodal) transport mass balance
*  IOBALGC                Global transport  mass balance (integrated in time)   
*  IOBALGH                Global flow mass balance (integrated in time)         
*  IOBALH                 Zonal flow mass balance at observation times          
*  IOCNSF                 Scheme for storage term in flow problem               
*  IOCNST                 Scheme for mass storage term in transport problem     
*  IOCRITRAP              Option of treatement on the direct problem            
*                         convergence criteria                                  
*  IODIM                  Maximum dimension of the problem                      
*  IOEQT                  Type of problem to be solved                          
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOFLSAT                Indicates the possibility that one part of the domain 
*                         reaches unsaturated state.                            
*  IOFMLF                 Flow Formulation number                               
*  IOFMLT                 Transport formulation number                          
*  IOINV_GS               Geostatistical inverse problem option (AAR PhD thesis)
*                         0- Conditional estimation (zonation approach, as usual)
*                         1- Some of zonal parameters are estimated geost. or
*                            are a linear combination of unknown parameters
*  IOINV                  Inverse problem option                                
*  IOPART                 If 1, the list of problem partitions is written       
*  IOPINITC               Option for the extrapolation of concentrations        
*                         at the next time step in the Newton process.          
*  IOPINITH               Option for the extrapolation of heads or pressures    
*                         at the next time step in the Newton process.          
*  IORDCH                 Option for chain reactions
*  IORTS                  Transport regime                                      
*  IOSMFL                 Option for simultaneous flow problems
*  IOSMTP                 Option for simultaneous tpt. problems    
*  IOSUCHUM               Indicates if measures are given in terms of pressure  
*                         or piezometric heads (set to 0), or in terms of       
*                         saturation degree                                     
*  IOTRLI                 Idem to IOFLLI, in the case of transport.             
*  IOTRS                  Flow regime                                           
*  IPROCESS               Controls if it is the first time that subroutine      
*                         LECDIM is called.                                     
*  ISOT                   Maximum hydraulic conductivity tensor anisotropy      
*                         degree in the problem.                                
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MAXNUMDIVFL            Maximum number of consecutives divergences for flow.
*  MAXNUMDIVTR            Maximum number of consecutives divergences for transport.
*  MXGRPZN                Maximum number of groups of zones. FIXED PARAMETER 
*                         used for dimension
*  MXLINCMB               Maximum number of unknowns defining a linear comb. of
*                         parameters
*  NBAND                  Half Bandwith (maximum difference between the         
*                         numbers of two nodes belonging to the same element)   
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NBLCVP                 Number of structures of a priori covariance matrix of 
*                         parameters to be read from COV file (deterministically
*                         estimated parameters). 
*  NDEVS                  Number of observation devices                       
*  NFLAGS                 Maximum number of allowed flags                       
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NFNT                   Number of time functions used for describing time     
*                         dependence of all transient parameters                
*  NGROUP_ZN              Number of groups of zones
*  NINT                   Number of observation times                           
*  NOPTS                  Used to dimension IOPTS
*  NPAR                   Total number of parameters to be estimated            
*  NPARALG                Maximum number of algorithm parameters, including     
*                         minimization and simulation (used for dimensioning)   
*  NPARF                  Number of transient parameters to be estimated        
*  NPARFPRG               Number of uncertain generic parameter zones involved  
*                         in the non-linear flow inverse problem                
*  NPARPRG                Total number of uncertain generic parameter zones     
*                         involved in the inverse problem                       
*  NPBFL                  Number of simultaneous flow problems                  
*  NPBTP                  Number of simultaneous transport problems             
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)               
*  NTDMT                  If diferent form zero, number of terms in matrix      
*                         diffusion(if zero, no diffusion)                      
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NUMTIT                 Total number of integration times                     
*  NUMTNOD                Total number of nodes used for calculating obs.       
*  NUMTOBS                Total number of observations                          
*  NWRITE                 Number of output options (used for dimensioning)      
*
* INTERNAL VARIABLES: SCALARS
*
*  I                                                                            
*  IERROR_AUX                                                                   
*  IINPUT_FILE_SYSTEM     Number of standard input                              
*  ILEN                                                                         
*  IOUTPUT_FILE_SYSTEM    Number of standard output                             
*  IUDIM                  Unit number of DIM file                               
*  J                                                                            
*  N_FIL                  Total number of files (input and output)              
*  N_FIL_DAT              Total number of input files                           
*  TOLER1                 Used to check XLAMHEAD  
*  TOLER2                 Used to check XLAMHEAD                                
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  FECHA                  Writes date and TITLE on MAIN FILE                    
*  FILEEX                 Cheks the existence of a file                         
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*  OPENFILES              Opens all input, output and auxiliary files           
*  SRC_NCARD              Loocks for the current card label on DIM FILE         
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  FECHA                  Writes date and TITLE on MAIN FILE                    
*  FILEEX                 Cheks the existence of a file                         
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*  OPENFILES              Opens all input, output and auxiliary files           
*  OPT_GROUPS_ZONES       Reads groups of zones estimation options and 
*                         dimensions
*  SRC_NCARD              Looks for the current card label on DIM FILE         
*
* HISTORY
*
*     AMS        1988     First coding
*     SCR      4-1997     Revision and verification
*     AMS      1-1998     Revision
*     AAR      7-2001     Revision. Inclusion of observation related variables
*                         and mass balance options and checking
*     AAR      4-2002     Inclusion of geostatistics related variables
*     AAR      5-2003     Inclusion of the concept of groups of zones
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER NCARD(15)*80,FILENAME(20)*20,FILETYPE(20)*3,ROOT*20,
     ;           TITLE*80,LEAUX*100,LEEL*100

       INTEGER*4 
     ;   NTYPAR,NPARALG,NWRITE,NFLAGS,NSTAT,
     ;   NZONE_PAR(NTYPAR),IOLG_PAR(NTYPAR,2),IPAR_INV(NPARALG),
     ;   IPAR_DIR(NPARALG),IOWRITE(NWRITE),IFLAGS(NFLAGS),IOPTS(NOPTS),
     ;   IOPT_GS(MXGRPZN,20),IO_KG_GS(MXGRPZN,16)
     ; !NUEVO
     ; ,LINMET(3,2)

       REAL*8 
     ;   PAR_WGT(NTYPAR),PAR_INV(NPARALG),PAR_DIR(NPARALG),
     ;   FOBJ_WGT(NSTAT)

       LOGICAL FILEEX

       DATA NCARD/'TITLE','PRT_OPTIONS','DEF_OPTIONS','DIMENSIONS',
     .            'ZONE_NUMBERS','LOG_ESTIMATION_OPTIONS',
     .            'OUTPUT_OPTIONS','MASS_BALANCE',
     .            'OPTIMIZATION_PARAMETERS','WEIGHTING_PARAMETERS',
     .            'DIRECT_PROB_CONV_PARAMETERS','ITERATIVE_SCHEME',
     ;            'WATSOLVE PARAMETERS',
     ;            'GEOSTATISTICAL INV. PROB. OPTIONS',
     ;            'IFLAGS'/


       DATA FILETYPE/'DIM','GRI','PAR','TIM','OBS','INI','GEO','COV',
     ;   'RES','PLT','MSH','MTH','MHH','MSC','MTC','MCC','SEC','CVM',
     &   'BLH','BLC'/
     
       DATA N_FIL_DAT/8/,N_FIL/20/,
     .      IINPUT_FILE_SYSTEM/5/,IOUTPUT_FILE_SYSTEM/6/

C------------------------- FIRST EXECUTABLE STATEMENT.
C------------------------- If it is not the first time we enter to this 
C------------------------- routine, jump the basic definition of the problem

       IUDIM=10
       MAINF=25
       IF (IPROCESS .NE. 0) GOTO 300

C------------------------- Opens dimensions and options file

       ROOT='                    '
       IF (FILEEX('ROOT')) THEN
          OPEN(UNIT=666,FILE='ROOT',STATUS='OLD')
          READ(666,'(A20)') ROOT
          CLOSE(666)
          IOPTS(50)=1    ! para saber que venimos del puto visual
       ELSE
          WRITE(IOUTPUT_FILE_SYSTEM,3000)
          READ(IINPUT_FILE_SYSTEM,1000) ROOT
 1000     FORMAT(A20)
 3000     FORMAT(1X,'NAME OF DIMENSIONS INPUT FILE ... ',$)
          IOPTS(50)=0    ! venimos de la forma normal
       END IF
        
C------------------------- Look for root system

       IF (INDEX(ROOT,'.') .EQ.  0) THEN
          ILEN=INDEX(ROOT,' ')
          IF (ILEN .GT. 20 ) THEN
             WRITE(IOUTPUT_FILE_SYSTEM,3100)
 3100        FORMAT(/,'ERROR: ROOT TOO LARGE, IS HAS TO HAVE LESS',
     ;           ' THAN 20 CHARACTERS')
             STOP
          ENDIF

C------------------------- Builds filenames

          DO I=1,N_FIL_DAT
             FILENAME(I)=ROOT(1:ILEN-1)//FILETYPE(I)//'.DAT'
          ENDDO
          DO I=N_FIL_DAT+1,N_FIL
             FILENAME(I)=ROOT(1:ILEN-1)//FILETYPE(I)//'.OUT'
          ENDDO
          IF (FILEEX(FILENAME(1))) THEN
             OPEN(UNIT=IUDIM,FILE=FILENAME(1),STATUS='OLD')

C------------------------- Reads TITLE, Card A3.1

             CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(1)
     ; ,FILENAME)
             LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
             READ(LEAUX,1100) TITLE 
          ELSE
             WRITE(IOUTPUT_FILE_SYSTEM,3200) FILENAME(1)
 3200        FORMAT(/2X,'DIMENSION INPUT FILE ',A20,3X,'NOT FOUND')
             STOP 'FILE NOT FOUND'
          ENDIF
       ELSE
          FILENAME(1)=ROOT
          IF (FILEEX(FILENAME(1))) THEN
             OPEN(UNIT=IUDIM,FILE=FILENAME(1),STATUS='OLD')
 1100        FORMAT(A80)

C------------------------- Reads input filenames from DIM FILE and TITLE
C------------------------- Cards groups A1,A2 and A3.1

             NROW=0
             READ(LEAUX,1100) (FILENAME(I),I=2,N_FIL)
       
             CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(1)
     ; ,FILENAME)
             LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
             READ(LEAUX,1100) TITLE
          ELSE
             WRITE(IOUTPUT_FILE_SYSTEM,3200) FILENAME(1)
             STOP 'FILE NOT FOUND'
          ENDIF
       ENDIF

       OPEN(MAINF,FILE=FILENAME(N_FIL_DAT+1),STATUS='UNKNOWN')

C------------------------- It writes filenames and date on MAIN FILE

       CALL FECHA (MAINF,TITLE)

       WRITE (MAINF,3300) (FILENAME(I),I=2,N_FIL_DAT)
 3300  FORMAT(////,
     .             1X,19('*'),' GENERAL INFORMATION ',19('*'),////,
     .             10X,'INPUT FILE NAMES',/,
     .             10X,17('='),//,
     .              5X,'GRID DEFINITION FILE.........: ',A20,/,
     .              5X,'PARAMETERS DEFINITION FILE...: ',A20,/,
     .              5X,'TIME FUNCTIONS...............: ',A20,/,
     .              5X,'OBSERVATIONS DATA............: ',A20,/,
     .              5X,'INITIAL CONDITIONS...........: ',A20,/,
     .              5X,'GEOSTATISTICS INFO...........: ',A20,/,
     .              5X,'A PRIORI COVARIANCE MATRIX...: ',A20,/)

C------------------------- Reads I/O options:
C------------------------- Card A3.2 Controling printing of input data

       CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(2)
     ; ,FILENAME)
       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ(LEAUX,1300) INPWR,IOWAR,IOPART


C------------------------- Writes in MAIN FILE if allowed

       IF (INPWR.NE.0) WRITE(MAINF,3400) INPWR,IOWAR,IOPART

 3400  FORMAT(////,10X,'OUTPUT OPTIONS',/,
     ;             10X,'==============',//,
     ; 5X,'INPWR   =',I5/,      
     ; 5X,'IOWAR   =',I5/,
     ; 5X,'IOPART  =',I5)

C------------------------- Reads I/O options
C------------------------- Card A3.3, Problem definitions options

       CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(3)
     ; ,FILENAME)
       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ (LEAUX,1300) IOEQT,IOINV,IOTRS,IOCNSF,IORTS,IOCNST,
     ;            IOVAR,IOFOD,IODIM,IOFLLI,IOTRLI,IOINV_GS

 1300  FORMAT(18I5)

C------------------------- Writes in MAIN FILE if allowed

       IF (INPWR.NE.0) WRITE(MAINF,3500) 
     ;  IOEQT,IOINV,IOTRS,IOCNSF,IORTS,IOCNST
     ; ,IOVAR,IOFOD,IODIM,IOFLLI,IOTRLI,IOINV_GS
   
 3500  FORMAT(////10X,'PROBLEM DEFINITION OPTIONS',/,
     ;            10X,'==========================',//,
     ; 5X,'EQUATION .......................................... =',I5,/,
     ; 5X,'INVERSE PROBLEM ................................... =',I5,/,
     ; 5X,'FLOW REGIME ....................................... =',I5,/,
     ; 5X,'SCHEME IN TIME (LUMPED OR CONSISTENT)(FLOW) ....... =',I5,/,
     ; 5X,'TRANSPORT REGIME .................................. =',I5,/,
     ; 5X,'SCHEME IN TIME (LUMPED OR CONSISTENT)(TRANSPORT) .. =',I5,/,
     ; 5X,'STATISTICAL ANALYSIS OPTION ....................... =',I5,/,
     ; 5X,'DESINTEGRATION OPTION (YES=T/NOT=F) ............... =',I5,/,
     ; 5X,'MAXIM DIMENSION OF THE PROBLEM..................... =',I5,/,
     ; 5X,'NON LINEAR OPTION IN FLOW EQUATION................. =',I5,/,
     ; 5X,'NON LINEAR OPTION IN TRANSPORT EQUATION............ =',I5,/,
     ; 5X,'GEOSTATISTICAL INVERSE PROBLEM OPTION ............. =',I5,/)

      IF((IOTRS.EQ.0.AND.IOFLLI.NE.0) .OR. (IOEQT.EQ.2.AND.IOFLLI.NE.0))
     ;CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;' INCONGRUENT NON LINEAR OPTION'
     ; //'IN FLOW EQUATION',NROW,0,IUDIM,1,1.04)

      IF((IORTS.EQ.0.AND.IOTRLI.NE.0) .OR. (IOEQT.EQ.1.AND.IOTRLI.NE.0))
     ;CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;' INCONGRUENT NON LINEAR OPTION IN'
     ;//'TRANSP. EQUATION',NROW,0,IUDIM,1,1.05)

C------------------------- Sets IORTS or IOTRS to zero if the related problem
C------------------------- is not solved

       IF (IOEQT.EQ.1) IORTS=0
       IF (IOEQT.EQ.2) IOTRS=0

C------------------------------------read iosparse, iodirect and IODENS_INI

       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
10000  FORMAT(5I5)
       READ (LEAUX,10000) IOSPARSE,IODIRECT,IODENS_INI, ITPTVAR,IOCONSRC

C-------------------------Checks solver scheme.

       IF((IODIRECT.EQ.1 .AND. IOSPARSE.EQ.1) .OR.
     ;         (IODIRECT.EQ.0 .AND. IOSPARSE.EQ.0)) THEN

          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;'SOLVER SCHEME AND STORAGE SCHEME INCONSISTENT'
     ;,NROW,2,IUDIM,2,1.04)

       END IF !IODIRECT.EQ.1 .AND. IOSPARSE.EQ.1) .OR. ...


C-------------------------Checks concentration source option.

      IF (IOCONSRC.EQ.1 .AND. (IODENS_INI.EQ.0 .OR. ITPTVAR.EQ.1)) THEN

          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     &'CONC. SOURCE OPT. INCONSISTENT WITH DENSITY OPT. OR TRANSP. VAR.
     & CHANGED TO 0.',NROW,5,IUDIM,0,0.0)

	    IOCONSRC = 0

	 END IF ! IOCONSRC.EQ.1 .AND. (IODENS_INI.EQ.0 .OR. ITPTVAR.EQ.1



      IF (INPWR.NE.0) WRITE(MAINF,10010) IOSPARSE,IODIRECT,IODENS_INI
     ;   ,ITPTVAR,IOCONSRC

10010 FORMAT(
     ; 5X,'SPARSE MATRIX STORAGE OPTION....................... =',I5,/,
     ; 5X,'DIRECT LINEAR SYSTEM SOLVER........................ =',I5,/,
     ; 5X,'VARIABLE DENSITY OPTION............................ =',I5,/,
     ; 5X,'TRANSPORT OF ENERGY OR MASS........................ =',I5,/
     & 5X,'CONCENTRATION SOURCES IN FLOW...................... =',I5,/)


C------------------------- Read dimensions
C------------------------- Card A4.1

       CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(4)
     ; ,FILENAME)
       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ(LEAUX,1300) NUMEL,NUMNP,LMXNDL,ISOT,NBAND,NPAR,NINT,
     ;  NFNT,NPARF,NTDMT,NFNL,NBLCVP,IOFLSAT,IOFMLF,IOFMLT

C-------------------------  Read observations related dimensions

       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ(LEAUX,1300) NDEVS,NUMTOBS,NUMTNOD,NUMTIT,NBANDCOV

C------------------------- Next variables are initialized for 
C------------------------- compatibility with older DIM files

       IOSMFL=0
       IOSMTP=0
       IORDCH=0
       NPBFL=0
       NPBTP=0
       IOPTS(31)=0      ! IOVRWC (variable water content)

       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ(LEAUX,1300,ERR=9100) 
     ;                     IOSMFL,IOSMTP,IORDCH,NPBFL,NPBTP,IOPTS(31)

C------------------------- If zero, the number of flow and transport problems
C------------------------- are set to one

       IF (NPBFL.LE.0 .AND. IOEQT.NE.2) NPBFL=1
       IF (NPBTP.LE.0 .AND. IOEQT.NE.1) NPBTP=1

C------------------------- Writes in MAIN FILE

 9100  IF (INPWR.NE.0) WRITE(MAINF,3600) NUMEL,NUMNP,LMXNDL,ISOT,NBAND,
     ;  NPAR,NINT,NFNT,NPARF,NTDMT,NFNL,NBLCVP,IOFLSAT,
     ;  IOFMLF,IOFMLT,NDEVS,NUMTOBS,NUMTNOD,NUMTIT,NBANDCOV,
     ;  IOSMFL,IOSMTP,IORDCH,NPBFL,NPBTP,IOPTS(31)

 3600  FORMAT(////10X,'DIMENSIONS',/,
     ;            10X,'==========',//,
     ; 5X,'NUMBER OF ELEMENTS .............................. =',I5,/,
     ; 5X,'NUMBER OF NODAL POINTS .......................... =',I5,/,
     ; 5X,'MAX. NUMBER OF NODES FOR ELEM. .................. =',I5,/,
     ; 5X,'ISOTROPY ........................................ =',I5,/,
     ; 5X,'BANDWIDTH ....................................... =',I5,/,
     ; 5X,'NUMBER OF UNCERTAIN PARAMETERS .................. =',I5,/,
     ; 5X,'NUMBER OF TIME INTERVALS ........................ =',I5,/,        
     ; 5X,'NUMBER OF TIME FUNCTIONS ........................ =',I5,/,
     ; 5X,'FLOW PARAMETERS FOR ESTIMATION .................. =',I5,/,
     ; 5X,'MATRIX DIFFUSION OPTION ......................... =',I5,/,
     ; 5X,'NON LINEAR FUNCTIONS NUMBER...................... =',I5,/,
     ; 5X,'STRUCTURES OF PARAMETERS COVARIANCE MATRIX....... =',I5,/,
     ; 5X,'UNSATURATED (1) OR SATURATED (0) FLOW ........... =',I5,/,
     ; 5X,'FORMULATION OF NON-LINEAR FLOW EQUATION ......... =',I5,/,
     ; 5X,'FORMULATION OF NON-LINEAR TRANSPORT EQUATION .... =',I5,/,
     ; 5X,'NUMBER OF DEVICES ............................... =',I5,/,
     ; 5X,'NUMBER OF OBSERVATIONS (ALL DEVICES) ............ =',I5,/,
     ; 5X,'NUMBER OF NODES USED TO DEFINE INTEGR. SPACES ... =',I5,/,
     ; 5X,'NUM. OF INTEGR. TIMES USED TO DESCRIBE ALL OBS .. =',I5,/,
     ; 5X,'BANDWIDTH OF THE COVAR. MATRIX .................. =',I5,/,
     ; 5X,'SIMULTANEOUS FLOW PROBLEMS ...................... =',I5,/,
     ; 5X,'SIMULTANEOUS TRANSPORT PROBLEMS ................. =',I5,/,
     ; 5X,'RADIOACTIVE CHAINS .............................. =',I5,/,
     ; 5X,'NUMBER OF FLOW PROBLEMS ......................... =',I5,/,
     ; 5X,'NUMBER OF TRANSPORT PROBLEMS .................... =',I5,/,
     ; 5X,'ACTUALIZATION OF WATER CONTENT .................. =',I5)

C------------------------- Check non linear functions

       IF(IOTRLI.EQ.0.AND.IOFLLI.EQ.0.AND.NFNL.NE.0)
     ; CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;' INCONGRUENT NON LINEAR FUNCTIONS' 
     ;//'NUMBER (NFNL)',NROW,0,IUDIM,1,1.06)

C------------------------- Check parameters cova. matrix structure

      IF (IOINV.GT.0 .AND. (NBLCVP.LT.0 .OR. NBLCVP.GT.NPAR)) 
     ;   CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'INVERSE PROBLEM REQUESTED AND NUMBER OF BLOCKS OF COVARIANCE' 
     ;//' MATRIX (COV FILE) IS OUT OF ORDER: < 0 OR >NPAR)',NROW,2,
     ;   IUDIM,2,1.07)

      IF (IOINV.LE.0 .AND. NBLCVP.GT.0) 
     ;   CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'DIRECT PROBLEM REQUESTED AND NUMBER OF BLOCKS OF COVARIANCE' 
     ;//' MATRIX (COV FILE) IS OUT OF ORDER: > 0',NROW,2,
     ;   IUDIM,1,1.07)

C------------------------- Check type of formulation choosen


C------------------------- Corriges values for dimension statements

       IF (NDEVS.EQ.0) NDEVS=1
       IF (NINT.EQ.0) NINT=1
       IDIMFNT=MAX(1,NFNT)
       IF (NFNL.EQ.0) NFNL=1
       IF (NPBFL.LE.0) NPBFL=1
       IF (NPBTP.LE.0) NPBTP=1


C------------------------- Reads numbers of zones of all parameters
C------------------------- Card A4.2 
       
       CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(5)
     ; ,FILENAME)
       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)

       READ(LEAUX,1300) NZTRA,NZSTG,NZARR,NZCHP,NZQQP,NZALF,NZDSP,
     ; NZDFM,NZPOR,NZFOD,NZCRD,NZCOE,NZDMT,NZPRG,NPARFPRG,NPARPRG,NZCLK,
     & NZFOF
       
       IF (INPWR.NE.0) WRITE(MAINF,3700) 
     ;    NZTRA,NZSTG,NZARR,NZCHP,NZQQP,NZALF,NZDSP,
     ;    NZDFM,NZPOR,NZFOD,NZCRD,NZCOE,NZDMT,NZPRG,NPARFPRG,NPARPRG,
     &    NZCLK,NZFOF

 3700  FORMAT(//,10X,'NUMBERS OF ZONES',/
     ;           10X,'================',//
     ; 5X,'TRANSMISIVITY ...... =',I5,/
     ; 5X,'STORAGE ............ =',I5,/
     ; 5X,'RECH. COEFF. ....... =',I5,/
     ; 5X,'PRESC. HEAD ........ =',I5,/
     ; 5X,'PRESC. FLOW ........ =',I5,/
     ; 5X,'LEAKAGE ............ =',I5,/
     ; 5X,'DISPERSIVITY ....... =',I5,/
     ; 5X,'MOLEC. DIFFUSION ... =',I5,/
     ; 5X,'POROSITY ........... =',I5,/
     ; 5X,'FIRST ORDER DECAY .. =',I5,/
     ; 5X,'RETARD. COEFF. ..... =',I5,/
     ; 5X,'EXTERNAL CONCENT.... =',I5,/
     ; 5X,'MATRIX DIFFUSION.... =',I5,/
     ; 5X,'GROUP PARAMETERS.... =',I5,/
     ; 5X,'FLOW GENERIC ESTIMATED PARAMETERS=',I5,/
     ; 5X,'TOTAL GENERIC ESTIMATED PARAMETERS=',I5,/
     ; 5X,'CONC. LEAKAGE....... =',I5,/
     ; 5X,'FORMATION FACTOR.... =',I5)

C------------------------- Corriges values for dimension statment relate to NPAR
       
       IF (NPAR.EQ.0 .AND. NPARPRG.EQ.0) NPAR=1

C------------------------- Check the number of zones

       IF (IOEQT.NE.2) THEN
          IF (NZTRA.EQ.0)
     ;      CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                'NUMBER TRANSMISSIVITY ZONES IS ZERO ',
     ;                 NROW,1,IUDIM,1,1.01)

          IF (IOTRS.NE.0) THEN
             IF (NZSTG.EQ.0)
     ;          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                     'NUMBER STORAGE COEFF. ZONES IS ZERO ',
     ;                      NROW,2,IUDIM,1,1.02)
          ELSE 
             IF (NZSTG.NE.0)
     ;         CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                   'STEADY STATE, BUT THE NUMBER OF STORAGE'
     ;                 //' COEFF. ZONES IS'
     ;                 //'NOT EQUAL ZERO ',NROW,2,IUDIM,0,0.00)
     ;      
          ENDIF 
       ENDIF

       IF (IOEQT.NE.1 .AND. NZPOR.EQ.0)
     ;   CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;             'NUMBER OF POROSITY ZONES IS ZERO ',
     ;              NROW,9,IUDIM,1,1.03)

C------------------------- Matrix diffusion checks

       IF (IOEQT.EQ.1. AND .NZDMT.NE.0)CALL ERROR 
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'INCORRECT MATRIX DIFUSSION ZONES ',NROW,12,IUDIM,1,1.10)
       IF (IOEQT.EQ.1. AND .NTDMT.NE.0)CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;      'INCORRECT MATRIX DIFUSSION OPTION ',NROW,11,IUDIM,1,1.11)
       IF (NTDMT.EQ.0. AND .NZDMT.NE.0)CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;      'INCORRECT MATRIX DIFUSSION ZONES ',NROW,12,IUDIM,1,1.12)
       IF (NTDMT.NE.0. AND .NZDMT.EQ.0)CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;      'INCORRECT MATRIX DIFUSSION ZONES ',NROW,12,IUDIM,1,1.13)
       IF (NTDMT.NE.0. AND .NZCRD.EQ.0)CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ; 'IF MATRIX DIF. IS CONSIDERED, USER MUST DEFINE A RETARD. ZONE'
     ;       ,NROW,10,IUDIM,1,1.14)
       IF (NTDMT.NE.0. AND .NZDFM.EQ.0)CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ; 'IF MATRIX DIF. IS CONSIDERED, USER MUST DEFINE A MOL. DIF. ZONE'
     ; ,NROW,8,IUDIM,1,1.15)

C------------------------- Non linear checks

       IF (IOTRLI.EQ.0.AND.IOFLLI.EQ.0.AND.NZPRG.NE.0)
     ;   CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;             'INCONGRUENT GROUP PARAMETERS ZONE NUMBER' 
     ;             ,NROW,0,IUDIM,1,1.16)

C------------------------- Generic parameters checks

       IF ((IOFLLI.EQ.0.AND.NPARFPRG.NE.0).OR.
     ;     ((IOFLLI+IOTRLI).EQ.0.AND.NPARPRG.NE.0).OR.
     ;     (NPARPRG.LT.NPARFPRG).OR.
     ;     (NFNL.EQ.0.AND.(NPARFPRG+NPARFPRG).NE.0).OR.
     ;     (NZPRG.EQ.0.AND.(NPARFPRG+NPARFPRG).NE.0))
     ;    CALL ERROR
     ;   (IERROR,IOWAR,MAINF,FILENAME,
     ;   ' NUMBER OF GENERIC ESTIMATED PARATETERS IS '
     ; //'INCOMPATIBLE WITH ANOTHER OPTIONS. CHECK IT',
     ;    NROW,1,IUDIM,2,1.17)

C------------------------- Updates number of total parameters for estimation

       NPAR=NPAR+NPARPRG
       NPARF=NPARF+NPARFPRG

C------------------------- Reads logarithmic estimation options
C------------------------- Card A5.1


 300   CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(6)
     ; ,FILENAME)
       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)

       READ(LEAUX,1300) IOLGTRA,IOLGSTG,IOLGARR,IOLGCHP,IOLGQQP,
     ;        IOLGALF,IOLGDSP,IOLGDFM,IOLGPOR,IOLGFOD,IOLGCRD,IOLGCOE,
     ;        IOLGPRG,IOSUCHUM,IOLGCLK,IOLGFOF

C------------------------- Writes in MAIN FILE if allowed


C------------------------- Checks option of type of measures

       IF (IOSUCHUM.NE.0.AND.
     ;    (IOFLSAT.EQ.0.OR.IOEQT.EQ.2.OR.IOINV.EQ.2))
     ;      CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;     'IT DOESNT HAVE SENSE MEASURES IN HUMIDITY THERMS. '
     ;   //'CHECK VALUES OF IOSUCHUM-IOEQT-IOINV-IOFLSAT',
     ;      NROW,1,IUDIM,2,1.18)
       
       IF (INPWR.NE.0) WRITE(MAINF,3800) 
     ;    IOLGTRA,IOLGSTG,IOLGARR,IOLGCHP,IOLGQQP,
     ;    IOLGALF,IOLGDSP,IOLGDFM,IOLGPOR,IOLGFOD,IOLGCRD,IOLGCOE,
     ;    IOLGPRG,IOSUCHUM,IOLGCLK,IOLGFOF

 3800 FORMAT(////,10X,' DERIVABILITY OPTIONS',/,
     .              9X,'======================',//,
     ; 5X,'LOG-DERIVATIVE OF TRANSMISIVITY (YES=T/NOT=F) .... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF STORAGE (YES=T/NOT=F) .......... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF RECHARGE (YES=T/NOT=F) ......... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF PRESC. HEAD (YES=T/NOT=F) ...... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF PRESC. FLOW (YES=T/NOT=F) ...... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF LEAKAGE (YES=T/NOT=F) .......... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF DISPERSIVITY (YES=T/NOT=F) ..... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF MOLEC. DIFF. (YES=T/NOT=F) ..... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF RAD. DESIN. (YES=T/NOT=F) ...... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF RET. COEFF. (YES=T/NOT=F) ...... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF POROSITY (YES=T/NOT=F) ......... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF EXT. CONC. (YES=T/NOT=F) ....... =',I5,/
     ; 5X,'LOG-DERIVATIVE OF GENERIC PARAMETERS (YES=T/NOT=F). =',I5,/
     ; 5X,'SUCTIONS OR HUMIDITY MEASURES FOR I.PROBLEM (0=SUCTIONS)=',I5
     &,/
     & 5X,'LOG-DERIVATIVE OF CONC. LEAKAGE (YES=T/NOT=F) .... =',I5,/
     & 5X,'LOG-DERIVATIVE OF FORMATION FACTOR (YES=T/NOT=F) . =',I5)
      
C------------------------- Read output options
C------------------------- Card A6.1

       CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(7)
     ; ,FILENAME)
       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ(LEAUX,1300) IOWRH,IOWRC,IOPLH,IOPLC,IOMHH,IOMHC,
     ;                  IOSEH,IOSEC,IOCMH,IOCMC,IOWPI,IOWIT,IOWVD

C------------------------- If flow equation is not solved, flow output options
C------------------------- are set to zero

       IF (IOEQT.EQ.2) THEN
          IOWRH=0
          IOPLH=0
          IOMHH=0
          IOSEH=0
          IOCMH=0
       ENDIF

C------------------------- If transport equation is not solved, transport
C------------------------- output options are set to zero

       IF (IOEQT.EQ.1) THEN
          IOWRC=0
          IOPLC=0
          IOMHC=0
          IOSEC=0
          IOCMC=0
       ENDIF

       IF (INPWR .NE. 0) WRITE(MAINF,3900) 
     ;    IOWRH,IOWRC,IOPLH,IOPLC,IOMHH,IOMHC,
     ;    IOSEH,IOSEC,IOCMH,IOCMC,IOWPI,IOWIT,IOWVD

 3900  FORMAT(////,10X,'OUTPUT OPTIONS',/,
     .             10X,'=============='//,
     .           5X,'HEAD RESIDUALS ...........................=',I5,/,
     .           5X,'CONCENTRATION RESIDUALS ..................=',I5,/,
     .           5X,'COMP. MEASURED HEAD VS TIME ..............=',I5,/,
     .           5X,'COMP. MEASURED CONCENTRATION VS TIME .....=',I5,/,
     .           5X,'HEAD CONTOUR MAPS ........................=',I5,/,
     .           5X,'CONCENTRATION CONTOUR MAPS ...............=',I5,/,
     .           5X,'NUM. HEAD CROSS SECTIONS .................=',I5,/,
     .           5X,'NUM OF CONCENTR. CROSS SECTIONS ..........=',I5,/,
     .           5X,'COMPUTED VS MEASURED HEAD ................=',I5,/,
     .           5X,'COMPUTED VS MEASURED CONCENTR.  ..........=',I5,/,
     .           5X,'INVERSE PROBLEM INFORMATION ..............=',I5,/,
     .           5X,'EVOLUTION OF PARAMETER VALUES ............=',I5,/,
     .           5X,'DARCY FLOW OUTPUT ........................=',I5)

C------------------------- Reads mass balance options
C------------------------- Card A7.1

       CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(8)
     ; ,FILENAME)
       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ(LEAUX,1300) IOBALH,IOBALC,IOBALGH,IOBALGC,IOBALDH,IOBALDC

C------------------------- Checks mass balance options


C------------------------- Flow problem and tpt. m.b. requested

       IF (IOEQT.EQ.1.AND. 
     ;    (IOBALC.NE.0.OR. IOBALGC.NE.0.OR.IOBALDC.NE.0)) THEN
         IF (IOWAR.NE.0) WRITE(MAINF,666) 
         IOBALC=0
         IOBALGC=0
         IOBALDC=0
       END IF

C------------------------- Transport problem and flow. m.b. requested

       IF (IOEQT.EQ.2.AND. 
     ;    (IOBALH.NE.0.OR. IOBALGH.NE.0.OR.IOBALDH.NE.0)) THEN
         IF (IOWAR.NE.0) WRITE(MAINF,667) 
         IOBALH=0
         IOBALGH=0
         IOBALDH=0
       END IF

C------------------------- Nodal f.m.b. requested without zonal balance 
C------------------------- computations

       IF (IOBALDH.NE.0.AND.(IOBALH+IOBALGH).EQ.0) THEN
         IOBALH=1
         IOBALGH=1
         IF (IOWAR.NE.0) WRITE(MAINF,668)
       END IF

C------------------------- Nodal tpt.m.b. requested without zonal balance 
C------------------------- computations

       IF (IOBALDC.NE.0.AND.(IOBALC+IOBALGC).EQ.0) THEN
         IOBALC=1
         IOBALGC=1
         IF (IOWAR.NE.0) WRITE(MAINF,669)
       END IF

C------------------------- Nodal f.m.b. option out of order (simul. times)

       IF (IOBALH.GT.NINT) THEN
         IOBALH=NINT
         IF (IOWAR.NE.0) WRITE(MAINF,670)
       END IF

C------------------------- Nodal t.m.b. option out of order (simul. times)

       IF (IOBALC.GT.NINT) THEN
         IOBALC=NINT
         IF (IOWAR.NE.0) WRITE(MAINF,671)
       END IF

       IF (IOBALDC.NE.0 .AND. IOBALDH.EQ.0) IOBALDH=IOBALDC

       IF (IOTRS.EQ.0.AND.IOBALGH.NE.0) THEN
         IOBALH=IOBALGH
         IOBALGH=0
       END IF

       IF (IORTS.EQ.0.AND.IOBALGC.NE.0) THEN
         IOBALC=IOBALGC
         IOBALGC=0
       END IF

C-------------------------  Writes mass balance options

       IF (INPWR.NE.0) THEN
          WRITE(MAINF,4000) IOBALH,IOBALC,IOBALGH,IOBALGC,IOBALDH
     ;                     ,IOBALDC


 4000     FORMAT(////,10X,'MASS BALANCE OPTIONS',/,
     ;             10X,'===================='//,
     ;       5X,'FLOW MASS BALANCE ....................... =',I5,/,
     ;       5X,'TRANSPORT MASS BALANCE .................. =',I5,/,
     ;       5X,'GLOBAL FLOW MASS BALANCE ................ =',I5,/,
     ;       5X,'GLOBAL TRANSPORT MASS BALANCE ........... =',I5,/,
     ;       5X,'NODAL FLOW MASS BALANCE ................. =',I5,/,
     ;       5X,'NODAL TRANSPORT MASS BALANCE ............ =',I5,/)

       END IF

 666   FORMAT (//,'WARNING:',
     & ' TRANSPORT NOT SOLVED BUT TRANSPORT MASS BALANCE OPTIONS ARE'
     & ' NON ZERO')

 667   FORMAT (//,'WARNING:',
     & ' FLOW NOT SOLVED BUT FLOW MASS BALANCE OPTIONS ARE NON ZERO')

 668   FORMAT (//,'WARNING:',
     & ' INCOHERENT DEFINITION OF FLOW MASS BALANCE OPTIONS.',/,
     & ' FLOW NODAL MASS BAL. REQUIRES ZONAL/GLOBAL MASS BAL.'
     & ' COMPUTATION',/,' WILL BE CHANGED')

 669   FORMAT (//,'WARNING:',
     & ' INCOHERENT DEFINITION OF TRANSPORT MASS BALANCE OPTIONS.',/,
     & ' TRANSPORT NODAL MASS BAL. REQUIRES ZONAL/GLOBAL MASS BAL.'
     & ' COMPUTATION',/,' WILL BE CHANGED')

 670   FORMAT(//,'WARNING:',
     ; ' IOBALH IS GREATER THAN MAXIMUM NUMBER OF OBSERVATION TIMES')
 671   FORMAT(//,'WARNING:',
     ; ' IOBALC IS GREATER THAN MAXIMUM NUMBER OF OBSERVATION TIMES')

C------------------------- Reads optimization parameters
C------------------------- Cards A8.1, A8.2

       CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(9)
     ; ,FILENAME)

       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ(LEAUX,1400)   !Card A8.1
     ;     XMARQ,NUMIN,NUMAX,PHIMIN,PHIMAX,GMNOR1,GMNOR,DMINF,COSMIN

       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ (LEAUX,1410)  !Card A8.2
     ;     MAXICOS,PERMX1,PERMX2,EPS,MAXITER,NMTERF1

1400   FORMAT (F10.0,2I5,6F10.0)
1410   FORMAT (I5,3F10.0,2I5)

C------------------------- Reads weighting coefficients
C------------------------- Cards A9.1 and A9.2

       CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(10)
     ; ,FILENAME)

       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ(LEAUX,1500)   !Card A9.1
     ;     XLAMTRA,XLAMSTG,XLAMARR,XLAMCHP,XLAMQQP,XLAMALF,XLAMPRGF

       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ(LEAUX,1500)   !Card A9.2
     ;     XLAMDSP,XLAMDFM,XLAMPOR,XLAMFOD,XLAMCRD,XLAMCOE,XLAMPRGT,
     &     XLAMCLK, XLAMFOF

       LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
       READ(LEAUX,1500)   !Card A9.2b
     ;     XLAMHEAD,XLAMCONC,XLAMHUMI,XLAMFLOW,XLAMGEO

C------------------------- Echoes a WARNING if XLAMHEAD is not 1.0

       TOLER=0.000001

       IF (XLAMHEAD.GT.1.D0+TOLER .OR. XLAMHEAD.LT.1.D0-TOLER)
     ;   CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'PIEZOMETRIC HEAD WEIGHT IS NOT 1D0',NROW,1,IUDIM,0,1.27)

 1500  FORMAT (9F10.0)
C------------------------- Writes mass balance options, optimization parameters
C------------------------- and weighting parameters in MAIN FILE

       IF (INPWR.NE.0) THEN

          WRITE(MAINF,4100) XMARQ,NUMIN,NUMAX,PHIMIN,PHIMAX,GMNOR1,
     ;        GMNOR,DMINF,COSMIN,MAXICOS,PERMX1,PERMX2,EPS,
     ;        MAXITER,NMTERF1

          WRITE(MAINF,7355) XLAMTRA,XLAMSTG,XLAMARR,XLAMCHP,
     ;        XLAMQQP,XLAMALF,XLAMPRGF,XLAMDSP,XLAMDFM,XLAMPOR,XLAMFOD,
     ;        XLAMCRD,XLAMCOE,XLAMCLK,XLAMFOF,XLAMPRGT,XLAMHEAD,XLAMCONC
     &       ,XLAMHUMI,XLAMFLOW,XLAMGEO

 4100     FORMAT(////,10X,'CONVERGENCE PARAMETERS',/,
     ;             10X,'===================='//,
     ;       5X,'MARQUARDT''S PARAMETER .... =',F10.3,/,
     ;       5X,'NUMIN ..................... =',I5,/,
     ;       5X,'NUMAX ..................... =',I5,/,
     ;       5X,'PHIMIN .................... =',F10.3,/,
     ;       5X,'PHIMAX .................... =',F10.3,/,
     ;       5X,'GMNOR1 .................... =',E10.3,/,
     ;       5X,'GMNOR ..................... =',E10.3,/,
     ;       5X,'DMINF ..................... =',E10.3,/,
     ;       5X,'MINIMUN COSINUS ........... =',E10.3,/,
     ;       5X,'MAXICOS ................... =',I5,/,
     ;       5X,'PERMX1 .................... =',F10.3,/,
     ;       5X,'PERMX2 .................... =',F10.3,/,
     ;       5X,'EPS ....................... =',E10.3,/,
     ;       5X,'MAX. NUMBER OF ITERATES ... =',I5,/,
     ;       5X,'NMTERF1 ................... =',I5,/)
 7355    FORMAT(////,10X,'PONDERATION PARAMETERS ' ,/,
     ;               10X,'=========== ========== ',//,
     ;         ' PONDERATION PARAMETERS FOR FLOW EQ.' ,//,
     ;        ' XLAMTRA .................... =',1P,G15.5,/,
     ;        ' XLAMSTG .................... =',1P,G15.5,/,
     ;        ' XLAMARR .................... =',1P,G15.5,/,
     ;        ' XLAMCHP .................... =',1P,G15.5,/,
     ;        ' XLAMQQP .................... =',1P,G15.5,/,
     ;        ' XLAMALF .................... =',1P,G15.5,/,
     ;        ' XLAMPRGF ................... =',1P,G15.5,//,
     ;        ' PONDERATION PARAMETERS FOR TRANSP. EQ.' ,//,
     ;        ' XLAMDSP .................... =',1P,G15.5,/,
     ;        ' XLAMDFM .................... =',1P,G15.5,/,
     ;        ' XLAMPOR .................... =',1P,G15.5,/,
     ;        ' XLAMLAM .................... =',1P,G15.5,/,
     ;        ' XLAMCRD .................... =',1P,G15.5,/,
     ;        ' XLAMCOE .................... =',1P,G15.5,/,
     &        ' XLAMCLK .................... =',1P,G15.5,/,
     &        ' XLAMFOF .................... =',1P,G15.5,/,         
     ;        ' XLAMPRGT ................... =',1P,G15.5,//,
     ;        ' PONDERATION PARAMETERS FOR STATE VAR.' ,//,
     ;        ' XLAMHEAD ................... =',1P,G15.5,/,
     ;        ' XLAMCONC ................... =',1P,G15.5,/,
     ;        ' XLAMHUMI ................... =',1P,G15.5,/,
     ;        ' XLAMFLOW ................... =',1P,G15.5,/,
     &        ' XLAMGEO ................... =',1P,G15.5 )
       ENDIF

C------------------------- Reads Newton-Raphson parameters process Cards A9.3 and
C------------------------- A9.4. Direct problem convergence parameters.

       DRELMXFL = 0D0
       DABSMXFL = 0D0
       RESIDMXF = 0D0
       ZEROF = 0D0
       DHITMX = 0D0
       MXNRTF = 0
       IOPINITH = 0
       IOWNRFL = 0D0
       DRELMXTR = 0D0
       DABSMXTR = 0D0
       RESIDMXT = 0D0
       ZEROT = 0D0
       DCITMX = 0D0
       MXNRTT = 0
       IOPINITC = 0
       IOWNRTR = 0D0
       OBJHED1 = 0D0
       RESIDMX1F = 0D0
       DABSMX1F = 0D0
       DRELMX1F = 0D0
       OBJHED2 = 0D0
       RESIDMX2F = 0D0
       DABSMX2F = 0D0
       DRELMX2F = 0D0
       OBJCON1 = 0D0
       RESIDMX1T = 0D0
       DABSMX1T = 0D0
       DRELMX1T = 0D0
       OBJCON2 = 0D0
       RESIDMX2T = 0D0
       DABSMX2T = 0D0
       DRELMX2T = 0D0

      IF (IOFLLI.NE.0 .OR. IOTRLI.NE.0 .OR. IODENS_INI .EQ.1) THEN
         
         CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(11)
     ; ,FILENAME)

C------------------------- Read general convergence parameter

         LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
         READ(LEAUX,1600) FCTNCV,FCTDEC,FCTINC,FCTDVNR,MINCAT,IOCRITRAP 
 1600    FORMAT(4F10.0,2I5)

C------------------------- write general convergence parameters to file 

         IF (INPWR.NE.0) WRITE(MAINF,4200) FCTNCV,FCTDEC,FCTINC,FCTDVNR
     ;      ,MINCAT,IOCRITRAP

 4200    FORMAT(///'    NEWTON  RHAPSON / PICARD PARAMETERS',/
     ;             '    ------  -------   ------ ----------',//
     ;' TIME INCREMENT FACTOR TYPE (FCTNCV)...................=',E10.4,/!fctncv CHECK
     ;' TIME INCREMENT FACTOR TYPE (FCTDEC)...................=',E10.4,/!fctdec  CHECK
     ;' TIME DECREASE FACTOR (FCTINC).........................=',E10.4,/!fctinc CHECK
     ;' MAXIM INCREMEMT"S RATIO AS DIVERGENCE CRITERION.......=',E10.4,/!fctdvnr
     ;' MINIM. NUMB. OF CONVERGENCES BEFORE INCREASING TIME...=',I5,/   !mincat NO CHECK
     ;' OPTION OF CHANGE FOR NEW-RAPS. CONVERGENCE CRITERIA ..=',I5,/)!iocritrap

C------------------------- error check
         
         IF(FCTNCV.LT.1)CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     &  'TIME INCREMENT FACTOR TYPE (FCTNCV) HAS TO BE GREATER THAN 0NE'
     &  ,NROW,0,IUDIM,1,1.19)

         IF(FCTINC.LT.1)CALL ERROR
     ;   (IERROR,IOWAR,MAINF,FILENAME,
     ;   'TIME INCREMENT FACTOR (FCTINC) HAS TO BE GREATER THAN ONE'
     ;   ,NROW,0,IUDIM,1,1.20)

         IF(FCTDEC.GT.1. OR .FCTDEC.LT.0)
     ;   CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     & ' TIME DECREMENT TYPE (FCTDEC) HAS TO BE LESS THAN ONE',NROW,0
     & ,IUDIM,1,1.21)

         IF (IOCRITRAP.NE.0 .AND. IOCRITRAP.NE.1) CALL ERROR
     ;   (IERROR,IOWAR,MAINF,FILENAME,
     ;   'IOCRITRAP HAS TO BE EITHER 0 OR 1'
     ;   ,NROW,0,IUDIM,1,1.22)


         IF (IOFLLI.NE.0 .OR. IODENS_INI.EQ.1) THEN

C------------------------- read flow convergence parameters 

            LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
            READ(LEAUX,1610)  DRELMXFL,DABSMXFL,RESIDMXF,ZEROF
     ;     ,DHITMX,MXNRTF,IOPINITH,IOWNRFL
1610        FORMAT(5E10.0,3I5)
   

C------------------------- Write flow convergence parameters to file

            IF (INPWR.NE.0) WRITE(MAINF,10040) '    FLOW EQUATION '

10040       FORMAT(/,'                     ',A18,/)
            IF (INPWR.NE.0) WRITE(MAINF,10030) DRELMXFL,DABSMXFL
     ;      ,RESIDMXF,ZEROF,MXNRTF,IOPINITH,IOWNRFL,DHITMX

10030        FORMAT(/
     ;' OF CONVERGENCE CRITERION RATIO (DRELMX)...............=',E10.4,/        !drelmxfl
     ;' OF CONVERGENCE CRITERION CHANGE (DABSMX)..............=',E10.4,/        !dabxmxfl
     ;' OF CONVERGENCE CRITERION RESIDUO (RESIDMX)............=',E10.4,/        !residmxf 
     ;' VALUE APROX. TO ZERO IN A CHANGE OF STATE EQUATION....=',E10.4,/        !zerof 
     ;' MAXIM NUMBER OF REDUCTION TIME........................=',I5,/           !mxnrtf 
     ;' OPTION OF INITIALIZATION FOR NEW-RAPS.................=',I5,/           !iopinith
     ;' NEWTON RAPSON"S PRINTOUT OPTION ......................=',I5,/           !iownrf
     ;' MAXIM INCREMEMT"S HEADS AS DIVERGENCE CRITERION.......='
     ;                                                  ,E10.4,/)!dhitmx

C------------------------- Newton data check for flow

             IF(RESIDMXF.LT.0.AND.IOFLLI.NE.0) CALL ERROR
     ;         (IERROR,IOWAR,MAINF,FILENAME,
     ;        'RESIDMX HAS TO BE GREATER THAN ZERO'
     ;        ,NROW,0,IUDIM,1,1.22)

             IF(ZEROF.LT.0.AND.IOFLLI.NE.0) CALL ERROR
     ;         (IERROR,IOWAR,MAINF,FILENAME,
     ;         'ZERO HAS TO BE GREATER THAN ZERO'
     ;         ,0,0,IUDIM,1,1.23)

             IF(MXNRTF.LT.1.AND.IOFLLI.NE.0) CALL ERROR
     ;         (IERROR,IOWAR,MAINF,FILENAME,
     ;         'MXNRT HAVE TO BE GREATER THAN ONE',
     ;          NROW,0,IUDIM,1,1.25)

         ENDIF !IF IOFLLI .NE. 0  .OR. IODENS_INI.EQ.1

	     IF (IOTRLI.NE.0 .OR. IODENS_INI.EQ.1) THEN
C------------------------- read transport convergence parameters

            LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
            READ(LEAUX,1610) DRELMXTR,DABSMXTR,RESIDMXT,ZEROT
     ;                   ,DCITMX,MXNRTT,IOPINITC,IOWNRTR
 
C------------------------- Writes on MAIN FILE if allowed

            IF (INPWR.NE.0) WRITE(MAINF,10040) 'TRANSPORT EQUATION'


            IF (INPWR.NE.0) WRITE(MAINF,10030) DRELMXTR,DABSMXTR
     ;      ,RESIDMXT,ZEROT,MXNRTT,IOPINITC,IOWNRTR,DCITMX

C------------------------- Newton data check for transport


            IF(RESIDMXT.LT.0.AND. IOTRLI.NE.0) CALL ERROR
     ;        (IERROR,IOWAR,MAINF,FILENAME,
     ;        'RESIDMX HAS TO BE GREATER THAN ZERO'
     ;        ,NROW,0,IUDIM,1,1.22)

             IF (ZEROT.LT.0.AND.IOTRLI.NE.0)CALL ERROR
     ;          (IERROR,IOWAR,MAINF,FILENAME,
     ;          'ZERO HAS TO BE GREATER THAN ZERO'
     ;          ,0,0,IUDIM,1,1.23)

             IF (MXNRTT.LT.1.AND.IOTRLI.NE.0) CALL ERROR
     ;         (IERROR,IOWAR,MAINF,FILENAME,
     ;         'MXNRT HAVE TO BE GREATER THAN ONE',
     ;          NROW,0,IUDIM,1,1.25)

         ENDIF ! IF (IOTRLI.NE.0 .OR. IODENS_INI.EQ.1) THEN

C------------------------- Checks inverse-problem and Newton-Raphson interactions

         IF ((IOPINITH.NE.0.AND.IOFLLI.EQ.0). OR .(IOPINITC.NE.0.AND.
     ;       IOTRLI.EQ.0). OR .
c     ;      (IOPINITH.NE.0.AND.(IOINV.NE.1.AND.IOINV.NE.3)).OR.
c     ;      (IOPINITC.NE.0.AND.(IOINV.NE.2.AND.IOINV.NE.3)).OR.
     ;       (IOCRITRAP.NE.0.AND.IOINV.LT.1)) THEN

             CALL ERROR
     ;        (IERROR,IOWAR,MAINF,FILENAME,
     ;        'INIT. OR CONVERG. CRITERIA OPTION FOR NEW-RAPSON IS BAD',
     ;        NROW,0,IUDIM,2,1.26)
         
         ENDIF

C------------------------- Reads parameters of convergence criteria variation
C------------------------- Card A9.5

         IF (IOCRITRAP.NE.0)THEN
           IF (IOFLLI.NE.0)THEN
             LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
             READ(LEAUX,'(8F10.0)',ERR=9000)
     ;       OBJHED1,RESIDMX1F,DABSMX1F,DRELMX1F,
     ;       OBJHED2,RESIDMX2F,DABSMX2F,DRELMX2F
           ENDIF !IOFLLI.NE.0
           IF (IOTRLI.NE.0)THEN
             LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
             READ(LEAUX,'(8F10.0)',ERR=9000)
     ;       OBJCON1,RESIDMX1T,DABSMX1T,DRELMX1T,
     ;       OBJCON2,RESIDMX2T,DABSMX2T,DRELMX2T
           ENDIF !IOTRLI.NE.0
           IF (INPWR.NE.0)THEN
             IF(IOFLLI.NE.0)WRITE(MAINF,4300)'   FLOW',
     ;       OBJHED1,RESIDMX1F,DABSMX1F,DRELMX1F,
     ;       OBJHED2,RESIDMX2F,DABSMX2F,DRELMX2F

             IF(IOTRLI.NE.0)WRITE(MAINF,4300)'   TRANSPORT',
     ;       OBJCON1,RESIDMX1T,DABSMX1T,DRELMX1T,
     ;       OBJCON2,RESIDMX2T,DABSMX2T,DRELMX2T

 4300      FORMAT(//' POINTS OF THE STRAIGHT LINE USED TO CHOOSE THE',/,
     ;     ' CONVERGENCE CRITERIA FOR NEWTON RAPHSON PROCESS AT A',/,
     ;     ' GIVEN MARQUARDT ITERATION.',A15,' EQUATION',/,
     ;     '            F.OBJ   ABS-CHANGE  REL-CHANGE  RESIDUO',/,
     ;     ' POINT 1',1X,2E10.4,3X,2E10.4,/,
     ;     ' POINT 2',1X,2E10.4,3X,2E10.4)
           ENDIF !INPWR.NE.0
         ENDIF   !IOCRITRAP.NE.0

       ENDIF !IOFLLI.NE.0 .OR. IOTRLI.NE.0 .OR. IODENS_INI .EQ.1

C------------------------- card a9.6 iterative schemes
      
       IF (IODENS_INI.EQ.1 .OR. IOTRLI.EQ.1 .OR. IOFLLI.EQ.1) THEN
          CALL SRC_NCARD
     ;    (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(12)
     ;    ,FILENAME)
       ENDIF

       IF (IOEQT.NE. 2 .AND. (IOFLLI.NE.0 .OR. IODENS_INI.EQ.1)) THEN !if we solve flow 

C------------------------- read flow iterative scheme

          LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
          READ(LEAUX,10050) ITERCHNGFL,ITERCONVFL ,NRITCTNFL,MAXNUMDIVFL
     &                     ,IDMET1FL  

10050     FORMAT(6I5)

C------------------------- if allowed, write to file

         IF (INPWR.NE.0) WRITE(MAINF,10060) ITERCHNGFL,ITERCONVFL 
     ;  ,NRITCTNFL,MAXNUMDIVFL,IDMET1FL

10060    FORMAT(/,
     ;'                    ITERATIVE SCHEMES',/,
     ;'                    -----------------',/,
     ;' FLOW EQUATION',/,
     ;' NUMBER OF ITERATIONS OF FIRST METHOD........=',I5,/,
     ;' SECOND METHOD STOPS AFTER ITERATION NUMBER..=',I5,/,
     ;' CHANGE TO FULL NEWTON SCHEME CONV.ITERATION.=',I5,/,
     &' MAXIMUN NUMBER OF CONSECUTIVE DIVERGENCES...=',I5,/,
     ;' IDENTIFIER OF FIRST METHOD..................=',I5 )

C------------------------- check errors flow iterative scheme
       

          IF (ITERCONVFL .LT.ITERCHNGFL) CALL ERROR
     ;       (IERROR,IOWAR,MAINF,FILENAME,
     ;       'ITERCONVFL  MUST BE AT LEAST ITERCHNGFL'
     ;       ,0,0,IUDIM,1,1.23)  

       ELSE  !if we have a linear flow problem

        ITERCONVFL = 1
        ITERCHNGFL = 1
        NRITCTNFL = 0
        MAXNUMDIVFL = 0
        IDMET1FL = 0

       ENDIF

       IF(IOEQT.NE.1 .AND. (IOTRLI.NE.0 .OR. IODENS_INI.EQ.1)) THEN !if we solve transport
C------------------------- read transport iterative scheme

          LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
          READ(LEAUX,10055) ITERCHNGTR,ITERCONVTR,NRITCTNTR,MAXNUMDIVTR
     &                      ,IDMET1TR  

10055     FORMAT(6I5)

C------------------------- if allowed, write to file

         IF (INPWR.NE.0) WRITE(MAINF,10070) ITERCHNGTR,ITERCONVTR
     ;  ,NRITCTNTR,MAXNUMDIVTR,IDMET1TR

10070    FORMAT(/,
     ;' TRANSPORT EQUATION',/,
     ;' NUMBER OF ITERATIONS OF FIRST METHOD........=',I5,/,
     ;' SECOND METHOD STOPS AFTER ITERATION NUMBER..=',I5,/,
     ;' CHANGE TO FULL NEWTON SCHEME CONV.ITERATION.=',I5,/,
     &' MAXIMUN NUMBER OF CONSECUTIVE DIVERGENCES...=',I5,/,
     ;' IDENTIFIER OF FIRST METHOD..................=',I5 )

C------------------------- check errors transport iterative scheme
       

          IF (ITERCONVTR.LT.ITERCHNGTR) CALL ERROR
     ;       (IERROR,IOWAR,MAINF,FILENAME,
     ;       'ITERCONVTR MUST BE AT LEAST ITERCHNGTR'
     ;       ,0,0,IUDIM,1,1.23)

       ELSE

	    ITERCONVTR = 1
	    ITERCHNGTR = 1
	    NRITCTNTR = 0
	    MAXNUMDIVTR = 0
	    IDMET1TR = 0
	   ENDIF

       IF(IODENS_INI.EQ.1) THEN !if we have variable density flow

C------------------------- read coupled problem iterative scheme

          LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
          READ(LEAUX,10065) ITERCHNGGL,NRITMET2GL,IDMET1GL
10065     FORMAT(3I5)

C------------------------- if allowed, write to file

          IF (INPWR.NE.0) WRITE(MAINF,10080) ITERCHNGGL,NRITMET2GL
     ;    ,IDMET1GL

10080     FORMAT(/,
     ;' COUPLED PROBLEM',/,
     ;' NUMBER OF ITERATIONS OF FIRST METHOD........=',I5,/,
     ;' SECOND METHOD STOPS AFTER ITERATION NUMBER..=',I5,/,
     ;' IDENTIFIER OF FIRST METHOD..................=',I5 )


          
C------------------------- check errors global iterative scheme
           
          IF (IDMET1GL.NE.1 .AND. IDMET1GL.NE.2) CALL ERROR
     ;       (IERROR,IOWAR,MAINF,FILENAME,
     ;       'CHOOSE EITHER NEWTON OR PICARD(IDMET1GL 1 OR 2)'
     ;       ,0,0,IUDIM,2,1.23)      

          IF (NRITMET2GL.LT.ITERCHNGGL) CALL ERROR
     ;       (IERROR,IOWAR,MAINF,FILENAME,
     ;       'NRITMET2GL MUST BE AT LEAST ITERCHNGGL'
     ;       ,0,0,IUDIM,1,1.23)

       ELSE

           ITERCHNGGL = 0
           NRITMET2GL = 0
           IDMET1GL = 0

       ENDIF !IF IODENS_INI .EQ.1

C------------------------- watsolv parameters

       IF (IODIRECT.EQ.0) THEN

         CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(13)
     ; ,FILENAME)

         LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
         READ(LEAUX,10090) IWALGO,IWNORTHSTART,IWNORTHMAX,IWPLAN_B,MAXNB
     ;                     ,IPRECOND,LEVEL,MAXNBF     

10090    FORMAT(8I5)

         IF (IPRECOND.GT.0) IPRECOND = 1

C------------------------- write to file if allowed

         IF (INPWR.NE.0) WRITE(MAINF,10100)IWALGO,IWNORTHSTART
     ;     ,IWNORTHMAX,IWPLAN_B,MAXNB,IPRECOND,LEVEL,MAXNBF
     ;   

10100    FORMAT(/,
     ;'             WATSOLV ITERATIVE SOLVER INFORMATION',/,
     ;'             ------- --------- ------ -----------',/,/,
     ;'IDENTIFIER STARTING ALGORITHM....................= ',I5,/,
     ;'INITIAL NR ORTHOGONAL VECTORS....................= ',I5,/,
     ;'MAXIMUM NR ORTHOGONAL VECTORS....................= ',I5,/,
     ;'ALTERNATIVE STRATEGY IDENTIFIER..................= ',I5,/,
     ;'NR OF CONNECTIONS OF MOST CONNECTED NODE.........= ',I5,/,
     ;'USE MATRIX PRECONDITIONING.......................= ',I5,/,
     ;'PRECONDITIONING LEVEL............................= ',I5,/,
     ;'NUMBER OF CONNECTIONS IN PRECONDITIONED MATRIX...= ',I5)

C------------------------- checking for errors


C------------------------- read watsolve convergence criteria

         LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
         READ(LEAUX,10110) RTWOTOL,RMAXTOL,SMAXTOL,NITMAX,IDETAIL
10110    FORMAT(3F10.0,2I5) 

C------------------------- write to file if allowed

         IF (INPWR.NE.0) WRITE(MAINF,10120) RTWOTOL
     ;     ,RMAXTOL,SMAXTOL,NITMAX,IDETAIL
      
10120 FORMAT(/,
     ;'RESIDUAL EUCLIDEAN TWO NORM CONVERGENCE ERROR CRITERIA..= '
     &,F10.5,/,
     ;'RESIDUAL INFINITY (MAXIMUM) NORM CONV.ERROR CRITERIA....= '
     &,F10.5,/,
     ;'SOLUTION UPDATE INFINITY (MAXIMUM) NORM SCALE CRITERIA..= '
     &,F10.5,/,
     ;'MAXIMUM NUMBER OF ITERATIONS TO PERFORM.................= '
     &,I5,/,
     ;'SOLVER PERFORMANCE DETAIL LEVEL.........................= '
     &,I5)

       ELSE

           IWALGO = 0
           IWNORTHSTART = 0
           IWNORTHMAX = 0
           IWPLAN_B = 0
           MAXNB = 0
           IPRECOND = 0
           LEVEL = 0
           MAXNBF = 0
           RTWOTOL = 0D0
           RMAXTOL = 0D0
           SMAXTOL = 0D0
           NITMAX = 0
           IDETAIL = 0

       ENDIF !IODIRECT.EQ.0


C------------------------- Reads IFLAGS Card.
C------------------------- This card is not necesari to be present in DIM file,
C------------------------- so user normally will omit it.

       IERROR_AUX=IERROR
       CALL SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUDIM    ,MAINF    ,NROW     ,NCARD(15)
     ; ,FILENAME)
       IF (IERROR.EQ. -99) THEN
         IERROR=IERROR_AUX
       ELSE
         LEAUX=LEEL(FILENAME,IUDIM,MAINF,NROW,INPWR)
         READ(LEAUX,1620) (IFLAGS(J),J=1,NFLAGS)
       END IF         
 1620  FORMAT(50I2)

C------------------------- Groups some variables in arrays for simplicity in 
C------------------------- arguments lists through the program
C------------------------- Commented lines are variables that are not read yet.

C------------------------- Number of zones
       NZONE_PAR(:) = 0
       NZONE_PAR (1)= NZTRA
       NZONE_PAR (2)= NZSTG
       NZONE_PAR (3)= NZARR
       NZONE_PAR (4)= NZCHP
       NZONE_PAR (5)= NZQQP
       NZONE_PAR (6)= NZALF
       NZONE_PAR (7)= NZDSP
       NZONE_PAR (8)= NZDSP
       NZONE_PAR (9)= NZDFM
       NZONE_PAR(10)= NZPOR
       NZONE_PAR(11)= NZFOD
       NZONE_PAR(12)= NZCRD
       NZONE_PAR(13)= NZCOE
       NZONE_PAR(14)= NZPRG
*       NZONE_PAR(15)= NZAGE
       NZONE_PAR(16)= NZDMT
*       NZONE_PAR(17)= Z.O.D.
       NZONE_PAR(18)= NZCLK
       NZONE_PAR(19)= NZFOF

*_______________________Sets to zero the number of zones of flow param. 
*_______________________if only tpt. is solved

        IF (IOEQT.EQ.2) THEN
           DO I=1,6
              NZONE_PAR(I)=0
           END DO
        END IF

*_______________________Sets to zero the number of zones of tpt. param. 
*_______________________if only FLOW is solved

        IF (IOEQT.EQ.1) THEN
           DO I=7,13
              NZONE_PAR(I)=0
           END DO
           IF (IOFLSAT.NE.0) NZONE_PAR(10)=NZPOR
           NZONE_PAR(16)=0
           NZONE_PAR(18)=0
           NZONE_PAR(19)=0
        END IF

C------------------------- Initializes some variables for mass balance comput.

       IF (IOPTS(17).NE.0.OR.IOPTS(19).NE.0) 

C------------------------- Total number of zones of flow parameters

     ;      NMAXF=NZONE_PAR(2)*(IODENS_INI+1)      ! Storavity
     ;           +NZONE_PAR(3)                     ! Areal recharge
     ;           +NZONE_PAR(4)                     ! Prescribed head
     ;           +NZONE_PAR(5)                     ! Prescribed flow
     ;           +NZONE_PAR(6)                     ! Mixed boun. cond.
     &           +NZONE_PAR(13)*IODENS_INI

       IF (IOPTS(18).NE.0.OR.IOPTS(20).NE.0) 

C------------------------- Total number of zones of transport parameters

     ;     NMAXT=NZONE_PAR(10)                    ! Porosity
     ;          +NZONE_PAR(11)                    ! First order decay reactions
     ;          +NZONE_PAR(17)                    ! Zero order reactions
     &          +NZONE_PAR(13)*6                  ! External concentration
     &          +NZONE_PAR(18)                    ! Concentration Leakage.
                                                  ! (expanded version)

C------------------------- NMAXF and NMAXT are used to dimension, so they
C------------------------- cannot be zero

       IF (NMAXF.EQ.0) NMAXF=1
       IF (NMAXT.EQ.0) NMAXT=1


C------------------------- Output options

       IOWRITE (1)= INPWR
       IOWRITE (2)= IOWAR
       IOWRITE (3)= IOWRH
       IOWRITE (4)= IOWRC
       IOWRITE (5)= IOPLH
       IOWRITE (6)= IOPLC
       IOWRITE (7)= IOMHH
       IOWRITE (8)= IOMHC
       IOWRITE (9)= IOSEH
       IOWRITE(10)= IOSEC
       IOWRITE(11)= IOCMH
       IOWRITE(12)= IOCMC
       IOWRITE(13)= IOWNRFL          
       IOWRITE(14)= IOWPI          
       IOWRITE(15)= IOWIT
       IOWRITE(16)= IOWVD  
       IOWRITE(17)= IOWNRTR
       
C------------------------- Weighting parameters

       PAR_WGT (1)= XLAMTRA
       PAR_WGT (2)= XLAMSTG
       PAR_WGT (3)= XLAMARR
       PAR_WGT (4)= XLAMCHP
       PAR_WGT (5)= XLAMQQP
       PAR_WGT (6)= XLAMALF
       PAR_WGT (7)= XLAMDSP
       PAR_WGT (8)= XLAMDSP
       PAR_WGT (9)= XLAMDFM
       PAR_WGT(10)= XLAMPOR
       PAR_WGT(11)= XLAMFOD
       PAR_WGT(12)= XLAMCRD
       PAR_WGT(13)= XLAMCOE
       PAR_WGT(14)= XLAMPRGF
*       PAR_WGT(15)= XLAMAGE
       PAR_WGT(16)= XLAMPRGT
*       PAR_WGT(17)= XLAMZOD???
       PAR_WGT(18)= XLAMCLK
       PAR_WGT(19)= XLAMFOF
C------------------------- Real optimization parameters 

       PAR_INV (1)= XMARQ
       PAR_INV (2)= PHIMIN
       PAR_INV (3)= PHIMAX
       PAR_INV (4)= GMNOR1
       PAR_INV (5)= GMNOR
       PAR_INV (6)= DMINF
       PAR_INV (7)= COSMIN
       PAR_INV (8)= PERMX1
       PAR_INV (9)= PERMX2
       PAR_INV(10)= EPS

C------------------------- Integer optimization parameters 

       IPAR_INV (1)= NUMIN
       IPAR_INV (2)= NUMAX
       IPAR_INV (3)= MAXICOS
       IPAR_INV (4)= MAXITER
       IPAR_INV (5)= NMTERF1

C------------------------- Logarithmic estimation option

       IOLG_PAR (1,1)= IOLGTRA
       IOLG_PAR (2,1)= IOLGSTG
       IOLG_PAR (3,1)= IOLGARR
       IOLG_PAR (4,1)= IOLGCHP
       IOLG_PAR (5,1)= IOLGQQP
       IOLG_PAR (6,1)= IOLGALF
       IOLG_PAR (7,1)= IOLGDSP
       IOLG_PAR (8,1)= IOLGDSP
       IOLG_PAR (9,1)= IOLGDFM
       IOLG_PAR(10,1)= IOLGPOR
       IOLG_PAR(11,1)= IOLGFOD
       IOLG_PAR(12,1)= IOLGCRD
       IOLG_PAR(13,1)= IOLGCOE
       IOLG_PAR(14,1)= IOLGPRG
*       IOLG_PAR(15,1)= IOLGAGE
*       IOLG_PAR(16,1)= IOLGMATDIF
       IOLG_PAR(18,1)= IOLGCLK
       IOLG_PAR(19,1)= IOLGFOF

C------------------------- Real direct problem parameters (simulation)

       PAR_DIR (1)= FCTNCV
       PAR_DIR (2)= FCTDEC
       PAR_DIR (3)= FCTINC
       PAR_DIR (4)= DRELMXFL
       PAR_DIR (5)= DABSMXFL
       PAR_DIR (6)= RESIDMXF
       PAR_DIR (7)= RESIDMXT
       PAR_DIR (8)= ZEROF
       PAR_DIR (9)= ZEROT
       PAR_DIR(10)= FCTDVNR
       PAR_DIR(11)= DHITMX !TEST FACTOR, =1 ALWAYS 
       PAR_DIR(12)= DCITMX !TEST FACTOR, =1 ALWAYS 
       PAR_DIR(13)= OBJHED1
       PAR_DIR(14)= RESIDMX1F
       PAR_DIR(15)= DABSMX1F
       PAR_DIR(16)= DRELMX1F
       PAR_DIR(17)= OBJHED2
       PAR_DIR(18)= RESIDMX2F
       PAR_DIR(19)= DABSMX2F
       PAR_DIR(20)= DRELMX2F
       PAR_DIR(21)= OBJCON1
       PAR_DIR(22)= RESIDMX1T
       PAR_DIR(23)= DABSMX1T
       PAR_DIR(24)= DRELMX1T
       PAR_DIR(25)= OBJCON2
       PAR_DIR(26)= RESIDMX2T
       PAR_DIR(27)= DABSMX2T
       PAR_DIR(28)= DRELMX2T
*       PAR_DIR(29)= THETAF
*       PAR_DIR(30)= THETAT
*       PAR_DIR(31)= EPSFLU
*       PAR_DIR(32)= EPSTRA
*       PAR_DIR(33)= ERRDMS
       PAR_DIR(34) = DRELMXTR
       PAR_DIR(35) = DABSMXTR
       PAR_DIR(36) = RTWOTOL 
       PAR_DIR(37)= RMAXTOL
       PAR_DIR(38)= SMAXTOL
 

C------------------------- Integer direct problem parameters (simulation)

       IPAR_DIR(1)= MXNRTF
       IPAR_DIR(2)= MXNRTT
       IPAR_DIR(3)= MINCAT
       IPAR_DIR(4)= IOPINITH
       IPAR_DIR(5)= IOPINITC
       IPAR_DIR(6)= IOCRITRAP
*       IPAR_DIR(7)= IOCAP
       IPAR_DIR(8)= ITERCHNGFL 
       IPAR_DIR(9)= ITERCONVFL
       IPAR_DIR(10)= NRITCTNFL
       IPAR_DIR(11)= ITERCHNGTR
       IPAR_DIR(12)= ITERCONVTR
       IPAR_DIR(13)= NRITCTNTR 
       IPAR_DIR(14)= ITERCHNGGL
       IPAR_DIR(15)= IWALGO  
       IPAR_DIR(16)= IWNORTHSTART
       IPAR_DIR(17)= IWNORTHMAX
       IPAR_DIR(18)= IWPLAN_B
       IPAR_DIR(19)= MAXNB
       IPAR_DIR(20)= NITMAX
       IPAR_DIR(21)= IDETAIL
       IPAR_DIR(22)= IPRECOND
       IPAR_DIR(23)= LEVEL
       IPAR_DIR(24)= MAXNBF
       IPAR_DIR(25)= NRITMET2GL
       IPAR_DIR(26)= MAXNUMDIVFL
       IPAR_DIR(27)= MAXNUMDIVTR

c-----------------------iterative scheme parameters
       IF(IODENS_INI.EQ.0 .AND. IOFLLI.EQ.0) IDMET1FL=1
       IF(IODENS_INI.EQ.0 .AND. IOTRLI.EQ.0) IDMET1TR=1

       LINMET(1,1)= IDMET1FL
       LINMET(2,1)= IDMET1TR
       LINMET(3,1)= MAX(IDMET1GL,1)
        



C------------------------- Weighting parameters for state variables 
C------------------------- (concentrations, heads, etc)
       FOBJ_WGT(:)= 0D0
       FOBJ_WGT(1)= XLAMHEAD
       FOBJ_WGT(2)= XLAMCONC
       FOBJ_WGT(3)= XLAMHUMI
       FOBJ_WGT(4)= XLAMFLOW
       FOBJ_WGT(6)= XLAMGEO

*_____________Opens input and output files

       CALL OPENFILES
     ; (IOBALC   ,IOBALGC  ,IOBALGH  ,IOBALH   ,IOCMC    ,IOCMH
     ; ,IOEQT    ,IOFLLI   ,IOINV    ,IOINV_GS ,IOMHC    ,IOMHH
     ; ,IPAR_DIR(5),IPAR_DIR(4),IOPLC,IOPLH    ,IORTS    ,IOSEC    
     ; ,IOSEH    ,IOTRLI   ,IOTRS    ,IOWIT    ,IOWPI    ,IPROCESS 
     ; ,MAINF    ,NBLCVP   ,NUMNP    ,FILENAME)

*_____________Reads geoestatistical/intepolation dimensions and options

       MXLINCMB = 0
       NGROUP_ZN = 0
       IF (IOINV_GS.GT.0) CALL OPT_GROUPS_ZONES
     ;(IERROR   ,INPWR    ,IOINV    ,IOWAR    ,ISOT      ,IUDIM     
     ;,16       ,MAINF    ,MXGRPZN  ,MXLINCMB ,NGROUP_ZN ,NTYPAR
     ;,FILENAME ,IOPT_GS  ,IO_KG_GS)

       RETURN

 9000  CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR WHEN READING NEWTON-RAPHSON ' 
     ;//'CONVERGENCE CRITERIA ARE NOT GOOD GIVEN',NROW,1,IUDIM,2,1.27)
 
       END
