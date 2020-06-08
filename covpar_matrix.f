      SUBROUTINE COVPAR_MATRIX
     ;(IDMBLDIM  ,IERROR     ,INPWR     ,IOWAR      ,IUCOV     ,MAINF
     ;,MXGRPZN    ,NBLCVP    ,NPARDET    ,NZPAR     ,COVPAR    ,FILENAME
     ;,IDMBLCVP  ,IGR_ZONE   ,IOPT_GS   ,IVEST     ,STPAR)

********************************************************************************
*
* PURPOSE Reads or assigns covariance matrix of deterministically estimated 
*         parameters. If parameters are estimated geostatistically, that part
*         of the covariance matrix will be calculated by kriging elsewhere 
*
* DESCRIPTION 
*
*  - Step 0: Declaration of variables
*  - Step 1: Assigns diagonal matrix
*    - Step 1.1: Assigns matrix shape
*    - Step 1.2: Assigns matrix components, if zonal param. is estimated determ.
*  - Step 2: Reads non diagonal matrix
*
*    LOOP OVER DETERMINISTICAL BLOCKS OF COVARIANCE MATRIX
*
*       - Step 2.1: Reads rows defining current block
*       - Step 2.2: Reads current block at covariance matrix
*         Loop over lines of NRECORDS records
*           - Step 2.2.A: Initializes local array
*           - Step 2.2.B: Reads line
*             Loop over records in this line 
*             - Step 2.2.C: Calculates local row and column
*             - Step 2.2.D: Calculates global row and column
*             - Step 2.2.E: Calculates position at COVPAR
*             - Step 2.2.F: Assigns COVPAR
*             - Step 2.2.G: Reads last line (GOTO STEP 2.2.B)
*             - Step 2.2.H: Updates counter of rows read
*  - Step 3: Echoes covariance matrix
*  - Step 4: Inverts covariance matrix
*  - Step 5: Echoes inverted covariance matrix
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVPAR                 Array containing inverse of a priori covariance matrix
*                         of parameters
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  IDMBLCVP               Array containing size of each covariance matrix block
*  IFLAGS                 Used mainly for code debugging
*  IGR_ZONE               Group to which a given zone belongs to
*  IOPT_GS                Array of geostatistical dimensions and options
*  IVEST                  Vector containing estimation index for all            
*                         parameters                                            
*  STPAR                  Vector containing standard deviation errors of        
*                         all parameters prioo information  
*  WORK                   Auxiliar array for inversion purposes                    
*
* INTERNAL VARIABLES: ARRAYS
*
*  AUXREAD                Auxiliar array used only to read from COV.DAT
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDMBLDIM               Used to dimension IDMBLCVP
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUCOV                  COV file unit number
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MXGRPZN                Maximum number of groups of zones
*  NBLCVP                 Number of covariance matrix blocks
*  NFLAGS                 Used to dimension IFLAGS
*  NPARDET                Number of deterministically estimated parameters
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Dummy counter of records (reading purposes)
*  IAUX                   Auxiliar. Position of diagonal term
*  IBGN                   Dummy pointer to initial positions
*  IBL                    Dummy counter for covariance blocks
*  ICONT                  Dummy counter of records
*  IDGT                   Precision at LIN1VP
*  IER                    Error flag of LIN1VP
*  IGROUP                 Group to which a given zone belongs to
*  ILIN                   Dummy counter for lines to be read
*  IOPESTGR               Estimation option of current group
*  IOPESTZON              Estimation option of current zone
*  IPOSCOV                Postion at array COVPAR                   
*  IREC                   Dummy counter of records
*  IROWGLOB               Global position (row) at COVPAR
*  IROWLOC                Local position of a record in this block (row)
*  IZPAR                  Dummy counter of zones                               
*  JCOLGLOB               Global position (column) at COVPAR
*  JCOLLOC                Local position of a record in this block (column)
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*  NFILPREV               Number of rows previously assigned at COVPAR
*  NLINES                 Number of complete lines to be read from COV.DAT
*  NROW                   Current record number                                 
*  NTOTALREC              Total number of records to read at current cov. block
*  NRECORDS               Number of records in a line
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*  LIN1VP                 IMSL routine to invert symmetric matrices
*
* HISTORY:   AAR (first coding), Jan 2003
*            AAR (Revision and inclusion of groups of zones), July 2003
*
********************************************************************************
       
C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NBLCVP,NPARDET,IDMBLDIM,NZPAR,MXGRPZN,IUCOV,MAINF
     ;         ,INPWR,IERROR,IOWAR
     ;         ,IDMBLCVP(IDMBLDIM),IVEST(NZPAR),IGR_ZONE(NZPAR)
     ;         ,IOPT_GS(MXGRPZN,20)
                                                                 ! Real external
      REAL*8 STPAR(NZPAR),COVPAR(NPARDET*(NPARDET+1)/2)
                                                              ! Integer internal
      INTEGER*4 IBGN,IZPAR,IPOSCOV,IGROUP,IOPESTGR,IOPESTZON,IBL
     ;         ,NROW,NTOTALREC,NLINES,NRECORDS,ICONT,I,ILIN,IREC
     ;         ,IROWLOC,JCOLLOC,IAUX,NFILPREV,IROWGLOB,JCOLGLOB
                                                                 ! Real internal
      REAL*8 AUXREAD(8)
                                                                    ! Characters
      CHARACTER FILENAME(20)*20,LEAUX*100,LEEL*100

C_______________________ Step 1: Assigns diagonal matrix

      IF (NBLCVP.EQ.0) THEN

C_______________________ Step 1.1: Assigns matrix shape

        DO IBL=1,IDMBLDIM
          IDMBLCVP(IBL)=1
        END DO

C_______________________ Step 1.2: Assigns matrix components, if zonal param. is
C_______________________           estimated deterministically
 
        IBGN=0
        DO IZPAR=1,NZPAR

          IGROUP=IGR_ZONE(IZPAR)      ! Group to which zone belongs to
          IOPESTZON=IVEST(IZPAR)      ! Zone estimation option

          IF (IOPESTZON.GT.0 .AND. IGROUP.GT.0) THEN 
            IOPESTGR=IOPT_GS(IGROUP,2)  ! Group estimation option
            IF (IOPESTGR.EQ.0) THEN ! Deterministical estim.
 
              IBGN=IBGN+1
              IPOSCOV=IBGN*(IBGN+1)/2


* Before:     COVPAR(IPOSCOV)=STPAR(IZPAR)*STPAR(IZPAR)

*** The inverse of the covariance matrix is directly assigned

              COVPAR(IPOSCOV)=1.0D0/STPAR(IZPAR)/STPAR(IZPAR)
            END IF
          END IF ! IOPESTZON.GT.0 .AND. IOPESTGR.EQ.0

        END DO ! IZPAR=1,NZPAR

      ELSE

C_______________________ Step 2: Reads non diagonal matrix (inverse)

        NFILPREV=0
        DO IBL=1,NBLCVP     ! Loop over deterministical blocks

C_______________________ Step 2.1: Reads rows defining current block

          LEAUX=LEEL(FILENAME,IUCOV,MAINF,NROW,INPWR)
          READ(LEAUX,1000,ERR=9000) IDMBLCVP(IBL)
 1000     FORMAT(I5)

C_______________________ Step 2.2: Reads current block at covariance matrix

          NTOTALREC=IDMBLCVP(IBL)*(IDMBLCVP(IBL)+1)/2   ! Records at block
          NLINES=INT(NTOTALREC/8)                       ! Complete lines
          NRECORDS=8                                    ! Records at line
          ICONT=0                              ! Init. counter of records
          
C_______________________ Loop over lines of NRECORDS records

 10       DO ILIN=1,NLINES

C_______________________ Step 2.2.A: Initializes local array

            CALL ZERO_ARRAY(AUXREAD,8)

C_______________________ Step 2.2.B: Reads line

            LEAUX=LEEL(FILENAME,IUCOV,MAINF,NROW,INPWR)
            READ(LEAUX,1100,ERR=9100) (AUXREAD(I),I=1,NRECORDS)
 1100       FORMAT(8F10.0)

            DO IREC=1,NRECORDS   ! Loop over records in this line 
              ICONT=ICONT+1      ! Updates number of assigned records 

C_______________________ Step 2.2.C: Calculates local row and column

              IROWLOC=
     ;           INT(0.5D0*(-1D0+DSQRT(1.D0+8.D0*FLOAT(ICONT)))+0.1D0)+1
              IAUX=IROWLOC*(IROWLOC-1)/2 
              IF (IAUX.EQ.ICONT) IROWLOC=IROWLOC-1

              JCOLLOC=ICONT-IAUX
              IF (IAUX.EQ.ICONT) JCOLLOC=JCOLLOC+IROWLOC
              
C_______________________ Step 2.2.D: Calculates global row and column

              IROWGLOB=IROWLOC+NFILPREV
              JCOLGLOB=JCOLLOC+NFILPREV

C_______________________ Step 2.2.E: Calculates position at COVPAR

              IPOSCOV=IROWGLOB*(IROWGLOB-1)/2+JCOLGLOB

C_______________________ Step 2.2.F: Assigns COVPAR

              COVPAR(IPOSCOV)=AUXREAD(IREC)

            END DO ! IREC=1,NRECORDS

          END DO ! ILIN=1,NLINES

C_______________________ Step 2.2.G: Reads last line

          IF (NRECORDS.EQ.8) THEN
            NRECORDS=NTOTALREC-8*NLINES
            NLINES=1
            GOTO 10
          END IF ! NLINES.GT.1

C_______________________ Step 2.2.H: Updates counter of rows read

          NFILPREV=NFILPREV+IDMBLCVP(IBL)

        END DO ! IBL=1,NBLCVP

      END IF ! NBLCVP.EQ.IPARDET

C_______________________ Step 3: Echoes covariance matrix

      IF (INPWR.NE.0) THEN

         WRITE(MAINF,2000)
 2000    FORMAT(//,6X,'INVERSE OF THE A PRIORI COVARIANCE MATRIX'
     ;                ' OF DETERMINISTICALLY ESTIMATED PARAMETERS',/
     ;            ,6X,'======= == === = ====== ========== ======'
     ;                ' == ================= ========= ==========',/)

         WRITE(MAINF,2100) COVPAR
 2100    FORMAT(7(1X,E10.3))

      END IF ! INPWR.NE.0


*** Not needed; old algorithm

C_______________________ Step 4: Inverts covariance matrix (ONLY IF GEOSTAT.
C_______________________         IS NOT DONE). If this is the case, inversion 
C_______________________         of full covariance matrix will be done elsewhere

*      IF (IOINV_GS.EQ.0) THEN
 
*        CALL EQUAL_ARRAY (WORK,COVPAR,NPARDET*(NPARDET+1)/2)
*        CALL LINV1P (WORK,NPARDET,COVPAR,IDGT,D1,D2,IER)
      
*        IF (IER.EQ.129) THEN
*           WRITE(MAINF,2200)
* 2200      FORMAT(//,' ERROR INVERTING A PRIORI COVARIANCE MATRIX OF'
*     ;             ' DETERMINISTICALLY ESTIMATED PARAMETERS.',/,' IT IS'
*     ;             ' NOT POSITIVE DEFINITE. FORCED STOP, SORRY',/)
*           STOP ' CRITICAL STOP. CHECK FILE RES.OUT'
*        END IF

C_______________________ Step 5: Echoes inverted covariance matrix

*        IF (IFLAGS(22).NE.0) THEN
*           WRITE(MAINF,2300)
* 2300      FORMAT(//,22X,'INVERTED A PRIORI COVARIANCE MATRIX',/
*     ;              ,22X,'======== = ====== ========== ======',/)

*           WRITE(MAINF,2100) COVPAR
*       END IF

*      END IF ! IOINV_GS.EQ.0

      RETURN
       
 9000 CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;      'GENERIC FORTRAN ERROR READING '//
     ;      'SHAPE OF COV. MATRIX',NROW,0,IUCOV,2,6.16)

      RETURN
       
 9100 CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;      'GENERIC FORTRAN ERROR READING '//
     ;      'PARAMETERS COVARIANCE MATRIX',NROW,0,IUCOV,2,6.16)

      RETURN

      END

