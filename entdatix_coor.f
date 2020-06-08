       SUBROUTINE ENTDATIX_COOR
     ; (IERROR   ,IDALF    ,IDALFT   ,IDCHP    ,IDCHPT   ,IDCON
     ; ,IDCONT   ,IDDMT    ,INPWR    ,IDQQP    ,IDQQPT   ,IOEQT    
     ; ,IOFLLI   ,IOTRLI   ,IOWAR    ,IUGRID   ,MAINF
     ; ,NROW     ,NTDMT    ,NUMNP    ,BTRA     ,CCAL     
     ; ,FILENAME ,HCAL     ,X        ,Y        ,Z        ,IDCLK)

*******************************************************************************
* PURPOSE 
*
*    Reads and checks data cards form GRI file corresponding to
*    nodal data: Cards groups B1 and B2. 
*
* DESCRIPTION
*     
*    This subroutine reads and checks nodal data form cards groups 
*    B1 and B2 on the following way:  
*
*         * Reads the current card trough character function LEEL,which
*           returns the LEAUX string
*         * Load the current card variables by reading LAUX string 
*         * Checks the read data and writes error messages on MAIN FILE
*           calling subroutine ERROR
*         * Write in MAIN FILE interpreted data, after each card is 
*           read.
*         * Data read from card 1.2: 
*            - Nodal coordinates are linearly interpolated for missing nodes
*            - BTRA,HCAL,CCAL, are defined for missing nodes through subroutine
*              ASS_NVAL for missing nodes before writing on MAIN file
*         * Data read from card 1.3 (Bounday conditions and parameter zones) 
*           is checked and missing nodes are assigned default values or those
*           read in the last node (if corresponging default value in card 1.1 
*           equals 0) through subroutine INTERP_NZNUM
*
* EXTERNAL VARIABLES: ARRAYS 
*     
*  BTRA                   Right hand side of transport discretized equation     
*  CCAL                   Used temporarily to store base concentration (CBASE)
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  HCAL                   Used temporarily to store base head (HBASE)
*  IBCOD                  Flow boundary condition index                         
*  IBTCO                  Transport boundary condition index                    
*  IXPARNP                Array containing zone number of each node j,          
*                         corresponding to INpar index parameter zone.          
*  X                      X-coord for a given node                              
*  Y                      Y-coord for a given node                              
*  Z                      Z-coord for a given node                              
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INALF                  Index for leackage                                    
*                         in array variables (IXPARNP and CFPARNP)              
*  INALFT                 Index for transient leackage                          
*                         in array variables (IXPARNP and CFPARNP)              
*  INCHP                  Index for prescribed head                             
*                         in array variables (IXPARNP and CFPARNP)              
*  INCHPT                 Index for transient prescribed head                   
*                         in array variables (IXPARNP and CFPARNP)              
*  INCON                  Index for external concentration                      
*                         in array variables (IXPARNP and CFPARNP)              
*  INCONT                 Index for transient external concentration            
*                         in array variables (IXPARNP and CFPARNP)              
*  INDMT                  Index for matrix diffusion                            
*                         in array variables (LXPAREL and CFPAREL)              
*  INPWR                  Allows writing on MAIN FILE                           
*  INQQP                  Index for prescribed flow                             
*                         in array variables (IXPARNP and CFPARNP)              
*  INQQPT                 Index for tramsient prescribed flow                   
*                         in array variables (IXPARNP and CFPARNP)              
*  IOEQT                  Type of problem to be solved                          
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOTRLI                 Idem to IOFLLI, in the case of transport.             
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUGRID                 Unit number of GRI file
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NTDMT                  If diferent form zero, number of terms in matrix      
*                         diffusion(if zero, no diffusion)                      
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS 
*
*  IN                     Current node number in read loop
*  C                      from GRI file      
*  AC                     Same as BTRANS(IN),but used when reading these 
*                         value from GRI file             
*  HB1                    Same as HCAL(IN),but used when reading these value
*                         from GRI file             
*  CB1                    Same as CCAL(IN),but used when reading these value
*                         from GRI file             
*  IDCHP                  Zone number default value for prescribed head in
*                         steady-state
*  IDCHPT                 Same as above for transient-state
*  IDQQP                  Zone number default value for prescribed flow in
*                         steady-state
*  IDQQPT                 Same as above for transient-state
*  IDALF                  Zone number default value for leackance in
*                         steady-state
*  IDALFT                 Same as above for transient-state
*  IDCON                  Zone number default value for external 
*                         concentration in steady-state
*  IDCONT                 Same as above for transient-state
*  IDDMT                  Zone number default value for matrix diffusion in
*                         steady-state
*  DFTACTH                Aquifer thickness default value 
*  DHBAS                  Default value for the bottom of the aquifer 
*  DCBAS                  Default value for concentration at the bottom of 
*                         the aquifer  
*  AUXAC                  Auxiliary var. containing BTRA value for missing 
*                         nodes to be written in MAIN file.       
*  AUXCCAL                Idem for HCAL.
*  AUXHCAL                Idem for HCAL.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ASS_NVAL               Assigns a generic node default value (def. value 
*                         non 0)or this corresponding to the last node read.
*  INTERP_NZNUM           Interpolates zone number for missing nodes       
*                       
*  ERROR                  Writes the current error message and error number 
*                         on MAIN FILE.
*  LEEL                   Returns a string value containing the current line
*                         of an input file, if no coment appears.
*
* HISTORY
*
*     AMS        1988     First coding
*     SCR      4-1997     Revision and verification
*     AMS      1-1998     Common elimination (and tabs!)
*     AMS      4-2000     Format correction and splitting in two routines
*                         (coordinates and zones)
*
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*100 LEEL,LEAUX,FILENAME(20)*20
         
       DIMENSION X(NUMNP),Y(NUMNP),Z(NUMNP)
     ;          ,BTRA(NUMNP),HCAL(NUMNP),CCAL(NUMNP)

C------------------------- FIRST EXECUTABLE STATEMENT.
C------------------------- Reads coordinates

*_______________________Reads defaults values for nodal parameters zone numbers
*_______________________Card B1.1

       LEAUX=LEEL(FILENAME,IUGRID,MAINF,NROW,INPWR)
       READ(LEAUX,1000,ERR=9000) IDCHP,IDCHPT,IDQQP,IDQQPT,IDALF,IDALFT,
     ;                    IDCON,IDCONT,IDDMT,DFTACTH,DHBAS,DCBAS,IDCLK
 1000  FORMAT(9I5,3F10.0,I5)

*_______________________Write nodal default values on MAIN file

       IF (INPWR.NE.0) WRITE(MAINF,3000) IDCHP,IDCHPT,
     ; IDQQP,IDQQPT,IDALF,IDALFT,IDCON,IDCONT,IDDMT,IDCLK,DFTACTH,DHBAS
     &,DCBAS

 3000   FORMAT(////,
     ;        1X,19('*'),' GRID  INFORMATION ',19('*'),////,
     ;        10X,'NODAL POINT INFORMATION'/,10X,23('='),//,
     ;  5X,'DEFAULT VALUES',/,5X,14('-'),/
     ;  5X,'PRESCRIBED HEAD (S.S.) [CHP]......=',I5,2X,'|',/,
     ;  5X,'PRESCRIBED HEAD (TRANSIENT) [CHPT]=',I5,2X,'|',/,
     ;  5X,'PRESCRIBED FLOW (S.S.) [QQP]......=',I5,2X,'|',/,
     ;  5X,'PRESCRIBED FLOW (TRANSIENT) [QQPT]=',I5,2X,'|',2X,
     ;  'ZONE NUMBERS',/,
     ;  5X,'LEAKAGE COEFFICIENT [ALF].........=',I5,2X,'|',/,
     ;  5X,'LEAKAGE (TRANSIENT) COEFF. [ALFT].=',I5,2X,'|',/,
     ;  5X,'EXTERNAL CONCENT.(S.S.) [CON].....=',I5,2X,'|',/,
     ;  5X,'EXTERNAL CONCENT.(TRANS.) [CONT]..=',I5,2X,'|',/,
     ;  5X,'MATRIX DIFUSSION.[DMT]............=',I5,2X,'|',/,
     &  5X,'CONC. LEAK. (TRANS.) COEFF. [CLK].=',I5,2X,'|',/,
     ;  5X,'AQUIFER THICKNESS ................=',F10.2,/,
     ;  5X,'MINIM HEAD ACCEPTED...............= ',E10.4,/,
     ;  5X,'MINIM CONC.ACCEPTED...............= ',E10.4,///,
     ;    ' NODES          COORDINATES           THICKNESS.',/,
     ;    1X,5('-'),10X,11('-'),12X,9('-'),/, 
     ;         '    N         X         Y         Z       AC.',
     ;              '        HBASE     CBASE')

*_______________________Check matrix diffusion option

       IF (NTDMT.EQ.0. AND .IDDMT.NE.0)
     ;        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,'INCORRECT DEFAULT' 
     ;      //'OF MATRIX DIFUSION',1,8,IUGRID,1,2.24)

       AUXHCAL=0.D0
       AUXCCAL=0.D0

*_______________________Starts nodal coordinates input loop
*_______________________Card B1.2

       I=0
       IN=0
       DO WHILE (IN.LT.NUMNP)
          IOLD=I
          I=I+1
          LEAUX=LEEL(FILENAME,IUGRID,MAINF,NROW,INPWR)      
          READ(LEAUX,1100,ERR=9100) IN,X1,Y1,Z1,AC,HB1,CB1

1100      FORMAT(I5,6F10.0)
          IF (IN.GT.NUMNP) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;         'NODE NUMBER IS TOO LARGE',NROW,1,IUGRID,1,2.01)

*_______________________Check increasing order of nodes

          IF (IOLD.GE.IN) CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;       'INCORRECT ORDER IN NODE NUMBERING ',NROW,1,IUGRID,1,2.02)

*_______________________Check the aquifer thickness

          IF (AC.LT.0) CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;       'NEGATIVE THICKNESS',NROW,13,IUGRID,1,2.23)

*_______________________Check base values

          IF (CB1.LT.0)
     ;       CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;            'CONC. BASE HAVE TO BE GREATER THAN ZERO',
     ;            0,0,IUGRID,1,2.25)

*_______________________Assigns values to coordinates vector

          X(IN)=X1
          Y(IN)=Y1
          Z(IN)=Z1

*_______________________Assigns defaults values to BTRA,HCAL,CCAL

          IF (IOEQT.NE.1)  CALL ASS_REAL_VAL(AC ,DFTACTH ,BTRA(IN))
          IF (IOFLLI.NE.0) CALL ASS_REAL_VAL(HB1,DHBAS   ,HCAL(IN))
          IF (IOTRLI.NE.0) CALL ASS_REAL_VAL(CB1,DCBAS   ,CCAL(IN))

*_______________________Interpolates nodal coordinates
      
          IF (I.NE.IN) THEN
             NN=IN-IOLD
             XINC=(X1-X(IOLD))/DFLOAT(NN)
             YINC=(Y1-Y(IOLD))/DFLOAT(NN)
             ZINC=(Z1-Z(IOLD))/DFLOAT(NN)
             DO II=1,NN-1
                X(I)=X(IOLD)+II*XINC
                Y(I)=Y(IOLD)+II*YINC
                Z(I)=Z(IOLD)+II*ZINC
                IF (IOEQT.NE.1) BTRA(I)=BTRA(IOLD)
                IF (IOFLLI.NE.0)THEN
                   HCAL(I)=HCAL(IOLD)
                   AUXHCAL=HCAL(I)
                ENDIF
                IF (IOTRLI.NE.0)THEN
                   CCAL(I)=CCAL(IOLD)
                   AUXCCAL=CCAL(I)
                ENDIF
                IF (INPWR.NE.0) THEN
                   IF (IOEQT.NE.1) AUXAC=BTRA(I)
                   WRITE(MAINF,3100) I,X(I),Y(I),Z(I),AUXAC,
     ;                               AUXHCAL,AUXCCAL
 3100              FORMAT(I5,6(1X,G10.3))
                ENDIF
                I=I+1
             ENDDO
          ENDIF
        
          IF (IOFLLI.NE.0) AUXHCAL=HCAL(I)
          IF (IOTRLI.NE.0) AUXCCAL=CCAL(I)

          IF (INPWR.NE.0) THEN
             IF (IOEQT.NE.1) AUXAC=BTRA(IN)
             WRITE(MAINF,3100) I,X(IN),Y(IN),Z(IN),AUXAC,
     ;                         AUXHCAL,AUXCCAL
          ENDIF

       ENDDO        ! NEXT NODE

       RETURN

C------------------------- Error messages
 
 9000       CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;          'GENERIC FORTRAN ERROR WHEN READING'
     ;    //' DEFAULT NODAL DATA (CARD B1.1) ',
     ;          NROW,0,IUGRID,2,3.15)

      RETURN

 9100       CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;          'GENERIC FORTRAN ERROR WHEN READING '
     ;    //'NODAL COORDINATES (CARD B1.2)',
     ;           NROW,0,IUGRID,2,3.15)

      RETURN
      END
