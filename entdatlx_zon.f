       SUBROUTINE  ENTDATLX_ZON
     ; (IERROR   ,INPWR    ,IOEQT    ,IOFLSAT  ,IOTRS    ,IOWAR
     ; ,IPROB    ,IUGRID   ,LDARR    ,LDARRT   ,LDCOE    ,LDCRD
     ; ,LDDFM    ,LDDSP    ,LDFOD    ,LDPOR    ,LDSTG    ,LDTRA
     ; ,MAINF    ,NROW     ,NUMEL    ,NZARR    ,NZCOE    ,NZCRD
     ; ,NZDFM    ,NZDSP    ,NZFOD    ,NZPOR    ,NZSTG    ,NZTRA
     ; ,FILENAME ,LDIM     ,LTYPE    ,LXARR    ,LXARRT   ,LXCOE
     ; ,LXCRD    ,LXDFM    ,LXDSP    ,LXFOD    ,LXPOR    ,LXSTG
     ; ,LXTRA)

********************************************************************************
* PURPOSE 
*
*    Reads and checks data cards form GRI file corresponding to zones of
*    element parameters: Card group B3.3
*
* DESCRIPTION
*
*    This subroutine reads and checks element data form card group B3.3
*    on the following way:  
*
*         * Reads the current card trough character function LEEL,which
*           returns the LEAUX string
*         * Load the current card variables by reading LAUX string 
*         * Checks the read data and writes error messages on MAIN FILE
*           calling subroutine ERROR
*         * Write in MAIN FILE interpreted data, after each card is 
*           read.
*
*         * Element parameter zones is checked and 
*           missing nodes are assigned default values (card B3.1) or those read
*           in the last node (if corresponging default value in card B3.1
*           equals 0) through subroutine INTERP_EZNUM.
* 
* EXTERNAL VARIABLES: ARRAYS 
*     
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  LDIM                   Vector containing fisical dimension of j-th element   
*  LTYPE                  Vector containing type for element j                  
*  LXARR                  Areal recharge (steady) zone number at a given element
*  LXARRT                 Areal recharge (trans.) zone number at a given element
*  LXCOE                  External concentration zone number at a given element 
*  LXCRD                  Retardation zone number at a given element            
*  LXDFM                  Molecular diffusion zone number at a given element    
*  LXDSP                  Dispersivity zone number at a given element           
*  LXFOD                  First order decay zone number at a given element      
*  LXPOR                  Porosity zone number at a given element               
*  LXSTG                  Storage coefficient zone number at a given element    
*  LXTRA                  Transmissivity zone number at a given element         
*
* EXTERNAL VARIABLES: SCALARS  
*
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOEQT                  Type of problem to be solved                          
*  IOFLSAT                Indicates the possibility that one part of the domain 
*                         reaches unsaturated state.                            
*  IOTRS                  Flow regime                                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUGRID                 Unit number of GRI file                               
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NROW                   Current record number                                 
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZARR                  Number of areal recharge zones                        
*  NZCOE                  Number of external concentration zones                
*  NZCRD                  Number of retardation Coefficient zones               
*  NZDFM                  Number of molecular difusion zones                    
*  NZDSP                  Number of dispersivity zones                          
*  NZFOD                  Number of zones of first order decay                  
*  NZPOR                  Number of porosity zones                              
*  NZSTG                  Number of storage Coefficient zones                   
*  NZTRA                  Number of transmissivity zones                        
*
* INTERNAL VARIABLES: ARRAYS 
*
*
* INTERNAL VARIABLES: SCALARS 
*
*  LDTRA                  Default value for transmisivity zone number
*  LDARR                  Default value for areal recharge
*  LDARRT                 Default value for transient areal recharge
*  LDSTG                  Default value for storage coeficient    
*  LDDSP                  Default value for dispersivity
*  LDDFM                  Default value for matrix diffusion
*  LDCOE                  Default value for external concentration
*  LDPOR                  Default value for porosity
*  LDCRD                  Default value for retardation coefficient    
*  LTRA                   Read value for transmisivity zone number
*  LARR                   Read value for areal recharge
*  LARRT                  Read value for transient areal recharge
*  LSTG                   Read value for storage coeficient    
*  LDSP                   Read value for dispersivity
*  LDFM                   Read value for matrix diffusion
*  LCOE                   Read value for external concentration
*  LPOR                   Read value for porosity
*  LCRD                   Read value for retardation coefficient    
*  LOLTRA                 Transmissivity zone number. Auxiliary var.
*  LOLARR                 Areal recharge zone number. Auxiliary var.
*  LOLARRT                Transient areal recharge zone number. Auxiliary var.
*  LOLSTG                 Storage coefficient zone number. Auxiliary var.
*  LOLDSP                 Dispersivity zone number. Auxiliary var.
*  LOLDFM                 Matrix diffusion zone number. Auxiliary var.
*  LOLCOE                 External concentration zone number. Auxiliary var.
*  LOLPOR                 Porosity zone number. Auxiliary var.
*  LOLCRD                 Retardation coeficient zone number. Auxiliary var.
*  NOLD                   Index containing last element read.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ASS_EVAL               Assigns default values for parameters zone number
*  ERROR                  Writes the current error message and error number     
*  INTERP_EZNUM           Interpolates zone numbers for missing elements
*  LDIMEN                 Computes the dimension of a given element             
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*  VERIFY_BW              Checks if the input bandwidth is coherent with the
*                         actual one
*
* HISTORY
*
*     AM         1988     Firs coding
*     SCR      4-1997     Revision and verification
*     AMS      1-1998     Common elimination
*     AMS      4-1999     Correction of the dimension of LXFOD array
*
********************************************************************************


       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*100 LEEL,LEAUX,FILENAME(20)*20  

       DIMENSION LTYPE(NUMEL),
     ;           LDIM(NUMEL),
     ;      LXTRA(NUMEL),LXARR(NUMEL),LXARRT(NUMEL),
     ;      LXSTG(NUMEL),LXDSP(NUMEL),LXDFM(NUMEL),LXCOE(NUMEL),
     ;      LXPOR(NUMEL),LXCRD(NUMEL),LXFOD(NUMEL)

C------------------------- FIRST EXECUTABLE STATEMENT.

C------------------------- Initializes auxiliary variables

       LOLTRA=0
       LOLSTG=0
       LOLARR=0
       LOLARRT=0
       LOLDSP=0
       LOLDFM=0
       LOLPOR=0
       LOLCRD=0
       LOLCOE=0
       LOLFOD=0

C------------------------- Prepare element zone labes for writting 
C------------------------- on MAIN file

       IF (INPWR.NE.0) THEN
          IF (IPROB.EQ.1) THEN
             WRITE(MAINF,3200)
 3200        FORMAT(//,5X,'Z O N E  N U M B E R S',/,5X,44('-'),//,
     ;          '  N.EL. TRA  STG  ARR  ART  ',
     ;             'DSP  DFM  POR  CRD  COE  FOD'/)
          ELSE
             WRITE(MAINF,3300) IPROB
 3300        FORMAT(3X,' PROBLEM NUMBER',I5)
          ENDIF
       ENDIF

       NOLD=0
       N=0               ! ANTES
       NE=0              !  ESTA NO EXISTIA
       DO WHILE (NE.LT.NUMEL)
          NOLD=NE           ! NOLD=NOLD+1
          N=NOLD+1          ! N=N+1

C------------------------- Reads definition of parameter zones       
C------------------------- Card B3.3

          LEAUX=LEEL(FILENAME,IUGRID,MAINF,NROW,INPWR)
          READ(LEAUX,1005,ERR=9200) NE,LTRA,LSTG,LARR,
     ;        LARRT,LDSP,LDFM,LPOR,LCRD,LCOE,LFOD
 1005     FORMAT(11I5)
 
          IF (NE.EQ.0) RETURN

C------------------------- Flow equation parameters

          IF (IOEQT.NE.2) THEN  

C------------------------- Checks transmissivity zone zumber


             IF (LTRA.GT.NZTRA .OR. LTRA.LE.0)
     ;        CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;       'INCORRECT TRANSMISSIVITY ZONE NUMBER (LXTRA) ',
     ;        NROW,5,IUGRID,1,3.06)

             CALL ASS_EVAL(NZTRA,LDTRA,LTRA,LXTRA(NE))

C------------------------- Assigns physical dimension to read transm. zone

             NZTR=LXTRA(NE)
             LDIM(NZTR)=LDIMEN(LTYPE(NE))  !Temporary use of LDIM variable

C------------------------- Checks zone number for storage and recharge
 
             IF (IOTRS.GE.1) THEN  

C------------------------- Checks storage zone number

                IF (LSTG.GT.NZSTG .OR. LSTG.LE.0) 
     ;            CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;            'INCORRECT STORAGE ZONE NUMBER (LXSTG) ',
     ;             NROW,6,IUGRID,1,3.07)
     
                CALL ASS_EVAL(NZSTG,LDSTG,LSTG,LXSTG(NE))            

             ELSE  

C------------------------- In Steady state, storage coefficient zone should be zero

                IF (LSTG.NE.0) 
     ;             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;                     'IOTRS IS ZERO AND LXSTG'
     ;                     //'IS NOT ZERO',NROW,6,IUGRID,0,3.19)

             ENDIF

C------------------------- Checks recharge zone number

             IF (LARR.GT.NZARR .OR. LARR.LT.0)
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;        'INCORRECT RECHARGE ZONE NUMBER (LXARR) '
     ;         ,NROW,8,IUGRID,1,3.09)

             CALL ASS_EVAL(NZARR,LDARR,LARR, LXARR(NE))

C------------------------- Checks transient recharge zone number

                IF ((NZARR.NE.0).AND.(LARRT.GT.NZARR.OR.LARRT.LT.0))
     ;           CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;           'INCORRECT TRANSIENT RECHARGE ZONE NUMBER (LXARRT) '
     ;            ,NROW,8,IUGRID,1,3.08)

                CALL ASS_EVAL(NZARR,LDARRT,LARRT,LXARRT(NE))

          END IF   ! IOEQT.NE.2

C------------------------- Transport equation parameters
 
          IF (IOEQT.NE.1.OR.IOFLSAT.NE.0) THEN
       
             IF (IOEQT.NE.1) THEN

C------------------------- Checks dispersivity zone number

                IF (NZDSP.NE.0 .AND. (LDSP.GT.NZDSP .OR. LDSP.LE.0) )
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;              'DISPERSIVITY ZONE NUMBER IS OUT OF BOUNDS '//
     ;              '(LXDSP)',NROW,9,IUGRID,1,3.10)

                CALL ASS_EVAL(NZDSP,LDDSP,LDSP,LXDSP(NE))

C------------------------- Checks diffusion zone number

                IF ((NZDFM.NE.0).AND.(LDFM.GT.NZDFM .OR. LDFM.LE.0)) 
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;         'DIFFUSION ZONE NUMBER IS OUT OF BOUNDS (LXDFM)',
     ;          NROW,10,IUGRID,1,3.11)

                CALL ASS_EVAL(NZDFM,LDDFM,LDFM,LXDFM(NE))

C------------------------- Checks porosity zone number

             ENDIF
        
             IF (LPOR.GT.NZPOR .OR. LPOR.LT.0)
     ;       CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;       'POROSITY ZONE NUMBER IS OUT OF BOUNDS (LXPOR) ',
     ;        NROW,11,IUGRID,1,3.12)


             CALL ASS_EVAL(NZPOR,LDPOR,LPOR,LXPOR(NE))
          
             IF (IOEQT.NE.1) THEN 

C------------------------- Checks retardation zone number

                IF (LCRD.GT.NZCRD .OR. LCRD.LT.0)
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;       'RETARDATION COEFF. ZONE NUMBER IS OUT OF BOUNDS (LXCRD)',
     ;        NROW,12,IUGRID,1,3.13)

                CALL ASS_EVAL(NZCRD,LDCRD,LCRD,LXCRD(NE))

C------------------------- Checks external concentration zone number

                IF (LCOE.GT.NZCOE.OR.LCOE.LT.0)
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;         'EXTERNAL CONCENTRATION ZONE NUMBER IS '//
     ;         'OUT OF BOUNDS (LXCOE)',
     ;          NROW,13,IUGRID,1,3.14)

                CALL ASS_EVAL(NZCOE,LDCOE,LCOE,LXCOE(NE))
          
C------------------------- Checks first order decay coef. zone number

                IF (LFOD.GT.NZFOD.OR.LFOD.LT.0)
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;         'FIRST ORDER DECAY COEF. ZONE NUMBER IS '//
     ;         'OUT OF BOUNDS (LXCOE)',
     ;          NROW,13,IUGRID,1,3.14)

                CALL ASS_EVAL(NZFOD,LDFOD,LFOD,LXFOD(NE))
             END IF
          END IF  ! IOEQT.NE.1.OR.IOFLSAT.NE.0

C------------------------- Zone number interpolation for 
C------------------------- missing elements

          IF (N.NE.NE) CALL INTERP_EZNUM
     ; (INPWR    ,IOEQT    ,IOFLSAT  ,IOTRS    ,LDARR    ,LDARRT
     ; ,LDCOE    ,LDCRD    ,LDDFM    ,LDDSP    ,LDFOD    ,LDPOR
     ; ,LDSTG    ,LDTRA    ,MAINF    ,N        ,NE       ,NOLD
     ; ,NUMEL    ,NZARR    ,NZCOE    ,NZCRD    ,NZDFM    ,NZDSP
     ; ,NZFOD    ,NZPOR    ,NZSTG    ,NZTRA    ,LDIM     ,LTYPE
     ; ,LXARR    ,LXARRT   ,LXCOE    ,LXCRD    ,LXDFM    ,LXDSP
     ; ,LXFOD    ,LXPOR    ,LXSTG    ,LXTRA)

C------------------------- Defines auxiliary variables to write
C------------------------- last read element zone numbers.

          IF (INPWR.NE.0) THEN
             IF (IOEQT.NE.2) THEN
                LOLTRA=LXTRA(NE)
                IF (NZARR.NE.0) LOLARR=LXARR(NE)
                IF (NZARR.NE.0) LOLARRT=LXARRT(NE)
                IF (IOTRS.NE.0) THEN
                   IF (NZSTG.NE.0) LOLSTG=LXSTG(NE)
                END IF
             END IF
             IF (IOEQT.NE.1) THEN
                IF (NZDSP.NE.0) LOLDSP=LXDSP(NE)
                IF (NZDFM.NE.0) LOLDFM=LXDFM(NE)
                IF (NZPOR.NE.0) LOLPOR=LXPOR(NE)
                IF (NZCRD.NE.0) LOLCRD=LXCRD(NE)
                IF (NZCOE.NE.0) LOLCOE=LXCOE(NE)
                IF (NZFOD.NE.0) LOLFOD=LXFOD(NE)
             ELSE IF(IOFLSAT.NE.0)THEN
                IF (NZPOR.NE.0) LOLPOR=LXPOR(NE)
             ENDIF             
             
C------------------------- Writes last read element zone numbers 
C------------------------- on MAIN file

             WRITE(MAINF,1000) NE,LOLTRA,LOLSTG,LOLARR,LOLARRT,
     ;             LOLDSP,LOLDFM,LOLPOR,LOLCRD,LOLCOE,LOLFOD
 1000        FORMAT(11I5)
          ENDIF

       ENDDO  !NEXT ELEMENT

       RETURN

 9200  CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;            'GENERIC FORTRAN ERROR WHEN READING ELEMENT '
     ;      //'ZONE NUMBER (CARD B3.3)',
     ;      NROW,0,IUGRID,1,3.17)

       RETURN
       END
