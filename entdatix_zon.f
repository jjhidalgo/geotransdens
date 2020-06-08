       SUBROUTINE ENTDATIX_ZON
     ; (IDALF    ,IDALFT   ,IDCHP    ,IDCHPT   ,IDCON    ,IDCONT   
     ; ,IDDMT    ,IDQQP    ,IDQQPT   ,IERROR   ,INALF    ,INALFT
     ; ,INCHP    ,INCHPT   ,INCON    ,INCONT   ,INDMT    ,INPWR
     ; ,INQQP    ,INQQPT   ,IOEQT    ,IORTS    ,IOTRS
     ; ,IOWAR    ,IPROB    ,IUGRID   ,MAINF    ,NPARNP   ,NROW
     ; ,NTDMT    ,NUMNP    ,NZALF    ,NZCHP    ,NZCOE    ,NZDMT
     ; ,NZQQP    ,FILENAME ,IBCOD    ,IBTCO    ,IXPARNP
     & ,IDCLK,NZCLK,INCLK)

*******************************************************************************
* PURPOSE 
*
*    Reads and checks data cards form GRI file corresponding to boundary 
*    and zones related to nodal parameters: Card group B2. 
*
* DESCRIPTION
*     
*    This subroutine reads and checks nodal data form cards group 
*    B2 on the following way:  
*
*         * Reads the current card trough character function LEEL,which
*           returns the LEAUX string
*         * Load the current card variables by reading LAUX string 
*         * Checks the read data and writes error messages on MAIN FILE
*           calling subroutine ERROR
*         * Write in MAIN FILE interpreted data, after each card is 
*           read.
*         * Data read from card 1.3 (Bounday conditions and parameter zones) 
*           is checked and missing nodes are assigned default values or those
*           read in the last node (if corresponging default value in card 1.1 
*           equals 0) through subroutine INTERP_NZNUM
*
* EXTERNAL VARIABLES: ARRAYS 
*     
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  IBCOD                  Flow boundary condition index                         
*  IBTCO                  Transport boundary condition index                    
*  IXPARNP                Array containing zone number of each node j,          
*                         corresponding to INpar index parameter zone.          
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
*  IORTS                  Transport regime                                      
*  IOTRS                  Flow regime                                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUGRID                 Unit number of GRI file
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NPARNP                 Number of nodal parameters in current problem         
*  NTDMT                  If diferent form zero, number of terms in matrix      
*                         diffusion(if zero, no diffusion)                      
*  NUMNP                  Number of nodes                                       
*  NZALF                  Number of leakance zones                              
*  NZCHP                  Number of prescribed head zones                       
*  NZCOE                  Number of external concentration zones                
*  NZDMT                  Number of matrix diffusion zones                      
*  NZQQP                  Number of prescribed flow zones                       
*
* INTERNAL VARIABLES: SCALARS 
*
*  IN                     Current node number in read loop
*  IBC                    Same as IBCOD(IN),but used when reading these value
*                         from GRI file      
*  IBT                    Same as IBTCOD(IN),but used when reading these 
*                         value from GRI file
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
*
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*100 LEEL,LEAUX,FILENAME(20)*20
         
       DIMENSION IXPARNP(NUMNP,NPARNP),IBCOD(NUMNP),IBTCO(NUMNP)

C------------------------- Initializes to zero some auxiliar variables

       IAUXIBD=  0
       IAUXIBT=  0
       IAUXALF=  0
       IAUXCHP=  0
       IAUXQQP=  0
       IAUXCHPT= 0
       IAUXQQPT= 0
       IAUXALFT= 0
       IAUXCON=  0
       IAUXCONT= 0
       IAUXDMT=  0
       IAUXCLK=  0

C------------------------- Writes on MAIN file 
  
       IF (INPWR.NE.0) THEN
          IF (IPROB.EQ.1) THEN            ! Only first time
             WRITE(MAINF,3200)
 3200        FORMAT(' NODES    BNDR.COND.             Z O N E S ',/,
     ;      1X,5('-'),4X,10('-'),13X,9('-'),/,'         N   FLW TRP',
     ;      '   CHP  CHPT QQP QQPT  ALF ALFT  CON  CONT DMT CLEAK ')
          ELSE
             WRITE(MAINF,3300) IPROB
 3300        FORMAT(3X,' PROBLEM NUMBER',I5)
          ENDIF
       ENDIF

       IN=0
       I=0

       DO WHILE (IN.LT.NUMNP)
          IOLD=I
          I=I+1

*_______________________Read equation parameters and nodal zone numbers
*_______________________Card B1.3
       
          LEAUX=LEEL(FILENAME,IUGRID,MAINF,NROW,INPWR)      
          READ(LEAUX,1200,ERR=9000) IN,IBC,IBT,ICHP,ICHPT,IQQP,
     ;                     IQQPT,IALF,IALFT,ICON,ICONT,IDMT,ICLK
 1200     FORMAT(13I5)

C------------------------- Flow parameters 

          IF (IOEQT.NE.2) THEN

C------------------------- Leakage zone number

             IF (IALF.GT.NZALF) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                    'LEAKAGE ZONE NUMBER IS TOO LARGE',
     ;                     NROW,8,IUGRID,1,2.03)

C------------------------- Prescribed head

             IF (ICHP.GT.NZCHP) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                    'PRESCRIBED HEAD ZONE NUMBER IS TOO LARGE',
     ;                    NROW,4,IUGRID,1,2.04)

C------------------------- Prescribed flow

             IF (IQQP.GT.NZQQP) CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                   'PRESCR. FLOW ZONE NUMBER IS TOO LARGE',
     ;                   NROW,6,IUGRID,1,2.05)            

C------------------------- Transient parameters

C------------------------- Prescribed head

                IF (ICHPT.GT.NZCHP) CALL ERROR
     ;             (IERROR,IOWAR,MAINF,FILENAME,
     ;             'PRESCR. TRANSIENT HEAD ZONE NUMBER IS TOO LARGE',
     ;             NROW,5,IUGRID,1,2.06)

C------------------------- Prescribed flow

                IF (IQQPT.GT.NZQQP) CALL ERROR
     ;             (IERROR,IOWAR,MAINF,FILENAME,
     ;             'PRESCR. TRANSIENT FLOW ZONE NUMBER IS TOO LARGE',
     ;             NROW,7,IUGRID,1,2.07)

C------------------------- Leakage

                IF (IALFT.GT.NZALF) CALL ERROR
     ;             (IERROR,IOWAR,MAINF,FILENAME,
     ;             'TRANSIENT LEAKAGE ZONE NUMBER IS TOO LARGE',
     ;             NROW,9,IUGRID,1,2.07)

C------------------------- Checks flow boundary conditions

             IF (IBC.LT.0 .OR. IBC.GT.4) CALL ERROR 
     ;          (IERROR,IOWAR,MAINF,FILENAME,
     ;          'FLOW BOUNDARY CONDITION INDEX IS OUT OF RANGE'
     ;          ,NROW,2,IUGRID,1,2.08)

             IF (IBC.LE.2 .AND. (IALF.NE.0 .OR. IALFT.NE.0)) 
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;          'INCORRECT LEAKAGE NODAL ZONE NUMBER (IXALF) '
     ;          ,NROW,9,IUGRID,1,2.09)

             IF (IBC.GT.2 .AND. ((IALF.EQ.0 .AND. IOTRS.NE.1)
     ;                .OR. (IALFT.EQ.0 .AND. IOTRS.NE.0) )) 
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;          'INCORRECT BOUND. COND. INDEX OR LEAKAGE ZONE NUMBER',
     ;          NROW,9,IUGRID,1,2.10)
  
             IF (IBC.NE.2 .AND. IBC.NE.4) THEN
                IF (IQQP.NE.0) 
     ;             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'WRONG BOUND. COND. OR PRESC. FLOW ZONE (STEADY)'
     ;             ,NROW,6,IUGRID,1,2.11)
                IF (IQQPT.NE.0) 
     ;             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'WRONG BOUND. COND. PRESC. FLOW ZONE (TRANS)'
     ;             ,NROW,7,IUGRID,1,2.12)
             ENDIF
  
             IF ( (IBC.EQ.2 .OR. IBC.EQ.4) .AND. 
     ;                              (IQQP.EQ.0 .AND. IQQPT.EQ.0) ) 
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;          'INCORRECT BOUND. COND. OR PRESC. FLOW ZONE NUMBER'
     ;          ,NROW,6,IUGRID,1,2.13)

             IF ( (IBC.EQ.1 .OR. IBC.EQ.3 .OR. IBC.EQ.4) .AND.
     ;                                  (ICHP.EQ.0 .AND. ICHPT.EQ.0) ) 
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;          'WRONG BOUND. COND. OR PRESC. HEAD ZONE NUMBER'
     ;          ,NROW,4,IUGRID,1,2.14)

             IF (IBC.EQ.0 .OR. IBC.EQ.2) THEN 
                IF (ICHP.NE.0) 
     ;             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'WRONG BOUND. COND. OR STEADY PRESC. HEAD ZONE'
     ;             ,NROW,6,IUGRID,1,2.15)
                IF (ICHPT.NE.0) 
     ;             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'WRONG BOUND. COND. OR TRANS. PRESC. HEAD ZONE'
     ;             ,NROW,7,IUGRID,1,2.16)
             ENDIF
          ENDIF                  !IOEQT.NE.2

*_______________________Checks transport parameters

          IF (IOEQT.NE.1) THEN 

C------------------------- Steady external concentration zone number

             IF (ICON.GT.NZCOE) 
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;          'STEADY EXT. CONCENTRATION ZONE NUMBER IS TOO LARGE'
     ;          ,NROW,10,IUGRID,1,2.17) 

C------------------------- Transient external concentration zone number

                IF (ICONT.GT.NZCOE)
     ;             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'TRANS. EXT. CONCENTRATION ZONE NUMBER IS TOO LARGE'
     ;             ,NROW,11,IUGRID,1,2.18)

C------------------------- Concentration leakage zone number

                IF (ICLK.GT.NZCLK)
     ;             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'CONCENTRATION LEAKAGE ZONE NUMBER IS TOO LARGE'
     ;             ,NROW,13,IUGRID,1,2.181)

C------------------------- Checks coherence between zones and 
C------------------------- transport boundary conditions 

             IF (IBT.LT.0 .OR. IBT.GT.5) 
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;          'TRANSP. BOUND. CONDITION INDEX IS OUT OF RANGE'
     ;          ,NROW,3,IUGRID,1,2.19)

             IF (IBT.EQ.0) THEN
                IF (ICON.NE.0) 
     ;             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'WRONG BOUND. COND. OR STEADY CONCENT. ZONE NUMBER'
     ;             ,NROW,3,IUGRID,1,2.20)
                IF (ICONT.NE.0) 
     ;             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'WRONG BOUND. COND. OR TRANS. CONCENT. ZONE NUMBER'
     ;             ,NROW,3,IUGRID,1,2.21)
             ELSE

                IF ( (ICON.EQ.0 .AND. IORTS.NE.1) .OR.
     ;                           (ICONT.EQ.0 .AND. IORTS.NE.0) )
     ;             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'WRONG BOUND. COND. OR CONCENT. ZONE NUMBER'
     ;             ,NROW,3,IUGRID,1,2.22)

	          IF ( (IBT.EQ.5 .AND. ICONT.EQ.0) )
     &             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'WRONG BOUND. COND. OR CONCENT. ZONE NUMBER'
     ;             ,NROW,3,IUGRID,1,2.23)

	          IF (IBT.LT.5 .AND. ICLK.GT.0)
     &             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'WRONG BOUND. COND. OR CONC. LEAK ZONE NUMBER'
     ;             ,NROW,3,IUGRID,1,2.24)
             ENDIF

*_______________________Checks matrix diffusion zones and terms

             IF (NTDMT.EQ.0 .AND. IDMT.NE.0)
     ;          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;          'NO MATRIX DIFF., BUT MATRIX DIFUSION ZONE NOT ZERO',
     ;          NROW,12,IUGRID,1,2.29)

             IF (IDMT.GT.NZDMT. OR .IDMT.LT.0)
     ;          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;         'MATRIX DIFFUSION ZONE OUT OF ORDER',
     ;         NROW,12,IUGRID,1,2.30)
          ENDIF

C------------------------- Assigns default values of zone numbers


C------------------------- Flow

          IF (IOEQT.NE.2) THEN
             IBCOD(IN)=IBC

             IF (NZCHP.NE.0)  
     ;          CALL ASS_INTEG_VAL(ICHP,IDCHP,IXPARNP(IN,INCHP))
             IF (NZQQP.NE.0)
     ;          CALL ASS_INTEG_VAL(IQQP,IDQQP,IXPARNP(IN,INQQP))
             IF (NZALF.NE.0)
     ;          CALL ASS_INTEG_VAL(IALF,IDALF,IXPARNP(IN,INALF))
                 
                IF (NZCHP.NE.0)
     ;             CALL ASS_INTEG_VAL(ICHPT,IDCHPT,IXPARNP(IN,INCHPT))
                IF (NZQQP.NE.0)
     ;             CALL ASS_INTEG_VAL(IQQPT,IDQQPT,IXPARNP(IN,INQQPT))
                IF (NZALF.NE.0)
     ;             CALL ASS_INTEG_VAL(IALFT,IDALFT,IXPARNP(IN,INALFT))

          ENDIF             ! (IOEQT.NE.2)

C------------------------- Transport 

          IF (IOEQT.NE.1) THEN
             IBTCO(IN)=IBT

             IF (NZCOE.NE.0) THEN
                CALL ASS_INTEG_VAL(ICON,IDCON,IXPARNP(IN,INCON))
                CALL ASS_INTEG_VAL(ICONT,IDCONT,IXPARNP(IN,INCONT))
             ENDIF

	       IF (NZCLK.NE.0) THEN
	          CALL ASS_INTEG_VAL(ICLK,IDCLK,IXPARNP(IN,INCLK))
	       END IF
  
             IF (NTDMT.NE.0 .AND. NZDMT.NE.0) 
     ;          CALL ASS_INTEG_VAL(IDMT,IDDMT,IXPARNP(IN,INDMT))
          ENDIF

C------------------------- Interpolates zone numbers

          IF (I.NE.IN)
     ;       CALL INTERP_NZNUM 
     ;(   I          ,IN         ,IOLD       ,INALF      ,INALFT
     ;   ,INCHP      ,INCHPT     ,INQQP      ,INQQPT     ,INCON
     ;   ,INCONT     ,INDMT      ,IOEQT      ,IOTRS      ,IORTS                      
     ;   ,NTDMT      ,NUMNP      ,INPWR      ,NPARNP     ,MAINF
     ;   ,IXPARNP    ,IBTCO      ,IBCOD      ,NZALF      ,NZCOE
     ;   ,NZCHP      ,NZQQP      ,INCLK      ,NZCLK)


          IF (INPWR.NE.0) THEN

             IF (IOEQT.NE.2) THEN    ! Flow

                IAUXIBD=IBCOD(IN)
                IF (NZALF.NE.0)IAUXALF=IXPARNP(IN,INALF)
                IF (NZCHP.NE.0)IAUXCHP=IXPARNP(IN,INCHP) 
                IF (NZQQP.NE.0)IAUXQQP=IXPARNP(IN,INQQP) 


                   IF (NZCHP.NE.0)IAUXCHPT=IXPARNP(IN,INCHPT)
                   IF (NZQQP.NE.0)IAUXQQPT=IXPARNP(IN,INQQPT)
                   IF (NZALF.NE.0)IAUXALFT=IXPARNP(IN,INALFT)

             ENDIF

             IF (IOEQT.NE.1) THEN                 !Transport 

                IAUXIBT=IBTCO(IN)
                IF (NZCOE.NE.0) IAUXCON=IXPARNP(IN,INCON)
                IF (NZCOE.NE.0) IAUXCONT=IXPARNP(IN,INCONT)
                IF (NTDMT.NE.0) IAUXDMT=IXPARNP(IN,INDMT)
	          IF (NZCLK.NE.0) IAUXCLK=IXPARNP(IN,INCLK)
             ENDIF
             WRITE(MAINF,3250) I,IAUXIBD,IAUXIBT,IAUXCHP,
     ;            IAUXCHPT,IAUXQQP,IAUXQQPT,IAUXALF,IAUXALFT,
     ;            IAUXCON,IAUXCONT,IAUXDMT,IAUXCLK
 3250        FORMAT(5X,13I5) 
          ENDIF                         !INPWR.NE.0
       ENDDO       !Next node

      RETURN

C______________________________ Error messages

 9000       CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;          'GENERIC FORTRAN ERROR WHEN READING '
     ;    //'BOUNDARY COND. AND PARAMETERS ZONES (CARD B1.3)',
     ;           NROW,0,IUGRID,2,3.15)

      RETURN
      END
