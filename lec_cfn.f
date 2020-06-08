      SUBROUTINE LEC_CFN
     ; (IERROR   ,INALF    ,INALFT   ,INCHP    ,INCHPT   ,INCON
     ; ,INCONT   ,INPWR    ,INQQP    ,INQQPT   ,IOEQT    ,IORTS
     ; ,IOTRS    ,IOWAR    ,IUPAR    ,MAINF    ,NPARNP   ,NROW
     ; ,NUMNP    ,NZALF    ,NZCHP    ,NZCOE    ,NZQQP    ,CFPARNP
     ; ,FILENAME ,IBCOD    ,IBTCO    ,INCLK    ,NZCLK)

*****************************************************************************
* PURPOSE
*     Reads nodal coefficients of all parmeters and boundary conditions of
*     flow and/or transport equations
*
* DESCRIPTION
*     Reads nodal coefficients of all parameters and boundary conditions of
*     flow and/or transport equations
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CFPARNP                Array containing node coefinient of node j            
*                         corresponding to INpar index parameter zone.          
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  IBCOD                  Flow boundary condition index                         
*  IBTCO                  Transport boundary condition index                    
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INALF                  Index for leakage                                    
*                         in array variables (IXPARNP and CFPARNP)              
*  INALFT                 Index for transient leakage                          
*                         in array variables (IXPARNP and CFPARNP)              
*  INCHP                  Index for prescribed head                             
*                         in array variables (IXPARNP and CFPARNP)              
*  INCHPT                 Index for transient prescribed head                   
*                         in array variables (IXPARNP and CFPARNP)              
*  INCON                  Index for external concentration                      
*                         in array variables (IXPARNP and CFPARNP)              
*  INCONT                 Index for transient external concentration            
*                         in array variables (IXPARNP and CFPARNP)              
*  INPWR                  Allows writing on MAIN FILE                           
*  INQQP                  Index for prescribed flow                             
*                         in array variables (IXPARNP and CFPARNP)              
*  INQQPT                 Index for tramsient prescribed flow                   
*                         in array variables (IXPARNP and CFPARNP)              
*  IOEQT                  Type of problem to be solved                          
*  IORTS                  Transport regime                                      
*  IOTRS                  Flow regime                                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUPAR                  Unit number of file PAR                               
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NPARNP                 Number of nodal parameters in current problem         
*  NROW                   Current record number                                 
*  NUMNP                  Number of nodes                                       
*  NZALF                  Number of leakage zones                              
*  NZCHP                  Number of prescribed head zones                       
*  NZCOE                  Number of external concentration zones                
*  NZQQP                  Number of prescribed flow zones                       
*
* INTERNAL VARIABLES: SCALARS
*
*  DALF                   Leakage default value.                               
*  DALFT                  Transient leakage default value.                               
*  DCHP                   Steady state prescribed head default value.           
*  DCHPT                  Transient prescribed head default value               
*  DCON                   External concentration default value.
*  DCONT                  Transient external concentration default value.
*  DQQP                   Steady state prescribed flow default value.           
*  DQQPT                  Transient prescribed flow default value.              
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
* HISTORY
*
*     AMS        1988     First coding
*     SCR      5-1997     Revision and verification
*     AMS      1-1998     Common elimination. 
*     AMS      4-2002     Correction of some little errors
*****************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*100 LEAUX,LEEL,FILENAME(20)*20
       DIMENSION IBCOD(NUMNP),IBTCO(NUMNP),CFPARNP(NUMNP,NPARNP)

C------------------------- FIRST EXECUTABLE STATEMENT.
       AUXCH = 0D0
       AUXQQ = 0D0
       AUXAL = 0D0
       AUXCHT  = 0D0
       AUXQQT  = 0D0
       AUXALT = 0D0
       AUXCO = 0D0
       AUXCOT = 0D0
       AUXCLK = 0D0
*_______________________Write title in MAIN file

       IF (INPWR.NE.0) WRITE(MAINF,3000)
 3000  FORMAT (////1X,19('*'),' COEFFICIENT INFORMATION ',19('*'))
       NROW=0

*_______________________Read default values

       LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)          
       READ(LEAUX,1000) DCHP,DCHPT,DQQP,DQQPT,DALF,DALFT,DCON,DCONT,DCLK
 1000  FORMAT(9F10.0)
       IF (INPWR.NE.0) 
     . WRITE(MAINF,3100) DCHP,DCHPT,DQQP,DQQPT,DALF,DALFT,DCON,DCONT
     &                  ,DCLK
 3100  FORMAT(////,
     . 10X,'NODAL COEFFICIENTS',/,10X,18('='),//,
     . 5X,'DEFAULT VALUES',/,5X,14('-'),/,
     . 5X,'PRES. HEAD [CHP]..................=',F10.3,/,
     . 5X,'PRES. HEAD (TRANS.- COND.)[CHPT]..=',F10.3,/,
     . 5X,'PRES. FLOW [QQP]..................=',F10.3,/,
     . 5X,'PRES. FLOW (TRANS.- COND.)[QQPT]..=',F10.3,/,
     . 5X,'LEAKAGE [ALF].....................=',F10.3,/,
     . 5X,'LEAKAGE (TRANS-COND)[ALFT]........=',F10.3,/,
     . 5X,'EXTRN. CONC. [CON]................=',F10.3,/,
     . 5X,'EXTRN. CONC.(TRANS.- CONC.) [CONT]=',F10.3,/,
     & 5X,'CONC.LEAKAGE [CLK]................=',F10.3)

*_______________________Initialize to defaults

*_______________________Flow parameters

       IF (IOEQT.NE.2) THEN

         IF (DCHP.NE.0.AND.NZCHP.NE.0) THEN
           DO  I=1,NUMNP
            CFPARNP(I,INCHP)=DCHP
           END DO
         END IF

         IF (DCHPT.NE.0.AND.NZCHP.NE.0) THEN
           DO I=1,NUMNP
            CFPARNP(I,INCHPT)=DCHPT
           END DO
         END IF

         IF (DQQP.NE.0.AND.NZQQP.NE.0) THEN
           DO I=1,NUMNP
            CFPARNP(I,INQQP)=DQQP
           END DO   
         END IF
       
         IF (DQQPT.NE.0 .AND. NZQQP.NE.0) THEN
           DO I=1,NUMNP
            CFPARNP(I,INQQPT)=DQQPT
           END DO 
         END IF

         IF (DALF.NE.0.AND.NZALF.NE.0) THEN
           DO I=1,NUMNP
            CFPARNP(I,INALF)=DALF
           END DO
         END IF

         IF (DALFT.NE.0.AND.NZALF.NE.0) THEN
           DO I=1,NUMNP
            CFPARNP(I,INALFT)=DALFT
           END DO
         END IF

       END IF

*_______________________Transport parameters

       IF (IOEQT.NE.1) THEN    

         IF (DCON.NE.0.AND.NZCOE.NE.0) THEN
           DO I=1,NUMNP
            CFPARNP(I,INCON)=DCON
           END DO
         END IF

         IF (DCONT.NE.0.AND.NZCOE.NE.0) THEN
           DO I=1,NUMNP
            CFPARNP(I,INCONT)=DCONT
           END DO
         END IF
         IF (DCLK.NE.0.AND.NZCLK.NE.0) THEN
           DO I=1,NUMNP
            CFPARNP(I,INCLK)=DCLK
           END DO
         END IF
        END IF

*_______________________Writes headER in MAIN file

        IF ( INPWR.NE.0)   WRITE(MAINF,3200) 
 3200   FORMAT(//,
     . ' NODE      CHP     CHPT      QQP     QQPT    ',
     . '  ALF     ALFT      CON     CONT       CLK')


*_______________________Starts loop reading nodal coeficients

        NOLD=0
        NN=0

        DO WHILE (NN.LT.NUMNP)
        
        LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)
        READ(LEAUX,1100) NN,CHP,CHPT,QQP,QQPT,ALF,ALFT,CON,CONT,CLK
 1100   FORMAT(I5,9F9.0)

*_______________________Checks order of nodes

        IF (NN.GT.NUMNP .OR. NN.LE.0) 
     ;    CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;               'NODE NUMBER  IS'//
     ;               'OUT OF RANGE ',NROW,1,IUPAR,1,4.01)

*_______________________Checks increasing order of nodes

        IF (NN.LE.NOLD) 
     ;    CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;               'INCORRECT ORDER IN NODE NUMBER ',
     ;                NROW,1,IUPAR,1,4.02)

*_______________________Assign flow parameters

        IF (IOEQT.NE.2) THEN   ! FLOW PARAMETERS
          IB=IBCOD(NN)      

*_______________________Leakage coeficient

          IF (IB.LE.2  .AND. (ALF.NE.0. OR .ALFT.NE.0)) THEN
            CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                'INCONSISTENT BOUND. COND. INDEX AND '//
     ;                'LEAKAGE NODAL COEFF. ',NROW,6,IUPAR,0,0.00)
          ENDIF

C------------------------- Coefficients must be assigned always, for 
C------------------------- using with several consecutive problems

          IF (CHP.NE.0.AND.NZCHP.NE.0)  CFPARNP(NN,INCHP)=CHP
          IF (CHPT.NE.0.AND.NZCHP.NE.0) CFPARNP(NN,INCHPT)=CHPT
          IF (ALF.NE.0.AND.NZALF.NE.0)  CFPARNP(NN,INALF)=ALF
          IF (ALFT.NE.0.AND.NZALF.NE.0) CFPARNP(NN,INALFT)=ALFT
          IF (QQP.NE.0.AND.NZQQP.NE.0)  CFPARNP(NN,INQQP)=QQP
          IF (QQPT.NE.0.AND.NZQQP.NE.0) CFPARNP(NN,INQQPT)=QQPT

C------------------------- Checks some possible inconsistencies

*_______________________Prescribed head

          IF (IB.EQ.2  .OR.  IB.EQ.0) THEN
            IF (CHP.NE.0  .OR. CHPT.NE.0) 
     ;        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                  'INCONSISTENT BOUND. COND. INDEX AND PRESC. '//
     ;                  'HEAD NODAL COEFF. ',NROW,2,IUPAR,0,0.00)

            IF (IOTRS.EQ.0 .AND. CHPT.NE.0.)
     ;        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                  'STEADY STATE REGIME BUT CFCHPT IS NOT EQUAL'//
     ;                  ' TO ZERO',NROW,3,IUPAR,0,0.00)
          
          END IF

*_______________________Prescribed flow

          IF (IB.NE.2 .AND. IB.NE.4) THEN
            IF (QQP.NE.0.OR.QQPT.NE.0.)
     ;        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                  'INCONSISTENT BOUND. COND. INDEX AND PRESC. '//
     ;                  'FLOW  NODAL COEFF. ',NROW,4,IUPAR,0,0.00)
            IF (IOTRS.EQ.0 .AND. QQPT.NE.0.) 
     ;        CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                  'STEADY STATE REGIME BUT CFQQPT IS NOT EQUAL TO'
     ;                   //' ZERO',NROW,5,IUPAR,0,0.00)
          END IF
        END IF   ! IOEQT.NE.2

*_______________________Assign transport parameters

        IF (IOEQT.NE.1) THEN  
          IB=IBTCO(NN)

C------------------------- Coefficients must be assigned always, for 
C------------------------- using with several consecutive problems

          IF (CON.NE.0.AND.NZCOE.NE.0)  CFPARNP(NN,INCON)=CON
          IF (CONT.NE.0.AND.NZCOE.NE.0) CFPARNP(NN,INCONT)=CONT
	    IF (CLK.NE.0.AND.NZCLK.NE.0) CFPARNP(NN,INCLK)=CLK

C------------------------- Checks some possible inconsistencies

*_______________________External concentration

          IF (IB.EQ.0 .AND. (CON.NE.0 .OR. CONT.NE.0))
     ;      CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                 'INCONSISTENT BOUND. COND. INDEX AND '//
     ;                 'EXT. CONCENTRATION NODAL COEFF.'
     ;                 ,NROW,8,IUPAR,0,0.00)

          IF (IORTS.EQ.0 .AND. CONT.NE.0.)
     ;      CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                'STEADY STATE REGIME BUT CFCONT IS NOT EQUAL TO '
     ;                 //'ZERO.',NROW,9,IUPAR,0,0.00)

*_______________________Concentration leakage

          IF (IB.LT.5 .AND. CLK.NE.0)
     ;      CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                 'INCONSISTENT BOUND. COND. INDEX AND '//
     ;                 'EXT. CONC. LEAKAGE NODAL COEFF.'
     ;                 ,NROW,10,IUPAR,0,0.00)


        END IF        ! IOEQT.NE.1

*_______________________Writes missing nodes coeficients in MAIN file

        IF (INPWR.GT.1) THEN
           I_INI=NOLD+1
        ELSE
           I_INI=NN
        ENDIF

         IF (INPWR.NE.0) THEN
          DO II=I_INI,NN      !FOR WRITE
            IF (IOEQT.NE.2) THEN
              IF (NZCHP.NE.0)AUXCH=CFPARNP(II,INCHP) 
              IF (NZQQP.NE.0)AUXQQ=CFPARNP(II,INQQP) 
              IF (NZALF.NE.0)AUXAL=CFPARNP(II,INALF) 
              IF (NZCHP.NE.0)AUXCHT=CFPARNP(II,INCHPT) 
              IF (NZQQP.NE.0)AUXQQT=CFPARNP(II,INQQPT) 
              IF (NZALF.NE.0)AUXALT=CFPARNP(II,INALFT)
            END IF
            IF (IOEQT.NE.1) THEN
              IF (NZCOE.NE.0)AUXCO=CFPARNP(II,INCON)
              IF (NZCOE.NE.0) AUXCOT=CFPARNP(II,INCONT)
	        IF (NZCLK.NE.0)AUXCLK=CFPARNP(II,INCLK)
            END IF
            WRITE(MAINF,3300)
     ;      II,AUXCH,AUXCHT,AUXQQ,AUXQQT,AUXAL,AUXALT,AUXCO,AUXCOT
     &      ,AUXCLK
 3300       FORMAT(I5,9F9.3)
          END DO
         ENDIF
         NOLD=NN
       END DO  ! SEARCH NEXT NODE

       RETURN

       END
