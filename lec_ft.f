       SUBROUTINE LEC_FT
     ; (EPSFLU   ,EPSTRA   ,IDIMFNT  ,IERROR   ,INPWR    ,IOCNSF
     ; ,IOCNST   ,IOEQT    ,IOFLLI   ,IOFMLF   ,IOFMLT    ,IORTS
     ; ,IOTRLI   ,IOTRS    ,IOWAR    ,IUTIM    ,MAINF     ,NFNT
     ; ,NINT     ,NROW     ,THETAF   ,THETAT   ,DTMXDS    ,FILENAME
     ; ,FNT      ,ISOLEQ   ,KINT     ,TIME)


*****************************************************************************
* PURPOSE
*     Reads discretization time data
*
* DESCRIPTION
*     Reads observation time, maximum desirable time increment (not used in 
*     linear problems), number of divisions between two sequential observation
*     times, time functions and weighting parameters. 
*
* EXTERNAL VARIABLES: ARRAYS
*
*  DTMXDS                 Maximum allowed time increment for nonlinear
*                         problems
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  FNT                    Array containing time functions                       
*  ISOLEQ                 Array containing the type of head/concentration 
*                         solution desired for the user at each obs. time
*  KINT                   Number of solution time increments, between           
*                         successive observation times.                   
*  TIME                   Observation times.                                    
*
* EXTERNAL VARIABLES: SCALARS
*
*  EPSFLU                 Time weighting parameter for nonlinear flow problems
*  EPSTRA                 Time weighting parameter for nonlinear transport
*                         problems
*  IDIMFNT                Used for dimensioning array FNT, it coincides with    
*                         NFNT if this latter is not zero                       
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOCNSF                 Scheme for storage term in flow problem               
*  IOCNST                 Scheme for mass storage term in transport problem     
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOFMLF                 Flow Formulation number                               
*  IOFMLT                 Transport formulation number                          
*  IORTS                  Transport regime                                      
*  IOTRLI                 Idem to IOFLLI, in the case of transport.             
*  IOTRS                  Flow regime                                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUTIM                  Unit number of TIM file                               
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFNT                   Number of time functions used for describing time     
*                         dependence of all transient parameters                
*  NINT                   Number of observation times                           
*  NROW                   Current record number                                 
*  THETAF                 Time weighting parameter for flow problems
*  THETAT                 Time weighting parameter for transport problems
*
* INTERNAL VARIABLES: SCALARS
*
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*
* HISTORY
*
*     AMS        1988     First coding
*     SCR      5-1997     Revision and verification
*     AMS      1-1998     Common elimination. addition of header
******************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*100 LEAUX,LEEL,FILENAME(20)*20
        
       DIMENSION TIME(NINT),FNT(IDIMFNT,NINT+1),KINT(NINT),DTMXDS(NINT)
     ;          ,ISOLEQ(NINT,4)

C-------------------------  Internal
      INTEGER*4::NI,NIOLD

C_______________________Reads times and division of time intervals

C_______________________Writes title in MAIN file

       IF (INPWR.NE.0) WRITE(MAINF,3000)
 3000  FORMAT(////,20X,'TIME INFORMATION',/,20X,16('*'),//,
     ; 10X,'TIME DISCRETIZATION',11X,'TYPE OF SOLUTION',
     ; 2X,'NUMBER OF PROBLEM',/,10X,19('-'),11X,16('-'),2X,
     ; 17('-'),/,' INTRV.   TIME    # DIV.  DES. T. INC.',
     ; 2('    FLOW  TRANSP.  '))

C------------------------- Starts loop through all time steps
C------------------------- Now filling gaps! (JHG-2011)

       NI = 0
       NIOLD = 0

       DO WHILE(NI.LT.NINT)

         LEAUX = LEEL(FILENAME,IUTIM,MAINF,NROW,INPWR)

         READ(LEAUX,1000,ERR=9000) NI,T,KT,DT1,IFLOW,ICONC,IPBFL,IPBTP
 1000    FORMAT (I5,F10.0,I5,F10.0,4I5)
C------------------------- Checks order of timestep
         IF (NI.GT.NINT .OR. NI.LE.0) THEN

           CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     &                'TIME STEP NUMBER  IS OUT OF RANGE ',
     &                NROW,1,IUTIM,1,0.00)
         END IF

C------------------------- Checks increasing order of time steps

        IF (NI.LT.NIOLD) THEN
          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     &               'INCORRECT ORDER IN TIME STEP ',
     &                NROW,1,IUTIM,1,0.00)
        END IF


C-------------------------Fills missing values

        IF (NI.GT.NIOLD+1) THEN
          DT = (T-TOLD)/(NI-NIOLD)
          NI_INI = NIOLD + 1

        ELSE
          DT = 0
          NI_INI = NI
        END IF

        DO I=NI_INI,NI

          TIME(I) = T - (NI-I)*DT
          KINT(I) = KT
          IF (IOFLLI.NE.0 .OR. IOTRLI.NE.0) DTMXDS(I)=DT1
          ISOLEQ(I,1)=IFLOW
          ISOLEQ(I,2)=ICONC
          ISOLEQ(I,3)=IPBFL
          ISOLEQ(I,4)=IPBTP

C------------------------- If no flow/transport problem number has be                    en
C------------------------- assigned, it is set equal to 1

         IF (IOEQT.NE.2 .AND. ISOLEQ(I,3).EQ.0) ISOLEQ(I,3)=1
         IF (IOEQT.NE.1 .AND. ISOLEQ(I,4).EQ.0) ISOLEQ(I,4)=1


         IF (INPWR.NE.0) THEN
            WRITE(MAINF,3100) I,TIME(I),KINT(I),DTMXDS(I)
     &     ,ISOLEQ(I,1),ISOLEQ(I,2),ISOLEQ(I,3),ISOLEQ(I,4)
         END IF

        END DO

 3100     FORMAT (I5,1X,G12.5,I5,3X,G12.5,2(2X,I5),5X,2(2X,I5))

         NIOLD = NI
         TOLD = T

       END DO !WHILE(NI.LT.NINT)

C------------------------- Last ISOLEQ (at NINT) is only used in mass balance to
C------------------------- identify the end of simulation time

       ISOLEQ(NINT,3)=0
       ISOLEQ(NINT,4)=0

C_______________________Reads time functions

       IF (INPWR.NE.0) WRITE(MAINF,3200)
 3200  FORMAT(////,10X,'TIME FUNCTIONS',/,
     ;             10X,'--------------')

       K=MOD(NINT,7)
       IF (K.NE.0) K=1

       DO NT=1,NFNT

C------------------------- DO loop through all lines containing time function 
C------------------------- values

          DO I=1,NINT/7+K

C------------------------- Location inside array FNT

             KPOS=7*(I-1)+1

             LEAUX=LEEL(FILENAME,IUTIM,MAINF,NROW,INPWR)          
             READ(LEAUX,1100,ERR=9100) (FNT(NT,J),
     ;                                       J=KPOS,MIN(KPOS+6,NINT))
 1100        FORMAT (7F10.0)
          END DO

          IF (INPWR.NE.0) WRITE(MAINF,3300) (FNT(NT,J),J=1,NINT)
 3300     FORMAT (7E10.3)
       
       END DO


C_______________________Reads parameters of weighting in resolution of
C_______________________linear and non-linear systems

       IF (IOTRS+IORTS.NE.0) THEN

          LEAUX=LEEL(FILENAME,IUTIM,MAINF,NROW,INPWR) 
          READ (LEAUX,1200) THETAF,THETAT,EPSFLU,EPSTRA

C_______________________Checks ponderation parameters

          IF (THETAF.LT.0 .OR. THETAF.GT.1.OR.
     ;        THETAT.LT.0 .OR. THETAT.GT.1.OR.
     ;        EPSFLU.LT.0 .OR. EPSFLU.GT.1.OR.
     ;        EPSTRA.LT.0 .OR. EPSTRA.GT.1)
     ;      CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;            'PONDERATION PARAMETER OUT OF ORDER',
     ;             NROW,0,IUTIM,1,9.02)


C_______________________Checks ponderations parameters acording to formulation

          IF ((IOFMLF.GT.3).AND.IOFLLI.NE.0.
     .       AND.(IOFMLF.NE.7.AND.IOFMLF.NE.8.AND.IOFMLF.NE.11))THEN

             IF(THETAF.NE.1D0.OR.EPSFLU.NE.1D0)
     ;         CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;            'ACORDING TO THE FORMULATION CHOOSEN, PARAMETERS '
     ;          //'EPSFLU AND THETAF HAVE TO BE EQUAL TO ONE.'
     ;          //' I WILL CHANGE THEM TO THIS VALUE.',
     ;             NROW,0,IUTIM,0,0.0)

             EPSFLU=1D0
     
             IF(IOFMLF.EQ.5.OR.IOFMLF.EQ.6.OR.IOFMLF.EQ.10)THETAF=1D0

          END IF
          
          IF ((IOFMLT.GT.3.).AND.IOTRLI.NE.0.
     .       AND.(IOFMLT.NE.7.AND.IOFMLT.NE.8.AND.IOFMLT.NE.11))THEN
             IF(THETAT.NE.1D0.OR.EPSTRA.NE.1D0)
     ;         CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;             'ACORDING TO THE FORMULATION CHOOSEN, PARAMETERS '
     ;           //'EPSTRA AND THETAT HAVE TO BE EQUAL TO ONE.'
     ;           //' I WILL CHANGE THEM TO THIS VALUE.',
     ;              NROW,0,IUTIM,0,0.0)

     
             EPSTRA=1D0

             IF(IOFMLT.EQ.5.OR.IOFMLT.EQ.6.OR.IOFMLT.EQ.10)THETAT=1D0
          END IF
              
          IF ((IOFMLF.EQ.2.OR.IOFMLF.EQ.3).AND.IOFLLI.NE.0)THEN
             IF(THETAF.NE.EPSFLU)
     ;         CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;             'ACORDING TO THE FORMULATION CHOOSEN, PARAMETERS '
     ;           //'THETAF AND EPSFLU HAVE TO BE EQUAL.'
     ;           //' I WILL CHANGE THEM ACORDING TO IT, TO THE GREATER'
     ;           //' VALUE.',NROW,0,IUTIM,0,0.0)
                  

            IF(THETAF.GT.EPSFLU)THEN
               EPSFLU=THETAF
            ELSE
               THETAF=EPSFLU
            ENDIF

          END IF

          IF ((IOFMLT.EQ.2.OR.IOFMLT.EQ.3).AND.IOTRLI.NE.0)THEN
             IF(THETAT.NE.EPSTRA)
     ;         CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;             'ACORDING TO THE FORMULATION CHOOSEN, PARAMETERS '
     ;           //'THETAT AND EPSTRA HAVE TO BE EQUAL.'
     ;           //' I WILL CHANGE THEM ACORDING TO IT, TO THE GREATER'
     ;           //'VALUE.',NROW,0,IUTIM,0,0.0)

             IF (THETAT.GT.EPSTRA)THEN
               EPSTRA=THETAT
             ELSE
               THETAT=EPSTRA
             ENDIF

          END IF


C_______________________Corriges scheme of storage term according to formulation

          IF (IOFMLF.GE.7.AND.IOCNSF.NE.1.AND.IOFMLF.NE.11)THEN

            CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;          'SCHEME OF STORAGE TERM IN FLOW EQUATION'
     ;        //'MUST BE LUMPED BECAUSE OF FORMULATION RESTRICTIONS.'
     ;        //' THE PROGRAM CHANGED IT.',NROW,0,IUTIM,0,0.0)


            IOCNSF=1
        
          ENDIF

C_______________________Corriges scheme of storage term according to formulation

          IF (IOFMLT.GE.7.AND.IOCNST.NE.1.AND.IOFMLT.NE.11)THEN

            CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;          'SCHEME OF STORAGE TERM IN TRANSPORT EQUATION'
     ;        //'MUST BE LUMPED BECAUSE OF FORMULATION RESTRICTIONS.'
     ;        //' THE PROGRAM CHANGED IT.',NROW,0,IUTIM,0,0.0)

            IOCNST=1

          END IF



 1200     FORMAT(4F10.0)
          IF (INPWR.NE.0) WRITE (MAINF,3400) THETAF,THETAT,EPSFLU,EPSTRA
 3400     FORMAT(//,'PARAMETERS OF WEIGHTING IN THE ',
     .           'RESOLUTION OF SISTEMS',/58('='),//
     .           'THETAF=......',F10.4,/
     .           'THETAT=......',F10.4,/
     .           'EPSFLU=......',F10.4,/
     .           'EPSTRA=......',F10.4,/)
 
       END IF
 
       RETURN 

 9000  CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;            'ABNORMAL READING TIME',NROW,0,IUTIM,1,9.03)

       RETURN
 9100  CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;            'ABNORMAL READING IN TIME FUNCTION',
     ;             NROW,J,IUTIM,1,9.04)

       RETURN
       END

