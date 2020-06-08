       SUBROUTINE ENDATINICOND 
     ; (IERROR   ,INPWR    ,IORTS    ,IOTRS    ,IOWAR    ,IUCAL
     ; ,MAINF    ,NFL_SIM  ,NTP_SIM  ,NUMNP    ,CCAL     ,FILENAME
     ; ,HCAL)

*****************************************************************************
*
* PURPOSE
*      Reads initial conditions (heads and/or concentrations)
*
* DESCRIPTION
*      Reads initial conditions (heads and/or concentrations)
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CCAL                   Computed concentrations.                              
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  HCAL                   Computed heads at node j                              
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IORTS                  Transport regime                                      
*  IOTRS                  Flow regime                                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUCAL                  Unit number of INI file
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFL_SIM                Number of simultaneous flow problems
*  NTP_SIM                Number of simultaneous transport problems
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  NROW                   Current record number                                 
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ENDATINICOND_AUX       Reads initial heads and/or concentrations and/or flows
*
* HISTORY
*
*     AMS        1988     First coding
*     SCR      5-1997     Revision 
*     AMS      1-1998     Revision to delete COMMONS
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*20 FILENAME (17)
       DIMENSION HCAL(NUMNP,NFL_SIM),CCAL(NUMNP,NTP_SIM)
 
C------------------------- FIRST EXECUTABLE STATEMENT.
C------------------------- Initializes the counter of record number

       NROW=0

C_______________________Writes read data on MAIN file

       IF (INPWR.NE.0) WRITE(MAINF,3000)
 3000  FORMAT(////10X,'INITIAL CONDITIONS',/10X,18('='))

C_______________________Flow problem

       IF (IOTRS.EQ.1) THEN
          DO IPROB=1,NFL_SIM

             CALL ENDATINICOND_AUX
     ; (  NUMNP    ,1          ,IUCAL    ,MAINF      ,IOWAR  
     ;   ,INPWR    ,IERROR     ,NROW     ,FILENAME   ,HCAL(1,IPROB) )

          ENDDO
       ENDIF

C_______________________Transport problem

       IF (IORTS.EQ.1) THEN
          DO IPROB=1,NTP_SIM

             CALL ENDATINICOND_AUX
     ; (  NUMNP    ,2          ,IUCAL    ,MAINF      ,IOWAR  
     ;   ,INPWR    ,IERROR     ,NROW     ,FILENAME   ,CCAL(1,IPROB) )

          ENDDO
       ENDIF

       RETURN
       END
