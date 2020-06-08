       SUBROUTINE ENTFLOW 
     ; (IERROR   ,INPWR    ,IOWAR    ,IUCAL    ,MAINF    ,NROW
     ; ,NUMNP    ,CAUDAL   ,FILENAME)

*****************************************************************************
* PURPOSE
*     Reads nodal flow (needed when IOEQT.EQ.2, ie, if do not solve flow
*     equation)
*
* DESCRIPTION
*     Reads nodal flow (needed when IOEQT.EQ.2, ie, if do not solve flow
*     equation)
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CAUDAL                 Input/output flow at every node.                      
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUCAL                  Unit number of INI file                               
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NROW                   Current record number                                 
*  NUMNP                  Number of nodes                                       
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ENDATINICOND_AUX       Reads NUMNP values of a variable (initial heads or
*                         concentrations or nodal flow
* HISTORY
*
*     AMS        1988     First coding
*     SCR      5-1997     Revision
*     AMS      1-1998     Common elimination and addition of header
*****************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*20 FILENAME(20)
       DIMENSION CAUDAL(NUMNP)

C------------------------- FIRST EXECUTABLE STATEMENT.
C------------------------- Reads flow data

       CALL ENDATINICOND_AUX
     ; (  NUMNP      ,3          ,IUCAL      ,MAINF      ,IOWAR  
     ;   ,INPWR      ,IERROR     ,NROW       ,FILENAME   ,CAUDAL      )

       RETURN
       END
