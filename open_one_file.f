      SUBROUTINE OPEN_ONE_FILE (FILE,MAINF,IUN) 

********************************************************************************
* PURPOSE 
*     Opens a specified file
*
* DESCRIPTION
*     Checks the existence of the FILE and if it exists, opens it. 
*     Otherwise stops the execution. 
* 
* EXTERNAL VARIABLES: SCALARS
*
*  FILE                   File name to be open
*  IUN                    Current FILE unit number    
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  FILEEX                 Checks the existence of a file
*
* HISTORY
*
*     AM         1988     First coding
*     SCR 04-abr-1997     Revision and verification
*
********************************************************************************

       CHARACTER*(*) FILE
       LOGICAL FILEEX
          
       IF (FILEEX(FILE)) THEN
          OPEN (UNIT=IUN,FILE=FILE,STATUS='OLD')
       ELSE
          WRITE(MAINF,3000) FILE
          STOP 'FILE NOT FOUND'
       END IF

       RETURN
3000   FORMAT(/2X,'INPUT FILE ..... ',A20,3X,'NOT FOUND')
       END
