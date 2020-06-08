       LOGICAL FUNCTION FILEEX (FILE)

*****************************************************************************
* PURPOSE
*     Checks the existence of a file
*
* EXTERNAL VARIABLES: STRINGS
*
*  FILE                   String containing the file name to check
*
* INTERNAL VARIABLES: LOGICALS
*
*       EX                True if FILE exists otherwise, false.
*
* HISTORY
*
*      AMS        1988     First coding
*      SCR      4-1997     Revision
*
****************************************************************************

       CHARACTER*20 FILE
       LOGICAL*4 EX
       INQUIRE(FILE=FILE,EXIST=EX)
       FILEEX=EX
       RETURN
       END
