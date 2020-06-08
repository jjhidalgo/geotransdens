       CHARACTER*100 FUNCTION LEEL(FILENAME,IUN,MAINF,NROW,INPWR)

********************************************************************************
* PURPOSE 
*
*    
*       Checks if the current line of an input data file is a comment line, and
*       set the comment character if '-->' appears. 
*
* DESCRIPTION
*
*       This function returns a string value coinciding with the current input 
*       file line. First of all reads the current line, ensuring that file pointer
*       is not at the end of file (in this case, an error message is processed through
*       subroutine ERROR), increments NROW (contains the current row number). 
*       If any comment char appears in this line:
*
*      -It searches two comment chars together. In this case,the line will
*       be written it in MAIN file.
*      -It searches substring '-->' appearing after a coment char. If so, the next 
*       char after '-->' is accepted as the new comment char.    
*      -Reads the sequent line.
*
* EXTERNAL VARIABLES: ARRAYS 
*
*  FILENAME                 Array containing all the input and output filenames
*                 
* EXTERNAL VARIABLES: SCALARS
*
*       NROW                Current file line number
*       IUN                 Unit number of current filename
*       INPWR               Allows writing on MAIN FILE
*
* EXTERNAL VARIABLES: STRINGS
*
*       LEEL                Value of the character function.It kepts the current
*                           of an input file.
*
* INTERNAL VARIABLES: STRINGS
*
*       COMM_CHAR           Contains the coment character (%, default)
*       DOSCOM_CHAR         Two character lenth string, containing two comment 
*                           character ('%%', default)
*
* INTERNAL VARIABLES: SCALARS 
*
*       KK                  It's a pointer value. Indicates the position in that 
*                           a comment char. followed by '-->' appears in LEEL. 
*
* SUBROUTINES REFERENCED
*
*     ERROR                 Writes the current error message and error number 
*                           on MAIN FILE.
* FUNCTIONS REFERENCED  
*
*
* HISTORY
*
*     AM         1988       First coding
*     SCR 04-abr-1997       Revision and verification
*
*
********************************************************************************

       CHARACTER*1 COMM_CHAR, DOSCOMM_CHAR*2, FILENAME(20)*20

       COMMON/FLAG/IERROR
       DATA COMM_CHAR, DOSCOMM_CHAR /'%','%%'/

*_______________________Reads the current line from IUN
   
 100   READ(IUN,1000,ERR=9000) LEEL
       NROW=NROW+1
1000   FORMAT(A100)

*_______________________Looks for commchar in LEEL

       IF (INDEX(LEEL,COMM_CHAR).NE.0) THEN

*_______________________Writes LEEL in MAIN file if DOSCOM_CHAR appears
*_______________________(for INPWR non zero)
 
          IF (INDEX(LEEL,DOSCOMM_CHAR).NE.0 .AND. INPWR.NE.0) THEN
             WRITE(MAINF,3000) LEEL
          END IF

*_______________________Looks for a change in COM_CHAR

          KK=INDEX(LEEL,COMM_CHAR//'-->')
          IF (KK.NE.0) THEN
             COMM_CHAR=LEEL(KK+4:KK+4)
             DOSCOMM_CHAR=COMM_CHAR//COMM_CHAR
          END IF
          KK=0
          GOTO 100
       END IF

 3000   FORMAT(//,A100)
       RETURN

*_______________________Error Call (if it tryes to read at the end of file)

 9000  CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;            ' END OF FILE ',NROW,1,IUN,1,1.28)
       
       END
