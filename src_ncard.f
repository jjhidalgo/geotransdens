      SUBROUTINE SRC_NCARD
     ; (IERROR   ,IOWAR    ,IUN      ,MAINF    ,NROW     ,N_CARD
     ; ,FILENAME)

****************************************************************************
* PURPOSE 
*    Searches in unit file IUN the card label N_CARD
* 
* DESCRIPTION
*     This subroutine rewinds the unit file IUN and reads it line by line 
*     until the read string equals to N_CARD. After this, the subroutine 
*     returns the NROW line number in which N_CARD appears, and the program 
*     is ready to read the next line. If N_CARD is not found, a message error 
*     is given by routine ERROR
* 
* EXTERNAL VARIABLES: ARRAYS 
*
*  FILENAME               Array containing all input and output filenames
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  IOWAR                  Allows writing warning messages
*  IUN                    Unit number of current filename
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NROW                   Current record number                                 
*  NCARD                  Contains the current card label
* 
* INTERNAL VARIABLES: STRINGS
*
*  SRCH                   The content of current record 
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number 
*                         on MAIN FILE.
*
* HISTORY
*
*     SCR      4-1997     First coding.
*     AMS      1-1998     Revision
*
****************************************************************************

       CHARACTER*80 SRCH,N_CARD,FILENAME(20)*20

C------------------------- FIRST EXECUTABLE STATEMENT.

       SRCH='     '
       NROW=0
       REWIND(IUN)

C------------------------- Loop looking for N_CARD since the begining of file

       DO WHILE (SRCH.NE.N_CARD)
          READ(IUN,1000,ERR=9000) SRCH
          NROW=NROW+1
       END DO

 1000  FORMAT (A80)
       RETURN

C------------------------- Call ERROR if NCARD is not found
C------------------------- Card 'IFLAG' is allowed to be omitted.

 9000  IF (N_CARD. NE .'IFLAGS                                  '
     ;               //'                                        ' )
     ;                  THEN
          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,' LABEL '//N_CARD//
     ;        ' NOT FOUND IN '//FILENAME(1),NROW,1,IUN,2,1.29)
       ELSE
          IERROR=-99
          RETURN
       ENDIF          

       END
