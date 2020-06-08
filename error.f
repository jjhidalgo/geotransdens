       SUBROUTINE ERROR (IERROR,IOWAR,MAINF,FILENAME,MSG,
     ;                   NROW,NCOL,IUN,II,ERNUM)
********************************************************************************
* PURPOSE 
*
*     This subroutine writes on MAIN file diferent error messages given from
*     different calls. It also counts the number of errors (IERROR), deppending of the
*     given error type (II external var.).  
* 
*
* DESCRIPTION
*
*     Error/warning messages, in MSG variable, are divided in two string vars
*     lower than 80 char. It's considered 4 error types, controled by II external
*     variable:
*
*               II=0: Warning. Message may be written in MAIN file, but ERROR var
*                     is not incremented  
*               II=1: Error. Message is written in MAIN file, and ERROR is
*                     incremented (ERROR=ERROR+1). If number of type 1 errors is
*                     greater than 10, program stops.
*               II=2: Fatal input error: Message is written in MAIN file, and programs 
*                     stops immediately.
*               II=3: Time execution error: Message is written in MAIN file and program
*                     is stopped. 
*
* 
* EXTERNAL VARIABLES: ARRAYS 
*
*     FILENAME            Array containing all the input and output filenames                 
*
*
* EXTERNAL VARIABLES: SCALARS
*
*     NROW                Current input file line where error is detected          
*     NCOL                Current input file column where error is detected
*     IUN                 Unit number of current filename
*     
* INTERNAL VARIABLES: SCALARS
*
*     IERROR              Counts number of error types 1 and 3
*     ERNUM               Error code number
*     II                  Error type (0,1,2 or 3)
*     LONG                Error message lenght       
*
* INTERNAL VARIABLES: STRINGS
*
*     MSG                 Error message
*     MSG1                First 80 characters of error message
*     MSG2                Error message part from 80th character
*     BLANK               80 character string lenght containing ' '
*
* HISTORY
*
*     AM         1988     First coding
*     SCR 04-abr-1997     Revision and verification
*
*
********************************************************************************
 
       REAL*4 ERNUM
       CHARACTER FILENAME(20)*20
       CHARACTER*(*) MSG
       CHARACTER*80 BLANK,MSGW1,MSGW2
       DATA BLANK/'                                                     
     ;                          '/
        
       
       MSGW1=BLANK
       MSGW2=BLANK
       LONG=LEN(MSG)

*_______________________Define MSG1 and MSG2 values

       IF (LONG.LE.80) THEN
          MSGW1=MSG
       ELSE
          I=80
*_______________________Looks for the first blank char. before 80th column
*_______________________to avoid breaking message words

          DO WHILE (MSG(I:I).NE.' ') 
            I=I-1
          END DO
          MSGW1=MSG(1:I)
          MSGW2=MSG((I+1):LONG)
       END IF       

*_______________________Checks error type_______________________________

*_______________________Warnings

       IF (II.EQ.0) THEN
         IF(IOWAR.GT.0) WRITE(MAINF,3000) MSGW1,MSGW2,NROW,NCOL,
     .                                    FILENAME(IUN-9)
 3000  FORMAT(/' WARNING: ',/,1X,A80,/,1X,A80,/'  IN ROW ',I5
     .     ,'  COLUMN' ,I5,'   IN FILE..........:',/3X,A20)

       END IF
C_______________________Error in data files

       IF (II.EQ.1) THEN
          IERROR=IERROR+1
          WRITE(MAINF,3100) ERNUM,MSGW1,MSGW2,NROW,NCOL,FILENAME(IUN-9) 
      
          IF (IERROR.GE.10) THEN
             WRITE(MAINF,3200)

 3200  FORMAT(////'  MORE THAT 10 ERRORS. STOP THE LECTURE OF DATES')

             STOP 'MORE THAN 10 ERRORS IN DATES READ'
          END IF
       END IF
C_______________________Fatal error

       IF (II.EQ.2) THEN
         WRITE(MAINF,3100) ERNUM,MSGW1,MSGW2,NROW,NCOL,FILENAME(IUN-9)

 3100 FORMAT(/,1X,' ERROR MESSAGE #',F4.2,/,1X,A80,/,1X,A80,
     ;           /'  IN ROW ',I5
     ;     ,'   COLUMN' ,I5,'   IN FILE........:',/3X,A20)

 

         STOP 'FATAL ERROR IN DATES READ'
       END IF
C_______________________Time excecution error

       IF (II.EQ.3) THEN
          IERROR=IERROR+1
          WRITE(MAINF,3300) MSGW1,MSGW2,IERROR
          STOP !EXECUTION TIME ERROR, I STOP
 3300     FORMAT(/' ERROR: ',/,1X,A80,/,1X,A79,/
     ;           ' ERROR NUMBER: ',I5)
      

       END IF
       RETURN 
       END
