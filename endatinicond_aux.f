       SUBROUTINE ENDATINICOND_AUX 
     ; (  NUMNP      ,INDEX      ,IUCAL      ,MAINF      ,IOWAR  
     ;   ,INPWR      ,IERROR     ,NROW       ,FILENAME   ,VAR         )

*****************************************************************************
*
* PURPOSE
*      Reads initial heads and/or concentrations and/or flows
*
* DESCRIPTION
*
*      Reads initial heads and/or concentrations and/or flows
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  VAR                    Array that depending on variable INDEX will 
*                         contain the heads, the concentrations or the flows
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INDEX                  If 1, initial heads are read
*                         If 2, initial concentrations are read                                                       
*                         If 3, nodal flows are read                                                       
*  INPWR                  Allows writing on MAIN FILE                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUCAL                  Unit number of INI file
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NROW                   Current record number                                 
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  IC                     Counter of the number of data read in array VAR, 
*                         including interpolated values
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*  VV                     Used to store the last read value of VAR
*  VVOLD                  Used to store the value of VAR prior to next read
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*  WRITE_ARRAY            Writes a vector in columns in file RES                
*
* HISTORY
*
*     AMS        1988     First coding
*     SCR      5-1997     Revision 
*     AMS      1-1998     Verification
*
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER FILENAME(20)*20,LEEL*100,LEAUX*100
       DIMENSION VAR(NUMNP)

       LEAUX=LEEL(FILENAME,IUCAL,MAINF,NROW,INPWR) !Skips title line
                                      
C------------------------- Writes header depending on INDEX

       IF (INPWR.NE.0) THEN
          IF (INDEX.EQ.1) WRITE(MAINF,3000)
          IF (INDEX.EQ.2) WRITE(MAINF,3100)
          IF (INDEX.EQ.3) WRITE(MAINF,3200)
 3000     FORMAT (/,10X,'INITIAL HEAD',/10X,11('-'))
 3100     FORMAT (/,10X,'INITIAL CONCENTRATION',/10X,20('-'))
 3200     FORMAT (/,10X,'FLOW',/10X,4('-'))
       END IF

C------------------------- Starts counter IC to 1, the current number of values
C------------------------- entered into VAR. 

       IC=1
       I=0
       DO WHILE (I.LT.NUMNP)

          LEAUX=LEEL(FILENAME,IUCAL,MAINF,NROW,INPWR)
          READ(LEAUX,1000,ERR=9000) I,VV
 1000     FORMAT(I5,F10.0)
          IF (IC.EQ.1) THEN
             IF (I.NE.IC) THEN
                CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;                 'ABSENT INITIAL CONDITION IN'//
     ;                 ' FIRST NODE ',NROW,1,IUCAL,1,8.01)
             ELSE
                VAR(I)=VV
             END IF
          ELSE

             IF (I.GT.NUMNP) THEN
                 CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     &               'NODE NUMBER GREATER THAN' //
     &               ' NUMBER OF NODES ' ,NROW,1,IUCAL,1,8.04)
             END IF
              
             IF (I.GE.IC) THEN
                IF (I.EQ.IC)THEN
                   VAR(IC)=VV
                ELSE
                   DO J=IC,I-1
                      VAR(J)=VVOLD
                   END DO 
                   VAR(I)=VV
                   IC=I
                ENDIF
             ELSE
                CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;                 ' INCORRECT ORDER IN NODE NUMBERING FOR'//
     ;                 ' INITIAL CONDITIONS ',NROW,1,IUCAL,1,8.02)
             END IF
          END IF
          VVOLD=VV
          IC=IC+1

       ENDDO !Next node

       IF (INPWR.NE.0) CALL WRITE_ARRAYN (MAINF,NUMNP,' ',VAR)

       RETURN

 9000  CALL ERROR 
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ; 'GENERIC FORTRAN ERROR WHEN READING INITIAL CONDITIONS'
     ; ,NROW,1,IUCAL,1,8.03)

       END
