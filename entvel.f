       SUBROUTINE ENTVEL 
     ; (IDIMQ    ,IERROR   ,INPWR    ,IODIM    ,IOWAR    ,IUCAL
     ; ,MAINF    ,NROW     ,NUMEL    ,ACTH     ,FILENAME ,LDIM
     ; ,LTYPE    ,QXYZ     ,VD       ,XNORVD)

*****************************************************************************
* PURPOSE
*     Reads darcy's velocity
*
* DESCRIPTION
*     Reads darcy's velocity
*
* EXTERNAL VARIABLES: ARRAYS
*
*  ACTH                   Aquifer thickness of every element. Cross section for 
*                         1-D elements, thickness for 2-D elements.             
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  LDIM                   Vector containing the dimension of each element       
*  LTYPE                  Vector containing the type of each element            
*  QXYZ                   Products between the different components of          
*                         Darcy's velocity divided by its norm                  
*  VD                     Darcy's velocity                                      
*  XNORVD                 Euclidean norm of Darcy's velocity                    
*
* INTERNAL VARIABLES: ARRAYS
*
*  V                      Auxiliary array to store darcy's velocity
*  VOLD                   Auxiliary array to store darcy's velocity 
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMQ                  Used to dimension array QXYZ                          
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IODIM                  Maximum dimension of any element included             
*                         in the problem                                        
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUCAL                  Unit number of INI file                               
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NROW                   Current record number                                 
*  NUMEL                  Number of elements                                    
*
* INTERNAL VARIABLES: SCALARS
*
*  ITIP                   Current element's type
*  LANI                   Current element's anisotropy
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
*     SCR      4-1997     Revision and verification
*     AMS      1-1998     Common elimination and addition of header
*****************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*100 LEEL,LEAUX,FILENAME(20)*20

       DIMENSION VD(IODIM,NUMEL),QXYZ(IDIMQ,NUMEL),LDIM(NUMEL),
     ;           XNORVD(NUMEL),ACTH(NUMEL),LTYPE(NUMEL),V(3),VOLD(3)

       IF (INPWR.NE.0) WRITE(MAINF,3000) 
 3000  FORMAT(//,' FLOW RATE PER UNIT WIDTH',//,' ELT.    ',
     ;        'VX      VY',/,' ')

       LEAUX=LEEL(FILENAME,IUCAL,MAINF,NROW,INPWR) !Skips title line

       
       IC=1
       DO WHILE (L.LT.NUMEL)

         LEAUX=LEEL(FILENAME,IUCAL,MAINF,NROW,INPWR)
       READ(LEAUX,1000,ERR=9000) L,(V(ID),ID=1,IODIM)
 1000  FORMAT(I5,3F10.0)

C_______________________Checks data

       IF (IC.EQ.1) THEN
          IF (L.NE.IC) THEN
             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;                  ' ABSENT INITIAL VELOCITY IN'//
     ;                  ' FIRST ELEMENT ',NROW,1,IUCAL,1,8.01)
          ELSE
             DO ID=1,IODIM
               VD(ID,L)=V(ID)
             END DO
          END IF
       ELSE
          IF (L.GE.IC) THEN
            IF (L.EQ.IC)THEN
              DO ID=1,IODIM
                VD(ID,L)=V(ID)
              END DO
            ELSE
               DO J=IC,L-1
                 DO ID=1,IODIM
                   VD(ID,J)=VOLD(ID)
                 END DO
               END DO

               DO ID=1,IODIM
                 VD(ID,L)=V(ID)
               END DO
               IC=L
            ENDIF
          ELSE
             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;            ' INCORRECT ORDER IN ELEMENT NUMBERING FOR'//
     ;            ' INITIAL VELOCITY ',NROW,1,IUCAL,1,8.02)
          END IF
       END IF

       DO ID=1,IODIM
         VOLD(ID)=V(ID)
       END DO
       IC=IC+1

       END DO

      

C_______________________Makes some auxiliary variables

         DO L=1,NUMEL
             IF (INPWR.NE.0) WRITE(MAINF,3100) L,(VD(ID,L),ID=1,IODIM)
 3100        FORMAT(I5,3G15.7)

             ITIP=LTYPE(L)               !Element's type
             LDIMAUX=LDIM(L)             !Element's dimension
             LANI=LDIMAUX*(LDIMAUX+1)/2
             IF (LDIMAUX.EQ.1) VD(1,L)=ACTH(L)*VD(1,L)
             XNOR=0.D0

             DO I=1,LDIMAUX
                XNOR=XNOR+VD(I,L)*VD(I,L)
             END DO

             XNOR=DSQRT(XNOR)

             IF (XNOR.LT.1.D-25) THEN
                XNORVD(L)=0.D+00

                DO I=1,IDIMQ
                   QXYZ(I,L)=0.D+00
                END DO

             ELSE
                
                XNORVD(L)=XNOR
                DO I=1,LDIMAUX
                   QXYZ(I,L)=VD(I,L)*VD(I,L)/XNOR            !Diagonal
                END DO

                DO I=LDIMAUX+1,LANI
                   QXYZ(I,L)=VD( MIN(I-LDIMAUX+1,3) ,L)*        !Non-Diagonal
     ;                       VD( MAX(I-4,1) ,L)/XNOR
                END DO
                 
             END IF

         END DO 

         RETURN

 9000  CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;            'GENERIC FORTRAN ERROR WHEN READING INITIAL VELOCITY'
     ;            ,NROW,1,IUCAL,1,8.03)


       END
