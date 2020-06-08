       SUBROUTINE  ENTDATLX_ELEM
     ; (IERROR   ,INPWR    ,IOEQT    ,IOWAR    ,IUGRID   ,LDARR
     ; ,LDARRT   ,LDCOE    ,LDCRD    ,LDDFM    ,LDDSP    ,LDFOD
     ; ,LDPOR    ,LDSTG    ,LDTRA    ,LMXNDL   ,MAINF    ,NBAND
     ; ,NROW     ,NUMEL    ,NUMNP    ,ACTH     ,FILENAME ,KXX
     ; ,LNNDEL   ,LTYPE    ,X        ,Y )

********************************************************************************
* PURPOSE 
*
*    Reads and checks data cards form GRI file corresponding to
*    element data: Card group B3. 
*
* DESCRIPTION
*
*    This subroutine reads and checks element data form card group B3
*    (B3.1,B3.2 and B3.3) on the following way:  
*
*         * Reads the current card trough character function LEEL,which
*           returns the LEAUX string
*         * Load the current card variables by reading LAUX string 
*         * Checks the read data and writes error messages on MAIN FILE
*           calling subroutine ERROR
*         * Write in MAIN FILE interpreted data, after each card is 
*           read.
*         * Data read from card B3.2: 
*             - Element corners (KXX(node,elem)) are lineary interpolated for
*               missing elements before writting on MAIN file
*             - LNNDEL and LTYPE are assumed to be the same for all missing 
*               elements and equal to the last element before missing ones. 
*
*         * Data read from card B3.3 (element parameter zones) is checked and 
*           missing nodes are assigned default values (card B3.1) or those read
*           in the last node (if corresponging default value in card B3.1
*           equals 0) through subroutine INTERP_EZNUM.
* 
* EXTERNAL VARIABLES: ARRAYS 
*     
*  ACTH                   Aquifer thickness of every element. Cross section for 
*                         1-D elements, thickness for 2-D elements.             
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing type for element j                  
*
* EXTERNAL VARIABLES: SCALARS  
*
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOEQT                  Type of problem to be solved                          
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUGRID                 Unit number of GRI file                               
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NBAND                  Half Bandwith (maximum difference between the
*                         numbers of two nodes belonging to the same element)
*  NROW                   Current record number                                 
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: ARRAYS 
*
*  NX(i)                  Vector containing node number of corner i for
*                         the current element. (i is in counterclockwise
*                         order)
*     
*
* INTERNAL VARIABLES: SCALARS 
*
*  LDTRA                  Default value for transmisivity zone number
*  LDARR                  Default value for areal recharge
*  LDARRT                 Default value for transient areal recharge
*  LDSTG                  Default value for storage coeficient    
*  LDDSP                  Default value for dispersivity
*  LDDFM                  Default value for matrix diffusion
*  LDCOE                  Default value for external concentration
*  LDPOR                  Default value for porosity
*  LDCRD                  Default value for retardation coefficient    
*  NOLD                   Index containing last element read
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ASS_EVAL               Assigns default values for parameters zone number
*  ERROR                  Writes the current error message and error number     
*  INTERP_EZNUM           Interpolates zone numbers for missing elements
*  LDIMEN                 Computes the dimension of a given element             
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*  VERIFY_BW              Checks if the input bandwidth is coherent with the
*                         actual one
*
* HISTORY
*
*     AM         1988     Firs coding
*     SCR      4-1997     Revision and verification
*     AMS      1-1998     Common elimination
*     AMS      4-1999     Splitting in two subroutines. One for grid elements
*                         and another for zones
*
********************************************************************************


       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*100 LEEL,LEAUX,FILENAME(20)*20  

       DIMENSION KXX(LMXNDL,NUMEL),LNNDEL(NUMEL),LTYPE(NUMEL),
     ;           NX(8),ACTH(NUMEL),X(NUMNP),Y(NUMNP)

C------------------------- FIRST EXECUTABLE STATEMENT.
C------------------------- Reads default values of parameter zones
C------------------------- Card 3.1

       LEAUX=LEEL(FILENAME,IUGRID,MAINF,NROW,INPWR)
       READ(LEAUX,1000,ERR=9000) LDTRA,LDSTG,LDARR,LDARRT,LDDSP,LDDFM,
     ;                           LDPOR,LDCRD,LDCOE,LDFOD
 1000  FORMAT(11I5)
 
C------------------------- Writes Card 3.1 in MAIN file

       IF (INPWR.NE.0) THEN
          WRITE(MAINF,3000) LDTRA,LDSTG,LDARR,
     ;         LDARRT,LDDSP,LDDFM,LDPOR,LDCRD,LDCOE,LDFOD
 3000     FORMAT(////,
     ;         10X,'ELEMENT INFORMATION'/10X,19('=')//
     ;          5X,'DEFAULT ZONE NUMBERS'/5X,20('-')/
     ;          5X,'TRANSMISSIVITY .......=',I5,2X,'|',/,
     ;          5X,'STORAGE COEFFICIENT ..=',I5,2X,'|',/,
     ;          5X,'RECHARGE (STEADY-ST) .=',I5,2X,'|',/,
     ;          5X,'RECHARGE (TRANSIENT) .=',I5,2X,'|',/,
     ;          5X,'DISPERSIVITY .........=',I5,2X,'|',2X,
     ;               'ZONE NUMBERS',/,
     ;          5X,'DIFFUSION COEFFICIENT =',I5,2X,'|',/,
     ;          5X,'POROSITY .............=',I5,2X,'|',/,
     ;          5X,'RETARDATION COEFFICNT =',I5,2X,'|',/,
     ;          5X,'EXTERNAL CONCENTRATION=',I5,2X,'|',//,
     ;          5X,'FIRST ORDER DEC. COEF.=',I5,2X,'|',//,
     ;         16X,'NODES',16X,/,11X,14('-'),//,
     ;          '  N.EL. TYPE N. NOD. N.1  N.2  N.3  N.4  ',
     ;             'N.5  N.6  N.7  N.8    AREA   ')
       END IF

C------------------------- Initializes element counter

       NE=0

C------------------------- Initializes ESPES (thickness) to zero

       ESPES=0.D0

C------------------------- Starts input element loop

       DO WHILE (NE.LT.NUMEL)
          NOLD=NE
          N=NOLD+1

C------------------------- Reads element type, nodes and tfickness
C------------------------- Card B3.2

          LEAUX=LEEL(FILENAME,IUGRID,MAINF,NROW,INPWR)
          READ(LEAUX,1100,ERR=9100) NE,ITI,NNE,NX(1),NX(2),NX(3),
     ;              NX(4),NX(5),NX(6),NX(7),NX(8),ESP
 1100     FORMAT(11I5,F10.0)
 1010     FORMAT(2I5,3X,9I5,2G11.4)

C------------------------- Checks element number

          IF (NE.GT.NUMEL) CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;       ' TOO LARGE ELEMENT NUMBER '
     ;        ,NROW,1,IUGRID,1,3.02)

C------------------------- Checks increasing order of elements numbers

          IF (NE.LT.N)
     ;      CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;      'INCORRECT ORDER IN ELEMENT'
     ;          //'NUMBERING',NROW,1,IUGRID,1,3.03)

C------------------------- Checks correct number of nodes in element

          IF (NNE.LT.2 .OR. NNE.GT.8) THEN
             WRITE(MAINF,3100) NNE
 3100        FORMAT('INCORRECT NUMBER OF NODES IN ELEMENT NUMBER ',I5)
             CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;     'INCORRECT NUMBER OF NODES IN ELEMENT',NROW,0,IUGRID,1,3.04)
          ENDIF   

C------------------------- Checks if node numbers of NE element are greater 
C------------------------- than NUMNP

          DO LX=1,NNE

             IF (NX(LX).GT.NUMNP .OR. NX(LX).LE.0) THEN
                CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;         'NODE NUMBER IS OUT OF RANGE ',NROW,2,IUGRID,1,3.05)
             ELSE
                KXX(LX,NE)=NX(LX)
             END IF
          END DO

          LTYPE(NE)=ITI
          LNNDEL(NE)=NNE
          IF (IOEQT.NE.1) ACTH(NE)=ESP

C------------------------- Checks the geometry of four nodes elements

       IF (NNE.EQ.4.AND.ITI.EQ.5) CALL CH_GEO_4N (LMXNDL, NUMEL, 
     ;                                     NUMNP, KXX, NX, X, Y, NE)

C------------------------- Element interpolation.
C------------------------- Checks element type. It must be the same for both
C------------------------- elements that are used to interpolate.


          IF (N.NE.NE) THEN

             IF (ITI.NE.LTYPE(NOLD)) THEN
                CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;          'WRONG TYPE ELEMENTS SEQUENCE',
     ;           NROW,0,IUGRID,1,3.18)
             ELSE

C------------------------- Interpolation loop

                DO JJ=1,NNE
                   NX(JJ)=(KXX(jj,NE)-KXX(jj,NOLD))/(NE-NOLD)
                END DO 

                NN=NOLD
                DO NINC=1,NE-N
                   NN=NN+1
                   DO JJ=1,NNE
                      KXX(JJ,NN)=KXX(JJ,NOLD)+NINC*NX(JJ)
                   END DO
                   LNNDEL(NN)=LNNDEL(NOLD)
                   LTYPE(NN)=LTYPE(NOLD)
                   IF (IOEQT.NE.1) ACTH(NN)=ACTH(NOLD)
                   IF (INPWR.NE.0) THEN  !WRITE THE INTERPOLATE
                      IF (IOEQT.NE.1) ESPES=ACTH(NN)
                         WRITE(MAINF,1010) NN,LTYPE(NN),LNNDEL(NN),
     ;                    (KXX(JJ,NN),JJ=1,NNE),(0,JJ=NNE+1,8),ESPES
                   ENDIF
                ENDDO  
             ENDIF
          ENDIF

          IF (IOEQT.NE.1) ESPES=ACTH(NE)
          IF (INPWR.NE.0) 
     ;       WRITE(MAINF,1010) NE,LTYPE(NE),LNNDEL(NE),
     ;         (KXX(JJ,NE),JJ=1,NNE),(0,JJ=NNE+1,8),ESPES
              
       ENDDO !Next element

C------------------------- Verify bandwidth

       CALL VERIFY_BW
     ;(   LNNDEL     ,KXX        ,NUMEL      ,NBAND       ,LMXNDL
     ;   ,IERROR     ,IOWAR      ,MAINF      ,IUGRID      ,FILENAME)

       RETURN

 9000  CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;            'GENERIC FORTRAN ERROR WHEN READING ELEMENT'
     ;      //' DEFAULT VALUES (CARD B3.1)',
     ;      NROW,0,IUGRID,1,3.15)

       RETURN

 9100  CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;            'GENERIC FORTRAN ERROR WHEN READING ELEMENT DATA '
     ;      //'(CARD B3.2)',NROW,0,IUGRID,1,3.16)

       RETURN

       END SUBROUTINE  ENTDATLX_ELEM
