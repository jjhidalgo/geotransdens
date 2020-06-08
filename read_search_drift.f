       SUBROUTINE READ_SEARCH_DRIFT
     ;(KTYPE_GS ,IERROR   ,INPWR    ,IOWAR        ,IUGEO
     ;,MAINF    ,NEXDR_GS ,NVAR_GS  ,SKMEAN    ,IPOLDRIFT_GS ,SEARCH_GS
     ;,TRIM_GS  ,FILENAME)

********************************************************************************
*
* PURPOSE Reads ellipsoid search parameters, kriging polynomial trend definition 
*         ,variable trimming limits, conditional simulation tolerance factor for
*         measurement assignment of a particular group of zones
*                                -- - ---------- ----- -- -----
*
* DESCRIPTION Reads sequentially cards G1 to G5, containing:
*
*             - Card G1: Primary variable (with sparse measurements) search 
*                        ellipsoid definition
*             - Card G2: Extensive variables search ellipsoid def. (common for 
*                        all extensive variables). Only if NVAR_GS>1
*             - Card G3: Definition of search ellipsoids angles (anisotropy 
*                        directions are common for all search ellips.)
*             - Assignation of all variables in array SEARCH_GS
*             - Card G4: Kriging polynomial trend definition. Only if SK or KED
*                        are done
*             - Card G5: Reads variable trimming limits (a value will not be 
*                        considered if VALUE.GE.TRIM_GS) for variables and
*                        external drift & conditional simulation tolerance
*             - Echoes all input variables to results file
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  IPOLDRIFT_GS            Array of 1's and 0's. A polynomial term is considered
*                         if value is 1. 
*                         - 1->3: Linear terms X-Y-Z
*                         - 4->6: Quadratic terms X-Y-Z
*                         - 7->9: Cross-quadratic terms: XY,XZ,YZ
*  SEARCH_GS              Array containing search constants 
*                         - 1: Search radius of primary variable
*                         - 2: Square of 1
*                         - 3: First anisotropy ratio of prim. var. ellip.
*                         - 4: Second anisotropy ratio of prim. var. ellip.
*                         - 5: Search radius of extensive variables
*                         - 6: Square of 5
*                         - 7: First anisotropy ratio of ext. var. ellip.
*                         - 8: Second anisotropy ratio of ext. var. ellip.
*                         - 9: First angle of anisotropy. See user's guide
*                         - 10: Second angle of anisotropy. See user's guide
*                         - 11: Third angle of anisotropy. See user's guide
*                         - 12: Tolerance factor for conditional simulations
*                         
*  TRIM_GS                Trimming limits for all variables+secondary attribute
*
* INTERNAL VARIABLES: ARRAYS
*
*  CHAR                   Character: 'True' or 'False'
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUGEO                  Unit number of GEO file
*  KTYPE_GS               Kriging/SC mode: 0) Simple kriging
*                                          1) Residual kriging
*                                          2) Kriging with locally varying mean
*                                          3) Kriging with external drift
*                                          4) Simple cokriging
*                                          5) Standardized ordinary cokriging
*                                          6) Traditional ordinary cokriging
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NVAR_GS                Number of variables (extensive+1 - the primary -)
*  SKMEAN                 Simple kriging mean
*
* INTERNAL VARIABLES: SCALARS
*
*  ITERM                  Dummy counter variable for polynomial drifts
*  LEAUX                  Auxiliar string where the last read line of the file 
*                         being read is stored
*  NROW                   Current record column number
*  RADIUSP                Dummy variable for reading purposes
*  RADIUSS                Dummy variable for reading purposes
*  RADSQDP                Dummy variable for reading purposes
*  RADSQDS                Dummy variable for reading purposes
*  RADIUS1                Dummy variable for reading purposes
*  RADIUS2                Dummy variable for reading purposes                   
*  SANG1                  Dummy variable for reading purposes
*  SANG2                  Dummy variable for reading purposes
*  SANG3                  Dummy variable for reading purposes
*  SANISP1                Dummy variable for reading purposes
*  SANISP2                Dummy variable for reading purposes
*  SANISS1                Dummy variable for reading purposes
*  SANISS2                Dummy variable for reading purposes
*  T1 -> T4               Dummy variable for reading purposes
*  TEX1 -> TEX4           Dummy variable for reading purposes
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*
* HISTORY: AAR    First coding (Dec-2001)
*          AAR    Revision (Feb-2002)
*          AAR    Inclusion of groups of zones (July-2003)
*
********************************************************************************

       IMPLICIT NONE

C______________________________________________ Step 0: Declaration of variables

       CHARACTER CHAR(9)*5,FILENAME(20)*20,LEAUX*100,LEEL*100

       INTEGER*4 IUGEO,MAINF,INPWR,IOWAR,NVAR_GS,IERROR
     ;          ,KTYPE_GS,NEXDR_GS
     ;          ,IPOLDRIFT_GS(9)

       REAL*8 SKMEAN,TRIM_GS(8),SEARCH_GS(11)

       REAL*8 RADIUSP,RADIUSS,SANG1,SANG2,SANG3,SANISP1,RADIUS1,RADIUS2
     ;       ,SANISP2,SANISS1,SANISS2,RADSQDP,RADSQDS,T1,T2,T3,T4,TEX1
     ;       ,TEX2,TEX3,TEX4

       INTEGER*4 ITERM,NROW

C_________________________ Card G1: Primary variable search ellipsoid definition

       LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
       READ (LEAUX,1000,ERR=9100) RADIUSP,RADIUS1,RADIUS2
 1000  FORMAT(7F10.0)

       RADSQDP=RADIUSP*RADIUSP
       IF (RADIUSP.LE.0D0) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'PRIM. VARIABLE SEARCH RADIUS MUST BE => 0',NROW,1,IUGEO,2,9.12)
                              
       SANISP1=RADIUS1/RADIUSP                              ! Anisotropy factors
       SANISP2=RADIUS2/RADIUSP

C___________________ Card G2: Secondary variables search ellipsoid def. (common)

       IF (NVAR_GS.GT.1) THEN
         LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
         READ (LEAUX,1000,ERR=9200) RADIUSS,RADIUS1,RADIUS2
         IF (RADIUSS.LE.0D0) CALL ERROR
     ;   (IERROR,IOWAR,MAINF,FILENAME,
     ;  'EXT. VARIABLE SEARCH RADIUS MUST BE => 0',NROW,1,IUGEO,2,9.12)

         RADSQDS=RADIUSS*RADIUSS
         SANISS1=RADIUS1/RADIUSS             ! Anisotropy factors
         SANISS2=RADIUS2/RADIUSS
       END IF

C______________________ Card G3: Definition of search ellipsoids (common) angles

       LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
       READ (LEAUX,1000,ERR=9300) SANG1,SANG2,SANG3

C______________ Assignation of all variables defining search ell.in array SEARCH

       SEARCH_GS(1)  = RADIUSP
       SEARCH_GS(2)  = RADIUSP*RADIUSP
       SEARCH_GS(3)  = SANISP1
       SEARCH_GS(4)  = SANISP2
       IF (NVAR_GS.GT.1) THEN
         SEARCH_GS(5)  = RADIUSS
         SEARCH_GS(6)  = RADIUSS*RADIUSS
         SEARCH_GS(7)  = SANISS1
         SEARCH_GS(8)  = SANISS2
       END IF
       SEARCH_GS(9)  = SANG1
       SEARCH_GS(10) = SANG2
       SEARCH_GS(11) = SANG3

C_______________________________Card G4: Simple kriging mean.

       IF (KTYPE_GS.EQ.0) THEN
         LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
         READ (LEAUX,1000,ERR=9400) SKMEAN
       END IF
C_______________________________Card G4BIS: Polynomial drift terms to be considered
C_______________________________         Only if resid. kriging or KED is done

       IF (KTYPE_GS.EQ.1 .OR. KTYPE_GS.EQ.3) THEN
         LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
         READ (LEAUX,1100,ERR=9400) (IPOLDRIFT_GS(ITERM),ITERM=1,9)
 1100    FORMAT(9I5)

C____________________________________________ Checks if IDRIF(I) is out of order

         DO ITERM=1,9
           IF (IPOLDRIFT_GS(ITERM).LT.0) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'POLYNOMIAL TREND INDEX OUT OF ORDER ',NROW,1,IUGEO,2,9.12)
           IF (IPOLDRIFT_GS(ITERM).GE.1) IPOLDRIFT_GS(ITERM)=1
         END DO

       END IF ! KTYPE_GS.EQ.1 .OR. KTYPE_GS.EQ.3

C____ Card G5: Reads variable trimming limits. A value will not be considered if
C____          VALUE.GE.TMAX

       LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
       READ (LEAUX,1200,ERR=9500) T1,T2,T3,T4,TEX1,TEX2,TEX3,TEX4

 1200    FORMAT(8F10.0)

       TRIM_GS(1)=T1
       IF(NVAR_GS.GE.2) TRIM_GS(2)=T2
       IF(NVAR_GS.GE.3) TRIM_GS(3)=T3
       IF(NVAR_GS.EQ.4) TRIM_GS(4)=T4

       IF (KTYPE_GS.EQ.2 .OR. KTYPE_GS.EQ.3) THEN
         TRIM_GS(5)=TEX1
         IF (NEXDR_GS.GE.2) TRIM_GS(6)=TEX2
         IF (NEXDR_GS.GE.3) TRIM_GS(7)=TEX3
         IF (NEXDR_GS.EQ.4) TRIM_GS(8)=TEX4
       END IF

C___________________________________________ Echoes all input variables to mainf

       IF (INPWR.NE.0) THEN

         WRITE(MAINF,2000) 
 2000    FORMAT(///,17X,'PRIMARY VARIABLE SEARCH ELLIP. DEFINITION',/
     ;             ,17X,'======= ======== ====== ====== =========='
     ;             ,//)

         WRITE(MAINF,2100) RADIUSP,SANISP1,SANISP2
 2100  FORMAT(5X,'SEARCH RADIUS (PRIM.+PILOT POINTS)... = ',F10.3,/
     ;       ,5X,'PRIM. VAR. HORIZ. ANISOTROPY RATIO... = ',F10.3,/
     ;       ,5X,'PRIM. VAR. VERTICAL ANISOTROPY RATIO. = ',F10.3)

         IF (NVAR_GS.GT.1) WRITE(MAINF,2200) 
 2200    FORMAT(///,17X,'SECONDARY VAR. SEARCH ELLIP. DEFINITION',/
     ;             ,17X,'========= ==== ====== ====== =========='
     ;             ,//)

         IF (NVAR_GS.GT.1) WRITE(MAINF,2300) RADIUSS,SANISS1,SANISS2
 2300  FORMAT(5X,'SEARCH RADIUS (SECONDARY)............ = ',F10.3,/
     ;       ,5X,'SEC. VAR. HORIZ. ANISOTROPY RATIO.... = ',F10.3,/
     ;       ,5X,'SEC. VAR. VERTICAL ANISOTROPY RATIO.. = ',F10.3,/)

         WRITE(MAINF,2400) 
 2400    FORMAT(///,17X,'SEARCH ELLIP. ANGLES',/
     ;             ,17X,'====== ====== ====== '
     ;             ,//)

         WRITE(MAINF,2500) SANG1,SANG2,SANG3
 2500    FORMAT(5X,'FIRST ANGLE OF SEARCH ELLIPSOID...... = ',F10.3,/
     ;         ,5X,'SECOND ANGLE OF SEARCH ELLIPSOID..... = ',F10.3,/
     ;         ,5X,'THIRD ANGLE OF SEARCH ELLIPSOID...... = ',F10.3)

         IF (KTYPE_GS.EQ.1 .OR. KTYPE_GS.EQ.3) THEN

            WRITE(MAINF,2600) 
 2600       FORMAT(///,17X,'POLYNOMIAL DRIFT TERMS',/
     ;                ,17X,'========== ===== ====='//)

            DO ITERM=1,9
              CHAR(ITERM)='FALSE'
              IF (IPOLDRIFT_GS(ITERM).EQ.0) THEN 
                CHAR(ITERM)='FALSE'
              ELSE IF (IPOLDRIFT_GS(ITERM).EQ.1) THEN 
                CHAR(ITERM)='TRUE '
              END IF
            END DO

            WRITE(MAINF,2700) (CHAR(ITERM),ITERM=1,9)
 2700       FORMAT( 5X,'LINEAR DRIFT IN X.............. = ',A5,/
     ;             ,5X,'LINEAR DRIFT IN Y.............. = ',A5,/          
     ;             ,5X,'LINEAR DRIFT IN Z.............. = ',A5,/          
     ;             ,5X,'QUADRATIC DRIFT IN X........... = ',A5,/          
     ;             ,5X,'QUADRATIC DRIFT IN Y........... = ',A5,/          
     ;             ,5X,'QUADRATIC DRIFT IN Z........... = ',A5,/          
     ;             ,5X,'CROSS-QUADRATIC DRIFT IN X-Y... = ',A5,/          
     ;             ,5X,'CROSS-QUADRATIC DRIFT IN X-Z... = ',A5,/          
     ;             ,5X,'CROSS-QUADRATIC DRIFT IN Y-Z... = ',A5)

         END IF ! KTYPE_GS.EQ.1 .OR. KTYPE_GS.EQ.3

         WRITE(MAINF,2800)
 2800    FORMAT(///,17X,'TRIMMING LIMITS',/
     ;             ,17X,'======== ======',//)

         WRITE(MAINF,2900) T1,T2,T3,T4
 2900    FORMAT(5X,'PRIM. VAR TRIMMING LIMIT................ =',E10.3,/
     ;         ,5X,'EXTENSIVE VAR. 1 TRIMMING LIMIT......... =',E10.3,/
     ;         ,5X,'EXTENSIVE VAR. 2 TRIMMING LIMIT......... =',E10.3,/
     ;         ,5X,'EXTENSIVE VAR. 3 TRIMMING LIMIT......... =',E10.3,/)

         IF (KTYPE_GS.EQ.3) THEN    ! Kriging with external drift
           WRITE(MAINF,3000) TEX1,TEX2,TEX3,TEX4
 3000    FORMAT(5X,'EXTERNAL DRIFT 1 TRIMMING LIMIT......... =',E10.3,/
     ;         ,5X,'EXTERNAL DRIFT 1 TRIMMING LIMIT......... =',E10.3,/
     ;         ,5X,'EXTERNAL DRIFT 1 TRIMMING LIMIT......... =',E10.3,/
     ;         ,5X,'EXTERNAL DRIFT 1 TRIMMING LIMIT......... =',E10.3,/)

         END IF

       END IF
C______________________________ Fatal error calls

       RETURN

 9100  CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING CARD G1 ',NROW,1,IUGEO,2,10.1)

       RETURN 

 9200  CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING CARD G2 ',NROW,1,IUGEO,2,10.2)

       RETURN 

 9300  CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING CARD G3 ',NROW,1,IUGEO,2,10.3)

       RETURN 

 9400  CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING CARD G4 ',NROW,1,IUGEO,2,10.4)

       RETURN 

 9500  CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING CARD G5.1 ',NROW,1,IUGEO,2,10.4)

       RETURN
       END

