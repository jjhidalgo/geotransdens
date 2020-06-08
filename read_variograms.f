       SUBROUTINE READ_VARIOGRAMS
     ;(IDIM     ,IERROR   ,INPWR    ,IOWAR    ,IUGEO    ,KRIGTYPE 
     ;,MAINF    ,MAXNST   ,NVAR     ,AA       ,ANG1
     ;,ANG2     ,ANG3     ,ANIS1    ,ANIS2    ,C0       ,CC       
     ;,IT       ,NST      ,FILENAME)

********************************************************************************
*
* PURPOSE Reads card G9, containing the definition of the corregionalization
*         model (variogram definition) of a given group of zones
*
* DESCRIPTION This subroutine works as summarized below:
*             Step 0: Declaration of variables
*             Step 1: Initialisation of NST value to -1 to flag all variograms
*                     (Only sense when cokriging is performed)
*             Step 2: Card G9.1. Number of nested structures and general nugget
*                     Reads as many nested structures of current variogram:
*                       Card G9.2. Structure type, parameters AA,AA1,AA2 and CC
*                                  Calculation of anisotropy ratios.Checks power
*                                  model requirements
*                       Card G9.3. Rotation angles of current structure
*                     Reads IFLAG: 0 if there are more variogram definitions.
*            Step 3: Fill in variograms j=i if they have not been explicitly
*                    entered. (Only sense when cokriging is performed)
*            Step 4: Checks linearity of corregionalization model(Matching
*                     number of nested structures and parameters, variogram type
*                     and parameters (except CC) for all variogram combinations
*            Step 5: If a linear model of corregionalization is being used,
*                    checks to ensure positive definiteness.
*            Step 6: Variogram parameters are echoed to MAINF (if allowed)
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AA                     Array containing AA parameter for all variograms
*  ANG1                   Array containing first angle of search ellipsoid for
*                         all variograms
*  ANG2                   Array containing first angle of search ellipsoid for
*                         all variograms
*  ANG3                   Array containing first angle of search ellipsoid for
*                         all variograms
*  ANIS1                  Array containing first anis. ratio of search ellipsoid
*                         for all variograms
*  ANIS2                  Array containing second ani. ratio of search ellipsoid
*                         for all variograms
*  C0                     Array containing general nugget for all variograms
*  CC                     Array containing CC parameter for all variograms
*  IT                     Array containing nested variog. models
*  NST                    Array containing number of nested struc. of all vari.
*                         data files                                            
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIM                   Used to dimension variogram arrays
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUGEO                  GEO file unit number
*  KRIGTYPE               Kriging/cokriging type
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MAXNST                 Maximum number of nested structures defining all vari.
*  NVAR                   Number of variables (prim. + all secondary)
*
* INTERNAL VARIABLES: SCALARS
*
*  AA1                    Dummy variable for reading purposes
*  AA2                    Dummy variable for reading purposes
*  IFLAG                  If different from 0, it writes storage partition
*  IND                    Variogram identifier
*  IND1                   Variogram identifier
*  IND2                   Variogram identifier
*  INDEX                  Variogram identifier
*  INDEX1                 Variogram identifier                                  
*  INDEX2                 Variogram identifier                                  
*  INEST                  Dummy counter variable of nested struc.
*  ISTART                 Variogram identifier                                  
*  ISTART1                Variogram identifier                                  
*  ISTART2                Variogram identifier                                  
*  IVAR                   First variable defining variogram
*  JVAR                   Second. var. defining variogram
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*  NROW                   Current record number                                 
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*
* HISTORY: AAR First coding (Dec-2001)
*          AAR Revision     (Jan-2002)
*
********************************************************************************


      IMPLICIT NONE
C______________________________________________ Step 0: Declaration of variables

      LOGICAL LINMOD,POSDEF

      INTEGER*4 KRIGTYPE,NVAR,MAXNST,IERROR,IOWAR,MAINF,IUGEO,NROW
     ;         ,INPWR,IDIM,IT(IDIM),NST(IDIM)

      REAL*8 C0(IDIM),AA(IDIM),CC(IDIM),ANG1(IDIM),ANG2(IDIM),ANG3(IDIM)
     ;      ,ANIS1(IDIM),ANIS2(IDIM)
 
      INTEGER*4 IVAR,JVAR,INEST,IND,IND1,IND2,IFLAG,ISTART,ISTART1,
     ;          ISTART2,INDEX,INDEX1,INDEX2,II,JJ,IJ,JI,IVAR2,JVAR2
     ;         ,ISTARTII,ISTARTJJ,ISTARTIJ,ISTARTJI,INDEXII
     ;         ,INDEXJJ,INDEXIJ,INDEXJI

      REAL*8 AA1,AA2,EPSILON

      CHARACTER VARIONAME(5)*17,VARNAME(5)*11, FILENAME(20)*20,LEAUX*100
     ;         ,LEEL*100
      DATA VARIONAME/'SPHERICAL MODEL  ','EXPONENTIAL MODEL',
     ;               'GAUSSIAN MODEL   ','POWER MODEL      ',
     ;               'HOLE EFFECT MODEL'/
      DATA VARNAME/'PRIMARY    ','EXTENSIVE 1',
     ;             'EXTENSIVE 2','EXTENSIVE 3','EXTENSIVE 4'/

      PARAMETER (EPSILON=1D-6)

C______________ Step 1: Initialisation of NST value to -1 to flag all variograms
C______________         (Only sense when cokriging is performed)

      IF (KRIGTYPE.GE.4) THEN
        DO IVAR=1,NVAR
          DO JVAR=1,NVAR
            IND=IVAR+(JVAR-1)*NVAR
            NST(IND)=-1
          END DO
        END DO
      END IF

C__________________ Step 2: Card G5:Reads as many variograms as are in GEO file.

      IFLAG=0
      DO WHILE (IFLAG.GE.0)
                                     ! Gamma (IVAR,JVAR) is going to be defined
        LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
        READ (LEAUX,1000,ERR=9500) IVAR,JVAR
 1000   FORMAT(2I5)

        IF (IVAR.GT.NVAR.OR.JVAR.GT.NVAR) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'ILLEGAL VALUE READING CARD G5 ',NROW,1,IUGEO,2,10.5)

        IND=IVAR+(JVAR-1)*NVAR
        ISTART=1+(IND-1)*MAXNST           ! Starting point of current GAMMA(I,J)

                     ! Card G5.1. Number of nested structures and general nugget

        LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
        READ (LEAUX,1100,ERR=9600) NST(IND),C0(IND)
 1100   FORMAT(I5,F10.0)

        IF(NST(IND).GT.MAXNST) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'NUMBER OF VARIOG. STRUCT. EXCEEDS MXNST CARD G5.1 '
     ; ,NROW,1,IUGEO,2,10.6)

                                            ! Loop reading all nested structures
        DO INEST=1,NST(IND)
          INDEX=ISTART+INEST-1
                               ! Card G5.2. Structure type, parameters AA and CC

          LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
          READ (LEAUX,1200,ERR=9700) IT(INDEX),AA(INDEX),AA1,AA2
     ;                              ,CC(INDEX)
 1200     FORMAT(I5,4F10.0)

                                                      ! Power model requirements

          IF (IT(INDEX).EQ.4.AND.KRIGTYPE.EQ.4) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'NO POWER MODEL WITH SIMPLE COKRIGING ',NROW,1,IUGEO,2,10.7)

          IF (IT(INDEX).EQ.4.AND.(AA(INDEX).LE.0D0.OR.
     ;                               AA(INDEX).GT.2D0)) CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'POWER MODEL. AA MUST BE >0 AND <=2',NROW,1,IUGEO,2,10.7)

          ANIS1(INDEX)=AA1/DMAX1(AA(INDEX),1.D-6)           ! Anisotropy factors
          ANIS2(INDEX)=AA2/DMAX1(AA(INDEX),1.D-6)

                               ! Card G5.3. Rotation angles of current structure
 
          LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
          READ (LEAUX,1300,ERR=9800) ANG1(INDEX),ANG2(INDEX)
     ;                              ,ANG3(INDEX)
 1300     FORMAT(5F10.0)           
          
        END DO                                           ! Next Nested structure

                                                                  ! Reads IFLAG
        LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
        READ (LEAUX,1400,ERR=9900) IFLAG
 1400   FORMAT(I5)                
      END DO                                         ! Next variogram definition

C_______________ Step 3: Fill in variograms j=i if they have not been explicitly
C_______________         entered. (Only sense when cokriging is performed)

      IF (KRIGTYPE.GE.4) THEN
        DO IVAR=1,NVAR
          DO JVAR=1,NVAR
            IND1=IVAR+(JVAR-1)*NVAR
            IND2=JVAR+(iVAR-1)*NVAR

                                                           ! Variogram not found

            IF (NST(IND1).EQ.-1.AND.NST(IND2).EQ.-1) THEN 
              WRITE(MAINF,2000) IVAR,JVAR
 2000         FORMAT(///,' NEEDED VARIOGRAM BETWEEN VARIABLES ',I5,
     ;                   'AND ',I5,'. FORCED STOP')
              STOP
            END IF
                                                          ! J defined, but not I
            IF (NST(IND1).EQ.-1) THEN
              NST(IND1)=NST(IND2)
              C0(IND1)=C0(IND2)
              ISTART1=1+(IND1-1)*MAXNST
              ISTART2=1+(IND2-1)*MAXNST
              DO INEST=1,NST(IND1)
                INDEX1=ISTART1+INEST-1
                INDEX2=ISTART2+INEST-1
                IT(INDEX1)=IT(INDEX2)
                CC(INDEX1)=CC(INDEX2)
                AA(INDEX1)=AA(INDEX2)
                ANG1(INDEX1)=ANG1(INDEX2)
                ANG2(INDEX1)=ANG2(INDEX2)
                ANG3(INDEX1)=ANG3(INDEX2)
                ANIS1(INDEX1)=ANIS1(INDEX2)
                ANIS2(INDEX1)=ANIS2(INDEX2)
              END DO
                                                          ! I defined, but not J
            ELSE IF (NST(IND2).EQ.-1) THEN
              NST(IND2)=NST(IND1)
              C0(IND2)=C0(IND1)
              ISTART1=1+(IND1-1)*MAXNST
              ISTART2=1+(IND2-1)*MAXNST
              DO INEST=1,NST(IND1)
                INDEX1=ISTART1+INEST-1
                INDEX2=ISTART2+INEST-1
                IT(INDEX2)=IT(INDEX1)
                CC(INDEX2)=CC(INDEX1)
                AA(INDEX2)=AA(INDEX1)
                ANG1(INDEX2)=ANG1(INDEX1)
                ANG2(INDEX2)=ANG2(INDEX1)
                ANG3(INDEX2)=ANG3(INDEX1)
                ANIS1(INDEX2)=ANIS1(INDEX1)
                ANIS2(INDEX2)=ANIS2(INDEX1)
              END DO
            ELSE
              CONTINUE
            END IF

          END DO                        ! IVAR
        END DO                          ! JVAR
      END IF                            ! Is cokriging being performed

C_____________ Step 4: Check for a linear model of corregionalization. Matching
C_____________         number of nested structures and parameters (except CC).
C_____________         All possible combinations of variograms are considered

      IF (NVAR.GT.1) THEN
        LINMOD=.TRUE.
        DO IVAR=1,NVAR
          DO JVAR=1,NVAR
                                                    ! Identifies first variogram
            IND1=IVAR+(JVAR-1)*NVAR

            DO IVAR2=1,NVAR
              DO JVAR2=1,NVAR
                                                   ! Identifies second variogram
                IND2=IVAR2+(JVAR2-1)*NVAR
                                                             ! Nested structures
                IF (NST(IND1).NE.NST(IND2)) THEN
                  LINMOD=.FALSE.
                  GOTO 10
                END IF

                ISTART1=1+(IND1-1)*MAXNST              
                ISTART2=1+(IND2-1)*MAXNST              

                DO INEST=1,NST(IND1)
                  INDEX1=ISTART1+INEST-1
                  INDEX2=ISTART2+INEST-1
                                                               ! Variogram types
                  IF (IT(INDEX1).NE.IT(INDEX2)) THEN
                    LINMOD=.FALSE.
                    GOTO 10
                  END IF
                                                                  ! Parameter AA
                  IF (DABS(AA(INDEX1)-AA(INDEX2)).GT.EPSILON) THEN
                    LINMOD=.FALSE.
                    GOTO 10
                  END IF
                                                                       ! Angles
                  IF (DABS(ANG1(INDEX1)-ANG1(INDEX2)).GT.EPSILON
     ;            .OR.DABS(ANG2(INDEX1)-ANG2(INDEX2)).GT.EPSILON 
     ;            .OR.DABS(ANG3(INDEX1)-ANG3(INDEX2)).GT.EPSILON) THEN
                    LINMOD=.FALSE.
                    GOTO 10
                  END IF
                                                            ! Anisotropy ratios
                  IF (DABS(ANIS1(INDEX1)-ANIS1(INDEX2)).GT.EPSILON
     ;           .OR.DABS(ANIS2(INDEX1)-ANIS2(INDEX2)).GT.EPSILON) THEN
                    LINMOD=.FALSE.
                    GOTO 10
                  END IF

                END DO                                         ! INEST
              END DO                                           ! JVAR1
            END DO                                             ! IVAR1
          END DO                                               ! JVAR
        END DO                                                 ! IVAR

 10     IF (LINMOD) THEN

C___________________ Step 5: A linear model of corregionalization is being used
C___________________         Checks to ensure positive definiteness.

          POSDEF=.TRUE.
          DO IVAR=1,NVAR
            DO JVAR=1,NVAR
              IF (IVAR.EQ.JVAR) GOTO 11
              II=IVAR+(IVAR-1)*NVAR                      ! Gamma (I,I) position
              JJ=JVAR+(JVAR-1)*NVAR                      ! Gamma (J,J) position
              IJ=IVAR+(JVAR-1)*NVAR                      ! Gamma (I,J) position
              JI=JVAR+(IVAR-1)*NVAR                      ! Gamma (J,I) position
              ISTARTII=1+(II-1)*MAXNST
              ISTARTJJ=1+(JJ-1)*MAXNST
              ISTARTIJ=1+(IJ-1)*MAXNST
              ISTARTJI=1+(JI-1)*MAXNST
                                                    ! Checks for nugget effects
              IF (C0(II).LE.0D0.OR.C0(JJ).LE.0D0.OR.
     ;            C0(II)*C0(JJ).LT.C0(IJ)*C0(JI)) THEN
                 POSDEF=.FALSE.
                 GOTO 20
              END IF
                                                      ! Checks variogram ranges
              DO INEST=1,NST(II)
                INDEXII=ISTARTII+INEST-1
                INDEXJJ=ISTARTJJ+INEST-1
                INDEXIJ=ISTARTIJ+INEST-1
                INDEXJI=ISTARTJI+INEST-1
  
                IF (CC(INDEXII).LE.0D0.OR.CC(INDEXJJ).LE.0D0.OR.
     ;              CC(INDEXII)*CC(INDEXJJ).LT.
     ;              CC(INDEXIJ)*CC(INDEXJI)) THEN
                   POSDEF=.FALSE.
                   GOTO 20
                END IF
              END DO                                           ! INEST
 11           CONTINUE
            END DO                                             ! JVAR
          END DO                                               ! IVAR

 20       IF (.NOT.POSDEF.AND.IOWAR.NE.0)

C_______________________ Step 5.1: The linear model of corregionalization is not
C_______________________           positive definite (warning).

     ;      CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'THE LINEAR MODEL OF CORREGIONALIZATION IS NOT POSITIVE '
     ;//'DEFINITE. THIS COULD LEAD TO SINGULAR KRIGING MATRICES AND '
     ;//'UNESTIMATED POINTS.',NROW,1,IUGEO,0,10.4)


        ELSE

C________________ Step 5: A non linear model of corregionalization is being used
C________________         (Warning)

           CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'A LINEAR MODEL OF CORREGIONALIZATION HAS NOT BEEN USED !!! '
     ;//'THIS COULD LEAD TO MANY SINGULAR KRIGING MATRICES AND '
     ;//'UNESTIMATED POINTS.',NROW,1,IUGEO,0,10.4)

        END IF

      END IF ! NVAR.GT.1

C__________________________ Step 6: Variogram parameters are now echoed to MAINF

      IF (INPWR.NE.0) THEN 
        WRITE(MAINF,2100)
 2100   FORMAT(///,27X,'CORREGIONALIZATION MODEL INFORMATION',/
     ;            ,27X,'================== ===== ===========',/)

        DO IVAR=1,NVAR
          DO JVAR=1,NVAR
            IND=IVAR+(JVAR-1)*NVAR
            WRITE(MAINF,2200) VARNAME(IVAR),VARNAME(JVAR),NST(IND)
     ;                       ,C0(IND)

 2200       FORMAT(/,10X,'VARIOGRAM FOR VARIABLES ',A11,' AND ',A11,/
     ;              ,10X,'========= === ========= ==========='
     ;                   ' === ===========',//
     ;              ,5X,'NUMBER OF STRUCTURES... = ',I5,/
     ;              ,5X,'GENERAL NUGGET......... = ',F10.3)

            ISTART=1+(IND-1)*MAXNST
            DO INEST=1,NST(IND)
              INDEX=ISTART+INEST-1
              WRITE(MAINF,2300) INEST,VARIONAME(IT(INDEX)),AA(INDEX)
     ;                         ,CC(INDEX),ANG1(INDEX)
     ;                         ,ANG2(INDEX),ANG3(INDEX)
     ;                         ,ANIS1(INDEX),ANIS2(INDEX)

 2300         FORMAT(//,5X,'NESTED STRUCTURE ',I5,/
     ;                 ,5X,'====== ========= =====',//
     ;                 ,5X,'STRUCTURE TYPE........ = ',A17,/
     ;                 ,5X,'PARAMETER AA.......... = ',F10.3,/ 
     ;                 ,5X,'PARAMETER CC.......... = ',F10.3,/
     ;                 ,5X,'FIRST ANIS. ANGLE..... = ',F10.3,/
     ;                 ,5X,'SECOND ANIS. ANGLE.... = ',F10.3,/
     ;                 ,5X,'THIRD ANIS. ANGLE..... = ',F10.3,/
     ;                 ,5X,'FIRST ANIS. RATIO..... = ',F10.3,/
     ;                 ,5X,'SECOND ANIS. RATIO.... = ',F10.3)
            END DO                                             ! INEST
          END DO                                               ! JVAR
        END DO                                                 ! IVAR
      END IF              
     
     
      RETURN 

 9500 CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING CARD G9 ',NROW,1,IUGEO,2,1.5)
      RETURN
      
 9600 CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING CARD G9.1 ',NROW,1,IUGEO,2,1.6)
      RETURN

 9700 CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING CARD G9.2 ',NROW,1,IUGEO,2,1.7)
      RETURN

 9800 CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING CARD G9.3 ',NROW,1,IUGEO,2,1.8)
      RETURN

 9900 CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING CARD G9.4 ',NROW,1,IUGEO,2,1.9)
      RETURN

      END
