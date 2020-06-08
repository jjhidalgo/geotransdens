      SUBROUTINE READ_LIN_COMB
     ;(IDIMWGT  ,IERROR     ,IGROUP     ,INPWR     ,IODIM     ,IOWAR
     ;,ISOT     ,IUGEO      ,MAINF      ,MXGRPZN   ,NLINCMB   ,NTYPAR
     ;,NZPAR    ,FILENAME   ,IGR_ZONE   ,INORPAR   ,IOPT_GS   ,IPNT_PAR
     ;,IVEND    ,IVEST      ,NZONE_PAR  ,PARZ      ,WGT_PAR)

********************************************************************************
*
* PURPOSE Reads a linear combination of zonal parameters. Assigns pointers and 
*         value of zonal parameters
*
* DESCRIPTION Summary:
*
*  - Step 0: Declaration of variables
*  - Step 1: Calculates number of complete lines. Writes main header
*  - Step 2: Reads parameter types, zones and weights defining the linear 
*            combination. Also. some quick checks are done
*    - Step 2.1: Initializes local arrays (just in case)
*    - Step 2.2: Reads param. types, zones and weights
*    - Step 2.3: LOOP OVER RECORDS
*      - Step 2.3.A: Identifies record number
*      - Step 2.3.B: Checks parameter type 'ortography' and finds type of param.
*      - Step 2.3.C: At this point, everything is checked. Checks coherency 
*                    between group type and zone type
*      - Step 2.3.D: Identifies estimation option of zone and group to which the
*                    zone belongs
*      - Step 2.3.E: If group is not estimated deterministically echoes an error 
*                    and stops
*      - Step 2.3.F: If zone is not estimated, echoes a warning
*      - Step 2.3.G: If zonal parameter is estimated, then assigns pointers and 
*                    weight to all zones belonging to actual group
*      - Step 2.3.H: Assigns PARZ to all components related to zones belonging 
*                    to actual group
*    - Step 2.4: Once checked, echoes the information
*  - Step 3: Updates number of lines and records and goes to step 2, reading 
*            uncomplete lines
*
* EXTERNAL VARIABLES: ARRAYS
*

*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  IGR_ZONE               Array containing group of zones to which a zone belongs
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR 
*  IOPT_GS                General options for inverse problem. Each row contains 
*                         information of a given group of zones
*  IPNT_PAR               Array containing pointers to DLT_PAR
*  IVEND                  Array of pointers to IPNT_PAR (last position)
*  IVEST                  Array of pointers to IPNT_PAR (initial position)
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PARZ                   Array containing values of zonal parameters
*  WGT_VAR                Weights defining the linear combination of unknowns
*
* INTERNAL VARIABLES: ARRAYS
*
*  IPINORPAR_AUX          Auxiliar data array with pointer to INORPAR
*  IPNZONE_PAR_AUX        Auxiliar data array with pointer to  NZONE_PAR
*  IZON                   Auxiliar array containing zone numbers defining the 
*                         linear combination, for reading purposes
*  WGT                    Auxiliar array containing weights for reading purpose
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMWGT                Used to dimension WGT_PARIERROR     
*  IGROUP                 Actual group of zones
*  INPWR                  Allows writing on MAIN FILE
*  IODIM                  Maximum dimension of the problem                      
*  IOWAR                  Allows writing warning messages
*  ISOT                   Maximum hydraulic conductivity tensor anisotropy      
*                         degree in the problem.                                
*  IUGEO                  GEO file unit number
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MXGRPZN                Maximum number of groups of zones. FIXED PARAMETER 
*                         used for dimension
*  NLINCMB                Number of terms defining the linear combination
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*
* INTERNAL VARIABLES: SCALARS
*
*  ICOMPLETE              Boolean (1:complete lines are read)
*  ICOMPO                 Actual component at arrays of pointers and weights
*  IDUM                   Dummy counter
*  IESTGROUPZONE          Estimation index of actual group of zones
*  IESTZONE               Estimation index of actual zone
*  IGROUPZONE             Group of zones to which a zone belongs to
*  ILIN                   Dummy counter
*  IPINORPAR              Value of INORPAR
*  IPNZONE_PAR            Value of NZONE_PAR 
*  IRECORD                Dummy counter
*  ISTART                 Initial position at arrays of pointers and weigths
*  ITERM                  Dummy counter
*  ITYP                   Dummy counter      
*  ITYPEGROUP             Type of parameter of actual group of zones
*  ITYPEZONE              Type of parameter of actual zone
*  IZCOMPDLT_PAR          Actual component at array DLT_PAR
*  IZONTYPE               Dummy counter
*  NLINES                 Number of lines to be read
*  NRECORDS               Number of records per line
*  NROW                   Row of actual record
*  NZONTOT                Number of zones previously defined
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  ZERO_ARRAY
*  ZERO_ARRAY_I
*
* HISTORY
*
*     AAR      7-2003     First coding
*
*****************************************************************************


C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NLINCMB,NZPAR,IUGEO,MAINF,INPWR,IOWAR,IERROR,NTYPAR
     ;         ,MXGRPZN,IGROUP,IDIMWGT,ISOT,IODIM
     ;         ,NZONE_PAR(NTYPAR),IOPT_GS(MXGRPZN,20),IGR_ZONE(NZPAR)
     ;         ,IVEST(NZPAR),INORPAR(NTYPAR),IPNT_PAR(NZPAR*IDIMWGT)
     ;         ,IVEND(NZPAR)
                                                                 ! Real external
      REAL*8 WGT_PAR(NZPAR*IDIMWGT),PARZ(NZPAR)
                                                              ! Integer internal
      INTEGER*4 ITERM,ICOMPLETE,NLINES,NRECORDS,ILIN,NROW,IRECORD
     ;         ,ITYP,IPINORPAR,IPNZONE_PAR,ITYPEGROUP,ITYPEZONE
     ;         ,IESTZONE,IGROUPZONE,IESTGROUPZONE,IZONTYPE,ICOMPO
     ;         ,IZCOMPDLT_PAR,IDUM,NZONTOT,ISTART
     ;         ,IZON(4),IPINORPAR_AUX(15),IPNZONE_PAR_AUX(15)
                                                                 ! Real internal
      REAL*8 WGT(4)
                                                                       ! Strings
      CHARACTER FILENAME(18)*20,LEAUX*100,LEEL*100,TYPE(4)*3
     ;         ,TYPEPAR(15)*3

      DATA TYPEPAR/'TRA','STG','ARR','CHP','QQP','ALF','DSL','DST','DFM'
     ;            ,'POR','FOD','CRD','CON','AGE','DMT'/

      DATA IPINORPAR_AUX/1,7,8,9,10,11,12,13,14,15,16,17,18,19,20/
      DATA IPNZONE_PAR_AUX/1,2,3,4,5,6,7,8,9,10,11,12,13,15,16/

C_______________________ Step 1: Calculates number of complete lines. Writes 
C_______________________         main header

      ICOMPLETE=1        ! Flag marking that complete lines are going to be read
      NLINES=INT(NLINCMB/4) ! Number of complete lines to be read
      NRECORDS=4            ! Number of groups of records on each line

      IF (INPWR.NE.0) WRITE(MAINF,2000) NLINCMB
 2000   FORMAT(5X,' NUMBER OF TERMS DEFINING LINEAR COMB.: ',I5
     ;      ,/,5X,' ====== == ===== ======== ====== ======',//
     ;        ,5X,' TYP ZONE    WEIGHT',/,5X,' === ====    ======',/)

C_______________________ Step 2: Reads parameter types, zones and weights 
C_______________________         defining the linear combination. Also. some 
C_______________________         quick checks are done

 10   DO ILIN=1,NLINES

C_______________________ Step 2.1: Initializes local arrays (just in case)

        DO ITERM=1,4
          TYPE(ITERM)='   '
        END DO
        CALL ZERO_ARRAY(WGT,4)
        CALL ZERO_ARRAY_I(IZON,4)

C_______________________ Step 2.2: Reads param. types, zones and weights

        LEAUX=LEEL(FILENAME,IUGEO,MAINF,NROW,INPWR)
        READ(LEAUX,1000,ERR=9000) 
     ;  (TYPE(ITERM),IZON(ITERM),WGT(ITERM),ITERM=1,NRECORDS)
 1000   FORMAT(4(1X,A3,I5,F10.0))

C_______________________ Step 2.3: LOOP OVER RECORDS

        DO ITERM=1,NRECORDS
   
C_______________________ Step 2.3.A: Identifies record number

          IF (ICOMPLETE.EQ.1) THEN
            IRECORD=(NLINES-1)*4+ITERM
          ELSE
            IRECORD=INT(NLINCMB/4)*4+ITERM
          END IF ! ICOMPLETE.EQ.1

C_______________________ Step 2.3.B: Checks parameter type 'ortography' and
C_______________________             finds type of parameter

          DO ITYP=1,15           ! Loop over types of parameters

C_______________________ Checks parameter type 'ortography'

            IF (TYPE(ITERM).EQ.TYPEPAR(ITYP)) THEN  ! Match found
               ITYPEZONE=ITYP
               IPINORPAR=INORPAR(IPINORPAR_AUX(ITYP))
               IPNZONE_PAR=NZONE_PAR(IPNZONE_PAR_AUX(ITYP))
               GOTO 100
            END IF ! TYPE(ITERM).EQ.TYPEPAR(ITYPE)

          END DO ! ITYP=1,NTYPAR_INTERP

C_______________________ Writes an error message due to wrong ortography. 
C_______________________ Forced stop

          WRITE(MAINF,2200) IRECORD
          WRITE(6,2200) IRECORD
 2200     FORMAT(/,' ERROR READING LINEAR COMB. OF PARAMETERS.'
     ;             ' PARAMETER TYPE NAME DOES NOT MATCH',/,' DATA BASE.'
     ;             ' RECORD: ',I5,/,' PLEASE, CHECK ORTOGRAPHY. FORCED'
     ;             ' STOP. BE CAREFUL WITH CAPITAL CHARACTERS.',/,
     ;             ' LIST OF AVAILALE PARAMETERS:',/)
          WRITE(MAINF,2300) (TYPEPAR(ITYP),ITYP=1,15)
          WRITE(6,2300) (TYPEPAR(ITYP),ITYP=1,15)
 2300     FORMAT(20(1X,A3))
          STOP

C_______________________ Ortography matched. Checks zone out of range

 100      IF (IZON(ITERM).LE.0 
     ;        .OR. IZON(ITERM).GT.IPNZONE_PAR) THEN
             WRITE(MAINF,2400) IRECORD
             WRITE(6,2400) IRECORD
 2400        FORMAT(/,' ERROR READING LINEAR COMB. OF PARAMETERS.'
     ;                ' PARAMETER ZONE OUT OF RANGE.',/,
     ;                ' RECORD: ',I5,/)
             STOP
          END IF

C_______________________ Step 2.3.C: At this point, everything is checked.
C_______________________             Checks coherency between group type
C_______________________             and zone type

          ITYPEGROUP=IOPT_GS(IGROUP,1)

          IF (IOWAR.NE.0 
     ;        .AND. TYPEPAR(ITYPEGROUP).NE.TYPEPAR(ITYPEZONE)) 
     ;       WRITE(MAINF,2500) TYPEPAR(ITYPEGROUP),IRECORD
     ;                        ,TYPEPAR(ITYPEZONE)
 2500        FORMAT(//,' WARNING: ACTUAL GROUP HAS A TYPE: ',A3,' WHILE'
     ;                 ' ZONAL PARAMETER ',I5,' IS OF TYPE: ',A3,/,
     ;                 ' CONCEPTUALLY, THIS IS NOT CORRECT',/)
             
C_______________________ Step 2.3.D: Identifies estimation option of zone
C_______________________             and group to which the zone belongs

           IESTZONE=IVEST(IPINORPAR+IZON(ITERM))
           IGROUPZONE=IGR_ZONE(IPINORPAR+IZON(ITERM))
           IESTGROUPZONE=IOPT_GS(IGROUPZONE,2)
      
C_______________________ Step 2.3.E: If group is not estimated deterministically
C_______________________             echoes an error and stops

           IF (IESTGROUPZONE.NE.0) THEN
              WRITE(MAINF,2600) IRECORD,IGROUPZONE
 2600         FORMAT(/,' RECORD: ',I5,' IS A ZONAL PARAMETER BELONGING'
     ;                 ' TO GROUP OF ZONES: ',I5,/,' THAT GROUP IS NOT'
     ;                 ' ESTIMATED DETERMINISTICALLY. FORCED STOP',/)
              STOP ' FORCED STOP. CHECK RES.OUT'
           END IF

C_______________________ Step 2.3.F: If zone is not estimated, echoes a warning
 
           IF (IESTZONE.EQ.0) THEN

              WRITE(MAINF,2700) IRECORD
 2700         FORMAT(' WARNING: RECORD ',I5,' IS A ZONAL PARAMETER'
     ;               ' THAT IS NOT ESTIMATED')

C_______________________ Step 2.3.G: If zonal parameter is estimated, then
C_______________________             assigns pointers and weight to all zones
C_______________________             belonging to actual group

           ELSE ! Zonal parameter is estimated

                                        ! Identifies actual component at DLT_PAR

              IZCOMPDLT_PAR=IPNT_PAR(IESTZONE)

   ! Identifies initial position of type of parameters of this group at IPNT_PAR

              NZONTOT=0
              DO IDUM=1,ITYPEGROUP-1
                IF (IDUM.EQ.1) THEN
                  NZONTOT=NZONE_PAR(1)*MAX(ISOT,IODIM)
                ELSE
                  NZONTOT=NZONTOT+NZONE_PAR(IDUM)
                END IF
              END DO ! IDUM=1,ITYPEGROUP
              ISTART=NZONTOT*IDIMWGT

              DO IZONTYPE=1,NZONE_PAR(ITYPEGROUP)
                IF (IGR_ZONE(IPINORPAR+IZONTYPE).EQ.IGROUP) THEN 

                  ICOMPO=ISTART+(IZONTYPE-1)*IDIMWGT+1
                  IF (ITYPEGROUP.EQ.1) 
     ;          ICOMPO=ISTART+(IZONTYPE-1)*MAX0(ISOT,IODIM)*IDIMWGT+1
 
                  IPNT_PAR(ICOMPO)=IZCOMPDLT_PAR
                  WGT_PAR(ICOMPO)=WGT(ITERM)

                  IF (IVEST(IPINORPAR+IZONTYPE).EQ.0 .OR. 
     ;                ICOMPO.LE.IVEST(IPINORPAR+IZONTYPE))
     ;                IVEST(IPINORPAR+IZONTYPE)=ICOMPO

                  IF (IVEND(IPINORPAR+IZONTYPE).EQ.0 .OR. 
     ;                ICOMPO.GE.IVEND(IPINORPAR+IZONTYPE))
     ;                IVEND(IPINORPAR+IZONTYPE)=ICOMPO

                END IF ! Zone belongs to group
              END DO ! IZONGROUP=1,NZONE_PAR(ITYPEGROUP)

           END IF ! Zonal parameter estimated

C_______________________ Step 2.3.H: Assigns PARZ to all components related
C_______________________             to zones belonging to actual group

            DO IZONTYPE=1,NZONE_PAR(ITYPEGROUP)
              IF (IGR_ZONE(IPINORPAR+IZONTYPE).EQ.IGROUP) THEN 

                 PARZ(IPINORPAR+IZONTYPE)=
     ;             PARZ(IPINORPAR+IZONTYPE)
     ;                 +WGT(ITERM)*PARZ(IPINORPAR+IZON(ITERM))
                 
              END IF ! Zone belongs to group
            END DO ! IZONGROUP=1,NZONE_PAR(ITYPEGROUP)

        END DO ! ITERM=1,NRECORDS

C_______________________ Step 2.4: Once checked, echoes the information

        IF (INPWR.NE.0) WRITE(MAINF,2100) 
     ;      (TYPE(ITERM),IZON(ITERM),WGT(ITERM),ITERM=1,NRECORDS)
 2100        FORMAT(6X,A3,I5,F10.3)

      END DO ! ILIN=1,NLINES

C_______________________ Step 3: Updates number of lines and records and goes
C_______________________         to step 2, reading uncomplete lines

      IF (ICOMPLETE.EQ.1) THEN   ! All complete lines have been read
        ICOMPLETE=0
        NRECORDS=NLINCMB-NLINES*4
        NLINES=1
        GOTO 10
      END IF ! ICOMPLETE.EQ.1

      RETURN      

 9000 CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;  'GENERIC FORTRAN ERROR READING LINEAR COMBINATION'
     ;  ,NROW,1,IUGEO,2,1.0)

      RETURN
      END
