      SUBROUTINE GET_ZONE_GROUP
     ;(IACTGROUP  ,IACTTYPE   ,IACTZONE   ,IERROR   ,IGROUP   ,IOWAR
     ;,IPLXPAREL  ,ITYPE      ,IUGEO      ,LMXNDL   ,MAINF    ,NPAREL   
     ;,NUMEL      ,NUMNP      ,NTYPAR     ,NZPAR    ,XPOINT   ,YPOINT   
     ;,ZPOINT     ,AREA       ,COORD      ,FILENAME ,IGR_ZONE ,INORPAR  
     ;,KXX        ,LTYPE      ,LXPAREL)

********************************************************************************
*
* PURPOSE Given the coordinates of a point (sampling location or pilot point)
*         actual parameter type, checks if point belongs to domain and if point
*         belongs to actual group of zones
*
* DESCRIPTION Flow chart:
*
*  - Step 0: Declaration of variables
*  - Step 1: Identifies pointers to TRANSIN arrays, given the type of parameter 
*            of actual group
*
*  MAIN LOOP OVER LOCATIONS. 
*     - Step 2: Identifies to which element point belongs to and checks if point 
*               belongs to the domain
*     - Step 3: ERROR:Point does not belong to the domain.
*     - Step 4: Finds out ZONE to which NEL belongs to. 
*     - Step 5: Echoes an error if zone is 0 (does not exist). No coherency 
*               between parameter type and definition of zones
*     - Step 6: WARNING:Point does not belong to the group.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  COORD                  Nodal coordinates                                     
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  IGR_ZONE               Array containing the relation zone->group (IVPAR(4))
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARZ, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LTYPE                  Vector containing the type of each element            
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*
* INTERNAL VARIABLES: ARRAYS
*
*  BF                     Values of basis functions
*
* EXTERNAL VARIABLES: SCALARS
*
*  IACTGROUP              Group to which point belongs to
*  IACTTYPE               Type of parameters of actual group
*  IACTZONE               Zone to which point belongs to
*  IERROR                 Current number of errors on input data                
*  IGROUP                 Actual group
*  IOWAR                  Program echoes warnings if IOWAR.NE.0
*  IPLXPAREL              Pointer to LXPAREL
*  ITYPE                  1: meas. locations; 2: pilot points locations
*  IUGEO                  GEO file unit number
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NPAREL                 Number of element parameters in current problem       
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  XPOINT,YPOINT,ZPOINT   Point coordinates
*
* INTERNAL VARIABLES: SCALARS
*
*  IPINORPAR              Pointer to array INORPAR
*  IPLXPAREL              Pointer to LXPAREL
*  IPOS                   Initial position at IGR_ZONE
*  NEL                    Element to which candidate p.point belongs to
*  NROW                   Current record number                                 
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  BASISFUNC_OBS          Determines to which element a pilot points belongs to
*
* HISTORY:  AAR   First coding (July-2003)
*
********************************************************************************
      
C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE

      INTEGER*4 IACTTYPE,ITYPE,IGROUP,LMXNDL,MAINF,NUMEL,NUMNP,IERROR
     ;         ,IOWAR,IUGEO,IACTZONE,NTYPAR,NZPAR,IACTGROUP,NPAREL
     ;         ,IPLXPAREL
     ;         ,KXX(LMXNDL,NUMEL),LTYPE(NUMEL),INORPAR(NTYPAR)
     ;         ,IGR_ZONE(NZPAR),LXPAREL(NUMEL,NPAREL)
      REAL*8 XPOINT,YPOINT,ZPOINT
     ;      ,AREA(NUMEL),COORD(NUMNP,3)
      INTEGER*4 IPINORPAR,NEL,NROW,IPOS
      REAL*8 BF(9)
      CHARACTER FILENAME(20)*20

C_______________________ Step 1: Identifies pointers to TRANSIN arrays, given 
C_______________________         type of parameter of actual group

      IPINORPAR=1     ! Pointer to INORPAR
      IPLXPAREL=1     ! Pointer to LXPAREL
      IF (IACTTYPE.EQ.2) THEN  ! Storage coeff.
         IPINORPAR=7
         IPLXPAREL=IACTTYPE
      ELSE IF (IACTTYPE.EQ.3 .OR. IACTTYPE.EQ.4) THEN  ! Areal recharge
         IPINORPAR=8
         IPLXPAREL=IACTTYPE
      ELSE IF (IACTTYPE.GT.4) THEN ! Others
         IPINORPAR=IACTTYPE+7
         IF (IACTTYPE.EQ.5 .OR. IACTTYPE.EQ.6) THEN
            IPLXPAREL=5
         ELSE
            IPLXPAREL=IACTTYPE-1
         END IF
      END IF

C________________________ Step 2: Finds out element to which point belongs to

      CALL BASISFUNC_OBS
     ;(LMXNDL     ,MAINF      ,NEL      ,NUMEL    ,NUMNP    ,XPOINT
     ;,YPOINT     ,ZPOINT     ,AREA     ,BF       ,KXX      ,LTYPE
     ;,COORD(1,1) ,COORD(1,2) ,COORD(1,3))       

C________________________ Step 3: ERROR:Point does not belong to the domain.

      IF (NEL.LE.0 .OR. NEL.GT.NUMEL) THEN

        IF (ITYPE.EQ.1) THEN   

          CALL ERROR           ! Location is a sampling point
     ; (IERROR,IOWAR,MAINF,FILENAME
     ; ,'MEASUREMENT LOCATION DOES NOT BELONG TO THE DOMAIN'
     ; ,NROW,1,IUGEO,2,9.8)
          
        ELSE                   ! Location is a pilot point
          CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME
     ; ,'PILOT POINT LOCATION DOES NOT BELONG TO THE DOMAIN'
     ; ,NROW,1,IUGEO,2,9.8)

        END IF ! ITYPE.EQ.1

      END IF ! NEL.LE.0 .OR. NEL.GT.NUMEL

C________________________ Step 4:Finds out zone 

      IACTZONE=LXPAREL(NEL,IPLXPAREL)                                  ! Zone
      IPOS=INORPAR(IPINORPAR)                  ! Initial position at IGR_ZONE

C________________________ Step 5: Echoes an error if zone is 0 (does not 
C________________________         exist). No coherency between parameter 
C________________________         type and definition of zones

      IF (IACTZONE.EQ.0) THEN
          CALL ERROR           ! Location is a sampling point
     ; (IERROR,IOWAR,MAINF,FILENAME
     ; ,'A MEASUREMENT-PILOT POINT OF ACTUAL PARAMETER TYPE IS DEFINED'
     ;//' INSIDE A NON-EXISTING ZONE',NROW,2,IUGEO,2,9.8)
      END IF        

C________________________ Step 6:Point does not belong to the group.

        IACTGROUP=IGR_ZONE(IPOS+IACTZONE)                      ! Actual group
        IF (IACTGROUP.NE.IGROUP) THEN

           IF (ITYPE.EQ.1) THEN   
             CALL ERROR           ! Location is a sampling point
     ; (IERROR,IOWAR,MAINF,FILENAME
     ; ,'MEASUREMENT LOCATION DOES NOT BELONG TO ACTUAL GROUP'
     ; ,NROW,1,IUGEO,0,9.8)
          
           ELSE                   ! Location is a pilot point
               CALL ERROR
     ; (IERROR,IOWAR,MAINF,FILENAME
     ; ,'PILOT POINT LOCATION DOES NOT BELONG TO ACTUAL GROUP'
     ; ,NROW,1,IUGEO,0,9.8)

           END IF ! ITYPE.EQ.1

           WRITE(MAINF,*)
           WRITE(MAINF,*) XPOINT,YPOINT,ZPOINT
           WRITE(MAINF,*)
           
        END IF ! IACTGROUP.NE.IGROUP


      RETURN
      END
