      SUBROUTINE RANDOM_PIPO
     ;(IACTTYPE    ,IGROUP    ,IORD_PP_GS    ,LMXNDL    ,MAINF
     ;,MXMEASPP_GS ,NFLAGS    ,NPAREL        ,NPP_GS    ,NTYPAR     
     ;,NUMEL       ,NUMNP     ,NZPAR         ,AREA      ,COORD     
     ;,COORDGR_GS  ,IFLAGS    ,IGR_ZONE      ,INORPAR   ,KXX       
     ;,LTYPE       ,LXPAREL   ,POSMEAS_GS)

********************************************************************************
*
* PURPOSE  Draws randomly the position of the pilot positions. They are drawn
*          absolutelly random (IORD_PP=1) or layed on a regular support mesh, 
*          given densities of pilot points (directions X-Y-Z) (IORD_PP=2;NOT YET)
*
* DESCRIPTION Algorithm can be summarized as follows:
*
*             - Step 0: Declaration of variables
*             - Step 1: Initialisation of auxiliar array for random number 
*                       generation. Pointers assignment
*             - Step 2: Identifies maximum and minimum and maximum coordinates
*                       of actual group of zones, calculated at ZONE_GEOMETRY
*             - MAIN LOOP OVER PILOT POINTS
*               - Step 3: Generates three random factors multiplying group sizes
*               - Step 4: Generates a random candidate
*               - Step 5: Checks that random candidate belongs to actual group
*                         If pilot point is outside the group, return to Step 3
*               - Step 6: Assigns pilot point coordinates
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  COORD                  Nodal coordinates                                     
*  COORDGR_GS             Array containing maximum and minimum coord. of all 
*                         groups of zones
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  IGR_ZONE               Array containing the relation zone->group (IVPAR(4))
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARZ, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LTYPE                  Vector containing the type of each element            
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  POSMEAS_GS             Array containing pilot points and sampling locations
*
* INTERNAL VARIABLES: ARRAYS
*
*  ISEED                  Auxiliar array for generate random numbers  
*  BF                     Values of basis functions
*  RANDOM                 Three random factors, multiplying grid sizes
*
* EXTERNAL VARIABLES: SCALARS
*
*  IACTTYPE               Type of parameters of actual group
*  IGROUP                 Actual group
*  IORD_PP_GS             Option for drawing randomly the pilot points
*                         - 0: fixed and read
*                         - 1: completelly random
*                         - 2: random but equally distributed, lating on a 
*                              regular support mesh
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MXMEASPP_GS            Used to dimension POSMEAS_GS
*  NFLAGS                 Used to dimension IFLAGS
*  NPAREL                 Number of element parameters in current problem       
*  NPP_GS                 Number of pilot points defining actual group
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*
* INTERNAL VARIABLES: SCALARS
*
*  IACTGROUP              Group to which candidate p.point belongs to
*  IACTZONE               Zone to which candidate p.point belongs to
*  IDUM                   Dummy counter
*  IPINORPAR              Pointer to array INORPAR
*  IPIPO                  Dummy counter of pilot points
*  IPLXPAREL              Pointer to LXPAREL
*  IPOS                   Initial position at IGR_ZONE
*  ITRIAL                 Counter of trials
*  NEL                    Element to which candidate p.point belongs to
*  XDUM                   Dummy for random number generation
*  XMAX_GR,YMAX_GR,ZMAX_GR  Maximum coordinates of actual group
*  XMIN_GR,YMIN_GR,ZMIN_GR  Minimum coordinates of actual group
*  XPIPO,YPIPO,ZPIPO      Candidate pilot point coordinates
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ACORNI                 Random number generator of twelveth order
*  BASISFUNC_OBS          Determines to which element a pilot points belongs to
*
* HISTORY:  AAR   First coding (Feb-2002)
*           AAR   Inclusion of groups of zones (July-2003)
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 IORD_PP_GS,NPP_GS,IACTTYPE,LMXNDL,MAINF,NUMEL,NUMNP
     ;         ,NPAREL,NTYPAR,NZPAR,IGROUP,MXMEASPP_GS,NFLAGS
     ;         ,KXX(LMXNDL,NUMEL),LTYPE(NUMEL),LXPAREL(NUMEL,NPAREL)
     ;         ,INORPAR(NTYPAR),IGR_ZONE(NZPAR),IFLAGS(NFLAGS)
                                                                 ! Real external
      REAL*8 COORDGR_GS(6),AREA(NUMEL),COORD(NUMNP,3)
     ;      ,POSMEAS_GS(MXMEASPP_GS,3)
                                                              ! Integer internal
      INTEGER*4 IPIPO,IDUM,IPINORPAR,IPLXPAREL,NEL,IACTZONE,IACTGROUP
     ;         ,IPOS,ITRIAL
     ;         ,ISEED(13)
                                                                 ! Real internal
      REAL*8 XMAX_GR,YMAX_GR,ZMAX_GR,XMIN_GR,YMIN_GR,ZMIN_GR,XDUM,ACORNI
     ;      ,XPIPO,YPIPO,ZPIPO
     ;      ,RANDOM(3),BF(9)

C_______________________ Step 1: Initialisation of auxiliar array for random 
C_______________________         number generation. Pointers assignment

      ISEED(1)=1073741825
      DO IDUM=1,1000
        XDUM=ACORNI(ISEED)
      END DO

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

C_______________________ Step 2: Identifies maximum and minimum and maximum 
C_______________________         coordinates of actual group of zones

      XMAX_GR=COORDGR_GS(1)
      YMAX_GR=COORDGR_GS(2)
      ZMAX_GR=COORDGR_GS(3)
      XMIN_GR=COORDGR_GS(4)
      YMIN_GR=COORDGR_GS(5)
      ZMIN_GR=COORDGR_GS(6)

C_______________________ MAIN LOOP OVER PILOT POINTS

      DO IPIPO=1,NPP_GS
        ITRIAL=0                               ! Counter of trials
        IF (IORD_PP_GS.EQ.1) THEN              ! Absolutelly random drawing

C_______________________ Step 3: Generates three random factors multiplying 
C_______________________         group sizes

 10     ITRIAL=ITRIAL+1
        IF (ITRIAL.EQ.50) WRITE(*,2000) IPIPO,IGROUP
 2000     FORMAT(/,' NUMBER OF TRIALS FOR DEFINING PILOT POINT: ',I5,
     ;             ' IS HUGE. MAYBE THE SIZE OF GROUP: ',I5,' IS TOO'
     ;             ' SMALL. I KEEP ON TRYING THAT. CUT THE PROGRAM'
     ;             ' AND FIX PILOT POINT POSITIONS FOR THIS GROUP',/)
        RANDOM(1)=ACORNI(ISEED)
        RANDOM(2)=ACORNI(ISEED)
        RANDOM(3)=ACORNI(ISEED)

C_______________________ Step 4: Generates a random candidate

        XPIPO=XMIN_GR+RANDOM(1)*(XMAX_GR-XMIN_GR)
        YPIPO=YMIN_GR+RANDOM(2)*(YMAX_GR-YMIN_GR)
        ZPIPO=ZMIN_GR+RANDOM(3)*(ZMAX_GR-ZMIN_GR)

C_______________________ Step 5: Checks that candidate belongs to actual group

        CALL BASISFUNC_OBS
     ;(LMXNDL     ,MAINF      ,NEL      ,NUMEL    ,NUMNP    ,XPIPO
     ;,YPIPO      ,ZPIPO      ,AREA     ,BF       ,KXX      ,LTYPE
     ;,COORD(1,1) ,COORD(1,2) ,COORD(1,3))       

                                 ! P. point does not belong to the domain. Retry
        IF (NEL.LE.0 .OR. NEL.GT.NUMEL) GOTO 10
                                        ! Finds out ZONE to which NEL belongs to
        IACTZONE=LXPAREL(NEL,IPLXPAREL)
        IF (IACTZONE.EQ.0) GOTO 10
                                      ! Finds out GROUP to which ZONE belongs to
        IPOS=INORPAR(IPINORPAR)                  ! Initial position at IGR_ZONE
        IACTGROUP=IGR_ZONE(IPOS+IACTZONE)                        ! Actual group
                                  ! P. point does not belong to the GROUP. Retry
        IF (IACTGROUP.NE.IGROUP) GOTO 10

C_____________________ Step 6: Assigns pilot point coordinates

        POSMEAS_GS(IPIPO,1)=XPIPO
        POSMEAS_GS(IPIPO,2)=YPIPO
        POSMEAS_GS(IPIPO,3)=ZPIPO

        ELSE                                   ! Regular mesh
           ! Code not yet available
        END IF ! IORD_PP_GS.EQ.1

      END DO ! IPIPO=1,NPP_GS

C_____________________ Step 7: Echoes pilot point coordinates

      IF (IFLAGS(21).NE.0) THEN

        WRITE(MAINF,2100)
 2100   FORMAT(/,5X,'PILOT POINT RANDOM COORDINATES',/
     ;          ,5X,'===== ===== ====== ===========',/)
        DO IPIPO=1,NPP_GS
          WRITE(MAINF,2200) IPIPO,(POSMEAS_GS(IPIPO,IDUM),IDUM=1,3)
 2200     FORMAT(5X,I5,3(5X,F10.3))
        END DO 
      END IF

      RETURN
      END
