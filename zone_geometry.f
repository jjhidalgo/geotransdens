      SUBROUTINE ZONE_GEOMETRY
     ;(IACTTYPE  ,IGROUP    ,INPWR     ,IODIM     ,LMXNDL    ,MAINF     
     ;,MAXSBX    ,MAXSBY    ,MAXSBZ    ,MXDISC_GS ,MXNZON_GS ,NDISCX    
     ;,NDISCY    ,NDISCZ    ,NPAREL    ,NTYPAR    ,NUMEL
     ;,NUMNP     ,NZON_GS   ,NZPAR     ,AREA      ,COGSB_GS  ,COORD     
     ;,COORDGR_GS,DIVZN_GS  ,IGR_ZONE  ,INORPAR   ,KXX
     ;,LNNDEL    ,LXPAREL   ,NZONE_PAR ,POSZN_GS  ,POSDIS_GS ,SIZSB_GS) 

********************************************************************************
*
* PURPOSE Calculates the position of the center of gravity of the zones of a 
*         particular group of zones and the offsets wrt zone cog of their discr.
*         points. Also, sets up the coarse super block search grid.
*
* DESCRIPTION Summary:
*
*  -  Step 0: Declaration of variables
*  -  Step 1: Init. max and min coordinates of group of zones
*  -  Step 2: Loop over zones of actual parameter type, working only with those
*             belonging to actual group of zones
*     -  Step 2.1: Init. max and min coord. of current zones
*     -  Step 2.2: Loop over mesh elements. Only elems. belonging to current 
*                 zone are considered
*        -  Step 2.2.1: Init. cog. of current element
*        -  Step 2.2.2: Loop over element nodes
*           -  Step 2.2.2a: Identifies node and coordinates
*           -  Step 2.2.2b: Update max and min of zone and group coordinates
*           -  Step 2.2.2c: Updates element cog
*        -  Step 2.2.3: Calculates element cog
*        -  Step 2.2.4: Adds contrib. of element cog. to zone cog.
*     -  Step 2.3: Calculates zone cog
*     -  Step 2.4: Calculates offsets of discretization points wrt to zone cog. 
*                  Discretization points are set fast cycling on Z, then on Y, 
*                  then on X. and equally distributed. There is a chance to set
*                  disc. points on the positions of the cog of elements 
*                  belonging to this zone. Only partition would be complicated
*  -  Step 3: Checks number of zones
*  -  Step 4: Sets up superblocks search grid (dimensions and cog. of first 
*             superblock)
*  -  Step 5: Stores maximum and minimum coordinates of group of zones. They 
*             will be used again if pilot points are drawn randomly
*  -  Step 6: Echoes position of zones cog
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  COGSB_GS               Center of gravity of first superblock of search grid
*  COORD                  Nodal coordinates                                     
*  COORDGR_GS             Maximum and minimum coodinates of all groups of zones
*  DIVZN_GS               Auxiliar array containing the sum of areas of elements
*                         drawing a given zone
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  IGR_ZONE               Group of zones to which a given zone belongs to. Part
*                         of array IVPAR on input
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  KXX                    Node numbers of every element (counterclockwise order)
*  LNNDEL                 Number of nodes at every element                      
*  LXPAREL                Zone number to which a given element belongs to        
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  POSZN_GS               Array containing the position of the center of gravity
*                         of the formation blocks                                   
*  POSDIS_GS              Array containing the position of block discretisation 
*                         points                                           
*  SIZSB_GS               Size of superblocks defining search grid
*
* INTERNAL VARIABLES: ARRAYS
*
*  COGELEM                Array containing the cent. of grav. of a given element
*
* EXTERNAL VARIABLES: SCALARS
*
*  IACTYTPE               Group of zones parameter type
*  IGROUP                 Actual group of zones
*  INPWR                  Allows writing on MAIN FILE                           
*  IODIM                  Maximum dimension of the problem                      
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MAXSBX                 Number of super blocks (X-direction)     
*  MAXSBY                 Number of super blocks (Y-direction)                  
*  MAXSBZ                 Number of super blocks (Z-direction)
*  MXDISC_GS              Number of discretisation points for a given block
*  MXNZN_GS               Number of zones of the most discretized formation
*  NBL_GS                 Number of zones of current formation
*  NDISCX                 Number of discretization points X-direction
*  NDISCY                 Number of discretization points Y-direction
*  NDISCZ                 Number of discretization points Z-direction
*  NFLAGS                 Used to dimension IFLAGS
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
*  DELTAX,DELTAY,DELTAZ   Maximum lenghts of current block (X-Y-Z)
*  I                      Dummy counter
*  ICONT                  Counter of zones of this group for check purposes
*  IZON                   Counter of zones
*  IDISC                  Counter of discretization points
*  INUD                   Counter of nodal points of current element L
*  IPINORPAR              Position at INORPAR of actual parameter type
*  IPLXPAREL              Position at LXPAREL of actual parameter type
*  IPNZPAR                Position at NZONE_PAR of actual parameter type
*  IX,IY,IZ               Counters of discretization points (X-Y-Z)
*  L                      Counter of elements
*  LX                     Zone of parameters to which element L belongs to
*  NODE                   Actual nodal point identifier
*  XDISC,YDISC,ZDISC      Offsets (updatable) of discr. points
*  XLOC,YLOC,ZLOC         Offsets of discr. points wrt zone cog
* XMAX_ZN,YMAX_ZN,ZMAX_ZN Maximum coordinates of a zone
* XMAX_GR,YMAX_GR,ZMAX_GR Maximum coordinates of actual group
* XMIN_BL,YMIN_BL,ZMIN_BL Minimum coordinates zone of a zone
* XMIN_GR,YMIN_GR,ZMIN_GR Minimum coordinates of actual group
*  XNODE,YNODE,ZNODE      Nodal point coordinates
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ZERO_ARRAY                                                                   
*
* HISTORY:  AAR  First coding (Dec-2001)
*           AAR  Revision and header (April-2002)
*           AAR  Revision and inclusion of groups of zones (July-2003)
*
********************************************************************************

C______________________ Step 0: Declaration/init. of variables

      IMPLICIT NONE 
                                                              ! Integer external
      INTEGER*4 NZON_GS,MXNZON_GS,NUMEL,LMXNDL,NUMNP,IODIM,NDISCX,NDISCY
     ;         ,NDISCZ,MXDISC_GS,MAXSBX,MAXSBY,MAXSBZ,NTYPAR,IACTTYPE 
     ;         ,NZPAR,IGROUP,MAINF,NPAREL,INPWR
     ;         ,NZONE_PAR(NTYPAR),LXPAREL(NUMEL,NPAREL),LNNDEL(NUMEL)
     ;         ,KXX(LMXNDL,NUMEL),INORPAR(NTYPAR),IGR_ZONE(NZPAR)
                                                                 ! Real external
      REAL*8 COORD(NUMNP,3),POSZN_GS(MXNZON_GS,3),DIVZN_GS(MXNZON_GS)
     ;      ,AREA(NUMEL),POSDIS_GS(MXDISC_GS,3,MXNZON_GS),SIZSB_GS(3)
     ;      ,COGSB_GS(3),COORDGR_GS(6)
                                                              ! Integer internal
      INTEGER*4 IZON,LX,L,INUD,NODE,IX,IY,IZ,IDISC,ICONT,IPOS,I,IPNZPAR
     ;         ,IPINORPAR,IPLXPAREL
                                                                 ! Real internal
      REAL*8 XMAX_GR,XMIN_GR,YMAX_GR,YMIN_GR,ZMAX_GR,ZMIN_GR
     ;      ,XMAX_ZN,XMIN_ZN,YMAX_ZN,YMIN_ZN,ZMAX_ZN,ZMIN_ZN
     ;      ,XNODE,YNODE,ZNODE,DELTAX,DELTAY,DELTAZ,XDISC,YDISC,ZDISC
     ;      ,XLOC,YLOC,ZLOC
     ;      ,COGELEM(3)

      CALL ZERO_ARRAY(DIVZN_GS,MXNZON_GS)

C______________________ Step 1: Init. max and min coordinates of formation

      XMAX_GR=-1.0D50
      XMIN_GR=1.0D50
      YMAX_GR=-1.0D50
      YMIN_GR=1.0D50
      ZMAX_GR=-1.0D50
      ZMIN_GR=1.0D50

C______________________ Step 2: Loop over zones of actual parameter type, 
C______________________         working only with those belonging to actual
C______________________         group of zones

      ICONT=0         ! Counter (internal) of zones belonging to actual group
      IPINORPAR=1     ! Pointer to INORPAR
      IPNZPAR=1       ! Pointer to NZONE_PAR
      IPLXPAREL=1     ! Pointer to LXPAREL
      IF (IACTTYPE.EQ.2) THEN  ! Storage coeff.
         IPINORPAR=7
         IPNZPAR=IACTTYPE
         IPLXPAREL=IACTTYPE
      ELSE IF (IACTTYPE.EQ.3 .OR. IACTTYPE.EQ.4) THEN  ! Areal recharge
         IPINORPAR=8
         IPNZPAR=3
         IPLXPAREL=IACTTYPE
      ELSE IF (IACTTYPE.GT.4) THEN ! Others
         IPINORPAR=IACTTYPE+7
         IPNZPAR=IACTTYPE+2
         IF (IACTTYPE.EQ.5 .OR. IACTTYPE.EQ.6) THEN
            IPLXPAREL=5
         ELSE
            IPLXPAREL=IACTTYPE-1
         END IF
      END IF

      DO IZON=1,NZONE_PAR(IPNZPAR)
        IPOS=INORPAR(IPINORPAR)                  ! Initial position at IGR_ZONE
        IF (IGR_ZONE(IPOS+IZON).EQ.IGROUP) THEN  ! Zone belongs to group

          ICONT=ICONT+1
          IF (ICONT.GT.NZON_GS) THEN
        WRITE(MAINF,1000) IGROUP,NZON_GS,ICONT
        WRITE(*,1000) IGROUP,NZON_GS,ICONT
        STOP
      END IF

C______________________ Step 2.1: Init. max and min coord. of current zone

          XMAX_ZN=-1.0D50
          XMIN_ZN=1.0D50
          YMAX_ZN=-1.0D50
          YMIN_ZN=1.0D50
          ZMAX_ZN=-1.0D50
          ZMIN_ZN=1.0D50

C______________________ Step 2.2: Loop over mesh elements. Only elems belonging
C______________________           to current zone are considered

          DO L=1,NUMEL
            LX=LXPAREL(L,IPLXPAREL)

            IF (LX.EQ.IZON) THEN

C______________________ Step 2.2.1: Init. cog. of current element

              CALL ZERO_ARRAY(COGELEM,3)

C______________________ Step 2.2.2: Loop over element nodes

              DO INUD=1,LNNDEL(L)
              
C______________________ Step 2.2.2a: Identifies node and coordinates

                 NODE=KXX(INUD,L)
                 XNODE=COORD(NODE,1)
                 YNODE=COORD(NODE,2)
                 ZNODE=COORD(NODE,3)

C______________________ Step 2.2.2b: Update max and min coord. of zone and group

                 IF (XNODE.GE.XMAX_ZN) XMAX_ZN=XNODE
                 IF (YNODE.GE.YMAX_ZN) YMAX_ZN=YNODE
                 IF (ZNODE.GE.ZMAX_ZN) ZMAX_ZN=ZNODE
                 IF (XNODE.LE.XMIN_ZN) XMIN_ZN=XNODE
                 IF (YNODE.LE.YMIN_ZN) YMIN_ZN=YNODE
                 IF (ZNODE.LE.ZMIN_ZN) ZMIN_ZN=ZNODE

                 IF (XMAX_ZN.GE.XMAX_GR) XMAX_GR=XMAX_ZN
                 IF (YMAX_ZN.GE.YMAX_GR) YMAX_GR=YMAX_ZN
                 IF (ZMAX_ZN.GE.ZMAX_GR) ZMAX_GR=ZMAX_ZN
                 IF (XMIN_ZN.LE.XMIN_GR) XMIN_GR=XMIN_ZN
                 IF (YMIN_ZN.LE.YMIN_GR) YMIN_GR=YMIN_ZN
                 IF (ZMIN_ZN.LE.ZMIN_GR) ZMIN_GR=ZMIN_ZN

C______________________ Step 2.2.2c: Updates element cog

                 COGELEM(1)=COGELEM(1)+XNODE
                 COGELEM(2)=COGELEM(2)+YNODE
                 COGELEM(3)=COGELEM(3)+ZNODE

              END DO ! INUD=1,NNUD

C______________________ Step 2.2.3: Calculates element cog

              COGELEM(1)=COGELEM(1)/DFLOAT(LNNDEL(L))
              COGELEM(2)=COGELEM(2)/DFLOAT(LNNDEL(L))
              COGELEM(3)=COGELEM(3)/DFLOAT(LNNDEL(L))

C______________________ Step 2.2.4: Adds contrib. of element cog. to block cog.
  
              POSZN_GS(ICONT,1)=POSZN_GS(ICONT,1)+COGELEM(1)*AREA(L)
              POSZN_GS(ICONT,2)=POSZN_GS(ICONT,2)+COGELEM(2)*AREA(L)
              POSZN_GS(ICONT,3)=POSZN_GS(ICONT,3)+COGELEM(3)*AREA(L)
              DIVZN_GS(ICONT)=DIVZN_GS(ICONT)+AREA(L)

            END IF ! LX.EQ.IZON

          END DO  ! L=1,NUMEL

C______________________ Step 2.3: Calculates zone cog

          POSZN_GS(ICONT,1)=POSZN_GS(ICONT,1)/DIVZN_GS(ICONT)
          POSZN_GS(ICONT,2)=POSZN_GS(ICONT,2)/DIVZN_GS(ICONT)
          POSZN_GS(ICONT,3)=POSZN_GS(ICONT,3)/DIVZN_GS(ICONT)

C______________________ Step 2.4: Calculates offsets of discretization points
C______________________           wrt to block cog. Discretization are set fast
C______________________           cycling on Z, then on Y, then on X. and equally
C______________________           distributed. There is a chance to set disc. 
C______________________           points on the positions of the cog of elements
C______________________           belonging to this block. Only partition would
C______________________           be complicated

          DELTAX=XMAX_ZN-XMIN_ZN
          DELTAY=YMAX_ZN-YMIN_ZN
          DELTAZ=ZMAX_ZN-ZMIN_ZN
          IF (IODIM.EQ.2) DELTAZ=1.D0

          XDISC=DELTAX/DMAX1(DFLOAT(NDISCX),1.D0)
          YDISC=DELTAY/DMAX1(DFLOAT(NDISCY),1.D0)
          ZDISC=DELTAZ/DMAX1(DFLOAT(NDISCZ),1.D0)

          XLOC=-1.D0*(DELTAX+XDISC)/2.D0        ! Initial X position
          IDISC=0                               ! Counter of discr. points
          DO IX=1,NDISCX
            XLOC=XLOC+XDISC
            YLOC=-1.D0*(DELTAY+YDISC)/2.D0      ! Initial Y position
            DO IY=1,NDISCY
              YLOC=YLOC+YDISC
              ZLOC=-1.D0*(DELTAZ+ZDISC)/2.D0
              DO IZ=1,NDISCZ
                ZLOC=ZLOC+ZDISC
                IDISC=IDISC+1
                POSDIS_GS(IDISC,1,ICONT)=XLOC
                POSDIS_GS(IDISC,2,ICONT)=YLOC
                POSDIS_GS(IDISC,3,ICONT)=ZLOC
              END DO ! IZ=1,NDISCZ
            END DO ! IY=1,NDISCY
          END DO ! IX=1,NDISCX

        END IF ! IGR_ZONE(IPOS+IZON).EQ.IGROUP
      END DO ! IZON=1,NZONE_PAR(IACTTYPE)

C______________________ Step 3: Checks number of zones

      IF (ICONT.LT.NZON_GS) THEN
        WRITE(MAINF,1000) IGROUP,NZON_GS,ICONT
        WRITE(*,1000) IGROUP,NZON_GS,ICONT
 1000   FORMAT(//,' ERROR CALCULATING GEOMETRY OF GROUP:',I5,/,
     ;            ' YOU DEFINED IT USING ',I5,' ZONES BUT I '
     ;            ' FOUND ',I5,' FORCED STOP. PLEASE, CHECK IT')
        STOP
      END IF

C______________________ Step 4: Sets up superblocks search grid (dimensions and
C______________________         cog. of first superblock)

      SIZSB_GS(1)=DMAX1(1.D0,(XMAX_GR-XMIN_GR)/DFLOAT(MAXSBX))
      SIZSB_GS(2)=DMAX1(1.D0,(YMAX_GR-YMIN_GR)/DFLOAT(MAXSBY))
      SIZSB_GS(3)=DMAX1(1.D0,(ZMAX_GR-ZMIN_GR)/DFLOAT(MAXSBZ))

      COGSB_GS(1)=XMIN_GR+0.5D0*SIZSB_GS(1)
      COGSB_GS(2)=YMIN_GR+0.5D0*SIZSB_GS(2)
      COGSB_GS(3)=ZMIN_GR+0.5D0*SIZSB_GS(3)
      IF (IODIM.EQ.2) COGSB_GS(3)=0.D0

C______________________ Step 5: Stores maximum and minimum coordinates of group
C______________________         of zones. They will be used again if pilot points
C______________________         are drawn randomly

      COORDGR_GS(1)=XMAX_GR
      COORDGR_GS(2)=YMAX_GR
      COORDGR_GS(3)=ZMAX_GR
      COORDGR_GS(4)=XMIN_GR
      COORDGR_GS(5)=YMIN_GR
      COORDGR_GS(6)=ZMIN_GR

C______________________ Step 6: Echoes positions of zones cog

      IF (INPWR.EQ.2) THEN

        WRITE(MAINF,2000) IGROUP
 2000   FORMAT(5X,'ZONES COG BELONGING TO GROUP: ',I5,/
     ;        ,5X,'===== === ========= == ======',//)
        DO I=1,NZON_GS
          WRITE(MAINF,2100) I,POSZN_GS(I,1),POSZN_GS(I,2),POSZN_GS(I,3)
 2100     FORMAT(I5,3(2X,F10.3))
        END DO
      END IF

      RETURN
      END
