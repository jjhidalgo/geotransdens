      SUBROUTINE ROTMAT_UNIVKRIG_COVMAX
     ;(COVMAX        ,KRIGTYPE      ,MAXROT        ,NESTED    
     ;,NRESTRIUNIV   ,NUGGET        ,RADIUS        ,SCALE_UNIV 
     ;,SEARCH_ANGLE1 ,SEARCH_ANGLE2 ,SEARCH_ANGLE3 ,SEARCH_ANIS1
     ;,SEARCH_ANIS2  ,ANGLE1        ,ANGLE2        ,ANGLE3
     ;,ANIS1         ,ANIS2         ,IDRIF         ,ITYPEVARIO
     ;,ROTMAT        ,SILL_NUGG)

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NESTED,MAXROT,KRIGTYPE,NRESTRIUNIV
     ;         ,ITYPEVARIO(NESTED),IDRIF(9)
                                                                 ! Real external
      REAL*8 COVMAX,NUGGET,SEARCH_ANGLE1,SEARCH_ANGLE2,SEARCH_ANGLE3
     ;      ,SEARCH_ANIS1,SEARCH_ANIS2,SCALE_UNIV,RADIUS
     ;      ,ROTMAT(MAXROT,3,3),ANGLE1(NESTED),ANGLE2(NESTED)
     ;      ,ANGLE3(NESTED),ANIS1(NESTED),ANIS2(NESTED)
     ;      ,SILL_NUGG(NESTED)
                                                              ! Integer internal
      INTEGER*4 ISTRU,ICOMPO
                                                                 ! Real internal
      REAL*8 MONOMIAL

C_______________________ Step 1: Calculation of nested structures of variogram

      MONOMIAL= 1.0D20     ! Parameter used as sill of monomial model variograms
      COVMAX = NUGGET      ! Initialization of maximum covariance: Nugget effect

      DO ISTRU=1,NESTED                    ! Loop over the nst nested structures

C__________________________ Step 1.1: Sets rotation matrix for current structure

        CALL SETROT
     ;(ANGLE1(ISTRU)   ,ANGLE2(ISTRU)   ,ANGLE3(ISTRU)   ,ANIS1(ISTRU)
     ;,ANIS2(ISTRU)    ,ISTRU           ,MAXROT          ,ROTMAT)

C___________________________ Step 1.2: Adds contribution of current structure to
C___________________________           maximum covariance

        IF (ITYPEVARIO(ISTRU).EQ.4) THEN              ! Monomial model variogram
           COVMAX= COVMAX + MONOMIAL                   ! COVMAX MUST BE INFINITE
         ELSE                                     ! Others: Adds SILL(IS)-NUGGET
           COVMAX = COVMAX + SILL_NUGG(ISTRU)      
         END IF

      END DO ! ISTRU=1,NESTED                                                   

C_______________________ Step 2: Calculates rotation matrix of search ellipsoid

      CALL SETROT
     ;(SEARCH_ANGLE1   ,SEARCH_ANGLE2   ,SEARCH_ANGLE3   ,SEARCH_ANIS1
     ;,SEARCH_ANIS2    ,NESTED+1        ,MAXROT          ,ROTMAT)

C_______________________ Step 3: Checks if universal kriging is done

      IF (KRIGTYPE.EQ.1 .OR. KRIGTYPE.EQ.3) THEN
        NRESTRIUNIV=0
        DO ICOMPO=1,9
          IF (IDRIF(ICOMPO).NE.0) NRESTRIUNIV=NRESTRIUNIV+1
        END DO ! IDRIF=1,9
      END IF ! KRIGTYPE.EQ.1 .OR. KRIGTYPE.EQ.3
      
C_______________________ Step 4: If universal kriging is done, calculates 
C_______________________         scale factor of universal kriging components

      IF (NRESTRIUNIV.NE.0) THEN

        IF (RADIUS*RADIUS.LT.1.0D0) THEN
           SCALE_UNIV = 2.0D0 * RADIUS / DMAX1(COVMAX,1.0D-4)
        ELSE
           SCALE_UNIV =(4.0D0 * RADIUS * RADIUS)/ DMAX1(COVMAX,1.0D-4)
        END IF

        SCALE_UNIV= 1.0D0 / SCALE_UNIV

      END IF ! NRESTRIUNIV.NE.0

      RETURN
      END
