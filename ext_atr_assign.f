      SUBROUTINE EXT_ATR_ASSIGN
     ;(COVMAX    ,IESTIM    ,MAINF    ,NEXDRIFTS    ,SCALE_EXDRIFT
     ;,SECEST1   ,SECEST2   ,SECEST3  ,SECEST4      ,EXTDRESTIM
     ;,TRIM_LIMIT)

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NEXDRIFTS,IESTIM,MAINF
                                                                 ! Real external
      REAL*8 SCALE_EXDRIFT,SECEST1,SECEST2,SECEST3,SECEST4,COVMAX
     ;      ,EXTDRESTIM(4),TRIM_LIMIT(4)
                                                              ! Integer external
      INTEGER*4 IEXT

C_______________________ Step 1: Loop over external drift terms, checking if
C_______________________         some of them are trimmed. In this case, echoes
C_______________________         a critical error and stops. Calculates scale 
C_______________________         factor of external drift terms

      DO IEXT=1,NEXDRIFTS

         IF (IEXT.EQ.1) EXTDRESTIM(IEXT) = SECEST1
         IF (IEXT.EQ.2) EXTDRESTIM(IEXT) = SECEST2
         IF (IEXT.EQ.3) EXTDRESTIM(IEXT) = SECEST3
         IF (IEXT.EQ.4) EXTDRESTIM(IEXT) = SECEST4

         IF (EXTDRESTIM(IEXT).GE.TRIM_LIMIT(IEXT)) THEN
            WRITE(6,1000) IEXT,IESTIM
            WRITE(MAINF,1000) IEXT,IESTIM
 1000       FORMAT(/,' ERROR: EXTERNAL DRIFT TERM NUMBER: ',I5,
     ;               ' IS TRIMMED AT ESTIMATION POINT NUMBER: ',I5,
     ;               ' CRITICAL STOP')
         END IF ! EXTDRESTIM(IEXT).GE.TRIM_LIMIT(IEXT)

         SCALE_EXDRIFT = DMAX1 ( COVMAX / DMAX1(EXTDRESTIM(IEXT),1D-4)
     ;                          ,SCALE_EXDRIFT)

      END DO ! IEXT=1,NEXDRIFTS

      RETURN
      END
