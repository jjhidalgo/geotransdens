      SUBROUTINE WRITE_DEBUG_INFO
     ;(IESTIM  ,LAGRANGE  ,NACCEPT_SAMPLES  ,XEST     ,YEST
     ;,ZEST    ,NSAMPLE   ,VARKRIG          ,WEIGHTS  ,XKRIG
     ;,YKRIG   ,ZKRIG)

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                             ! Integer external
      INTEGER*4 IESTIM,NACCEPT_SAMPLES
     ;         ,NSAMPLE(NACCEPT_SAMPLES)
                                                                ! Real external
      REAL*8 XEST,YEST,ZEST,LAGRANGE
     ;      ,XKRIG(NACCEPT_SAMPLES),YKRIG(NACCEPT_SAMPLES)
     ;      ,ZKRIG(NACCEPT_SAMPLES),VARKRIG(NACCEPT_SAMPLES)
     ;      ,WEIGHTS(NACCEPT_SAMPLES)
                                                             ! Integer internal
      INTEGER*4 ISAMPLE
                                                                ! Real internal
      WRITE(666,1000) IESTIM,XEST,YEST,ZEST
 1000 FORMAT(//,' ESTIMATION POINT ',I5,' AT ',3E10.3,':')

      WRITE(666,1100) LAGRANGE
 1100 FORMAT(/,' LAGRANGE MULT. ....=',E10.3)

      WRITE(666,1200) 
 1200 FORMAT(/,9X,'X',9X,'Y',9X,'Z',5X,'VALUE',4X,'WEIGHT',' SAMPLE_ID'
     ;       /,9X,'=',9X,'=',9X,'=',5X,'=====',4X,'======',' =========')

      WRITE(666,'(5E10.3,I5)') ( XKRIG(ISAMPLE)+XEST
     ;                           ,YKRIG(ISAMPLE)+YEST
     ;                           ,ZKRIG(ISAMPLE)+ZEST
     ;                           ,VARKRIG(ISAMPLE)
     ;                           ,WEIGHTS(ISAMPLE)
     ;                           ,NSAMPLE(ISAMPLE) 
     ;                           ,ISAMPLE=1,NACCEPT_SAMPLES)

      RETURN
      END
