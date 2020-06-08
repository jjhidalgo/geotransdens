      SUBROUTINE CALC_MEAN_POLDRIFT
     ;(NDISC    ,SCALE_UNIV    ,IDRIF    ,MEAN_DRIFT    ,XOFF
     ;,YOFF     ,ZOFF)

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NDISC
     ;         ,IDRIF(9)
                                                                 ! Real external
      REAL*8 SCALE_UNIV
     ;      ,MEAN_DRIFT(9),XOFF(NDISC),YOFF(NDISC),ZOFF(NDISC)
                                                              ! Integer external
      INTEGER*4 IDISC,ICOMPO

C_______________________ Step 1: Initialization of mean values

      CALL ZERO_ARRAY (MEAN_DRIFT,9)

C_______________________ Step 2: Calculation of mean values over discretization
C_______________________         points

      DO IDISC=1,NDISC
                                                                 ! Linear terms
        IF (IDRIF(1).NE.0) MEAN_DRIFT(1) = MEAN_DRIFT(1) + XOFF(IDISC)
        IF (IDRIF(2).NE.0) MEAN_DRIFT(2) = MEAN_DRIFT(2) + YOFF(IDISC)
        IF (IDRIF(3).NE.0) MEAN_DRIFT(3) = MEAN_DRIFT(3) + ZOFF(IDISC)
                                                              ! Quadratic terms
        IF (IDRIF(4).NE.0) 
     ;     MEAN_DRIFT(4) = MEAN_DRIFT(4) + XOFF(IDISC) * XOFF(IDISC)
        IF (IDRIF(5).NE.0) 
     ;     MEAN_DRIFT(5) = MEAN_DRIFT(5) + YOFF(IDISC) * YOFF(IDISC)
        IF (IDRIF(6).NE.0) 
     ;     MEAN_DRIFT(6) = MEAN_DRIFT(6) + ZOFF(IDISC) * ZOFF(IDISC)
                                                        ! Cross-Quadratic terms
        IF (IDRIF(7).NE.0) 
     ;     MEAN_DRIFT(7) = MEAN_DRIFT(7) + XOFF(IDISC) * YOFF(IDISC)
        IF (IDRIF(8).NE.0) 
     ;     MEAN_DRIFT(8) = MEAN_DRIFT(8) + XOFF(IDISC) * ZOFF(IDISC)
        IF (IDRIF(9).NE.0) 
     ;     MEAN_DRIFT(9) = MEAN_DRIFT(9) + YOFF(IDISC) * ZOFF(IDISC)
      END DO ! IDISC=1,NDISC

      DO ICOMPO=1,9
        MEAN_DRIFT(9) = MEAN_DRIFT(9) * SCALE_UNIV / FLOAT(NDISC)
      END DO ! ICOMPO=1,9

      RETURN
      END
