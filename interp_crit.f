       SUBROUTINE INTERP_CRIT (CRIT,CRIT1,CRIT2,OBJ,OBJ1,OBJ2)
********************************************************************************
***                       INTERPOLES LINEARY                                 ***
********************************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)

       IF (OBJ.LE.OBJ1)THEN
         CRIT=CRIT1
       ELSE IF(OBJ.GE.OBJ2)THEN
         CRIT=CRIT2
       ELSE
         SLOPE=(CRIT2-CRIT1)/(OBJ2-OBJ1)
         RHS=CRIT1-SLOPE*OBJ1
         CRIT=SLOPE*OBJ+RHS
       ENDIF
       RETURN
       END
