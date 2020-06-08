       SUBROUTINE MOD_RAPSON_CRITERIA
     ;(  DABSMX      ,DABSMX1       ,DABSMX2       ,DRELMX
     ;  ,DRELMX1     ,DRELMX2       ,OBJ           ,OBJ1
     ;  ,OBJ2        ,RESIDMX       ,RESIDMX1      ,RESIDMX2 )
********************************************************************************
***    MODIFIES THE NEWTON-RAPSON'S CONVERGENCE CRITERIA ACCORDING TO THE    ***
*****                ACTUAL VALUE OF THE OBJECTIVE FUNCTION                *****
********************************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)

C_______________________Computes the three differents criteria as a funtion of OBJ

       CALL INTERP_CRIT (DRELMX,DRELMX1,DRELMX2,OBJ,OBJ1,OBJ2)
       CALL INTERP_CRIT (DABSMX,DABSMX1,DABSMX2,OBJ,OBJ1,OBJ2)
       CALL INTERP_CRIT (RESIDMX,RESIDMX1,RESIDMX2,OBJ,OBJ1,OBJ2)

       RETURN
       END

