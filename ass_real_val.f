       SUBROUTINE ASS_REAL_VAL(V_READ,V_DEF,V_ASS)
********************************************************************************
* PURPOSE 
*
*       This subroutine V_ASSigns V_DEF value to V_ASS if V_READ equals 0
*       and V_DEF non equals zero; otherwise V_ASSigns V_READ value to V_ASS 
*
* EXTERNAL VARIABLES: SCALARS 
*
*       V_DEF
*       V_ASS
*       V_READ
*
* HISTORY
*
*     SCR      11-04-1997     First coding
*
*
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       IF (V_READ.EQ.0.AND.V_DEF.NE.0) THEN
         V_ASS=V_DEF
       ELSE
         V_ASS=V_READ
       END IF
       RETURN
       END
