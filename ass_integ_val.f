       SUBROUTINE ASS_INTEG_VAL(I_READ,I_DEF,I_ASS)
********************************************************************************
* PURPOSE 
*
*       This subroutine I_ASSigns I_DEF value to I_ASS if I_READ equals 0
*       and I_DEF non equals zero; otherwise I_ASSigns I_READ value to I_ASS 
*
* EXTERNAL VARIABLES: SCALARS 
*
*       I_DEF
*       I_ASS
*       I_READ
*
* HISTORY
*
*     SCR      11-04-1997     First coding
*
*
********************************************************************************

       IF (I_READ.EQ.0.AND.I_DEF.NE.0) THEN
         I_ASS=I_DEF
       ELSE
         I_ASS=I_READ
       END IF
       RETURN
       END
