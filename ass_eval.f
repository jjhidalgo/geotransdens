       SUBROUTINE ASS_EVAL(IOPAR ,LDPAR ,LPAR  ,LXPAR)
       
********************************************************************************
* PURPOSE 
*
*     This subroutine assigns an specific element parameter zone number  
*     [LXPAREL(INPAR,NE), the last variable in argument]
* 
* DESCRIPTION
*
*     If IOPAR non equals zero:  LXPAR is assigned defaulat value LDPAR
*     value if this last is not zero, otherwise LXPAR is assigned LPAR value. 
*
*
* EXTERNAL VARIABLES: SCALARS
*
*     IOPAR               If eq. zero, current parameter is not considered
*     LDPAR               Default vaule
*     LPAR                Read or alternative value
*     LXPAR               Element parameter zone number to define.
*                     
* HISTORY
*
*     SCR 15-abr-1997     Firt coding
*
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       
       IF (IOPAR.NE.0) THEN
       
          IF (LDPAR.NE.0) THEN
              LXPAR=LDPAR
          ELSE
              LXPAR=LPAR
          END IF
       END IF
       
       
       RETURN
       END 
       
