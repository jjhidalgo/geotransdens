        INTEGER*4 FUNCTION LDIMEN (ITIP)

*****************************************************************************
* PURPOSE
*     Computes the dimension of a element from its type
*
* DESCRIPTION
*     Computes the dimension of a element from its type
*
* EXTERNAL VARIABLES: SCALARS
*
*     ITIP                Type of the current element
*
*
* HISTORY
*
*     SCR      5-1997     First coding
*     AMS      4-1988     Full modification
*****************************************************************************

        IF (ITIP.EQ.1 .OR. ITIP.EQ.2) THEN
           LDIMEN=1
        ELSE IF ( (ITIP.GE.3 .AND. ITIP.LE.8) .OR. ITIP.EQ.10) THEN
           LDIMEN=2
        ELSE IF (ITIP.EQ.9 .OR. ITIP.EQ.11) THEN
           LDIMEN=3
        ENDIF

        RETURN
        END
