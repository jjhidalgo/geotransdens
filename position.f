      INTEGER FUNCTION POSITION (KOL,NUMTOBS,ROW)

********************************************************************************
*
* PURPOSE to calculate the array-index in the array covinv or any
* other array in the same storage mode, starting from the row and
* column position of the element and from the nr of rows
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan. 2003): Revision and formatting
*
********************************************************************************

      IMPLICIT NONE
      INTEGER*4 ROW,KOL,BAND,J,NUMTOBS


      BAND = ROW - KOL 
      POSITION=0
      IF (ABS(BAND) .GE. 1) THEN     !IF NONDIAGONAL ELEMENT
        DO J=0,BAND-1     
          POSITION=POSITION+NUMTOBS-J
        ENDDO
      ENDIF                  !ENDIF NONDIAGONAL ELEMENT
      POSITION=POSITION+KOL

      RETURN
      END 
