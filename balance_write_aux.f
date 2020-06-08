      SUBROUTINE BALANCE_WRITE_AUX 
     ;(BM,CZONE,NZ,LINE,FORMATINI,FORMATEND,TOTAL)

*******************************************************************************
*
* PURPOSE Writes flow or transport mass balance information for a paticular 
*         zonal/nodal process.
*
* DESCRIPTION It is called from BALANCE_WRITE, once per zone (zonal mass balance)
*             or once per nodal point (nodal mass balance)
*
* EXTERNAL VARIABLES: ARRAYS
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  BM                     Array containing mass balance information
*  CZONE                  Current zone of process type
*  END                    Final position of writing at LINE
*  INI                    Initial position of writing at LINE
*  NZ                     Number of zones of current process
*  TOTAL                  Mass balance contribution of current process
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY: AAR: First coding (Dec-2000)
*
********************************************************************************

      IMPLICIT NONE
      INTEGER*4 NZ,FORMATEND,FORMATINI,CZONE
      REAL*8 BM,TOTAL
      CHARACTER*156 LINE

C______________________________ When zone is not defined, writes ------

         IF (CZONE.GT.NZ) THEN
            LINE(FORMATINI:FORMATEND)='    ----    '
         ELSE
            TOTAL=TOTAL+BM
            WRITE(LINE(FORMATINI:FORMATEND),1000) BM
         END IF

 1000 FORMAT(G12.6)
      RETURN
      END
