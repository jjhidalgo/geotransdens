       SUBROUTINE CHECK_CONVERGENCE
     &           (DABSMX   ,DELVMAX  ,DRELMX   ,DRELVMX  ,INDCONVV
     &           ,INDENDDT ,IOWNR    ,NCONVI   ,RESIDMX  ,RESVMAX
     &           ,IOVLI    ,IONEWT)

*********************************************************
*
* PURPOSE  Check convergences criterion and writes hopeful information
*
* DESCRIPTION The criterions are: - fitting relative maximum change in state 
*                                   variable 
*                                 - fitting maximum change in state variable 
*                                   controlled
*                                 - fitting maximum residuo controlled
*             Process must verify both first and second or the third one.
*
* EXTERNAL VARIABLES: SCALARS
*
*  DABSMX                 Absolute convergence criterion                        
*  DELVMAX                Maximum change in state variable                                                                      
*  DRELMX                 Relative convergence criterion                        
*  DRELVMX                Maximum relative change in state variable                                                      
*  INDCONVV               Indicator variable (=1 --> convergence reached)                                                      
*  INDENDDT               Controls the coincidence of the next computing        
*                         time with the end of the current time observation     
*                         interval                                              
*  IOWNR                  Option of printing information about the direct       
*                         problem iterative process evolution                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NCONVI                 Number of convergences reached during the direct 
*                         problem resolution
*  RESIDMX                Maximum residuo permitted by user                                                                      
*  RESVMAX                Maximum resiuo obtained during the solution                                                      
*  IOVLI                  Indicates wether the problem considered is linear
*
* INTERNAL VARIABLES: SCALARS
*
*  EQUATION               String containing the tyoe of equation solved.
*
*
* HISTORY: First coding : German Galarza
*          Revision: Andres Alcolea (also in August-98)
*
********************************************************************************

        IMPLICIT NONE

      INTEGER*4::INDCONVV ,INDENDDT ,IOWNR    ,NCONVI   ,IOVLI
     &          ,IONEWT

      REAL*8::DRELVMX,DRELMX,DELVMAX,DABSMX,RESVMAX,RESIDMX


C-------------------- Convergence reached.

      IF ((DRELVMX.LT.DRELMX .AND. DELVMAX.LT.DABSMX)
     &    .OR. RESVMAX.LT.RESIDMX) THEN

          INDCONVV=1

          NCONVI=NCONVI+1

          IF (IOWNR.NE.0) THEN

              IF (MOD(NCONVI,IOWNR).EQ.0 .OR. INDENDDT.EQ.1) THEN
c-parche
c                  WRITE(MAINF,1000) NCONVI

c1000              FORMAT
c    &        (/,' TOTAL NUMBER OF CONV. (PER ITERATION OF I.P.)=',I5,/)
c-finfparce
              END IF !MOD(NCONVI,IOWNR).EQ.0 .OR. INDENDDT.EQ.1

          END IF !IOWNR.NE.0

      ELSE

C-------------------- Convergence NOT reached.

          INDCONVV=0

      ENDIF

C-------------------- If the problem is linear,
C-------------------- then it has converged if it was solved with Picard.
C-------------------- If it was solved with Newton, the maximum update
C-------------------- might have prevented the update to reach the solution
      IF (IOVLI .EQ.0 .AND. IONEWT.EQ.0) INDCONVV = 1

      RETURN

      END SUBROUTINE CHECK_CONVERGENCE
