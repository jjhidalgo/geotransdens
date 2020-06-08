       SUBROUTINE WRITE_NONLIN_INFO
     ;(DELVMAX  ,EQUATION ,IDELVMAX ,INDENDDT ,IOREDTIMV,IOWNR
     ;,IRESVMAX ,ITRAP    ,MAINF    ,NCONVI   ,RESVMAX  ,TABSOLUT ,TINC)

*********************************************************
*
* PURPOSE Writes some useful information about the Newton-Raphson process
*         evolution, for both flow and transport equation.
*
* DESCRIPTION Writes the solution time and the current time increment. 
*             Next, it writes the iteration number, the maximum residuo
*             and the corresponding node, the maximum change in the state 
*             variable and the corresponding node. 
*
* EXTERNAL VARIABLES: SCALARS
*
*  DELVMAX                Maximum change in the state variable (h or c)                                                                      
*  EQUATION               String containing the equation name 
*                         (FLOW or TRANSPORT)                                                      
*  IDELVMAX               Node associated with the max. chane in the state 
*                         variable                                                       
*  INDENDDT               Controls the coincidence of the next computing        
*                         time with the end of the current time observation     
*                         interval                                              
*  IOWNR                  Option of printing information about the direct       
*                         problem iterative process evolution                   
*  IRESVMAX               Node associated with the maximum residuo.                                                       
*  ITRAP                  Newton-Raphson iteration number. This counter variable is 
*                         initialized after convergence is reached (hopefully)
*                         or after reduction time has been done.                                                      
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NCONVI                 Number of convergences reached during the N-R process                                                      
*  RESVMAX                Maximum residuo.                                                      
*  TABSOLUT               Current absolut computation time                      
*  TINC                   Current time increment                                
*
* HISTORY: First coding: Germán Galarza
*          Revision:     Andres Alcolea (also in August-98)
*
********************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*10 EQUATION

C_____________________________ If user needs the results...

      IF (IOWNR.NE.0) THEN

C_____________________________ Information will be write each IOWNR convergences
C_____________________________ or if solution time is an observation time.

		IF (MOD(NCONVI,IOWNR).EQ.0.OR.INDENDDT.EQ.1.OR.NCONVI.EQ.0) THEN

C_____________________________ Writes the 'title'

			IF (ITRAP.EQ.1) THEN

				WRITE(MAINF,1001)EQUATION,TABSOLUT,TINC


 1001				FORMAT(/,'PROCESS EVOLUTION SOLVING',A10,'EQUATION',
     ; /,'SOLUTION TIME=',E13.7,'  CURRENT TIME INCREMENT=',E13.7,//,
     ; ' ITERATION      MAXIMUM        NODE MAX.     MAXIMUM      ',
     ; 'NODE OF MAX.',/,
     ; 16X,'RESIDUE        RESIDUE    VAR. CHANGE    VAR. CHANGE')

			ENDIF

		WRITE(MAINF,1002)ITRAP,RESVMAX,IRESVMAX,DELVMAX,IDELVMAX

 1002		FORMAT(I10,4X,G12.5,3X,I8,3X,G12.5,3X,I12)

		ENDIF ! MOD(NCONVI...

		IF(IOREDTIMV.EQ.1) THEN

			WRITE(MAINF,*)
     &       (' TIME STEP REDUCTION BECAUSE OF CONVERGENCE PROBLEMS')
		END IF

      ENDIF                 !IOWNR  <> 0

      RETURN

      END SUBROUTINE WRITE_NONLIN_INFO
