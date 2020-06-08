       SUBROUTINE COMPUTE_RESID_MAX
     &           (IAD      ,IADN      ,INDFLTR  ,IOCOUPLED,IRESVMAX 
     &           ,NUMNP     ,NUMNV
     &           ,RESVMAX  ,RESVMAXOLD         ,IBVAR     ,RHS
     &           ,ADSC     ,KXX      ,LNNDEL   ,NUMEL     ,IADSC_COLS
     &           ,IADSC_ROWS         ,ITYPADSC ,VCALIT    ,SOLUTION
     &           ,LMXNDL   ,IONEWT)     

*********************************************************
*
* PURPOSE
*
*
* DESCRIPTION
*
* EXTERNAL VARIABLES: ARRAYS
*
*  RHS                                                                          
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IRESVMAX                                                                     
*  NUMNP                  Number of nodes                                       
*  RESVMAX                                                                      
*  RESVMAXOLD                                                                   
*
* INTERNAL VARIABLES: SCALARS
*
*  I                                                                            
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
********************************************************************************

       IMPLICIT NONE

	INTEGER*4::I,INDFLTR,IOCOUPLED,IPOS,IRESVMAX,ISTEP,NUMNP,NUMNV
     &          ,IADSC_COLS,IADSC_ROWS,IONEWT,ITYPADSC,LMXNDL,NUMEL
	INTEGER*4::IBVAR(NUMNP),KXX(LMXNDL,NUMEL),LNNDEL(NUMEL),IAD(*)
     &          ,IADN(*)

      REAL*8::RESVMAX  ,RESVMAXOLD
	REAL*8::RHS(NUMNV),ADSC(IADSC_ROWS,IADSC_COLS),VCALIT(NUMNP)
     &       ,SOLUTION(NUMNP)

       RESVMAXOLD=RESVMAX   !update the previous maximum residuo

       RESVMAX=0D0

      ISTEP = 0

C------------ This way the subroutine can be used for
C------------ both coupled and non-coupled case.
	IF (IONEWT.EQ.1) THEN
C------------ Residue if Newton-Raphson method is used.

C------------ This way the subroutine can be used for
C------------ both coupled and non-coupled case.

		IF (IOCOUPLED.EQ.1 .AND. INDFLTR.EQ.0) THEN

			ISTEP = 1

		END IF

		 DO I=1,NUMNP

		   IF (IBVAR(I).NE.1) THEN
          
			 IPOS = (1 + IOCOUPLED)*I - ISTEP
			 IF (DABS(RHS(IPOS)).GT.RESVMAX) THEN

			   IRESVMAX=I

			   RESVMAX=DABS(RHS(IPOS))      !Updates the current maximum residuo

			 END IF
            
		   ENDIF

		 END DO

	ELSE
C------------ Residue if Picard method is used.

          SOLUTION=0.D0

		CALL PROD_MAT_VEC
     ;        (1.D0    ,IAD         ,IADN
     ;        ,IADSC_COLS  ,IADSC_ROWS  ,NUMNP
     ;        ,1       ,ITYPADSC    ,1           ,LMXNDL         
     ;        ,NUMEL   ,NUMNP       ,KXX         ,LNNDEL         
     ;        ,ADSC    ,SOLUTION    ,VCALIT)
		
	    SOLUTION=SOLUTION-RHS
		 DO I=1,NUMNP

		   IF (IBVAR(I).NE.1) THEN
          
			 IF (DABS(SOLUTION(I)).GT.RESVMAX) THEN

			   IRESVMAX=I

			   RESVMAX=DABS(SOLUTION(I))      !Updates the current maximum residuo

			 END IF
            
		   ENDIF

		 END DO
		 SOLUTION=0.D0
	 END IF

       RETURN

       END
