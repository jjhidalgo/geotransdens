      SUBROUTINE WRITE_TEMP_STATE_VARS
     &          (CCALIT   ,HCALIT   ,INDENDDT ,INDSSTR  ,INTI
     &          ,IOEQT    ,IOMHC    ,IOMHH    ,ISOLFL   ,ISOLTR
     &          ,NINT     ,NPBFL    ,NPBTP     ,NUMNP)

	IMPLICIT NONE

C------------------- External

	INTEGER*4::INDENDDT, INDSSTR, INTI, IOEQT, IOMHC, IOMHH,ISOLFL
     &          ,ISOLTR,NINT,NPBFL,NPBTP,NUMNP

	REAL*8::CCALIT(NUMNP,NPBTP), HCALIT(NUMNP,NPBFL)

C------------------- Internal

      Integer*4::I

C------------------- If solving flow...

	IF (IOEQT.NE.2 .AND. ISOLFL.GT.0) THEN

C------------------- ...heads are writen for the first and last times
C------------------- and every IOMHH times.

		IF (IOMHH.NE.0 .AND. INDENDDT.NE.0) THEN !IOWRITE(7)=IOMHH

			IF (INDSSTR.EQ.0 .OR. MOD(INTI+1,IOMHH).EQ.0 
     &           .OR. INTI.EQ.NINT-1 ) THEN

                  DO I=1,NPBFL
			        WRITE(94) HCALIT(:,I)
	            END DO !I=1,NPBFL

			ENDIF
		ENDIF
	ENDIF ! IOEQT.NE.2

C------------------- If solving transport...

	IF (IOEQT.NE. 1 .AND. ISOLTR.GT.0) THEN

C------------------- ...heads are writen for the first and last times
C------------------- and every IOMHH times.
		IF (IOMHC.NE.0 .AND. INDENDDT.NE.0) THEN !IOWRITE(8)=IOMHC
			IF (INDSSTR.EQ.0 .OR. MOD(INTI+1,IOMHC).EQ.0 
     &           .OR. INTI.EQ.NINT-1 ) THEN
			
                  DO I=1,NPBTP
                      WRITE(93) CCALIT(:,I)
	            END DO !I=1,NPBTP

			ENDIF
		ENDIF
      ENDIF !IOEQT.NE.1



      END SUBROUTINE WRITE_TEMP_STATE_VARS