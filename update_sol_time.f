        SUBROUTINE UPDATE_SOL_TIME
     ;(DTIMEF    ,DTIMET   ,DTINITIAL
     ;,IENTRY   ,INDENDDT  ,INTI     ,IOEQT    ,IREDTIMC
     ;,IREDTIMH ,NINT     ,TABSOLUT  ,THETAF   ,THETAT   ,TICAL
     ;,TICALAN  ,TINC     ,TINTERVOBS,TIME     ,ISOLEQ)


C PURPOSE This routine update both the relative and the absolute solution
C          time. Then, the next solution of the flow and or transport
C          equations will be done at time computed here.
C ACRONYM: UPDATE SOLution TIME
C DESCRIPTION: Serie of arithmetic operations attending to update the time
C              solution.
C REFERENCES: author GALARZA G.
C EXTERNALS
C ientry=indicates how many times TRANSIN has acceed to the routine
C        during the current observation intervale
C ticalan=diference between the last successful solution time and the initial
C         one of the current observation intervale.
C tical=diference between the updated solution time and the initial time
C         of the current observation intervale.
C dtinitial=initial time step size adopted for this observation intervale
C dtimef=realtive (to the observation intervale length) value of (t(k+theta)-tobs(inti)),
C         where the first term is the time where the flow equation will be  solved, 
C         and the second is the initial time of the current observation intervale.
C dtimet= idem for transport
C tabsolut=absolute current solution time
C time=vector containing the observation times
C inti= current observation intervale index.
C epsflu=weighted relative flow solution time
C epstra=weighted relative transport solution time
C tintervobs=length of the current observation intervale
C iredtimh=reduction time step in flow solution indicator
C iredtimc=reduction time step in transport solution indicator
C nint=total observation times number
C HISTORY programed in may of 1997
C INTERNALS:

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       DIMENSION TIME(NINT),ISOLEQ(NINT,4)


       IF (IENTRY.EQ.1) THEN  !first step of the intervale

         TICALAN=0D0
         TICAL=DTINITIAL

       ELSE IF ((IREDTIMH+IREDTIMC).NE.0) THEN !reduction of the time step

         TICAL=TICALAN+TINC

       ELSE                   !generig step

         TICALAN=TICAL
         TICAL=TICALAN+TINC

       ENDIF

C_________________________Solving only flow or coupled flow and transport

       IF (IOEQT.EQ.1 .OR. IOEQT.EQ.3) THEN
          IF (ISOLEQ(INTI,1).EQ.0) THEN
             DTIMEF=(TICALAN+THETAF*TINC)/TINTERVOBS ! only solving trans. flow
          ELSE
             DTIMEF=0.D0
          ENDIF
       ENDIF

C________________________Solving only transport or coupled flow and transport

       IF (IOEQT.EQ.2 .OR. IOEQT.EQ.3) THEN
          IF (ISOLEQ(INTI,2).EQ.0) THEN
             DTIMET=(TICALAN+THETAT*TINC)/TINTERVOBS  ! only solving trans. tpt
          ELSE
             DTIMET=0.D0
          ENDIF
       ENDIF

C_________________________Both cases

       TABSOLUT=TIME(INTI)+TICAL

       IF (DABS(TICAL-TINTERVOBS).LT.1.D-5*DABS(TINC)) INDENDDT=1

       RETURN
       END
