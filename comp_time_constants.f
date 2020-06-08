       SUBROUTINE COMP_TIME_CONSTANTS
     ;    (KINT       ,IOINV      ,INTI       ,NUMITER    ,IOPINVDT
     ;    ,NINT       ,TIME       ,DTINITIAL  ,DTMXDS     ,DTMAXIO
     ;    ,DTPREVINV  ,DTAVGIO    ,TINTERVOBS ,IOFLLI     ,IOTRLI)

C PURPOSE: This routine computes some parameters which controle the variable
C time in the simulation during a particular INTI observation time intervale
C ACRONYM: COMPute TIME parameters CONSTANTS during this observation
C          intervale
C DESCRIPTION: serie of simple computations
C REFERENCES: author GALARZA G.
C EXTERNALS
C INTERNALS
C HISTORY programed in may of 1997 By G.G.


       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       DIMENSION
     ;    DTPREVINV(NINT)   ,TIME(NINT)   ,DTMXDS(NINT)   ,KINT(NINT)

       TINTERVOBS=TIME(INTI+1)-TIME(INTI)
       DTAVGIO=TINTERVOBS/KINT(INTI)

       IF (IOPINVDT.NE.0.AND.NUMITER.GT.1.AND.IOINV.GT.0) THEN

         DTINITIAL=DTPREVINV(INTI)

       ELSE

         DTINITIAL=DTAVGIO

       ENDIF

       IF (IOFLLI+IOTRLI.NE.0) DTMAXIO=DTMXDS(INTI)

       RETURN
       END
