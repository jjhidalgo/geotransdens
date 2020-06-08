       SUBROUTINE UPDATE_FILES_SCRATCH
     ;(IOFLLI   ,IOINV    ,IOPINITC ,IOPINITH ,IOTRLI   ,ISUMFO
     ;,MAXITER  ,NUMITER)


C PURPOSE: Updating or opening files containing solutions (of a total simulation) 
C of the flow or transport equations, which will be used for initializing the nonlinear 
C process in the current inverse problem iteration
C ACRONYM: STORing SOLutions within SCRATCH files
C DESCRIPTION: During the first inverse problem iteration we open files, during later
C iterations, solutions are stored in the corresponding non-provisional files, obtained
C from provisional files, to be used in the initialization. Nevertheless, the updating
C is conditioned to the successful of the previous inverse problem iteration (ISUMFO=0)
C EXTERNALS
C HYSTORY: Programed by G.Galarza at december of 1997

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C ------------------------------------------------ files used for the flow problem
       IF (IOFLLI.NE.0.AND.MAXITER.GT.1.AND.IOPINITH.NE.0)THEN

         IF((ISUMFO.EQ.0.AND.IOINV.NE.2).OR.    !updates data for initializing during 
     ;     (IOINV.EQ.2.AND.NUMITER.EQ.2))THEN   !the current inverse iteration
           REWIND(56)
           REWIND(57)
 8         READ(57,END=18)TSOL,NORDEN
           WRITE(56)TSOL,NORDEN
           GOTO 8
 18        REWIND(57)
           DO KREG=1,NORDEN
             READ(51,REC=KREG)HCALIT
             WRITE(50,REC=KREG)HCALIT
           END DO

         ELSE IF(ISUMFO.NE.0.AND.IOINV.NE.2)THEN

           REWIND(57)

         ENDIF
       ENDIF
         
C ------------------------------------------------ files used for the transport problem

       IF (IOTRLI.NE.0.AND.MAXITER.GT.1.AND.IOPINITC.NE.0)THEN

         IF(ISUMFO.EQ.0)THEN      !updates data for initial.

           REWIND(54)
           REWIND(55)
 9         READ(55,END=19)TSOL,NORDEN
           WRITE(54)TSOL,NORDEN
           GOTO 9
 19        REWIND(55)
           DO KREG=1,NORDEN
             READ(53,REC=KREG)CCALIT
             WRITE(52,REC=KREG)CCALIT
           END DO

         ELSE IF(ISUMFO.NE.0)THEN

           REWIND(55)

         ENDIF

       ENDIF
        
       RETURN
       END
