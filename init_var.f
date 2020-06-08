          SUBROUTINE INIT_VAR 
     ;(IOINV    ,IOPINITVAR ,IOREGIM  ,ISUMFO   ,MAXITER  ,NUMITER
     ;,NUMNP    ,VARCALIT )                                            

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VARCALIT(NUMNP)

       IF (IOREGIM.NE.0.AND.MAXITER.GT.1.AND.IOPINITVAR.NE.0)THEN
         IF (NUMITER.EQ.1)THEN                
             
C_______________________Open aulixiary files at first iteration

           OPEN(UNIT=50, !FILE='INITH.DAT',
     ;     STATUS='SCRATCH',FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=2*NUMNP)
           OPEN(UNIT=51, !FILE='INITPROVH.DAT',
     ;     STATUS='SCRATCH',FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=2*NUMNP)
           OPEN(UNIT=56, !FILE='INDICEH.DAT',
     ;     STATUS='SCRATCH',FORM='UNFORMATTED')
           OPEN(UNIT=57, !FILE='INDICEPROVH.DAT',
     ;     STATUS='SCRATCH',FORM='UNFORMATTED')

         ELSE 

           IF((ISUMFO.EQ.0.AND.IOINV.NE.2).OR.    
     ;       (IOINV.EQ.2.AND.NUMITER.EQ.2))THEN
    
C_______________________Updates data for initial.

             REWIND(56)
             REWIND(57)
             IEND=0
            
             DO WHILE (IEND.EQ.0)
               READ(57,END=20)TSOL,NORDEN
               WRITE(56)TSOL,NORDEN
             END DO
            
   20        REWIND(57)
            
             DO KREG=1,NORDEN
               READ(51,REC=KREG)VARCALIT
               WRITE(50,REC=KREG)VARCALIT
             END DO

           ELSE

             IF(ISUMFO.NE.0.AND.IOINV.NE.2)THEN
               REWIND(57)
             END IF
           END IF
         END IF
       END IF
       RETURN
       END

