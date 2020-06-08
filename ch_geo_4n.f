       SUBROUTINE CH_GEO_4N (LMXNDL, NUMEL, NUMNP, KXX, NX, X, Y, NE)

*************************************************************************
***  CHECK GEOMETRICS CHARACTERISTICS OF FOUR NODES ELEMENTS
*************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION KXX(LMXNDL,NUMEL),NX(8),X(NUMNP),Y(NUMNP)


       YMIN=DMIN1(Y(NX(1)),Y(NX(2)),Y(NX(3)),Y(NX(4))) 

       J2=0
       ICONT=0
       DO I=1,4
          IF (Y(NX(I)).EQ.YMIN) THEN
             IF (ICONT.EQ.0) THEN
                YMIN=Y(NX(I))
                J1=I
                ICONT=1
             ELSE
                J2=I
             ENDIF
          ENDIF
       ENDDO

       IF (J2.EQ.0) THEN
          WRITE(6,*) ' NO ES UN RECTANGULO EL ELTO',NE
          RETURN
       ENDIF

       IF (X(NX(J1)).GT.X(NX(J2))) THEN
          IAUX=J1
          J1=J2
          J2=IAUX
       ENDIF
         
       DO I=1,4
          IF (I.NE.J1 .AND. I.NE.J2) THEN
             IF (X(NX(I)).EQ.X(NX(J1))) THEN
                J4=I
             ELSE IF (X(NX(I)).EQ.X(NX(J2))) THEN
                J3=I
             ELSE
                WRITE(6,*) ' NO ES UN RECTANGULO EL ELTO',NE
                RETURN
             ENDIF
          ENDIF
       ENDDO

C------------------------- Assigns the new node numbering

       KXX(1,NE)=NX(J1)
       KXX(2,NE)=NX(J2)
       KXX(3,NE)=NX(J3)
       KXX(4,NE)=NX(J4)

       RETURN
       END
