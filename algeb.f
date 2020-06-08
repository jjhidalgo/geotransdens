       SUBROUTINE ZERO_ARRAY_I (INP,N)
       IMPLICIT REAL*8 (A-H,O-Z)
       DATA INTZERO/0/
       DIMENSION INP(N)
       DO 10 I=1,N
 10       INP(I)=INTZERO
       RETURN
       END

********************************************************************************
********************************************************************************

       SUBROUTINE ZERO_ARRAY (A,N)
       IMPLICIT REAL*8 (A-H,O-Z)
       DATA ZERO/0.D0/
       DIMENSION A(N)
       DO 10 I=1,N
 10       A(I)=ZERO
       RETURN
       END

************************************************************************
************************************************************************

       SUBROUTINE ARRAY_DIFF (V1,V2,N,DEL)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION V1(N),V2(N)
       DO 10 I=1,N
 10       V1(I)=(V2(I)-V1(I))/DEL
       RETURN
       END

************************************************************************
************************************************************************

       SUBROUTINE ARRAY_POND (V1,V2,N,PON)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION V1(N),V2(N)
       CPON=1.D0-PON
       DO 10 I=1,N
 10       V1(I)=PON*V1(I)+CPON*V2(I)
       RETURN
       END

************************************************************************
************************************************************************

       SUBROUTINE EQUAL_ARRAY (A,B,N)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION A(N),B(N)
       DO 10 I=1,N
 10       A(I)=B(I)
       RETURN
       END


************************************************************************
************************************************************************

       SUBROUTINE COMP_AUX (VAUX1,VAUX2,VCALIT,VCALAN,THETA,TINC,NUMNP)

******************************************************************
*****  COMPUTES HAUX1, HAUX2, CAUX1 & CAUX2
******************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION VAUX1(NUMNP),VAUX2(NUMNP),VCALIT(NUMNP),VCALAN(NUMNP)

       CH=1.D0-THETA
       DO 10 I=1,NUMNP
          VAUX1(I)=THETA*VCALIT(I)-CH*VCALAN(I)
 10       VAUX2(I)=(VCALIT(I)-VCALAN(I))/TINC

       RETURN
       END

************************************************************************
************************************************************************

       SUBROUTINE VALMAX (ARRAY,VALMX,NUMNP,IRES)

******************************************************************
****  COMPUTES ARRAY'S MAXIMUM VALUE 
******************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION ARRAY(NUMNP)

       VALMX=0D0
       DO I=1,NUMNP
         IF (ARRAY(I).GT.VALMX) THEN
            VALMX=ARRAY(I)
            IRES=I
         END IF
       END DO

       RETURN
       END
******************************************************************

       SUBROUTINE ENS_IND_T9 (AFLU,DFLU,BFLU,HCAL,IBCOD)

******************************************************************
*** COMPUTES CONTRIBUTIONS OF PREVIOUS TIME STEP TO BFLU AND
*** CORRIGES BOUNDARY CONDITIONS OF TYPE 1.
*** (DFLU IS IN CONSISTENT FORM)
******************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       INCLUDE 'COMMON.FOR'
       DIMENSION AFLU(NUMNP,NBAND1),DFLU(NUMNP,IDIMDFLU),BFLU(NUMNP),
     ;           HCAL(NUMNP),IBCOD(NUMNP)

       TH=THETAF-1
       DD=1/TINC

*** MULTIPLY BY LOWER TRIANGLE

       IN=NBAND1
       DO 10 K=1,NBAND1
          JI=1
          DO 20 I=IN,NUMNP
             IF(IBCOD(I).EQ.1) GOTO 20
             BFLU(I)=BFLU(I)+(TH*AFLU(I,K)+DD*DFLU(I,K))*HCAL(JI)
 20          JI=JI+1
 10       IN=IN-1

*** MULTIPLY BY UPPER TRIANGLE

       IN=NBAND1
       DO 40 K=1,NBAND
          JI=1
          DO 30 I=IN,NUMNP
             IF( IBCOD(JI).EQ.1 ) GOTO 30
             BFLU(JI)=BFLU(JI)+(TH*AFLU(I,K)+DD*DFLU(I,K))*HCAL(I)
  30         JI=JI+1
  40      IN=IN-1

       RETURN
       END

******************************************************************

       SUBROUTINE ENS_IND_T10 (AFLU,DFLU,BFLU,HCAL,IBCOD,TINC
     ;,THETAF,NUMNP,NBAND,IDIMDFLU)      !prov


******************************************************************
*** COMPUTES CONTRIBUTIONS OF PREVIOUS TIME STEP TO BFLU AND
*** CORRIGES BOUNDARY CONDITIONS OF TYPE 1.
*** (DFLU IS IN LUMPED FORM)
******************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       DIMENSION AFLU(NUMNP,NBAND+1),DFLU(NUMNP,IDIMDFLU),BFLU(NUMNP),
     ;           HCAL(NUMNP),IBCOD(NUMNP)
       TH=THETAF-1
       DD=1/TINC

       NBAND1 =NBAND+1
 
*** MULTIPLY BY LOWER TRIANGLE

       IN=NBAND1
       DO 10 K=1,NBAND
          JI=1
          DO 20 I=IN,NUMNP
             IF(IBCOD(I).EQ.1) GOTO 20
             BFLU(I)=BFLU(I)+TH*AFLU(I,K)*HCAL(JI)
 20          JI=JI+1
 10       IN=IN-1

*** MULTIPLY BY UPPER TRIANGLE

       IN=NBAND1
       DO 40 K=1,NBAND
          JI=1
          DO 30 I=IN,NUMNP
             IF( IBCOD(JI).EQ.1 ) GOTO 30
             BFLU(JI)=BFLU(JI)+TH*AFLU(I,K)*HCAL(I)
  30         JI=JI+1
  40      IN=IN-1

*** MULTIPLY BY DIAGONAL


       DO 50 I=1,NUMNP                                              
 50       IF(IBCOD(I).NE.1) BFLU(I)=BFLU(I)+
     ;                      (TH*AFLU(I,NBAND1)+DD*DFLU(I,1))*HCAL(I)

       RETURN
       END

********************************************************************************
********************************************************************************
       SUBROUTINE WRIVECTOR(NUND,VAR,IDIM)
********************************************************************************
***    WRITE A VECTOR IN A SCRATCH FILE, COMPONENT BY COMPONENT               ***
********************************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION VAR(IDIM)
       DO I=1,IDIM
         WRITE(NUND)VAR(I)
       END DO
       RETURN
       END
********************************************************************************
       SUBROUTINE READVECTOR(NUND,VAR,IDIM)
********************************************************************************
***    READ A VECTOR IN A SCRATCH FILE, COMPONENT BY COMPONENT               ***
********************************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION VAR(IDIM)
       DO I=1,IDIM
         READ(NUND)VAR(I)
       END DO
       RETURN
       END

********************************************************************************
      SUBROUTINE MUL_SYMMAT_VEC(IDIMMAT,MATRIX,RESULT,VECTOR)
********************************************************************************
***    MULTIPLIES A SYMMETRIC MATRIX (VECTOR) TIMES A VECTOR                 ***
********************************************************************************

      IMPLICIT NONE
      INTEGER*4 IDIMMAT,IFIL,JCOL,IPOSMAT
      REAL*8 MATRIX(IDIMMAT*(IDIMMAT+1)/2),RESULT(IDIMMAT)
     ;      ,VECTOR(IDIMMAT)

      CALL ZERO_ARRAY(RESULT,IDIMMAT)
      DO IFIL=1,IDIMMAT
        DO JCOL=1,IDIMMAT

          IF (IFIL.GT.JCOL) THEN
            IPOSMAT=IFIL*(IFIL-1)/2+JCOL
          ELSE
            IPOSMAT=JCOL*(JCOL-1)/2+IFIL
          END IF

          RESULT(IFIL)=RESULT(IFIL)+MATRIX(IPOSMAT)*VECTOR(JCOL)
        END DO
      END DO

      RETURN
      END

************************************************************************
************************************************************************

       SUBROUTINE CHECK_ATRA (NBAND2,NUMNP,ATRA)

       IMPLICIT REAL*8 (A-H,O-Z)

       DIMENSION ATRA(NUMNP,NBAND2)

       DO I=1,NUMNP
          SUM=0.D0
          DO J=1,NBAND2
             SUM=SUM+ATRA(I,J)
          ENDDO
          WRITE(90,*) I,SUM
       ENDDO

       RETURN
       END

**************************************************************************
**************************************************************************

      SUBROUTINE IO_SUB (NAME_SUBR,INDEX)

********************************************************************************
*
* PURPOSE Echoes a message with the name of actual subroutine. Depending on
*         INDEX variable, message is written at the very beginning of the routine
*         (INDEX=0) or at the very end (INDEX=1)
*
********************************************************************************

      IMPLICIT NONE
      INTEGER*4 INDEX
      CHARACTER*(*) NAME_SUBR

      IF (INDEX.EQ.0) THEN
        WRITE (*,10) NAME_SUBR
      ELSE IF (INDEX.EQ.1) THEN
        WRITE (*,20) NAME_SUBR
      ENDIF

 10   FORMAT(//,' Coming in: ',A25)
 20   FORMAT(//,' Coming out: ',A25)

      RETURN
      END
