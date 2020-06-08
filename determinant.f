      SUBROUTINE DETERMINANT(DET,LENGTH,N,NP,HESS)
      IMPLICIT NONE
      INTEGER*4 N,NP,INDX(N),NMAX,LENGTH,POS,I,J
      REAL*8 D,A(NP,NP),TINY,HESS(LENGTH),DET
      PARAMETER (NMAX=500,TINY=1.0E-20) 
      INTEGER*4 IMAX,K,CODE
      REAL*8  AMAX,DUM, SUM,VV(NMAX)
 
*
*subroutine based on one from
*NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
*
* Copyright (C) 1986-1992 by Cambridge University Press.
* Programs Copyright (C) 1986-1992 by Numerical Recipes Software.
* Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
* the LU decomposition of a rowwise permutation of itself. a and n are input. a is output,
* arranged as in equation (2.3.14) above; indx(1:n) is an output vector that records the
* row permutation eected by the partial pivoting; d is output as .1 depending on whether
* the number of row interchanges was even or odd, respectively.
*
* comments by Luit Jan:
* subroutine arguments have been changed so that ony the determinant is 
* given as output. The lu decomposition  is not stored.
* input has been changed to work with symmetric storage mode
*
      !vv stores the implicit scaling of each row.

* 1) write the input matrix b , which is in simetric storage mode,
*    into a , which is in full storage mode




      POS=0
      DO I=1,NP
         DO J=1,I
            POS=POS+1
            A(I,J)=HESS(POS)
            A(J,I)=HESS(POS)
         ENDDO
      ENDDO
      IF (NP .EQ. 1) THEN
         DET=A(1,1)
         GOTO 100
      ENDIF
* 2) use the NR-routine to caluclate the lu decompostion of A.
*    This is the numerical recipies routine

      D=1
      CODE=0

      DO I=1,N
       AMAX=0.d0
       DO J=1,N
        IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
       END DO ! j loop
       IF(AMAX.LT.TINY) THEN
        CODE = 1
        RETURN
       END IF
       VV(I) = 1.d0 / AMAX
      END DO ! i loop

      DO J=1,N
        DO I=1,J-1
          SUM = A(I,J)
          DO K=1,I-1
            SUM = SUM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUM
        END DO ! i loop
        AMAX = 0.d0
        DO I=J,N
          SUM = A(I,J)
          DO K=1,J-1
            SUM = SUM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUM
          DUM = VV(I)*DABS(SUM)
          IF(DUM.GE.AMAX) THEN
            IMAX = I
            AMAX = DUM
          END IF
        END DO ! i loop  
   
        IF(J.NE.IMAX) THEN
          DO K=1,N
            DUM = A(IMAX,K)
            A(IMAX,K) = A(J,K)
            A(J,K) = DUM
          END DO ! k loop
          D = -D
          VV(IMAX) = VV(J)
        END IF

        INDX(J) = IMAX
        IF(DABS(A(J,J)) .LT. TINY) A(J,J) = TINY

        IF(J.NE.N) THEN
          DUM = 1.d0 / A(J,J)
          DO I=J+1,N
            A(I,J) = A(I,J)*DUM
          END DO ! i loop
        END IF 
      END DO ! j loop

 
* 3) calculate the determinant of the matrix

        DO J=1,N
           D=D*A(J,J)
        ENDDO 
        DET=D

100   CONTINUE
      RETURN
      END
