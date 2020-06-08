      SUBROUTINE WRITE_SQR_MATRIX
     &          (A_MAT    ,IACOLS   ,IAROWS   ,IAD_S    ,IADN_S
     &          ,ISPARSE  ,ITYPE    ,MAXNB    ,NUMNP    ,NBAND
     &          ,UNIT)


      IMPLICIT NONE

C-------------------- External

      INTEGER*4::UNIT,NUMNP,IACOLS,IAROWS,ISPARSE,ITYPE,MAXNB,NBAND

      INTEGER*4::IAD_S(MAXNB*ISPARSE,NUMNP*ISPARSE)
     &          ,IADN_S(NUMNP*ISPARSE)

      REAL*8::A_MAT(IAROWS,IACOLS)

C-------------------- Internal

      INTEGER*4::I,ICOL,IDIAG,INODE,IROW,J,JNODE,NCOLS

      CHARACTER::strFmt1*20

      REAL*8,ALLOCATABLE::B_MAT(:,:)

      INTEGER*4,PARAMETER::
     &                     IBANDSYM   = 1
     &                    ,IBANDNOSYM = 0

C------------------------- First Executable Statement

      ALLOCATE(B_MAT(NUMNP,NUMNP))
      
      B_MAT = 0D0

C------------------------- Auxiliar string for format.

      strFmt1 = ''
      WRITE (strFmt1,*) NUMNP
      strFmt1 = '(I5,'//Trim(AdjustL(strFmt1))//'G22.15)'

      IF (ISPARSE.EQ.0) THEN

C------------------------- Posición de la diagonal en la matriz en banda no
C------------------------- simétrica grande.

          IF (ITYPE.EQ.IBANDSYM) THEN

              IDIAG = NBAND + 1

              NCOLS = NBAND+1

          ELSE IF (ITYPE.EQ.IBANDNOSYM) THEN

              IDIAG = NBAND + 1

              NCOLS = 2*NBAND+1

          END IF !ITYPE.EQ. ...


C------------------------- Se recorren las filas de LHS, equivale a
C------------------------- recorrer los nudos.

          DO INODE=1,NUMNP

C------------------------- Se recorren las columnas

              DO ICOL = 1,NCOLS

C------------------------- Se calcula la posición global del elemento
C------------------------- (INODE,ICOL) en una matriz de (N,N)

                  JNODE = ICOL + INODE - IDIAG

              
C------------------------- Almacenamiento en la matriz auxiliar para la
C------------------------- escritura posterior. No hay es necesario actualizar
C------------------------- las posiciones que son siempre nulas (y para las que
C------------------------- el algoritmo falla).

                   IF(JNODE.GT.0 .AND. JNODE.LE.(NUMNP)) THEN

C------------------------- Triangulo inferior (caso simétrico) y matriz
C------------------------- completa (caso no simétrico).

                      B_MAT(INODE,JNODE) =  A_MAT(INODE,ICOL)

C------------------------- Triangulo superior (sólo caso simétrico).

	                IF (ITYPE.EQ.IBANDSYM) THEN
	                    B_MAT(JNODE,INODE) = A_MAT(INODE,ICOL)
	                END IF

                   END IF

              END DO !ICOL = 1,NCOLS

          END DO !INODE=1,NUMNP


      ELSE

          DO I=1,NUMNP

              DO J=1,NUMNP
               
                  CALL FIND(I,J,IROW,IAD_S,IADN_S,MAXNB,NUMNP)

                  IF (IROW.GT.0) THEN

                      B_MAT(I,J) = A_MAT(IROW,I)

                  END IF !IROW.GT.0

              END DO !J=1,NUMNP

          END DO !I=1,NUMNP

      END IF !ISPARSE.EQ.0

C------------------------- Escritura de las matrices auxiliares.

      DO I=1,NUMNP


          WRITE(UNIT,strFmt1) I,(B_MAT(I,J),J=1,NUMNP)

c    1     FORMAT(I5,<NUMNP>G22.15)
    
      END DO !I=1,NUMNP


      DEALLOCATE(B_MAT)


      END SUBROUTINE WRITE_SQR_MATRIX
