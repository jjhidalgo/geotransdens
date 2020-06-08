      SUBROUTINE UNSHUFFLE_LHS
     &          (INTI     ,LHS      ,NUMNP    ,NBAND)

      IMPLICIT NONE
C------------------------- External

      INTEGER*4::INTI,NUMNP,NBAND

     
      REAL*8::LHS(2*NUMNP,4*NBAND+3)

C------------------------- Internal

      INTEGER*4::I,ICOL,IDIAG,INODE,IROWF,IROWT,J,JNODE,NCOLS

      CHARACTER::strFmt1*20

      REAL*8,ALLOCATABLE::DER(:,:)

C------------------------- First Executable Statement

      ALLOCATE(DER(2*NUMNP,2*NUMNP))
      
      DER = 0D0

C------------------------- Auxiliar string for format.

      strFmt1 = ''
      WRITE (strFmt1,*) NUMNP
      strFmt1 = '(I5,'//Trim(AdjustL(strFmt1))//'G22.15)'

C------------------------- Posición de la diagonal en la matriz en banda no
C------------------------- simétrica grande.

      IDIAG = 2*NBAND + 2

      NCOLS = 4*NBAND+3

C------------------------- Número de columnas en una matriz en banda no
C------------------------- simétrica normal.



C------------------------- Se recorren las filas de LHS, equivale a
C------------------------- recorrer los nudos.

      DO INODE=1,2*NUMNP

C------------------------- Se rrecorren las columnas

          DO ICOL = 1,NCOLS

C------------------------- Se calcula la posición global del elemento
C------------------------- (INODE,ICOL) en una matriz de (2*N,2*N)

              JNODE = ICOL + INODE - IDIAG

              
C------------------------- Almacenamiento en la matriz auxiliar para la
C------------------------- escritura posterior. No hay es necesario actualizar
C------------------------- las posiciones que son siempre nulas (y para las que
C------------------------- el algoritmo falla).

               IF(JNODE.GT.0 .AND. JNODE.LE.(2*NUMNP)) THEN
                  DER(INODE,JNODE) =  LHS(INODE,ICOL)
	         END IF

          END DO !ICOL = 1,NCOLS

      END DO !INODE=1,NUMNP


C------------------------- Escritura de las matrices auxiliares.

      WRITE(611,*) 'dFdF en INTI = ',INTI
      WRITE(612,*) 'dFdT en INTI = ',INTI
      WRITE(622,*) 'dTdT en INTI = ',INTI
      WRITE(621,*) 'dTdF en INTI = ',INTI
	
      DO I=1,NUMNP

          IROWF = 2*I -1
          IROWT = 2*I

	  !ICOLF = 2*J -1
          !ICOLT = 2*J

          WRITE(611,strFmt1) I,(DER(IROWF,2*J -1),J=1,NUMNP)
	  WRITE(612,strFmt1) I,(DER(IROWF,2*J),J=1,NUMNP)
	  WRITE(621,strFmt1) I,(DER(IROWT,2*J -1),J=1,NUMNP)
	  WRITE(622,strFmt1) I,(DER(IROWT,2*J),J=1,NUMNP)

c    1     FORMAT(I5,<NUMNP>G22.15)
    
      END DO !I=1,NUMNP

      WRITE(611,*) ''
      WRITE(612,*) ''
      WRITE(622,*) ''
      WRITE(621,*) ''

      DEALLOCATE(DER)

      END SUBROUTINE UNSHUFFLE_LHS
