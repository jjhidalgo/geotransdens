      SUBROUTINE COMP_FLOWGRAV
     &          (AREA       ,BFLU       ,DFLUDFLU   ,DPARELDH   ,GRAVEL
     &          ,GRDFF      ,IOCALCDEVF ,IODIM      ,KXX        ,LDIM
     &          ,LMXNDL     ,LNNDEL     ,NUMNP      ,NPPEL      ,NUMEL
     &          ,PAREL)


      IMPLICIT NONE

C-------------------- External

      INTEGER*4::IOCALCDEVF,IODIM,LMXNDL,NUMNP,NPPEL,NUMEL
	INTEGER*4::KXX(LMXNDL,NUMEL),LDIM(NUMEL),LNNDEL(NUMEL)
	REAL*8::AREA(NUMEL),BFLU(NUMNP),DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &,DPARELDH(NPPEL,NUMEL),GRAVEL(NUMEL,3)
     &,GRDFF(IODIM,LMXNDL,NUMEL)
     &      ,PAREL(NUMEL,NPPEL)

C-------------------- Internal

      INTEGER*4::I,IK_POS,J,J1,J2,K,L,LD,NNUD,NODE
	REAL*8::AREAL,SUMGRAVI

	INTEGER*4::IND(3, 3, 3)

C------------------------- Array IND is used to simplify the coding of
C------------------------- conductivity matrix times the gradient of FEM
C------------------------- interpolation functions and gravity

      DATA ((IND(I,J,3),I=1,3),J=1,3)/1,4,5,4,2,6,5,6,3/
      DATA ((IND(I,J,2),I=1,2),J=1,2)/1,3,3,2/,IND(1,1,1)/1/

C-------------------- First executable statement.

      DO L=1,NUMEL

	    AREAL = AREA(L)
	    LD = LDIM(L)
	    NNUD = LNNDEL(L)

          SUMGRAVI = 0D0

          DO I=1,NNUD

              NODE = KXX(I,L)
              SUMGRAVI = 0D0

              DO J1=1,LD

                  DO J2=1,LD

                      SUMGRAVI = SUMGRAVI - PAREL(L,IND(J1,J2,LD))*
     &                                  GRAVEL(L,J2)*GRDFF(J1,I,L)*AREAL

                      IF (IOCALCDEVF.NE.0) THEN

                          DO K=1,NNUD !IDIMDTRH ? - 1?

	                        IK_POS = (I - 1)*NNUD + K
                              DFLUDFLU(L,IK_POS) = DFLUDFLU(L,IK_POS)
     &                       + DPARELDH(IND(J1,J2,LD),L)*GRAVEL(L,J2)*
     &                           GRDFF(J1,I,L)*AREAL

                          END DO !K=1,IDIMDTRH

                      END IF !IOCALCDEVF.NE.0

                  END DO !J2=1,LD

              END DO !J1=1,LD

              BFLU(NODE) = BFLU(NODE)+SUMGRAVI

          END DO !I=1,NNUD

	END DO !L=1,NUMEL

	END SUBROUTINE COMP_FLOWGRAV