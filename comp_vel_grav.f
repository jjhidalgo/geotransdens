      SUBROUTINE COMP_VEL_GRAV
     &          (AREA       ,DVDH       ,GRAVEL     ,IOCALCDEVF
     &          ,IODIM      ,ISOZ       ,L          ,LDIM
     &          ,LMXNDL     ,LNNDEL     ,LXPAREL    ,NPAREL
     &          ,NPPEL      ,NUMEL      ,NZTRA      ,VEL)
********************************************************************************
*
* PURPOSE
*
*    Computes gravitational term of Darcy's velocity and its derivatives for
*    constant density.
*
* DESCRIPTION
*
*    Computes gravitational term of Darcy's velocity and its derivatives for
*    constant density.
*
*
********************************************************************************

      IMPLICIT NONE

C---------------------------  External.

      INTEGER*4::IOCALCDEVF,IODIM,L,LMXNDL,NPAREL,NPPEL,NUMEL,NZTRA

	REAL*8::AREA(NUMEL),DPARELDH(NPPEL,NUMEL),DVDH(LMXNDL,IODIM)
     &       ,GRAVEL(NUMEL,3),PAREL(NUMEL,NPPEL),VEL(IODIM)

	INTEGER*4::ISOZ(NZTRA),LDIM(NUMEL),LNNDEL(NUMEL)
     &,LXPAREL(NUMEL,NPAREL)

C---------------------------  Internal.

      INTEGER*4::I,IDIM,ISMAX,ISZ,J,K,KDIM,LD,NNUD,NZONE

      REAL*8::AREAL

	REAL*8::DTRANS(6),GRAV(3),SUMA(3),TRACT(6)
	INTEGER*4::IND(3, 3, 3)

C------------------------- Array IND is used to simplify the coding of
C------------------------- conductivity matrix times the gradient of FEM
C------------------------- interpolation functions and gravity

      DATA ((IND(I,J,3),I=1,3),J=1,3)/1,4,5,4,2,6,5,6,3/
      DATA ((IND(I,J,2),I=1,2),J=1,2)/1,3,3,2/,IND(1,1,1)/1/

C-------------------- First executable statement.

      DATA ((IND(I,J,3),I=1,3),J=1,3)/1,4,5,4,2,6,5,6,3/
      DATA ((IND(I,J,2),I=1,2),J=1,2)/1,3,3,2/,IND(1,1,1)/1/


C--------------------------- For each element...

C--------------------------- sets element properties

      LD = LDIM(L)
	NNUD = LNNDEL(L)
	AREAL = AREA(L)
      NZONE = LXPAREL(L,1)      ! T zone of current element
      ISZ = ISOZ(NZONE)         ! Anisotropy degree of the zone
      ISMAX = MAX(ISZ,LD)       ! Number of T-tensor components

C--------------------------- Identifies the value of the T components

      TRACT(1:ISMAX) = PAREL(L,1:ISMAX)
      DTRANS(1:ISMAX) =  DPARELDH(1:ISMAX,L)
      GRAV(1:3) = GRAVEL(L,1:3)

      SUMA = 0D0
              
C--------------------------- Finally the product of GRAV*K

      DO IDIM=1,LD

          DO KDIM=1,LD 

              SUMA(IDIM) = SUMA(IDIM)
     &                   + TRACT(IND(IDIM,KDIM,LD))*GRAV(IDIM)

              IF (IOCALCDEVF.NE.0) THEN

                  DO K=1,NNUD
                      DVDH(K,IDIM) = DVDH(K,IDIM)
     &              + DTRANS(IND(IDIM,KDIM,LD))*GRAV(IDIM)

                  END DO !K=1,IDIMDTRH

              END IF !IOCALCDEVF.NE.0
          END DO !KDIM=1,LD
      END DO !IDIM=1,LD

C--------------------------- The computed values are stored.
C--------------------------- The minus sign comes from Darcy's law.

      VEL(1:LD) = VEL(1:LD) + SUMA(1:LD)/AREAL

	END SUBROUTINE COMP_VEL_GRAV