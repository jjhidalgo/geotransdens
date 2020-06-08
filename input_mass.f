      SUBROUTINE INPUT_MASS
     &          (ATRADSC  ,IAD_S    ,IADD_S   ,IADN_S   ,IATRADSC_COLS
     &          ,IATRADSC_ROWS      ,IBTCO    ,ITYPTRADSC
     &          ,KXX      ,LMXNDL   ,LNNDEL   ,MAXNB    ,NBAND
     &          ,NPPNP    ,NUMEL    ,NUMNP    ,PARNP    ,THETAT)
********************************************************************************
*
* PURPOSE
*
*      Adds contribution of input mass boundary condition to  diagonal of
*      ATRADSC.
*
* DESCRIPTION
*
*
* EXTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
* INTERNAL VARIABLES: SCALARS
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY
*
*     JHG     02-2006     First coding.
*
********************************************************************************


      IMPLICIT NONE


C------------------------- External

      INTEGER*4::IATRADSC_COLS      ,IATRADSC_ROWS      ,ITYPTRADSC
     &          ,LMXNDL   ,MAXNB    ,NBAND    ,NPPNP    ,NUMEL
     &          ,NUMNP

      REAL*8::THETAT

      INTEGER*4::IADN_S(NUMNP),IAD_S(MAXNB, NUMNP),IADD_S(NUMNP)
     &          ,IBTCO(NUMNP) ,KXX(LMXNDL,NUMEL)  ,LNNDEL(NUMEL)

      REAL*8::ATRADSC(IATRADSC_ROWS,IATRADSC_COLS)
     &       ,PARNP(NUMNP,NPPNP)

C------------------------- Internal

      INTEGER*4::I        ,IBT

      REAL*8,ALLOCATABLE::MASS(:)


C------------------------- First executable statment

      ALLOCATE(MASS(NUMNP))


      DO I=1,NUMNP

          IBT = IBTCO(I)

          IF (IBT.EQ.4) THEN

              MASS(I) = PARNP(I,4)

          ELSE

              MASS(I) = 0D0

          END IF !IBT.EQ.4

      END DO !I=1,NUMNP


      CALL ASSEMBLE
     &   (THETAT   ,1        ,NUMNP    ,IATRADSC_COLS
     &   ,IATRADSC_ROWS      ,1        ,ITYPTRADSC         ,LMXNDL
     &   ,MAXNB    ,NBAND    ,NUMNP    ,NUMEL    ,MASS     ,ATRADSC
     &   ,IAD_S    ,IADD_S   ,IADN_S   ,KXX      ,LNNDEL)

      DEALLOCATE(MASS)

      END SUBROUTINE INPUT_MASS