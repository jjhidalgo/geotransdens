      SUBROUTINE COMP_DIV_Q
     &          (AFLU     ,AREA     ,BUOYANCY ,COORD    ,DENSITY
     &          ,DIVQ     ,GP_COORD ,GRADLOC  ,HAUX1    ,IAD_S
     &          ,IADN_S   ,IDIMAFLU ,IODENS   ,IODIM    ,ISOZ
     &          ,ITYPAFLU ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL
     &          ,LTYPE    ,LXPAREL  ,MAXNB    ,MAXPG    ,NPAREL
     &          ,NPPEL    ,NTYPEL   ,NUMEL    ,NUMNP    ,NZTRA
     &          ,PAREL    ,POINTWEIGHT        ,THETAT)
********************************************************************************
*
* PURPOSE Computes Div·q
*
*  Div·q = AFLU+h(k+thetaf) - Buoyancy
*
********************************************************************************

       
      IMPLICIT NONE

C--------------------------- External

      INTEGER*4::IDIMAFLU ,IODENS   ,IODIM    ,ITYPAFLU ,LMXNDL
     &          ,MAXNB    ,MAXPG    ,NPAREL   ,NPPEL    ,NTYPEL
     &          ,NUMEL    ,NUMNP    ,NZTRA

      REAL*8::BETAC    ,THETAT


      INTEGER*4::IAD_S(MAXNB,NUMNP) ,IADN_S(NUMNP)        ,ISOZ(NZTRA)
     &          ,KXX(LMXNDL,NUMEL)  ,LDIM(NUMEL)          ,LNNDEL(NUMEL)
     &          ,LTYPE(NUMEL)       ,LXPAREL(NUMEL,NPAREL)


      REAL*8::AFLU(NUMEL,IDIMAFLU)         ,AREA(NUMEL)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL) ,COORD(NUMNP,3)
     &       ,DIVQ(NUMNP)                  ,DENSITY(NUMEL)
     &       ,GP_COORD(6,8,IODIM)          ,GRADLOC(IODIM,LMXNDL,MAXPG)
     &       ,HAUX1(NUMNP)                 ,PAREL(NUMEL,NPPEL)
     &       ,POINTWEIGHT(MAXPG,NTYPEL)

C--------------------------- Internal

      Real*8::DUMMY(1)

C---------------------------  First executable statement.

C-------------------- Buoyancy correction if density dependent flow.

      IF (IODENS.EQ.1) THEN

          CALL COMP_BFLU_BUOYANCY
     &        (AREA     ,BETAC    ,DIVQ     ,BUOYANCY ,COORD
     &        ,DUMMY    ,DENSITY  ,DUMMY    ,GP_COORD ,GRADLOC
     &        ,0        ,IODIM    ,ISOZ     ,KXX      ,LDIM
     &        ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPAREL  ,MAXPG
     &        ,NPAREL   ,NPPEL    ,NTYPEL   ,NUMEL    ,NUMNP
     &        ,NZTRA    ,PAREL    ,POINTWEIGHT        ,THETAT)

C-------------------- Sign for BFLU is opposite to the one
C-------------------- needed for div·q

          DIVQ(:) = -1D0*DIVQ(:)

      END IF !IODENS.EQ.1

C--------------------  Computes AFLU*h(k+thetaf)

      CALL PROD_MAT_VEC
     &    (1D0      ,IAD_S    ,IADN_S   ,IDIMAFLU ,NUMEL    ,NUMNP
     &    ,1        ,ITYPAFLU ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &    ,KXX      ,LNNDEL   ,AFLU     ,DIVQ  ,HAUX1)


      END SUBROUTINE COMP_DIV_Q