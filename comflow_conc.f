      SUBROUTINE COMFLOW_CONC
     &          (AFLU     ,AREA     ,ARRC     ,ATRA     ,BETAC
     &          ,BUOYANCY ,CAUDAL   ,CAUX1    ,CAUX2    ,CFLU
     &          ,COECEL   ,CONCFLOW ,COORD    ,CREF     ,DENSITY
     &          ,DENSREF  ,DFLU     ,DTRA     ,GP_COORD ,GRADLOC
     &          ,HAUX1    ,HAUX2    ,IAD_S    ,IADN_S   ,IBCOD
     &          ,IBTCO    ,IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU
     &          ,IDIMDTRA ,IODENS   ,IODIM    ,IORECATRA,IREGFL
     &          ,IREGTR   ,ISOZ     ,ITYPAFLU ,ITYPATRA ,ITYPCFLU
     &          ,ITYPDFLU ,ITYPDTRA ,KXX      ,LDIM     ,LMXNDL
     &          ,LNNDEL   ,LTYPE    ,LXARR    ,LXPAREL  ,MAXNB
     &          ,MAXPG    ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL
     &          ,NUMEL    ,NUMNP    ,NZARR    ,NZTRA    ,PAREL
     &          ,PARNP    ,POINTWEIGHT        ,THETAT)
********************************************************************************
*
* PURPOSE
*
*      Computes concentration sources at nodes with prescribed concentration
*      input mass or concentration leakage boundary condition.
*
* DESCRIPTION
*
*      Computes concentration sources at nodes with prescribed concentration
*      input mass or concentration leakage boundary condition.
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

C------------------------- External

      INTEGER*4::IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU ,IDIMDTRA
     &          ,IODENS   ,IODIM    ,IORECATRA,IREGFL   ,IREGTR
     &          ,ITYPAFLU ,ITYPATRA ,ITYPCFLU ,ITYPDFLU ,ITYPDTRA
     &          ,LMXNDL   ,MAXNB    ,MAXPG    ,NPAREL   ,NPPEL
     &          ,NPPNP    ,NTYPEL   ,NUMEL    ,NUMNP    ,NZARR
     &          ,NZTRA

      REAL*8::BETAC    ,CREF     ,DENSREF  ,THETAT


      INTEGER*4::IAD_S(MAXNB,NUMNP)  ,IADN_S(NUMNP)
     &          ,IBCOD(NUMNP)        ,IBTCO(NUMNP)
     &          ,ISOZ(NZTRA)         ,KXX(LMXNDL,NUMEL)
     &          ,LDIM(NUMEL)         ,LNNDEL(NUMEL)
     &          ,LTYPE(NUMEL)        ,LXARR(NUMNP)
     &          ,LXPAREL(NUMEL,NPAREL)

      REAL*8::AFLU(NUMEL,IDIMAFLU)          ,ATRA(NUMEL,IDIMATRA)
     &       ,AREA(NUMEL)                   ,ARRC(NUMEL)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL)  ,CAUDAL(NUMNP)
     &       ,CAUX1(NUMNP)                  ,CAUX2(NUMNP)
     &       ,CFLU(NUMEL,IDIMCFLU)          ,COECEL(NUMEL)
     &       ,COORD(NUMNP,3)                ,CONCFLOW(NUMNP)
     &       ,DFLU(NUMEL,IDIMDFLU)          ,DENSITY(NUMEL)
     &       ,DTRA(NUMEL,IDIMDTRA)          ,GP_COORD(6,8,IODIM)
     &       ,GRADLOC(IODIM,LMXNDL,MAXPG)   ,HAUX1(NUMNP)
     &       ,HAUX2(NUMNP)                  ,PAREL(NUMEL,NPPEL)
     &       ,PARNP(NUMNP,NPPNP)            ,POINTWEIGHT(MAXPG,NTYPEL)

C------------------------- Internal

      INTEGER*4::I,IBT

C------------------------- First executable statment

      CONCFLOW(1:NUMNP) = 0D0

C------------------------- Prescribed concentration nodes.

      IF (ANY(IBTCO(1:NUMNP).EQ.1)) THEN

          CALL COMFLOW_PRESC_CONC
     &        (AFLU     ,AREA     ,ARRC     ,ATRA     ,BETAC
     &        ,BUOYANCY ,CAUDAL   ,CAUX1    ,CAUX2    ,CFLU
     &        ,COECEL   ,COORD    ,CREF     ,DFLU     ,DENSITY
     &        ,DENSREF  ,DTRA     ,CONCFLOW ,GP_COORD ,GRADLOC
     &        ,HAUX1    ,HAUX2    ,IAD_S    ,IADN_S   ,IBCOD
     &        ,IBTCO    ,IDIMAFLU ,IDIMATRA ,IDIMCFLU ,IDIMDFLU
     &        ,IDIMDTRA ,IODENS   ,IODIM    ,IORECATRA,IREGFL
     &        ,IREGTR   ,ISOZ     ,ITYPAFLU ,ITYPATRA ,ITYPCFLU
     &        ,ITYPDFLU ,ITYPDTRA ,KXX      ,LDIM     ,LMXNDL
     &        ,LNNDEL   ,LTYPE    ,LXARR    ,LXPAREL  ,MAXNB
     &        ,MAXPG    ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL
     &        ,NUMEL    ,NUMNP    ,NZARR    ,NZTRA    ,PAREL
     &        ,PARNP    ,POINTWEIGHT        ,THETAT)

      END IF ! ANY(IBTCO(1:,NUMNP)).EQ.1


      DO I=1,NUMNP

          IBT = IBTCO(I)

          SELECT CASE (IBT)

C------------------------- Input mass nodes

              CASE(4)

                  CONCFLOW(I) = PARNP(I,4)

C------------------------- Concentration leakage nodes.

              CASE(5) ! Concentration leakage

                  CONCFLOW(I) = PARNP(I,6)*(PARNP(I,4) - CAUX1(I))

          END SELECT !IBT

      END DO !I=1,NUMNP

      END SUBROUTINE COMFLOW_CONC
