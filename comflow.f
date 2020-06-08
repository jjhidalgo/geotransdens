      SUBROUTINE COMFLOW
     &          (AFLU     ,AREA     ,BETAC    ,BUOYANCY ,CAUDAL
     &          ,CAUX1    ,CAUX2    ,CFLU     ,COORD    ,CREF
     &          ,DBUOYANCY,DENSITY  ,DENSREF  ,DFLU     ,DFLUDFLU
     &          ,DFLUDTRA ,DPARELDH ,GP_COORD ,GRADLOC  ,GRAVEL
     &          ,GRDFF    ,HAUX1    ,HAUX2    ,IBCOD    ,IDIMAFLU
     &          ,IDIMDFLU ,IFLAGS   ,IODENS   ,IODIM    ,IOFLLI
     &          ,ISOLFL   ,ISOZ     ,KXX      ,LDIM     ,LMXNDL
     &          ,LNNDEL   ,LTYPE    ,LXARR    ,LXPAREL  ,MAXPG
     &          ,NFLAGS   ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL
     &          ,NUMEL    ,NUMNP    ,NZTRA    ,PAREL    ,PARNP
     &          ,POINTWEIGHT        ,THETAT,CONCFLOW,IOCONSRC,IBTCO)

********************************************************************************
*
* PURPOSE
*
*      Computes nodal flow at every node and at the current time
*      (probably at k+theta)
*
* DESCRIPTION
*
*      Computes nodal flow at every node and at the current time
*      (probably at k+theta)
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AFLU                   Matrix of finite elements equations for flow problem
*                         No boundary conditions are included on it.
*  CAUDAL                 Input/output flow at every node.
*  DFLU                   Matrix of finite elements equations for flow
*                         problem related to storage term.
*  HAUX1                  Array containing heads, ponderated by THETAF time
*                         weight. Is equal to THETAF*HCAL+(1-THETAF)*HCALAN
*  HAUX2                  Array containing difference of heads in two
*                         consecutives times. Is equal to HCAL-HCALAN/TIME STEP
*  IBCOD                  Flow boundary condition index
*  PARNP                  Parameter values at every node and current time for
*                         all nodal parameters (each value is computed as the
*                         product of up to four terms:
*                            coeff*zonal value*time funct.*nonl. funct. )
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMDFLU               Used to dimension array DFLU
*  ISOLFL                 If 1, steady flow has been solved. If 2, transient.
*  IOCNSF                 Scheme for storage term in flow problem
*                                 1.Lumped
*                                 2.Consistent.
*                         numbers of two nodes belonging to the same element)
*  NPPNP                  Number of nodal parameters
*  NUMNP                  Number of nodes
*  THETAF                 Time weighting parameter for flow problems
*  TINC                   Current time increment
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Nodal counter
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  MUL_SS
*  MUL_TT
*
* HISTORY
*
*     JCA      4-1998     First coding
*     AMS     12-1998     Revision and modification to include both steady and
*                         transient in one subroutine
*     AMS      1-2002     Revision to simplify call to MUL_SS and MUL_TT 
*     AMS      4-2002     Correction of recharge
*     JHG     10-2003     Modification to include elementwise computation of
*                         AFLU, DFLU and CFLU.
*                         Inclusion of density-dependent flow contribution.
********************************************************************************


      IMPLICIT NONE

      INTEGER*4::IDIMAFLU ,IDIMDFLU ,IOCONSRC ,IODENS   ,IODIM
     &          ,IOFLLI   ,ISOLFL   ,LMXNDL   ,MAXPG    ,NFLAGS
     &          ,NPAREL   ,NPPEL    ,NPPNP    ,NTYPEL   ,NUMEL
     &          ,NUMNP    ,NZTRA
     
      REAL*8::BETAC    ,CREF     ,DENSREF  ,THETAT

      REAL*8::DENS

      INTEGER*4::IBCOD(NUMNP)     ,IBTCO(NUMNP)      ,IFLAGS(NFLAGS)
     &          ,ISOZ(NZTRA)      ,KXX(LMXNDL,NUMEL) ,LDIM(NUMEL)
     &          ,LNNDEL(NUMEL)    ,LTYPE(NUMEL)      ,LXARR(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL)

      REAL*8::AFLU(NUMEL,IDIMAFLU)                ,AREA(NUMEL)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL)        ,CAUDAL(NUMNP)
     &       ,CAUX1(NUMNP)                        ,CAUX2(NUMNP)
     &       ,CFLU(NUMEL,IDIMDFLU)                ,COORD(NUMNP,3)
     &       ,DBUOYANCY(IODIM,LMXNDL*LMXNDL,NUMEL)
     &       ,DENSITY(NUMEL)                      ,DFLU(NUMEL,IDIMDFLU)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL)       ,DPARELDH(NPPEL,NUMEL)
     &       ,GP_COORD(6,8,IODIM)
     &       ,GRADLOC(IODIM,LMXNDL,MAXPG)         ,GRAVEL(NUMEL,3)
     &       ,GRDFF(IODIM,LMXNDL,NUMEL)           ,HAUX1(NUMNP)
     &       ,HAUX2(NUMNP)                        ,PAREL(NUMEL,NPPEL)
     &       ,PARNP(NUMNP,NPPNP)
     &       ,POINTWEIGHT(MAXPG,NTYPEL),concflow(numnp)


C------------------------ Internal Scalars
      INTEGER*4::I

      REAL*8::CNODE    ,DENSNODE ,LEAKFLOW

C------------------------ Initialitation.

      CAUDAL = 0D0

C----------------------------------------------
C------------------------ Prescribed head -----
C----------------------------------------------

      CALL COMFLOW_PRESC_HEAD
     &    (AFLU     ,AREA     ,BETAC    ,BUOYANCY ,CAUDAL
     &    ,CAUX2    ,CFLU     ,COORD    ,CREF     ,DBUOYANCY
     &    ,DENSITY  ,DENSREF  ,DFLU     ,DFLUDFLU ,DFLUDTRA
     &    ,DPARELDH ,GP_COORD ,GRADLOC  ,GRAVEL   ,GRDFF
     &    ,HAUX1    ,HAUX2    ,IBCOD    ,IDIMAFLU ,IDIMDFLU
     &    ,IODENS   ,IODIM    ,IOFLLI   ,ISOLFL   ,ISOZ
     &    ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &    ,LXARR    ,LXPAREL  ,MAXPG    ,NPAREL   ,NPPEL
     &    ,NTYPEL   ,NUMEL    ,NUMNP    ,NZTRA    ,PAREL
     &    ,POINTWEIGHT        ,THETAT)

C------------------------- Now the rest of the boundary conditions.

      DO I=1,NUMNP

C------------------------- First the flow is computed according to
C------------------------- the boundary condition.

          SELECT CASE (IBCOD(I))

C---------------------------------------------
C------------------------ Prescribed head-----
C---------------------------------------------

C------------------------- Contribution of external density is removed.

              CASE (1)

                  IF (IODENS.EQ.1) THEN

C------------------------- Correction for nodes where also concentration
C------------------------- is prescribed when solute sources are
C------------------------- included in the flow equation.
C------------------------- This correction has not to be done again in
C------------------------- the mass balance.

                      IF (IOCONSRC.EQ.1 .AND. IBTCO(I).EQ.1) THEN

                          CAUDAL(I) = CAUDAL(I) - CONCFLOW(I)

	                END IF !IOCONSRC.EQ.1 .AND. IBTCO(I).EQ.1

                      IF (CAUDAL(I).GT.0) THEN

                          CNODE = PARNP(I,4)

                    ELSE

                          CNODE = CAUX1(I)

                    END IF !CAUDAL(I).GT.0

                      DENSNODE = DENS(DENSREF,BETAC,CNODE,CREF)

                      CAUDAL(I) = CAUDAL(I) / DENSNODE
            
                  END IF !IODENS.EQ.1 .AND. CAUDAL(I).GT.0

C----------------------------------------------
C------------------------ Prescribed flow -----
C----------------------------------------------

              CASE(2) 

                  CAUDAL(I) = PARNP(I,2)

C----------------------------------------------
C------------------------ Leakage -------------
C----------------------------------------------

              CASE(3) 

                  CAUDAL(I) = PARNP(I,3)*(PARNP(I,1)-HAUX1(I))

C----------------------------------------------
C------------------------ Leakage and flow ----
C----------------------------------------------

              CASE(4)

                CAUDAL(I) = PARNP(I,2)
     &                     + PARNP(I,3)*(PARNP(I,1)-HAUX1(I))

          END SELECT !IBCOD(INODE)

      END DO !I=1,NNUD


C------------------------- Writing of flow to output file

      IF (IFLAGS(21).EQ.1) THEN

          WRITE(7,*) ' CAUDAL ...'

          DO I=1,NUMNP

              WRITE(7,*) I,CAUDAL(I)

          ENDDO

      END IF !IFLAGS(21).EQ.1

      END SUBROUTINE COMFLOW
