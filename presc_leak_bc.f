      SUBROUTINE PRESC_LEAK_BC
     &          (A_DSC      ,ALFA     ,BETAC     ,CAUX1    ,CREF
     &          ,DENSREF    ,IAD_S    ,IADD_S    ,IADN_S   ,IBC
     &          ,IDSC_COLS  ,IDSC_ROWS,INALF     ,INDCHANGES
     &          ,INDFLTR    ,INPRESV  ,IOCONSRC  ,IODENS   ,IONEWT
     &          ,ITPTVAR    ,ITYPADSC ,KXX       ,LMXNDL   ,LNNDEL
     &          ,MAXNB      ,NBAND    ,NPPNP     ,NUMEL    ,NUMNP
     &          ,PARNP      ,THETA    ,VAUX1     ,WSPECHEAT)
   
c     ******************************************************************
c     Adds leakeage contribution to LHS.
c
c     ******************************************************************
      IMPLICIT NONE 

C--------------------------- EXTERNAL VARIABLES: INTEGERS

      INTEGER*4::IOCONSRC,IONEWT,NUMNP,NPPNP,INDCHANGES ,LMXNDL,NUMEL
     &          ,IDSC_COLS,IDSC_ROWS,MAXNB  ,NBAND,ITYPADSC,IODENS
     &          ,INALF,INPRESV,INDFLTR,ITPTVAR

      REAL*8::THETA,CREF,BETAC,DENS,DENSREF,WSPECHEAT

C--------------------------- EXTERNAL VARIABLES, ARRAYS

      INTEGER*4::IBC(NUMNP)          ,KXX(LMXNDL,NUMEL) ,IADN_S(NUMNP)
     &          ,IAD_S(MAXNB, NUMNP) ,IADD_S(NUMNP)     ,LNNDEL(NUMEL)

      REAL*8::PARNP(NUMNP,NPPNP)          ,ALFA(NUMNP)  ,VAUX1(NUMNP),
     &        A_DSC(IDSC_ROWS,IDSC_COLS)  ,CAUX1(NUMNP)

      REAL*8,ALLOCATABLE::LEAK_CONT(:)

C--------------------------- INTERNAL VARIABLES, SCALARS

      INTEGER*4::I,IB,IOLEAK
      REAL*8::ALFX     ,FACTOR   ,LEAK_DENS,LEAK_FLOW,THT1     ,VPRES


C--------------------------- Step 1. Initialises array of derivatives
C--------------------------- with the contribution the diagonal of the
C--------------------------- system matrix.
 
      ALLOCATE(LEAK_CONT(NUMNP))

      LEAK_CONT = 0D0
      THT1 = THETA - 1D0

C------------------------- Whe solving heat transport the heat flow
C------------------------- has to be divided by the water specific heat
C------------------------- to be consistent with the rest of the equation.

      IF (INDFLTR.EQ.1 .AND. ITPTVAR.EQ.1) THEN

          IF (IODENS.EQ.1) THEN

              FACTOR = 1D0/WSPECHEAT

          ELSE

C------------------------- If density is constant, water density appears
C------------------------- in this term (see TRANSDENS Guia rapida, 3.7).

              FACTOR = 1D0/(DENSREF*WSPECHEAT)

          END IF

      ELSE

          FACTOR = 1D0

      END IF !ITPTVAR.EQ.1


C--------------------------- Step 2: Begins main loop over nodal points

      DO I=1,NUMNP

          IB = IBC(I)     ! Boundary condition.
          IOLEAK = 0

C--------------------------- Checks if node has leakage boundary condition
C--------------------------- (3,4 if flow, 5 if transport)

          IF ( (INDFLTR.EQ.0 .AND. IB.GE.3) .OR.
     &         (INDFLTR.EQ.1 .AND. IB.EQ.5)) THEN

              IOLEAK = 1

          END IF !INDFLTR.EQ.0 .AND. IB.GE.3 ...

          IF (IOLEAK.EQ.1) THEN

C--------------------------- Step 2.1: Identifies values calculated
C--------------------------- at COMP_PARAM_FLOW.

              ALFX = PARNP(I,INALF)*FACTOR ! Leakage coefficient value
              VPRES = PARNP(I,INPRESV)     ! Prescribed head/conc. value
              LEAK_FLOW = ALFX*(VPRES-VAUX1(I))

C--------------------------- Density of leakage flow is computed (only flow).

              IF (IODENS.EQ.1 .AND. INDFLTR.EQ.0) THEN

C--------------------------- If mixed flow boundary conditions,
C--------------------------- prescribed flow is taken into account
C--------------------------- to compute total flow.

                  IF (IBC(I) .EQ. 4) THEN

                      LEAK_FLOW = LEAK_FLOW + PARNP(I,2)

                  END IF !IBCOD(I). EQ. 4

                  IF (LEAK_FLOW.GT.0D0) THEN

                      LEAK_DENS = DENS(DENSREF,BETAC,PARNP(I,4),CREF)

                  ELSE

                      LEAK_DENS = DENS(DENSREF,BETAC,CAUX1(I),CREF)

                  END IF !LEAK_FLOW.GT.0D0

              ELSE

                  LEAK_DENS = 1.D0

              END IF !IODENS.EQ.1

C--------------------------- Step 2.2.
C--------------------------- Computation of leakage contribution
C--------------------------- to the diagonal of the matrix (LEAK_CONT)
C--------------------------- and to RHS.

              IF (IONEWT.EQ.0) THEN

C--------------------------- Step 2.2.2. Contribution to the
C--------------------------- diagonal of the system matrix.

                  IF (INDCHANGES.EQ.0) THEN 

C--------------------------- If the system matrix does not change,
C--------------------------- (i.e. it is not computed again)
C--------------------------- the contribution of the leakage in the
C--------------------------- previous time step (ALFA) has to be removed.
C--------------------------- So, the contribution to the diagonal
C--------------------------- is Leakage(K+1)-Leakage(K).

                      LEAK_CONT(I) = LEAK_DENS*THETA*ALFX - ALFA(I)

                  ELSE

                      IF (INDFLTR.EQ.0) THEN !FLOW

                          LEAK_CONT(I) = LEAK_DENS*THETA*ALFX

                      ELSE !TRANSPORT

C--------------------------- Constant density or heat transport or
C--------------------------- solute transport with variable density and
C--------------------------- no concentration sources [alf*(w*-w)].

                          IF (IODENS.EQ.0 .OR. ITPTVAR.EQ.1 .OR.
     &           ((IODENS.EQ.1 .AND. ITPTVAR.EQ.0 .AND. IOCONSRC.EQ.0)
     &           )) THEN

                              LEAK_CONT(I) = LEAK_DENS*THETA*ALFX

C--------------------------- Solute transport with variable density and
C--------------------------- concentration sources [alf*(w*-w)*(1-w)].

	                    ELSE IF (IODENS.EQ.1 .AND. ITPTVAR.EQ.0
     &                            .AND. IOCONSRC.EQ.1) THEN

                              LEAK_CONT(I) = LEAK_DENS*THETA
     &                          *ALFX*(1D0 - VAUX1(I))

                          END IF!IODENS.EQ.0 ...

                      END IF !INDFLTR.EQ.0

                  END IF !INDCHANGES


              END IF !IONEWT.EQ.0

C--------------------------- Alfa stored for next time step

              IF (INDFLTR.EQ.0) THEN

                  ALFA(I) = LEAK_DENS*THETA*ALFX

              END IF !INDFLTR.EQ.0

          END IF !IB.GE.3

      END DO !I=1,NUMNP


C--------------------------- Adds leakage contribution to the 
C--------------------------- diagonal of system matrix (A_DSC)

      IF (IONEWT.EQ.0) THEN


          CALL ASSEMBLE      
     &      (1D0       ,1              ,NUMNP     ,IDSC_COLS  ,IDSC_ROWS
     &      ,1         ,ITYPADSC       ,LMXNDL    ,MAXNB      
     &      ,NBAND     ,NUMNP
     &      ,NUMEL     ,LEAK_CONT      ,A_DSC     ,IAD_S        ,IADD_S
     &      ,IADN_S    ,KXX            ,LNNDEL)


      END IF !IONEWT.EQ.0

      DEALLOCATE(LEAK_CONT)

      END SUBROUTINE PRESC_LEAK_BC
