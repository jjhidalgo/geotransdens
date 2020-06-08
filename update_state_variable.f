       SUBROUTINE UPDATE_STATE_VARIABLE
     &           (DELVMAX  ,DELTAVGL  ,DELVMAXOLD,DRELVMX  ,DVITMX   
     &           ,INDFLTR  ,IDELVGL   ,IDELVMAX  ,IOVLIN   ,IOCOUPLED
     &           ,IODENS   ,IONEWT    ,NUMNP     ,ZEROF    ,DELTAITER
     &           ,SOLUTION ,VCALAN    ,VCALIT)

*********************************************************
*
* PURPOSE  Updating the state variable once solved a non-linear system.
*          Calculates maximun (absolute and relatives) changes used when
*          cheking convergence.
*
* ACRONYM UPDATEs the STATE VARIABLE vector
*
* DESCRIPTION 
*     The state variable (head, concentration or the vector containing both variable shuffled)
*  at each node of the finite element grid is increased by a delta, which was obtined after
*  solution of the Newton-Raphson system.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  SOLUTION           Increment of state variable, if Newton method or 
*                     solution in k+1,l+1 if Picard method.
*  VCALAN             Solution in previous time (k).
*  VCALIT             Solution in last iteration (k+1,l).
*
* INTERNAL VARIABLES: ARRAYS
*
*  DELTA_V            Increment of state variable.
*
* EXTERNAL VARIABLES: SCALARS
*
*  DELVMAX            Maximum increment of state variable in current iteration.                                                                      
*  DELVMAXOLD         Maximum increment of state variable in last iteration.
*  DRELVMX            Maximum relative increment of state variable in current iteration.
*  DVITMX             Maximum allowed  increment of state variable per iteration.
*  IONEWT             Lienarization method used to solve non-linear system.
*                       0. Picard.
*                       1. Newton - Raphson.
*  NUMV               Number of nodes.
*
* INTERNAL VARIABLES: SCALARS
*
*  DELTOT             Total increment between last time step and current iteration.
*  DREL               Auxiliar variable to calculate DRELVMX.
*  FACTOR             Multiplicative factor to limitate increment.
*  I                  Counter in DO... END DO statement.
*  X1                 Auxiliar to calculate DELVMAX.
*
* HYSTORY
*
*      G.Galarza  10-1997     First coding
*      JHG        11-2003     Adapted to use any lenght of vector VCAL.
*                             Modificated to be used with solutions 
*                             obtained with Newton - Raphson or Picard
*                             linearization methods.
*
********************************************************************************




      IMPLICIT NONE

      INTEGER*4:: I,INDFLTR,IDELVMAX,IOCOUPLED,IONEWT,IPOS_INI,ISTEP
     &           ,NUMNP,IODENS,IOVLIN,IDELVGL,INDEX

      REAL*8:: DELVMAX,DELVMAXOLD,DELTOT,DREL,DRELVMX,DVITMX,FACTOR,X1,
     &         ZEROF,DELTAVGL

      REAL*8:: SOLUTION((IOCOUPLED+1)*NUMNP), VCALIT(NUMNP)
     &        ,VCALAN(NUMNP)
c     &        ,DELTA_V(NUMNP)
     &        ,DELTAITER((IOCOUPLED+1)*NUMNP)

      real*8,allocatable::delta_v(:)
      allocate(delta_v(numnp))
c----------------------------- initialize the maximum updates	
      DELVMAXOLD=DELVMAX
      DELVMAX=0D0
      DRELVMX=0D0



C----------------------------- In order to localize the relevant elements of the solution
C----------------------------- vector, we check if we are solving a flow or transport system
C----------------------------- or a coupled system. In the last case, the heads and 
C----------------------------- concentrations are mixed, and therefore  separated by one
C----------------------------- position; in the first case, they are not. 

      IF (IOCOUPLED.EQ.1) THEN

          ISTEP=2

          IF (INDFLTR.EQ.0) THEN

              IPOS_INI = 1

          ELSE

              IPOS_INI = 2

          END IF

      ELSE

          IPOS_INI = 1
          ISTEP = 1

      END IF


C------------------------- Computes state variable increment between iterations.


      DO I=1,NUMNP

          IF (IOCOUPLED .EQ.1) THEN

              INDEX = (I-1)*2  + IPOS_INI

          ELSE

              INDEX = I

          ENDIF

          IF (IONEWT.NE.1) THEN

              DELTA_V(I) = SOLUTION(I) - VCALIT(I)

         ELSE

             DELTA_V(I) = SOLUTION(INDEX)

         END IF

      ENDDO

C------------------------- Now we compute : Maximum absolute change in current iteration.
C-------------------------                  Maximum relative change

      DO I = 1,NUMNP

C------------------------- Maximum absolute change in current iteration.

          X1 = DABS(DELTA_V(I))

          IF (X1.GT.DELVMAX)THEN

              DELVMAX=X1
              IDELVMAX=I

          END IF !X1.GT.DELVMAX

C------------------------- Maximum relative change

          DELTOT = VCALIT(I) + DELTA_V(I) - VCALAN(I)

          IF (DELTOT.GT.ZEROF) THEN

              DREL = DABS(DELTA_V(I)/DELTOT)
              IF(DREL.GT.DRELVMX) DRELVMX = DREL

          END IF !DELTOT.GT.ZEROF

      END DO !I=1,NUMV


C------------------------- Maximum permitted change 

      IF (DELVMAX.GT.DVITMX .AND.IONEWT.EQ.1) THEN

          FACTOR=DVITMX/DELVMAX
          DELVMAX=DVITMX

      ELSE

          FACTOR=1D0

      ENDIF

C------------------------- Update the state variable.

      DO I = 1,NUMNP

C------------------------- for newton's method (including for coupled newton)

          IF (IONEWT.GT.0) THEN

              VCALIT(I) = VCALIT(I) + DELTA_V(I)*FACTOR

C------------------------- for picard

          ELSEIF (IONEWT.EQ.0 .AND. IOCOUPLED.EQ.0)THEN

              VCALIT(I) = SOLUTION(I)

          END IF

      ENDDO

C------------------------If we have variable density...


      IF (IODENS.EQ.1) THEN

C------------------------ ...and a nonlinear problem

          IF (IOVLIN .EQ. 1) THEN

              DELTAVGL = 0D0

              DO I=1,NUMNP

c------------------------  ...update the total state variable increments

                  IF (IOCOUPLED.EQ.1) THEN

                      INDEX = (I-1)*2+IPOS_INI

                  ELSE

                      INDEX = I

                  ENDIF

                  DELTAITER(INDEX)=DELTAITER(INDEX)+DELTA_V(I)

c------------------------ ... and find the largest one

                 IF (ABS(DELTAITER(I)).GT.DELTAVGL) THEN

                     DELTAVGL = ABS(DELTAITER(I))
                     IDELVGL = I

                 ENDIF

              ENDDO
C------------------------ ...if the proble is linear, keeps the previously
C------------------------ ...computed increment.

          ELSE
	
              DELTAVGL = DELVMAX

          END IF !IOVLIN .EQ. 1

      END IF !IODENS.EQ.1

      deallocate(delta_v)

      END SUBROUTINE UPDATE_STATE_VARIABLE
