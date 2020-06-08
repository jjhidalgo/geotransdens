      SUBROUTINE CALC_DENS
     &          (BETAC    ,CCALIT   ,CREF    ,DENSREF  ,DENSITY
     &          ,IDIMDENS ,IOVRWC   ,KXX     ,LMXNDL   ,LNNDEL
     &          ,NUMEL    ,NUMNP)
**************************************************************
* PURPOSE:
*
* Manages the computation of the density elementwise or nodewise.
*  
*
* DESCRIPTION
*
* Manages the computation of the density elementwise or nodewise by calling
* the routines that use density by nodes or elements as needed.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CCALIT                 Computed concentrations in the current iteration.    
*  DENSITY                Array containing the density of every element.
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  BETAC                Derivative of density w.r.t. concentration divided by density.
*  CREF                 Reference concentration.
*  DENSREF              Reference density.
*  IOVRWC                   Equal to IOPTS(31).
*                           1.WATVOL calculated elementwise
*                           2.WATVOL calculated nodewise.
*                           DENSITY is calculate in the same way that WATVOL.
*  LMXNDL                 Maximum number of nodes per element                   
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  CALC_DENS_I            Computes DENSITY by nodes.
*  CALC_DENS_L            Computes DENSITY by elements.
*
* HISTORY
*
*     JHG         10-2003    First coding.
*
*******************************************************************************


      IMPLICIT NONE
      
      INTEGER*4::IDIMDENS,IOVRWC, LMXNDL, NUMEL, NUMNP
      INTEGER*4::KXX(LMXNDL,NUMEL), LNNDEL(NUMEL)

      REAL*8::BETAC, CREF, DENSREF

      REAL*8::CCALIT(NUMNP)
      REAL*8::DENSITY(IDIMDENS)

c-parche-provisional-hasta-que-la-densidad-vaya-por-nudos
c      IF (IOVRWC.LE.1) THEN

          CALL CALC_DENS_L
     &        (BETAC    ,CCALIT   ,CREF    ,DENSREF
     &        ,DENSITY  ,KXX      ,LMXNDL  ,LNNDEL
     &        ,NUMEL    ,NUMNP)

c      ELSE

c          CALL CALC_DENS_I
c    &         (BETAC    ,CCALIT   ,CREF    ,DENSREF
c    &         ,DENSITY  ,NUMNP)

c      END IF

      END SUBROUTINE CALC_DENS
**************************************************************
**************************************************************




      SUBROUTINE CALC_DENS_L
     &          (BETAC    ,CCALIT   ,CREF    ,DENSREF
     &          ,DENSITY  ,KXX      ,LMXNDL  ,LNNDEL
     &          ,NUMEL    ,NUMNP)

**************************************************************
* PURPOSE:
*
* To calculate the density elementwise.
*
* DESCRIPTION
*
* Density is calculated with the mean concentration 
* of the element.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CCALIT                 Computed concentrations in the previous time step.    
*  DENSITY                Array containing the density of every element.
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  BETAC                  Derivative of density w.r.t. concentration divided by density.
*  CREF                   Reference concentration.
*  DENSREF                Reference density.
*  LMXNDL                 Maximum number of nodes per element
*  NUMEL                  Number of elements
*  NUMNP                  Number of nodes
*
*
* INTERNAL VARIABLES: SCALARS
*
*  INODE                    Node index.
*  L                        Element index.
*  MEANCONC                 Mean concentration in the element.
*  NNUD                     Number of nodes in current element.
*  NODE                     Current node.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  DENS                    Function to compute density.
*
* HISTORY
*
*    JHG         10-2003    First coding.

**************************************************************
      IMPLICIT NONE

      INTEGER*4::INODE, L, LMXNDL, NUMEL, NUMNP, NNUD, NODE
      INTEGER*4::LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)

      REAL*8::BETAC, CREF, DENSREF, MEANCONC
      REAL*8::CCALIT(NUMNP),DENSITY(NUMEL)

      REAL*8::DENS

C------------------------- For each element
      DO L=1,NUMEL
          NNUD=LNNDEL(L)
          MEANCONC=0D0
C------------------------- the mean concentration is calculated
          DO INODE=1,NNUD
              NODE=KXX(INODE,L)
              MEANCONC=MEANCONC+CCALIT(NODE)
          ENDDO

          MEANCONC=MEANCONC/NNUD

C------------------------- and the density evaluated with the mean concetration value.

          DENSITY(L) = DENS(DENSREF,BETAC,MEANCONC,CREF)


      ENDDO !Next element.


      END SUBROUTINE CALC_DENS_L





**************************************************************
**************************************************************


      SUBROUTINE CALC_DENS_I
     &          (BETAC    ,CCALIT   ,CREF    ,DENSREF  ,DENSITY  ,NUMNP)


**************************************************************
* PURPOSE
*
* To calculate the density nodewise.
*
* DESCRIPTION
*
* Density is calculated as the density of the mean concentration 
* of the element.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CCALIT                 Computed concentrations in the previous time step.    
*  DENSITY                Array containing the density of every element.
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  BETAC                    Derivative of density w.r.t. concentration divided by density.
*  CREF                    Reference concentration.
*  DENSREF                Reference density.
*  NUMNP                  Number of nodes                                       
*
*
* INTERNAL VARIABLES: SCALARS
*
*  INODE                    Node index.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  DENS                    Function to compute density.
*
* HISTORY
*
*    JHG         10-2003    First coding.

**************************************************************
      IMPLICIT NONE

      INTEGER*4::INODE, NUMNP

      REAL*8::BETAC, CREF, DENSREF
      REAL*8::CCALIT(NUMNP),DENSITY(NUMNP)
      REAL*8::DENS

C------------------------- For each node
      DO INODE=1,NUMNP

C------------------------- Density is calculated with node's concentration.         

          DENSITY(INODE) = DENS(DENSREF,BETAC,CCALIT(INODE),CREF)

      ENDDO !Next node.

      RETURN
      END SUBROUTINE CALC_DENS_I



**************************************************************
**************************************************************


      FUNCTION DENS(DEN_REF,BETAC,C,C_REF) RESULT (DEN)

**************************************************************
* PURPOSE
*
* To calculate the density.
*  
* DESCRIPTION
*
* Density is calculated according to the formula
*
*            D = D_0 * EXP (B_c* (C_-C_0))
*
* ARGUMENTS
*
*  BETAC                    Derivative of density w.r.t. concentration divided by density.
*  C                        Cocncentration.
*  C_REF                    Reference concentration.
*  DEN_REF                  Referemce density.
*
*
* HISTORY
*
*    JHG         10-2003    First coding.
*
**************************************************************
      IMPLICIT NONE

      REAL*8::DEN, DEN_REF, BETAC, C, C_REF

      DEN=DEN_REF*DEXP(BETAC*(C-C_REF))

      END FUNCTION DENS
