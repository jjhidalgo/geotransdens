      SUBROUTINE UPDATE_VECTORS_AUX
     ;(NUMV     ,TINC    ,THETA     ,VCALAN   
     ;,VCALIT   ,VAUX1  ,VAUX2  )

********************************************************************************
*
* PURPOSE
* Computes auxiliar arrays VAUX1 and VAUX2 after
* system has been solved.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  VAUX1              Array containing heads/concnetration ponderated by THETA.
*  VAUX2              Array containing difference of heads/concentrations in two
*                     consecutives times. Is equal to VCALIT-HCALAN/TIME STEP
*  VCALAN             Solution in previous time (k).
*  VCALIT             Solution in last iteration (k+1,l).
*
* INTERNAL VARIABLES: ARRAYS
*
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  NUMV               Number of nodes.
*  THETA              Time weighting parameter
*  TINC               Current time increment
*
*
* INTERNAL VARIABLES: SCALARS
*
*  I                  Counter in DO... END DO statement.
*
*
* HISTORY
*      L. J. Slooten  09-2003     First coding
*      JHG            11-2003     Adapted to use any length of vectors VACALAN, etc.
*                                 Replacement of VACAL with VACALIT (due to changes
*                                 in UPDATE_STATE_VARIABLE).
*                                 Elimination of loop through problems.
********************************************************************************

      IMPLICIT NONE
      
      INTEGER*4::I, NUMV

      REAL*8 ::TINC, THETA

      REAL*8 VCALAN(NUMV),VCALIT(NUMV),VAUX1(NUMV),
     &       VAUX2(NUMV)

           
      DO I=1,NUMV

          VAUX1(I) = (1.D0 - THETA)*VCALAN(I) + THETA*VCALIT(I)

          VAUX2(I)=(VCALIT(I) - VCALAN(I)) / TINC

      END DO
      
      RETURN
      END SUBROUTINE UPDATE_VECTORS_AUX
      