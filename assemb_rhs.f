      SUBROUTINE ASSEMB_RHS
     &          (A_MAT    ,B_VEC    ,CAUX2    ,CFLU    ,D_MAT
     &          ,IAD      ,IADN     ,IDIMA_MAT
     &          ,IDIMCFLU ,IDIMD_MAT
     &          ,INDFLTR  ,INDSSTR  ,IODENS   ,IONEWT
     &          ,ITYP_A   ,ITYP_D   ,ITYP_C   ,KXX
     &          ,LNNDEL   ,LMXNDL   ,NUMEL    ,NUMNP   ,THETA   ,TINC
     &          ,VAUX1    ,VAUX2    ,VCALAN)

****************************************************************
*  PURPOSE
*  To calculate the right hand side of the flow and transport equations
*  both for linear and nonlinear systems
*
*  DESCRIPTION
*
*
****************************************************************

      IMPLICIT NONE 

C  EXTERNAL VARIABLES: SCALARS
      INTEGER*4::IDIMA_MAT,IDIMCFLU ,IDIMD_MAT
     &          ,INDFLTR
     &          ,INDSSTR  ,IODENS   ,IONEWT   
     &          ,ITYP_A   ,ITYP_D   ,ITYP_C
     &          ,LMXNDL   ,NUMEL    ,NUMNP

      REAL*8     THETA    ,TINC

C  EXTERNAL VARIABLES: ARRAYS

      INTEGER*4::LNNDEL(NUMEL), KXX(LMXNDL,NUMEL),IAD(*),IADN(*)

      REAL*8::A_MAT(NUMEL,IDIMA_MAT) ,B_VEC(NUMNP)
     &       ,CAUX2(NUMNP) ,CFLU(NUMEL,IDIMCFLU)
     &       ,D_MAT(NUMEL,IDIMD_MAT)
     &       ,VAUX1(NUMNP),VAUX2(NUMNP),VCALAN(NUMNP) 


C------------------------- Assembles stiffness matrix and, if solving a 
C------------------------- transient problem, storage matrix.

C------------------------- If linearization is done by Picard method,
C------------------------- the vector that multiplies A matrix is the
C------------------------- solution in the previous time step and the
C------------------------- factor is 1 minus the time-weighting
C------------------------- parameter.
C------------------------- If Newton's method is used, the vector that
C------------------------- multiplies A matrix is the time-weighted
C------------------------- solution and the factor is equal to 1.

C------------------------- Something similar happens to storage matrix.
C------------------------- If Picard iterations are being used, the
C------------------------- vector that multiplies it is the solution in
C------------------------- the previous time step and the factor the
C------------------------- inverse of the time increment.
C------------------------- If Newton's method is used, the vector that
C------------------------- multiplies A matrix is the solution increment
C------------------------- divided by the time increment and the factor
C------------------------- is equal to 1.


      IF (IONEWT.EQ.0) THEN

C------------------------------Picard's iterations------------------

          CALL PROD_MAT_VEC
     ;(THETA-1D0     ,IAD   ,IADN
     ;,IDIMA_MAT      ,NUMEL          ,NUMNP   
     ;,1             ,ITYP_A         ,1              ,LMXNDL         
     ;,NUMEL         ,NUMNP          ,KXX            ,LNNDEL         
     ;,A_MAT         ,B_VEC          ,VCALAN)

          IF (INDSSTR.EQ.1) THEN

              CALL PROD_MAT_VEC
     ;(1D0/TINC      ,IAD          ,IADN
     ;,IDIMD_MAT     ,NUMEL        ,NUMNP   
     ;,1             ,ITYP_D       ,1                ,LMXNDL         
     ;,NUMEL         ,NUMNP        ,KXX              ,LNNDEL         
     ;,D_MAT         ,B_VEC        ,VCALAN)

          END IF !INDSSTR.EQ.0

      ELSE 

C------------------------------Newton's method------------------

          CALL PROD_MAT_VEC
     ;(-1D0          ,IAD            ,IADN
     ;,IDIMA_MAT     ,NUMEL          ,NUMNP   
     ;,1             ,ITYP_A         ,1              ,LMXNDL         
     ;,NUMEL         ,NUMNP          ,KXX            ,LNNDEL         
     ;,A_MAT         ,B_VEC          ,VAUX1)
     

          IF (INDSSTR.EQ.1) THEN

              CALL PROD_MAT_VEC
     ;(-1D0          ,IAD            ,IADN
     ;,IDIMD_MAT     ,NUMEL          ,NUMNP   
     ;,1             ,ITYP_D         ,1              ,LMXNDL         
     ;,NUMEL         ,NUMNP          ,KXX            ,LNNDEL         
     ;,D_MAT         ,B_VEC          ,VAUX2)


          END IF !INDSSTR.EQ.1

      END IF !IONEWT.EQ.0


C------------------------- If solving a density-dependent, flow
C------------------------- and transient problem, CFLU matrix
C------------------------- has to be assembled
C------------------------- (for both Picard and Newton's method).

      IF (IODENS.EQ.1 .AND. INDFLTR.EQ.0 .AND. INDSSTR.EQ.1) THEN


                  CALL PROD_MAT_VEC 
     ;(-1D0         ,IAD            ,IADN     
     ;,IDIMCFLU     ,NUMEL          ,NUMNP   
     ;,1            ,ITYP_C          ,1              ,LMXNDL         
     ;,NUMEL        ,NUMNP           ,KXX            ,LNNDEL         
     ;,CFLU         ,B_VEC           ,CAUX2)


      END IF !IODENS.EQ.1 .AND. INDFLTR.EQ.0 .AND. INDSSTR.EQ.1

      END SUBROUTINE ASSEMB_RHS
