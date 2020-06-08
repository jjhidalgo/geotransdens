      SUBROUTINE PRESC_BC_COUPLED
     &(IDIMA1   ,IDIMA2   ,ISPARSE   ,MAXNN   ,NBAND
     &,NPPNP    ,NUMNP    ,A_DSC     ,B       ,CCALIT  ,HCALIT  
     &,IADD     ,IADN     ,IBCOD     ,IBTCO   ,PARNP) 
***************************************************************
*
*  PURPOSE
*  To implement boundary conditions of fixed head or concentration
*  to the coupled system of flow and transport with newton´s method
*
***************************************************************

      IMPLICIT NONE

C------------------------- EXTERNAL VARIABLES. SCALARS

      INTEGER*4::IDIMA1, IDIMA2,ISPARSE,NUMNP,MAXNN,NBAND,LMDIAG
     &          ,NPPNP

C------------------------- EXTENAL VARIABLES: ARRAYS

      INTEGER*4::IBCOD(NUMNP), IBTCO(NUMNP), IADD(MAXNN),IADN(MAXNN)

      REAL*8::B(NUMNP*2), A_DSC(IDIMA1,IDIMA2),HCALIT(NUMNP)
     &       ,CCALIT(NUMNP),PARNP(NUMNP,NPPNP)

C------------------------- INTERNAL VARIABLES: SCALARS

      INTEGER*4::I, IROW


      
      LMDIAG = 2*NBAND + 2
 
C------------------------- Loop over nodes

      DO I=1,NUMNP

C------------------------- if fixed head boundary

          IF (IBCOD(I).EQ.1) THEN

              IROW = I*2 - 1

C------------------------- set row to zero

              CALL ZERO_ROW(IDIMA1,IDIMA2,ISPARSE,MAXNN,IROW,IADN,A_DSC)
              
              IF (ISPARSE.EQ.0) A_DSC(IROW,LMDIAG) = 1D0
              IF (ISPARSE.EQ.1) A_DSC(IADD(IROW),IROW) = 1D0

              B(IROW) = PARNP(I,1) - HCALIT(I)
     
          END IF !IBCOD(I).EQ.1
           
C------------------------- if fixed concentration boundary

          IF (IBTCO(I).EQ.1) THEN
             
              IROW = I*2

C------------------------- set row to zero

              CALL ZERO_ROW(IDIMA1,IDIMA2,ISPARSE,MAXNN,IROW,IADN,A_DSC)

              IF (ISPARSE.EQ.0) A_DSC(IROW,LMDIAG) = 1D0
              IF (ISPARSE.EQ.1) A_DSC(IADD(IROW),IROW) = 1D0 

              B(IROW) = PARNP(I,4) - CCALIT(I)


          END IF !IBTCO(I).EQ.1

      END DO !I=1,NUMNP

      END SUBROUTINE PRESC_BC_COUPLED



